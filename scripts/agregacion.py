#!/usr/bin/env python3
"""Numero de agregacion de CTAB a partir de una trayectoria GROMACS.

Dos moleculas pertenecen al mismo agregado si alguna pareja de sus cuentas de cola
(T2-T5) esta a una distancia <= --cutoff (nm), con condiciones periodicas. Los
agregados son las componentes conexas de ese grafo (equivale a Hoshen-Kopelman).

Criterio de la referencia: Perez-Sanchez et al. 2016 usan colas a <= 0.76 nm
(Martini 2.2). Stavert et al. 2025 (Martini 3) no se verifico en el texto leido;
por eso se ofrece --sweep para ver la sensibilidad al corte. Validar el corte con
la g(r) cola-cola antes de dar el numero como definitivo.

Definiciones (todas se reportan, el articulo no dice cual usa):
  Nn   = media numerica del tamano de agregado (incluye monomeros libres)
  Nw   = media ponderada en peso, sum(s^2)/sum(s)
  *_mic = lo mismo, pero solo agregados con tamano >= --min-micela
  f_libre = fraccion de CTAB en agregados de tamano < --min-micela

Uso:
  python3 scripts/agregacion.py --top npt.gro --traj md.xtc --begin 200 --out agg.csv
"""
import argparse
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree

TAIL = ("T2", "T3", "T4", "T5")


def tamanos(pos, box, mol, n_mol, cutoff):
    """Tamanos de agregado (array de enteros) para un fotograma. pos y box en nm."""
    pos = np.mod(pos, box)
    tree = cKDTree(pos, boxsize=box)
    pares = tree.query_pairs(cutoff, output_type="ndarray")
    if len(pares) == 0:
        return np.ones(n_mol, dtype=int)
    a, b = mol[pares[:, 0]], mol[pares[:, 1]]
    k = a != b
    a, b = a[k], b[k]
    g = coo_matrix((np.ones(len(a), dtype=np.int8), (a, b)), shape=(n_mol, n_mol))
    _, lab = connected_components(g, directed=False)
    return np.bincount(lab)


def metricas(s, min_mic):
    n_tot = s.sum()
    mic = s[s >= min_mic]
    out = {
        "n_agregados": len(s),
        "Nn": s.mean(),
        "Nw": (s.astype(float) ** 2).sum() / n_tot,
        "n_micelas": len(mic),
        "Nn_mic": mic.mean() if len(mic) else np.nan,
        "Nw_mic": (mic.astype(float) ** 2).sum() / mic.sum() if len(mic) else np.nan,
        "f_libre": s[s < min_mic].sum() / n_tot,
        "max": s.max(),
    }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--top", required=True, help="topologia con nombres (npt.gro)")
    ap.add_argument("--traj", help="trayectoria (md.xtc); si falta, usa solo --top")
    ap.add_argument("--cutoff", type=float, default=0.76, help="corte cola-cola en nm (0.76)")
    ap.add_argument("--min-micela", type=int, default=10, help="tamano minimo para contar como micela (10)")
    ap.add_argument("--begin", type=float, default=0.0, help="ns desde donde analizar (0)")
    ap.add_argument("--end", type=float, default=None, help="ns hasta donde analizar")
    ap.add_argument("--step", type=int, default=1, help="usar 1 de cada k fotogramas (1)")
    ap.add_argument("--blocks", type=int, default=5, help="bloques para el error (5)")
    ap.add_argument("--sweep", action="store_true", help="repetir con varios cortes y mostrar sensibilidad")
    ap.add_argument("--out", default="agregacion.csv")
    a = ap.parse_args()

    import MDAnalysis as mda

    u = mda.Universe(a.top, a.traj) if a.traj else mda.Universe(a.top)
    sel = u.select_atoms("resname CTAB and name " + " ".join(TAIL))
    if len(sel) == 0:
        sys.exit("no se encontraron cuentas de cola de CTAB (resname CTAB, nombres T2-T5)")
    _, mol = np.unique(sel.resindices, return_inverse=True)
    n_mol = mol.max() + 1
    cortes = [0.6, 0.7, 0.76, 0.85, 1.0] if a.sweep else [a.cutoff]
    if a.sweep and a.cutoff not in cortes:
        cortes.append(a.cutoff)

    filas = {c: [] for c in cortes}
    tiempos = []
    for ts in u.trajectory[::a.step]:
        t_ns = ts.time / 1000.0
        if t_ns < a.begin or (a.end is not None and t_ns > a.end):
            continue
        box = ts.dimensions[:3] / 10.0
        pos = sel.positions / 10.0
        tiempos.append(t_ns)
        for c in cortes:
            filas[c].append(metricas(tamanos(pos, box, mol, n_mol, c), a.min_micela))
    if not tiempos:
        sys.exit("ningun fotograma en el rango pedido")

    cols = ["n_agregados", "Nn", "Nw", "n_micelas", "Nn_mic", "Nw_mic", "f_libre", "max"]
    with open(a.out, "w") as f:
        f.write("cutoff_nm,t_ns," + ",".join(cols) + "\n")
        for c in cortes:
            for t, m in zip(tiempos, filas[c]):
                f.write(f"{c},{t:.3f}," + ",".join(f"{m[k]:.4f}" for k in cols) + "\n")

    print(f"{n_mol} moleculas CTAB, {len(tiempos)} fotogramas ({tiempos[0]:.1f}-{tiempos[-1]:.1f} ns), "
          f"min-micela={a.min_micela}")
    for c in cortes:
        print(f"\ncorte {c} nm" + ("  <- principal" if c == a.cutoff else ""))
        for k in ("Nn", "Nw", "Nn_mic", "Nw_mic", "n_micelas", "f_libre"):
            x = np.array([m[k] for m in filas[c]], dtype=float)
            bl = [b for b in np.array_split(x, min(a.blocks, len(x))) if len(b)]
            medias = [np.nanmean(b) for b in bl]
            err = np.std(medias, ddof=1) / np.sqrt(len(medias)) if len(medias) > 1 else float("nan")
            print(f"  {k:10s} {np.nanmean(x):9.2f} +- {err:6.2f}")
    print(f"\nCSV por fotograma: {a.out}")


if __name__ == "__main__":
    main()
