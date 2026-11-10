#!/usr/bin/env python3
"""Construye un sistema CTAB/agua en Martini 3 siguiendo Stavert et al., Mol. Simul. 51 (2025) 188.

Parametros del tensioactivo: Tablas 1 y 2 del articulo. Supuestos propios (no dados en el texto
leido) estan marcados con SUPUESTO en comentarios y en el README.
"""
import argparse
import math
import os
import shutil
import sys

import numpy as np

MW_CTAB = 364.45          # g/mol, incluye Br
MW_WBEAD = 72.06          # 4 H2O por bead W de Martini 3
RHO_W = 8.357             # beads W por nm3 a 1 g/cm3
A_LAT = RHO_W ** (-1 / 3)  # 0.4929 nm

CTAB_ITP = """; CTA+ Martini 3, Stavert et al. Mol. Simul. 51 (2025) 188, Tablas 1-2
; T2-T5 = cola (C1), H1 = cabeza (Q2)
; SUPUESTO: la Tabla 1 del texto extraido lista T3-T5 (0.483); aqui se usa T4-T5 (cadena lineal).
; SUPUESTO: angulos funct 2 (coseno armonico, estandar de Martini); el articulo da solo theta y k.
[ moleculetype ]
CTAB 1

[ atoms ]
1 Q2 1 CTAB H1 1  1.0 72.0
2 C1 1 CTAB T2 2  0.0 72.0
3 C1 1 CTAB T3 3  0.0 72.0
4 C1 1 CTAB T4 4  0.0 72.0
5 C1 1 CTAB T5 5  0.0 72.0

[ bonds ]
1 2 1 0.375 39000
2 3 1 0.478 3450
3 4 1 0.478 3450
4 5 1 0.483 3450

[ angles ]
1 2 3 2 156 180
2 3 4 2 150 100
3 4 5 2 150 100
"""

MDP_COMMON = """cutoff-scheme            = Verlet
coulombtype              = reaction-field
rcoulomb                 = 1.1
epsilon_r                = 15
epsilon_rf               = 0
vdw_type                 = cutoff
vdw-modifier             = Potential-shift-verlet
rvdw                     = 1.1
verlet-buffer-tolerance  = 0.005
tcoupl                   = v-rescale
tc-grps                  = System
tau_t                    = 1.0
ref_t                    = {T}
nstcomm                  = 100
comm-mode                = Linear
"""

MDP_EM = """integrator = steep
emtol      = 100
emstep     = 0.01
nsteps     = 20000
""" + MDP_COMMON.replace("tcoupl                   = v-rescale\n", "").replace("tc-grps                  = System\ntau_t                    = 1.0\nref_t                    = {T}\n", "")

MDP_NVT = """integrator = md
dt         = 0.002
nsteps     = {nsteps}
nstxout-compressed = 0
nstenergy  = 500
nstlog     = 500
gen_vel    = yes
gen_temp   = {T}
gen_seed   = {seed}
pcoupl     = no
""" + MDP_COMMON

MDP_NPT = """integrator = md
dt         = 0.002
nsteps     = {nsteps}
nstenergy  = 500
nstlog     = 500
continuation = yes
pcoupl     = c-rescale
pcoupltype = isotropic
tau_p      = 4.0
compressibility = 3e-4
ref_p      = 1.0
""" + MDP_COMMON

MDP_MD = """integrator = md
dt         = 0.030
nsteps     = {nsteps}
nstxout-compressed = {nstx}
compressed-x-grps  = System
nstenergy  = 5000
nstlog     = 5000
continuation = yes
pcoupl     = parrinello-rahman
pcoupltype = isotropic
tau_p      = 12.0
compressibility = 3e-4
ref_p      = 1.0
""" + MDP_COMMON


def build(wt, n_ctab, seed, outdir, T, prod_us):
    rng = np.random.default_rng(seed)
    n_w_target = n_ctab * MW_CTAB * (1 / (wt / 100) - 1) / MW_WBEAD
    n_sites_target = n_w_target + 6 * n_ctab
    n = round(n_sites_target ** (1 / 3))
    sites = n ** 3
    n_w = sites - 6 * n_ctab
    wt_real = 100 * n_ctab * MW_CTAB / (n_ctab * MW_CTAB + n_w * MW_WBEAD)
    box = n * A_LAT
    occ = np.zeros((n, n, n), dtype=bool)
    chains = []
    dirs = [np.array(v) for v in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    for _ in range(n_ctab):
        for _try in range(10000):
            i = rng.integers(0, n, 3)
            d = dirs[rng.integers(0, 3)] * rng.choice([-1, 1])
            cells = [tuple((i + k * d) % n) for k in range(5)]
            if not any(occ[c] for c in cells):
                for c in cells:
                    occ[c] = True
                chains.append(cells)
                break
        else:
            sys.exit("no se pudo colocar una cadena; baje la concentracion")
    free = np.argwhere(~occ)
    rng.shuffle(free)
    br = free[:n_ctab]
    water = free[n_ctab:]
    assert len(water) == n_w, (len(water), n_w)
    os.makedirs(outdir, exist_ok=True)
    lines = []
    atom = 0
    res = 0

    def add(resname, name, idx):
        nonlocal atom
        atom += 1
        x = (np.array(idx) + 0.5) * A_LAT
        lines.append(f"{res % 100000:5d}{resname:<5s}{name:>5s}{atom % 100000:5d}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}")

    names = ["H1", "T2", "T3", "T4", "T5"]
    for cells in chains:
        res += 1
        for nm, c in zip(names, cells):
            add("CTAB", nm, c)
    for c in br:
        res += 1
        add("BR", "BR", c)
    for c in water:
        res += 1
        add("W", "W", c)
    with open(os.path.join(outdir, "conf.gro"), "w") as f:
        f.write(f"CTAB {wt_real:.2f} wt% N={n_ctab} seed={seed}\n{atom}\n")
        f.write("\n".join(lines) + "\n")
        f.write(f"{box:10.5f}{box:10.5f}{box:10.5f}\n")
    with open(os.path.join(outdir, "ctab.itp"), "w") as f:
        f.write(CTAB_ITP)
    with open(os.path.join(outdir, "topol.top"), "w") as f:
        f.write('#include "../../ff/martini_v3.0.0.itp"\n#include "../../ff/martini_v3.0.0_ions_v1.itp"\n'
                '#include "../../ff/martini_v3.0.0_solvents_v1.itp"\n#include "ctab.itp"\n\n'
                f"[ system ]\nCTAB {wt_real:.2f} wt%\n\n[ molecules ]\nCTAB {n_ctab}\nBR {n_ctab}\nW {n_w}\n")
    nsteps_md = round(prod_us * 1e6 / 0.030)  # us -> ps, dt = 0.030 ps
    fmt = dict(T=T, seed=seed)
    for name, tpl, extra in (
        ("em", MDP_EM, {}),
        ("nvt", MDP_NVT, {"nsteps": 50000}),      # SUPUESTO: 100 ps a 2 fs
        ("npt", MDP_NPT, {"nsteps": 100000}),     # SUPUESTO: 200 ps a 2 fs
        ("md", MDP_MD, {"nsteps": nsteps_md, "nstx": 33333}),  # frame cada ~1 ns
    ):
        with open(os.path.join(outdir, name + ".mdp"), "w") as f:
            f.write(tpl.format(**fmt, **extra))
    with open(os.path.join(outdir, "run.sh"), "w") as f:
        f.write("""#!/bin/sh
set -e
gmx grompp -f em.mdp -c conf.gro -p topol.top -o em.tpr -maxwarn 0
gmx mdrun -deffnm em
gmx grompp -f nvt.mdp -c em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt
gmx grompp -f npt.mdp -c nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md
""")
    os.chmod(os.path.join(outdir, "run.sh"), 0o755)
    print(f"{outdir}: {atom} beads, caja {box:.2f} nm, {wt_real:.2f} wt% (objetivo {wt}), agua {n_w}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--wt", type=float, required=True)
    ap.add_argument("--n", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", required=True)
    ap.add_argument("--T", type=float, default=298.0)   # SUPUESTO: Stavert da 298 K para AA; T del CG no aparece en el texto leido
    ap.add_argument("--us", type=float, default=3.0)
    a = ap.parse_args()
    build(a.wt, a.n, a.seed, a.out, a.T, a.us)
