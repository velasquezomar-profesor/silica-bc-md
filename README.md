# silica-bc-md: piloto CG (Martini 3) siguiendo a Stavert et al.

Etapa S0a: validación del flujo con CTAB. Objetivo (Stavert, Mol. Simul. 51 (2025) 188): número de agregación medio ≈ 125 a 3 % p/p y ≈ 167 a 6 % p/p.

## Uso

    ./scripts/get_martini3.sh
    ./scripts/make_S0a.sh
    cd S0a/ctab_6wt_r1 && ./run.sh

`get_martini3.sh` descarga los 3 archivos de Martini 3.0.0 (16 MB, commit fijado y verificado con sha256; no se versionan). `make_S0a.sh` regenera los `conf.gro` (no se versionan; son deterministas con la semilla). Para réplicas, añadir semillas en el bucle del script.

`run.sh` corre minimización, NVT (2 fs), NPT (2 fs, c-rescale) y producción (30 fs).

## Del artículo (Tablas 1-2 y métodos)

- CTA⁺: cabeza H1 = Q2, cola T2–T5 = C1. Enlaces H1–T2 0.375 nm/39 000; T2–T3 y T3–T4 0.478/3450. Ángulos H1–T2–T3 156°/180; T2–T3–T4 y T3–T4–T5 150°/100.
- Producción: dt 30 fs, leap-frog, v-rescale, Parrinello-Rahman, corte Verlet 1.1 nm, reaction-field con εr = 15.
- Equilibración: minimización steepest-descent, NVT con v-rescale a 2 fs, NpT a 2 fs.

## Supuestos míos (no vienen en el texto leído; revisar)

1. Enlace T4–T5 = 0.483 nm / 3450: la Tabla 1 extraída del PDF dice "T3–T5", que probablemente es una errata (cadena lineal).
2. Ángulos con funct 2 (coseno armónico), el estándar de Martini.
3. Bromuro: bead `SQ4` del archivo de iones de Martini 3 (`martini_v3.0.0_ions_v1.itp`). El artículo no indica qué modelo usa.
4. Sin bead anticongelante (en Martini 3 no se usa por defecto).
5. Temperatura 298 K. El artículo da 298 K para los atomísticos; la del CG no aparece en el texto leído.
6. Duraciones de equilibración (NVT 100 ps, NPT 200 ps) y barostato de equilibración C-rescale (el Mol. Simul. dice Berendsen; el RSC dice C-rescale).
7. Configuración inicial: retícula cúbica con cadenas de CTAB colocadas al azar sobre ejes de la retícula; la minimización la relaja. Si el empaquetado produce solapes, cambiar de semilla.
8. "% p/p" = masa de CTAB (364.45 g/mol, con Br) sobre masa total, con un bead W = 4 H₂O (72.06 g/mol). El artículo no define la base exacta.

## Pendiente

- Medir el número de agregación (Hoshen-Kopelman; en Pérez-Sánchez 2016 las colas a ≤ 0.76 nm) y compararlo con 125/167.
- Parametrizar BC (etapa A) siguiendo el mismo procedimiento.

## Estado de validación (2026-11-10)

Corrección: la primera versión calculaba mal los pasos de producción (100 000 pasos = 3 ns en vez de 3 µs). Corregido: 100 000 000 pasos a 30 fs. Con la versión errónea, 6 % p/p (85 184 beads) dio 1822 ns/día en una RTX 2000 Ada con GROMACS 2026.3 (medida corta de 3 ns, posible contención con otro trabajo en la misma GPU al inicio; provisional).


Probado con GROMACS 2025.6-dev (solo CPU, sin GPU): `grompp` acepta los dos sistemas; la minimización del sistema de 6 % p/p llegó a Fmax ≈ 298 en 5000 pasos (el límite ahora es 20 000); NVT de 2000 pasos a 2 fs sin errores (T ≈ 259 K a los 4 ps, aún relajando). No se ha corrido ninguna producción ni medido el número de agregación.
