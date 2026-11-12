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

## Análisis: número de agregación

    python3 -m pip install MDAnalysis scipy
    python3 scripts/agregacion.py --top S0a/ctab_6wt_r1/npt.gro --traj S0a/ctab_6wt_r1/md.xtc --begin 200 --out agg_6wt.csv
    python3 scripts/agregacion.py --top S0a/ctab_3wt_r1/npt.gro --traj S0a/ctab_3wt_r1/md.xtc --begin 200 --out agg_3wt.csv --sweep

Dos moléculas están en el mismo agregado si alguna pareja de cuentas de cola (T2-T5) está a ≤ 0.76 nm (criterio de Pérez-Sánchez 2016, Martini 2.2), con condiciones periódicas; los agregados son las componentes conexas. Reporta Nn (numérica), Nw (en peso), las mismas solo para agregados ≥ 10 (`--min-micela`) y la fracción libre, con error por bloques. `--sweep` repite con cortes de 0.6 a 1.0 nm para ver la sensibilidad.

Límites: el corte de 0.76 nm viene de Martini 2.2 y no se verificó en Stavert (Martini 3); el artículo no dice si su 125/167 es Nn o Nw. Comparar con ambas y mirar el `--sweep`. Probado solo con un sistema sintético (8 micelas de 60 y 20 monómeros, una micela cruzando la caja): recupera 60 y 8. No se ha corrido sobre las trayectorias reales.

## Réplicas con semillas adicionales (2026-11-11)

En la corrida original (semilla 1) el tamaño de las micelas se congela hacia los 600-800 ns: 11 micelas en 3 % (Nn 90.9, Nw 104.4) y 9 en 6 % (Nn 111.1, Nw 131.5), sin cambios después. Con 1000 cadenas, N = 1000/(número de micelas), así que una sola corrida es una realización de un proceso estocástico (125 = 1000/8 y 167 ≈ 1000/6 en el artículo). Para ver la distribución se añaden réplicas de 1.5 µs:

    git pull
    WTS="3 6" SEEDS="2 3" US=1.5 ./scripts/make_S0a.sh
    ./scripts/lanzar_replicas.sh 6 2 3     # en la máquina de 6 %: 2 réplicas a la vez en la GPU
    ./scripts/lanzar_replicas.sh 3 2 3     # en la máquina de 3 %

`lanzar_replicas.sh` pone `MDRUN_FLAGS="-ntmpi 1 -ntomp 10"` (con GPU, GROMACS exige `-ntmpi` si se fija `-ntomp`) para repartir los 20 hilos entre las dos réplicas (`run.sh` lo usa si existe). No relanza si ya hay dos `mdrun` activos (usar `FORCE=1` para saltarlo). Analizar con `scripts/agregacion.py --begin 800` (el reparto se congela antes).

## Pendiente

- Correr `scripts/agregacion.py` sobre las trayectorias de S0a y compararlo con 125/167.
- Parametrizar BC (etapa A) siguiendo el mismo procedimiento.

## Estado de validación (2026-11-10)

Corrección: la primera versión calculaba mal los pasos de producción (100 000 pasos = 3 ns en vez de 3 µs). Corregido: 100 000 000 pasos a 30 fs. Con la versión errónea, 6 % p/p (85 184 beads) dio 1822 ns/día en una RTX 2000 Ada con GROMACS 2026.3 (medida corta de 3 ns, posible contención con otro trabajo en la misma GPU al inicio; provisional).


Probado con GROMACS 2025.6-dev (solo CPU, sin GPU): `grompp` acepta los dos sistemas; la minimización del sistema de 6 % p/p llegó a Fmax ≈ 298 en 5000 pasos (el límite ahora es 20 000); NVT de 2000 pasos a 2 fs sin errores (T ≈ 259 K a los 4 ps, aún relajando). No se ha corrido ninguna producción ni medido el número de agregación.
