#!/bin/sh
# Espera a que no quede ningún "gmx" corriendo y prolonga (o reanuda) réplicas desde su punto de control.
# Uso: scripts/extender.sh <wt> <extender_ps|0> <semilla> [<semilla> ...]
#   <extender_ps> = ps a sumar al md.tpr (1500000 = 1.5 us); 0 = reanudar con el md.tpr tal cual
#   (caso de la semilla 1, cuyo md.tpr ya es de 3 us).
# Ejemplo: nohup ./scripts/extender.sh 3 1500000 2 3 > extender.out 2>&1 &
# WAIT=0 omite la espera. Cada trabajo usa -ntmpi 1 -ntomp 10 (dos por GPU).
set -e
cd "$(dirname "$0")/.."
wt=$1; ext=$2
[ -n "$wt" ] && [ -n "$ext" ] && [ $# -ge 3 ] || { echo "uso: $0 <wt> <extender_ps|0> <semilla> [...]"; exit 1; }
shift 2
if [ "${WAIT:-1}" != "0" ]; then
  echo "$(date) esperando a que terminen los gmx activos"
  while pgrep -x gmx >/dev/null 2>&1; do sleep 60; done
fi
for seed in "$@"; do
  d=S0a/ctab_${wt}wt_r${seed}
  [ -f "$d/md.cpt" ] || { echo "falta $d/md.cpt"; exit 1; }
  if [ "$ext" != "0" ]; then
    [ -f "$d/md_ext.tpr" ] && { echo "$d ya tiene md_ext.tpr; no se extiende dos veces"; exit 1; }
    gmx convert-tpr -s "$d/md.tpr" -extend "$ext" -o "$d/md_ext.tpr" > "$d/convert.out" 2>&1
    tpr=md_ext.tpr
  else
    tpr=md.tpr
  fi
  (cd "$d" && nohup gmx mdrun -ntmpi 1 -ntomp 10 -s "$tpr" -cpi md.cpt -deffnm md >> md_ext.out 2>&1 &)
  echo "$(date) lanzada $d con $tpr"
  sleep 5
done
