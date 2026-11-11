#!/bin/sh
# Lanza varias réplicas a la vez en la misma GPU, 10 hilos OpenMP cada una.
# Uso: scripts/lanzar_replicas.sh 6 2 3      (6 % p/p, semillas 2 y 3)
# Antes: WTS=6 SEEDS="2 3" US=1.5 ./scripts/make_S0a.sh
set -e
cd "$(dirname "$0")/.."
wt=$1; shift
[ -n "$wt" ] && [ $# -gt 0 ] || { echo "uso: $0 <wt> <semilla> [<semilla> ...]"; exit 1; }
for seed in "$@"; do
  d=S0a/ctab_${wt}wt_r${seed}
  [ -f "$d/conf.gro" ] || { echo "falta $d/conf.gro (correr make_S0a.sh con SEEDS)"; exit 1; }
  if pgrep -f "mdrun.*-deffnm" >/dev/null && [ "$(pgrep -fc 'mdrun.*-deffnm')" -ge 2 ] && [ -z "$FORCE" ]; then
    echo "ya hay 2 o más mdrun activos; use FORCE=1 para lanzar de todos modos"; exit 1
  fi
  (cd "$d" && MDRUN_FLAGS="-ntomp 10" nohup ./run.sh > run.out 2>&1 &)
  echo "lanzada $d"
  sleep 2
done
