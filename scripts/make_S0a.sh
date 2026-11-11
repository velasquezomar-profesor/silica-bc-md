#!/bin/sh
# Genera los sistemas de S0a. Por defecto: 3 y 6 % p/p, semilla 1, 3 us (corrida original).
# Réplicas: WTS="3 6" SEEDS="2 3" US=1.5 ./scripts/make_S0a.sh
set -e
cd "$(dirname "$0")/.."
WTS=${WTS:-"3 6"}
SEEDS=${SEEDS:-"1"}
US=${US:-3}
for wt in $WTS; do
  for seed in $SEEDS; do
    python3 scripts/build_ctab.py --wt $wt --n 1000 --seed $seed --out S0a/ctab_${wt}wt_r${seed} --us $US
  done
done
