#!/bin/sh
set -e
cd "$(dirname "$0")/.."
for wt in 3 6; do
  for seed in 1; do
    python3 scripts/build_ctab.py --wt $wt --n 1000 --seed $seed --out S0a/ctab_${wt}wt_r${seed} --us 3
  done
done
