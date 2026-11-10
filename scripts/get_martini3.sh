#!/bin/sh
set -e
cd "$(dirname "$0")/../ff"
C=784591ebdc91d762ed4df986c4650546c938f776
B=https://raw.githubusercontent.com/marrink-lab/martini-forcefields/$C/martini_forcefields/regular/v3.0.0/gmx_files
for f in martini_v3.0.0.itp martini_v3.0.0_ions_v1.itp martini_v3.0.0_solvents_v1.itp; do
  curl -fsSL -o "$f" "$B/$f"
done
shasum -a 256 -c <<SUMS
0ecbf0ff334b2ec21dc5e3151aca42167b977bc2aed41eb2a44bd8c997754894  martini_v3.0.0.itp
699aba132755b0ef13a0d7b41e6230a27777781bdc8075507413b89499c2be4f  martini_v3.0.0_ions_v1.itp
e6a78414b317c38a1a816eded4aee6082a5da1ed229a126359f0cc92314dcb4e  martini_v3.0.0_solvents_v1.itp
SUMS
