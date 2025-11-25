#!/usr/bin/env bash

set -euo pipefail

OUT_ROOT="rg_apt_only_results"

MODELS=(3dRNA_DNA_apt)
#MODELS=(af3_apt chai_apt boltz_apt rf2na_apt)

for M in "${MODELS[@]}"; do

  find "$M" -type d -name gromacs -print0 | while IFS= read -r -d '' G; do

    C=$(basename "$(dirname "$(dirname "$G")")") 

    OUT_DIR="$(pwd)/${OUT_ROOT}/${M}"           

    mkdir -p "$OUT_DIR"

    (

      cd "$G"

      if [[ -f md.tpr && -f md_nopbc.xtc && -f index.ndx ]]; then

        echo "Processing ${M} / ${C} ..."

        echo SOLU | gmx gyrate -s md.tpr -f md_nopbc.xtc -n index.ndx -b 1 -tu ns -o rg_apt_only.xvg

        awk 'BEGIN{OFS=",";print "Time(ns)","Rg(nm)"} /^[^@#]/ {print $1,$2}' rg_apt_only.xvg > "${OUT_DIR}/rg_${C}.csv"

      else

        echo "Skip ${G} (missing files)"

      fi

    )

  done

done

echo "All Rg calculations completed. Results in ${OUT_ROOT}/"