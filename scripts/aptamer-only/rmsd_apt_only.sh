#!/usr/bin/env bash

set -euo pipefail

OUT_ROOT="rmsd_apt_only_results"

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

        # both choose 0（SOLU）

        printf "0\n0\n" | gmx rms -s md.tpr -f md_nopbc.xtc -n index.ndx -tu ns -o rmsd_apt_only.xvg

        awk 'BEGIN{OFS=",";print "Time(ns)","RMSD(nm)"} /^[^@#]/{print $1,$2}' rmsd_apt_only.xvg > "${OUT_DIR}/rmsd_${C}.csv"

      else

        echo "Skip ${G} (missing files: md.tpr / md_nopbc.xtc / index.ndx)"

      fi

    )

  done

done

echo "All RMSD calculations completed. Results in ${OUT_ROOT}/"