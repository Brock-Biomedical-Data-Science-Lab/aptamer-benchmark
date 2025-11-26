#!/usr/bin/env bash
set -Eeuo pipefail

# MODELS list. Add model directories here if needed.
# MODELS=(af3_apt chai_apt boltz_apt rf2na_apt)
MODELS=(3dRNA_DNA_apt)

# Group ID in index.ndx used for centering and output (0 = SOLU).
# Modify this number if a different group should be used.
GROUP=0

LOG="$(pwd)/trjconv_pbc.log"
: > "$LOG"

for M in "${MODELS[@]}"; do

  while IFS= read -r -d '' G; do

    echo -e "\nProcessing directory: $G" | tee -a "$LOG"

    if [[ -f "$G/md.tpr" && -f "$G/md.xtc" && -f "$G/index.ndx" ]]; then

      (
        cd "$G"
        rm -f md_nopbc.xtc

        # Run gmx trjconv in one line, with real-time output and logging
        echo -e "${GROUP}\n${GROUP}" | gmx trjconv \
          -s md.tpr -f md.xtc -n index.ndx \
          -o md_nopbc.xtc -pbc mol -center \
          2>&1 | tee -a "$LOG"
      )

      if [[ -s "$G/md_nopbc.xtc" ]]; then
        echo "Completed: $G" | tee -a "$LOG"
      else
        echo "Failed: $G (output file missing or empty)" | tee -a "$LOG"
      fi

    else
      missing=()
      [[ -f "$G/md.tpr" ]]   || missing+=("md.tpr")
      [[ -f "$G/md.xtc" ]]   || missing+=("md.xtc")
      [[ -f "$G/index.ndx" ]] || missing+=("index.ndx")

      echo "Missing required files: ${missing[*]}" | tee -a "$LOG"
    fi

  done < <(find "$M" -type d -name gromacs -print0)

done

echo -e "\nAll processing completed. Log saved to: $LOG"