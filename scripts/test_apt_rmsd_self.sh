#!/usr/bin/env bash
set -euo pipefail

# ======= Config =======
COMPLEX_LIST=(7lri 7szu 7v5n 7zko 7zqs 8bw5 8d29 8tfd 8tqs 8zbf 9gxh)
BEST_CSV="best_mmpbsa.csv"

# POS 
POS_MODELS=(af3 boltz chai rf2na gt)
declare -A POS_BASE_PATH=(
  [af3]="af3_complex_md"
  [boltz]="boltz_complex_md"
  [chai]="chai_complex_md"
  [rf2na]="rf2na_complex_md"
  [gt]="gt_2.0water"
)

# NEG 
NEG_MODELS=(af3 boltz chai rf2na)
declare -A NEG_BASE_PATH=(
  [af3]="af3_neg_md"
  [boltz]="boltz_neg_md"
  [chai]="chai_neg_md"
  [rf2na]="rf2na_neg_md"
)

# Output dirs
mkdir -p rmsd_apt_results/pos/{af3,boltz,chai,rf2na,gt} rmsd_apt_results/neg/{af3,boltz,chai,rf2na}

# ======= Functions =======
run_rmsd_in_dir() {
  local gmx_dir="$1"
  local tag_raw="$2"
  local tag="${tag_raw,,}"

  [[ -d "$gmx_dir" ]] || return 1
  pushd "$gmx_dir" >/dev/null || return 1

  [[ -f md.tpr && -f md_nopbc.xtc && -f index_mmpbsa.ndx ]] || { popd >/dev/null; return 1; }

  local last_id
  last_id=$(awk 'BEGIN{i=-1} /^\[/{i++} END{print i}' index_mmpbsa.ndx) || last_id=0
  local fit_id=$last_id

  gmx rms -s md.tpr -f md_nopbc.xtc -n index_mmpbsa.ndx -tu ns -o "rmsd_${tag}_apt.xvg" <<EOF
$fit_id
$last_id
EOF

  popd >/dev/null
}

xvg_to_csv() {
  local xvg_file="$1"
  local out_file="$2"
  mkdir -p "$(dirname "$out_file")"
  awk 'BEGIN{OFS=",";print "Time(ns)","RMSD(nm)"} 
       !/^[@#]/ && NF>=2 {print $1,$2}' "$xvg_file" > "$out_file"
}

convert_pos_to_csv() {
  while IFS=',' read -r complex af3_id af3_delta chai_id chai_delta boltz_id boltz_delta rf2na_id rf2na_delta gt_delta; do
    [[ "$complex" == "complex_name" ]] && continue

    # af3
    if [[ -n "$af3_id" ]]; then
      g=$(find "${POS_BASE_PATH[af3]}/$complex/$af3_id" -type d -name gromacs | head -n1 || true)
      if [[ -n "$g" ]]; then
        run_rmsd_in_dir "$g" "${complex}_${af3_id}"
        out_file="rmsd_apt_results/pos/af3/${complex}_${af3_id}_apt.csv"
        xvg_to_csv "$g/rmsd_${complex}_${af3_id}_apt.xvg" "$out_file"
        echo "POS-CSV: af3/$complex/$af3_id -> $out_file"
      fi
    fi

    # chai
    if [[ -n "$chai_id" ]]; then
      g=$(find "${POS_BASE_PATH[chai]}/$complex/$chai_id" -type d -name gromacs | head -n1 || true)
      if [[ -n "$g" ]]; then
        run_rmsd_in_dir "$g" "${complex}_${chai_id}"
        out_file="rmsd_apt_results/pos/chai/${complex}_${chai_id}_apt.csv"
        xvg_to_csv "$g/rmsd_${complex}_${chai_id}_apt.xvg" "$out_file"
        echo "POS-CSV: chai/$complex/$chai_id -> $out_file"
      fi
    fi

    # boltz
    if [[ -n "$boltz_id" ]]; then
      g=$(find "${POS_BASE_PATH[boltz]}/$complex/$boltz_id" -type d -name gromacs | head -n1 || true)
      if [[ -n "$g" ]]; then
        run_rmsd_in_dir "$g" "${complex}_${boltz_id}"
        out_file="rmsd_apt_results/pos/boltz/${complex}_${boltz_id}_apt.csv"
        xvg_to_csv "$g/rmsd_${complex}_${boltz_id}_apt.xvg" "$out_file"
        echo "POS-CSV: boltz/$complex/$boltz_id -> $out_file"
      fi
    fi

    # rf2na
    if [[ -n "$rf2na_id" ]]; then
      g=$(find "${POS_BASE_PATH[rf2na]}/$complex/$rf2na_id" -type d -name gromacs | head -n1 || true)
      if [[ -n "$g" ]]; then
        run_rmsd_in_dir "$g" "${complex}_${rf2na_id}"
        out_file="rmsd_apt_results/pos/rf2na/${complex}_${rf2na_id}_apt.csv"
        xvg_to_csv "$g/rmsd_${complex}_${rf2na_id}_apt.xvg" "$out_file"
        echo "POS-CSV: rf2na/$complex/$rf2na_id -> $out_file"
      fi
    fi

    # gt
    g=$(find "${POS_BASE_PATH[gt]}/$complex" -type d -name gromacs | head -n1 || true)
    if [[ -n "$g" ]]; then
      run_rmsd_in_dir "$g" "$complex"
      out_file="rmsd_apt_results/pos/gt/${complex}_apt.csv"
      xvg_to_csv "$g/rmsd_${complex}_apt.xvg" "$out_file"
      echo "POS-CSV: gt/$complex -> $out_file"
    fi

  done < "$BEST_CSV"
}

convert_neg_to_csv() {
  for model in "${NEG_MODELS[@]}"; do
    local base="${NEG_BASE_PATH[$model]}"

    shopt -s nullglob
    local dirs=("$base"/*/charmm-gui-*/gromacs)
    shopt -u nullglob
    [[ ${#dirs[@]} -gt 0 ]] || continue

    for d in "${dirs[@]}"; do
      local complex_dir
      complex_dir=$(basename "$(dirname "$(dirname "$d")")")
      local lc_complex="${complex_dir,,}"
      local uc_complex="${complex_dir^^}"
      run_rmsd_in_dir "$d" "$lc_complex" 
     
      local xvg="$d/rmsd_${lc_complex}_apt.xvg"
      if [[ -f "$xvg" ]]; then
        local out_file="rmsd_apt_results/neg/$model/${uc_complex}.csv"
        xvg_to_csv "$xvg" "$out_file"
        echo "NEG-CSV: $(basename "$xvg") -> $out_file"
      else
        echo "NEG SKIP: no $xvg"
      fi
    done
  done
}

# ======= Run All =======
#convert_pos_to_csv
convert_neg_to_csv
echo "Done. CSVs under rmsd_apt_results/."