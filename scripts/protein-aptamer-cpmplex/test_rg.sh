#!/usr/bin/env bash
set -euo pipefail

MODELS=(af3 boltz chai rf2na gt)

declare -A BASE_PATH=(
  [af3]="af3_complex_md"
  [boltz]="boltz_complex_md"
  [chai]="chai_complex_md"
  [rf2na]="rf2na_complex_md"
  [gt]="gt_2.0water"
)

BEST_CSV="best_mmpbsa.csv"

declare -A AF3_IDX BOLTZ_IDX CHAI_IDX RF2NA_IDX
if [[ -f "$BEST_CSV" ]]; then
  while IFS=' ' read -r complex af3_id chai_id boltz_id rf2na_id; do
    [[ -z "${complex:-}" || "${complex,,}" == "complex_name" ]] && continue
    AF3_IDX["$complex"]="${af3_id:-}"
    CHAI_IDX["$complex"]="${chai_id:-}"
    BOLTZ_IDX["$complex"]="${boltz_id:-}"
    RF2NA_IDX["$complex"]="${rf2na_id:-}"
  done < <(
    awk -F',' '
      BEGIN{OFS=" "}
      NR==1{ for(i=1;i<=NF;i++){ h=$i; gsub(/^[ \t]+|[ \t]+$/,"",h); idx[h]=i } ; next }
      {
        cn = (idx["complex_name"]   ? $(idx["complex_name"])   : "");
        a  = (idx["af3_model_id"]   ? $(idx["af3_model_id"])   : "");
        c  = (idx["chai_model_id"]  ? $(idx["chai_model_id"])  : "");
        b  = (idx["boltz_model_id"] ? $(idx["boltz_model_id"]) : "");
        r  = (idx["rf2na_model_id"] ? $(idx["rf2na_model_id"]) : "");
        gsub(/^[ \t]+|[ \t]+$/,"",cn); gsub(/^[ \t]+|[ \t]+$/,"",a);
        gsub(/^[ \t]+|[ \t]+$/,"",c);  gsub(/^[ \t]+|[ \t]+$/,"",b);
        gsub(/^[ \t]+|[ \t]+$/,"",r);
        if (cn!="") print cn, a, c, b, r;
      }
    ' "$BEST_CSV"
  )
else
  echo "BEST_CSV not found: $BEST_CSV"; exit 1
fi

COMPLEX_LIST=(7lri 7szu 7v5n 7zko 7zqs 8bw5 8d29 8tfd 8tqs 8zbf 9gxh)

CSV_ROOT=~/rg_results
for m in "${MODELS[@]}"; do
  mkdir -p "$CSV_ROOT/${m}_apt" "$CSV_ROOT/${m}_complex"
done

LOG=rg_workflow.log
: > "$LOG"

have_required_files() {
  [[ -f index_mmpbsa.ndx && -f md_nopbc.xtc && -f md.tpr ]]
}

get_group_idx_by_name() {
  local name="$1"
  awk -v tgt="$name" '
    BEGIN{idx=-1;sec=-1}
    /^\[/{sec++}
    $0 ~ "^[[:space:]]*\\[[[:space:]]*"tgt"[[:space:]]*\\]" {idx=sec}
    END{if(idx>=0) print idx}
  ' index_mmpbsa.ndx
}

make_protapt_and_get_idx() {
  # create index_rg.ndx，combine Protein | Aptamer to ProtApt，and return last group index
  local prot="$1" apt="$2"
  gmx make_ndx -f md.gro -n index_mmpbsa.ndx -o index_rg.ndx <<EOF >/dev/null
$prot | $apt
name X ProtApt
q
EOF
  grep '^\[' index_rg.ndx | nl -ba -v 0 | tail -n1 | awk '{print $1}'
}

write_xvg_to_csv() {
  awk 'BEGIN{OFS=",";print "Time(ns)","Rg(nm)"} /^[[:space:]]*[0-9]/{print $1,$2}' "$1" > "$2"
}

for model in "${MODELS[@]}"; do
  echo -e "\n Model: $model" | tee -a "$LOG"
  base="${BASE_PATH[$model]}"

  for complex in "${COMPLEX_LIST[@]}"; do
    gmx_dir=""
    idx=""
    case "$model" in
      af3)   idx="${AF3_IDX[$complex]:-}" ;;
      boltz) idx="${BOLTZ_IDX[$complex]:-}" ;;
      chai)  idx="${CHAI_IDX[$complex]:-}" ;;
      rf2na) idx="${RF2NA_IDX[$complex]:-}" ;;  
    esac

    if [[ "$model" == "gt" ]]; then
      gmx_dir=$(find "$base/$complex" -type d -name "gromacs" | head -n1 || true)
    else
      gmx_dir=$(find "$base/$complex/$idx" -type d -path "*/charmm-gui-*/gromacs" | head -n1 || true)
    fi

    if [[ -z "${gmx_dir:-}" ]]; then
      echo "Skip $model/$complex: no gromacs dir (idx=${idx:-none})" | tee -a "$LOG"
      continue
    fi

    echo "$model / $complex @ $gmx_dir" | tee -a "$LOG"
    pushd "$gmx_dir" >/dev/null

    if ! have_required_files; then
      echo "Missing index_mmpbsa.ndx or md_nopbc.xtc or md.tpr" | tee -a "$LOG"
      popd >/dev/null
      continue
    fi

    apt_idx="$(get_group_idx_by_name 'Aptamer' || true)"
    prot_idx="$(get_group_idx_by_name 'Protein' || true)"
    if [[ -z "${apt_idx:-}" || -z "${prot_idx:-}" ]]; then
      echo "Missing required groups [Aptamer] or [ProtenH]" | tee -a "$LOG"
      popd >/dev/null
      continue
    fi

    protapt_grp="$(make_protapt_and_get_idx "$prot_idx" "$apt_idx")"

    idx_tag="" 

    # 1) Aptamer Rg
    echo "$apt_idx" | gmx gyrate -s md.tpr -f md_nopbc.xtc -n index_mmpbsa.ndx \
      -o rg_aptamer.xvg -b 1 -tu ns
    out_csv_apt="$CSV_ROOT/${model}_apt/rg_${complex}${idx_tag}.csv"
    write_xvg_to_csv rg_aptamer.xvg "$out_csv_apt"

    # 2) ProtApt Rg
    echo "$protapt_grp" | gmx gyrate -s md.tpr -f md_nopbc.xtc -n index_rg.ndx \
      -o rg_complex.xvg -b 1 -tu ns
    out_csv_cplx="$CSV_ROOT/${model}_complex/rg_${complex}${idx_tag}.csv"
    write_xvg_to_csv rg_complex.xvg "$out_csv_cplx"

    popd >/dev/null
    echo "Done $model/$complex" | tee -a "$LOG"
  done
done

echo -e "\n Finished.
• All CSVs under $CSV_ROOT
• .xvg files remain in each gromacs folder." | tee -a "$LOG"