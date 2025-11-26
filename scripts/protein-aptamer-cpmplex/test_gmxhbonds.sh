#!/usr/bin/env bash
set -euo pipefail

BEST_CSV="best_mmpbsa.csv"   # contain columns：complex_name, af3_model_id, chai_model_id, boltz_model_id, rf2na_model_id

declare -A paths=(
  [af3_complex_md]="af3"
  [chai_complex_md]="chai"
  [boltz_complex_md]="boltz"
  [rf2na_complex_md]="rf2na"   
  [gt_2.0water]="gt"
)

declare -A af3_models
declare -A chai_models
declare -A boltz_models
declare -A rf2na_models   

if [[ ! -f "$BEST_CSV" ]]; then
  echo "Did not find CSV: $BEST_CSV"
  exit 1
fi

mapfile -t CSV_ROWS < <(
  awk -F',' '
    NR==1{
      for(i=1;i<=NF;i++){
        h=tolower($i)
        gsub(/^[ \t]+|[ \t]+$/,"",h)
        if(h=="complex_name") c=i
        else if(h=="af3_model_id") a=i
        else if(h=="chai_model_id") ch=i
        else if(h=="boltz_model_id") b=i
        else if(h=="rf2na_model_id") r2=i
      }
      next
    }
    NR>=2{
      cn=(c? $c:""); af=(a? $a:""); chv=(ch? $ch:""); bo=(b? $b:""); rfv=(r2? $r2:"");
      gsub(/^[ \t]+|[ \t]+$/,"",cn)
      gsub(/^[ \t]+|[ \t]+$/,"",af)
      gsub(/^[ \t]+|[ \t]+$/,"",chv)
      gsub(/^[ \t]+|[ \t]+$/,"",bo)
      gsub(/^[ \t]+|[ \t]+$/,"",rfv)
      if(cn!=""){ print cn "\t" af "\t" chv "\t" bo "\t" rfv }
    }
  ' "$BEST_CSV"
)

COMPLEXES=()
for row in "${CSV_ROWS[@]}"; do
  IFS=$'\t' read -r comp af chai boltz rf2na <<<"$row"
  [[ -n "${af}"    && "${af}" =~ ^[01]$    ]] && af3_models["$comp"]="${af}"
  [[ -n "${chai}"  && "${chai}" =~ ^[01]$  ]] && chai_models["$comp"]="${chai}"
  [[ -n "${boltz}" && "${boltz}" =~ ^[01]$ ]] && boltz_models["$comp"]="${boltz}"
  [[ -n "${rf2na}" && "${rf2na}" =~ ^[01]$ ]] && rf2na_models["$comp"]="${rf2na}"
  COMPLEXES+=("$comp")
done

get_group_indices() {
  local ndxfile="$1"
  mapfile -t _grp < <(grep -n '^\[' "$ndxfile" || true)
  local protein_idx="" aptamer_idx="" protein_group="" aptamer_group=""
  for i in "${!_grp[@]}"; do
    local line="${_grp[$i]}"
    local name="${line#*\[}"; name="${name%\]*}"; name="$(echo "$name" | xargs)"
    if [[ "$name" == "Protein" ]]; then
      protein_idx="$i"; protein_group="$name"
    elif [[ "$name" == "Aptamer" ]]; then
      aptamer_idx="$i"; aptamer_group="$name"
    fi
  done
  echo "$protein_idx" "$aptamer_idx" "$protein_group" "$aptamer_group"
}

run_hbond_in_gromacs_dir() {
  local model="$1" complex="$2" gmxp="$3"
  echo "Processing $model / $complex → $gmxp"
  if [[ ! -d "$gmxp" ]]; then
    echo "gromacs directory not exists：$gmxp"; return 0
  fi
  cd "$gmxp" || { echo "Cannot enter $gmxp"; return 0; }

  local xtc=""
  [[ -f md.xtc ]] && xtc="md.xtc" || { [[ -f md_nopbc.xtc ]] && xtc="md_nopbc.xtc"; }
  if [[ -z "$xtc" || ! -f md.tpr ]]; then
    echo "Missing md(.nopbc).xtc or md.tpr"
    cd - >/dev/null; return 0
  fi
  if [[ ! -f index_mmpbsa.ndx ]]; then
    echo "Missing index_mmpbsa.ndx"
    cd - >/dev/null; return 0
  fi

  read -r protein_idx aptamer_idx protein_group aptamer_group < <(get_group_indices "index_mmpbsa.ndx")
  if [[ -z "$protein_idx" || -z "$aptamer_idx" ]]; then
    echo "Missing Protein or Aptamer group in index_mmpbsa.ndx"
    cd - >/dev/null; return 0
  fi
  echo "   → Protein-H idx=$protein_idx ($protein_group)；Aptamer idx=$aptamer_idx ($aptamer_group)"

  # gmx hbond
  gmx hbond \
    -f "$xtc" -s md.tpr -n index_mmpbsa.ndx \
    -r "group \"$protein_group\"" \
    -t "group \"$aptamer_group\"" \
    -num hbnum.xvg \
    -dist hbdist.xvg \
    -ang  hbang.xvg \
    -dan  hbdan.xvg \
    -o    hbond.ndx \
    -b 1000 -dt 10 -tu ps \
    -xvg xmgrace

  echo "H-bond done for $complex ($model)."
  cd - >/dev/null
}

for base in "${!paths[@]}"; do
  model="${paths[$base]}"

  for complex in "${COMPLEXES[@]}"; do
    if [[ "$model" == "gt" ]]; then
      charmm_dir=$(find "$base/$complex" -maxdepth 1 -type d -name "charmm-gui-*" | head -n1 || true)
      if [[ -z "${charmm_dir:-}" ]]; then
        echo "Missing charmm-gui folder under [$model/$complex], skip"; continue
      fi
      run_hbond_in_gromacs_dir "$model" "$complex" "$charmm_dir/gromacs"
      continue
    fi

    subdir=""
    case "$model" in
      af3)   subdir="${af3_models[$complex]:-}";;
      chai)  subdir="${chai_models[$complex]:-}";;
      boltz) subdir="${boltz_models[$complex]:-}";;
      rf2na) subdir="${rf2na_models[$complex]:-}";;
    esac

    if [[ -n "${subdir:-}" ]]; then
      charmm_dir=$(find "$base/$complex/$subdir" -maxdepth 1 -type d -name "charmm-gui-*" | head -n1 || true)
      if [[ -z "${charmm_dir:-}" ]]; then
        echo "Missing charmm-gui-* in $base/$complex/$subdir [$model/$complex], skip"
        continue
      fi
      run_hbond_in_gromacs_dir "$model" "$complex" "$charmm_dir/gromacs"
    else
      ran_any=0
      for try_idx in 0 1; do
        charmm_dir=$(find "$base/$complex/$try_idx" -maxdepth 1 -type d -name "charmm-gui-*" | head -n1 || true)
        if [[ -n "${charmm_dir:-}" ]]; then
          run_hbond_in_gromacs_dir "$model" "$complex" "$charmm_dir/gromacs"
          ran_any=1
        fi
      done
      if [[ $ran_any -eq 0 ]]; then
        echo Missing 0/1 in [$model/$complex] CSV, skip"
      fi
    fi
  done
done