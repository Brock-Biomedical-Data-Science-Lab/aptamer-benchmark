#!/bin/bash

root_dir=~/hbplus

find "$root_dir" -type f -name "*.pdb" | while read -r pdb_file; do
  echo "Running hbplus for: $pdb_file"
  (
    cd "$(dirname "$pdb_file")" || exit
    hbplus "$(basename "$pdb_file")"
  )
done
