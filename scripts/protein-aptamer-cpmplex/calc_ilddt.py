#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, csv, sys

BASE_ROOT = os.path.expanduser("~/tmscore_ilddt_in")
OUTPUT_CSV = os.path.join(BASE_ROOT, "ilddt.csv")
MODELS = ["af3", "chai", "boltz", "rf2na"]

RE_SEED_SAMPLE = re.compile(r"^(?P<cpx>[A-Za-z0-9]+)_seed(?P<seed>\d+)_sample(?P<sample>\d+)\.pdb$")
RE_CHAI = re.compile(r"^(?P<cpx>[A-Za-z0-9]+)_seed(?P<seed>\d+)_model(?P<model>\d+)\.pdb$", re.IGNORECASE)

try:
    from ost import io
    from ost.mol.alg import scoring
except Exception:
    sys.stderr.write("Cannot import OpenStructure (ost). Ensure it is installed.\n")
    raise

ALL_IDS = [(s, sm) for s in range(1, 6) for sm in range(0, 5)]

def find_complexes(base_root: str, models: list[str]) -> list[str]:
    found = set()
    for m in models:
        mdir = os.path.join(base_root, m)
        if not os.path.isdir(mdir):
            continue
        for name in os.listdir(mdir):
            if os.path.isdir(os.path.join(mdir, name)):
                found.add(name.lower())
    return sorted(found)

def find_reference_pdb(base_root: str, complex_name: str, models_priority: list[str]) -> str | None:
    for m in models_priority:
        mdir = os.path.join(base_root, m)
        if not os.path.isdir(mdir):
            continue
        for sub in os.listdir(mdir):
            if sub.lower() != complex_name:
                continue
            cdir = os.path.join(mdir, sub)
            for fn in os.listdir(cdir):
                if fn.lower() == f"{complex_name}.pdb":
                    return os.path.join(cdir, fn)
    return None

def scan_seed_sample_dir(base_root, model, complex_name):
    out = {}
    cdir = os.path.join(base_root, model, complex_name)
    if not os.path.isdir(cdir):
        return out
    for fn in os.listdir(cdir):
        m = RE_SEED_SAMPLE.match(fn)
        if not m:
            continue
        if m.group("cpx").lower() != complex_name:
            continue
        seed = int(m.group("seed"))
        sample = int(m.group("sample"))
        out[(seed, sample)] = os.path.join(cdir, fn)
    return out

def scan_chai(base_root, complex_name):
    out = {}
    cdir = os.path.join(base_root, "chai", complex_name)
    if not os.path.isdir(cdir):
        return out
    for fn in os.listdir(cdir):
        m = RE_CHAI.match(fn)
        if not m:
            continue
        if m.group("cpx").lower() != complex_name:
            continue
        seed = int(m.group("seed"))
        sample = int(m.group("model"))
        out[(seed, sample)] = os.path.join(cdir, fn)
    return out

def remap_to_grid_5x5(pred_map: dict[tuple[int,int], str]) -> dict[tuple[int,int], str]:
    items = sorted(pred_map.items(), key=lambda kv: (kv[0][0], kv[0][1], kv[1]))
    paths = [p for (_, p) in items]
    new_map = {}
    for (seed, sample), path in zip(ALL_IDS, paths[:25]):
        new_map[(seed, sample)] = path
    return new_map

def calc_ilddt(ref_ent, ref_pdb_path: str, pred_path: str) -> float | None:
    try:
        ent_pred = io.LoadPDB(pred_path).CreateFullView()
        s = scoring.Scorer(ref_ent, ent_pred)
        return float(s.ilddt)
    except Exception as e:
        print(f"iLDDT failed: {pred_path} → {e}")
        return None

def main():
    complexes = find_complexes(BASE_ROOT, MODELS)
    if not complexes:
        print("No complexes found")
        return

    header = ["Complex", "ID"] + [f"{m}_iLDDT" for m in MODELS]
    rows_written = 0

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)

        for complex_name in complexes:
            ref_pdb = find_reference_pdb(BASE_ROOT, complex_name, MODELS)
            if not ref_pdb:
                print(f"Missing reference {complex_name}.pdb, skip")
                continue

            print(f"\nProcessing {complex_name} (ref: {ref_pdb})")
            try:
                ref_ent = io.LoadPDB(ref_pdb).CreateFullView()
            except Exception as e:
                print(f"Failed to load ref {ref_pdb}: {e}")
                continue

            af3_map   = scan_seed_sample_dir(BASE_ROOT, "af3",   complex_name)
            boltz_map = scan_seed_sample_dir(BASE_ROOT, "boltz", complex_name)
            chai_map  = scan_chai(BASE_ROOT, complex_name)
            rf2na_raw = scan_seed_sample_dir(BASE_ROOT, "rf2na", complex_name)

            rf2na_map = remap_to_grid_5x5(rf2na_raw)

            model_maps = {
                "af3": af3_map,
                "chai": chai_map,
                "boltz": boltz_map,
                "rf2na": rf2na_map
            }

            for seed, sample in ALL_IDS:
                id_str = f"seed{seed}_sample{sample}"
                row = [complex_name, id_str]
                for m in MODELS:
                    pred = model_maps[m].get((seed, sample))
                    if not pred:
                        row.append("NA")
                        continue
                    val = calc_ilddt(ref_ent, ref_pdb, pred)
                    row.append(f"{val:.4f}" if isinstance(val, float) else "NA")
                writer.writerow(row)
                rows_written += 1

    print(f"Done. Wrote {rows_written} rows → {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
