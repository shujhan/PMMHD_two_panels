#!/usr/bin/env python3
import os
import json
import subprocess
import time


def parse_species(deck_path):
    """Read species names from deck.json."""
    with open(deck_path, "r") as f:
        deck = json.load(f)

    species = []
    if "initial_list" in deck:
        for sp in deck["initial_list"]:
            if "name" in sp:
                species.append(sp["name"])
    return species


def ensure_dirs(base_dir, species_names):
    """Create all output directories for each species."""
    for species in species_names:
        base = os.path.join(base_dir, "simulation_output", species)
        subdirs = ["xs", "ys", "w0s", "j0s", "q0s", "u1s", "u2s", "b1s", "b2s", "panels"]
        for sd in subdirs:
            path = os.path.join(base, sd)
            os.makedirs(path, exist_ok=True)


def main():
    # project root
    script_dir = os.path.abspath(os.path.dirname(__file__))
    deck_path = os.path.join(script_dir, "deck.json")

    if not os.path.isfile(deck_path):
        raise FileNotFoundError(f"Deck file not found: {deck_path}")

    # Parse species and make directories
    species_names = parse_species(deck_path)
    if not species_names:
        raise ValueError("No species found in deck.json under 'initial_list'.")

    ensure_dirs(script_dir, species_names)

    # Run simulation
    print(f"\n[RUN] Executing simulation in: {script_dir}")

    use_gpu = False
    if use_gpu:
        PMMHD_exe = os.path.join(script_dir, "PMMHD_gpu")
    else:
        PMMHD_exe = os.path.join(script_dir, "PMMHD_cpu")

    if not os.path.isfile(PMMHD_exe):
        raise FileNotFoundError(f"Executable not found: {PMMHD_exe}")

    PMMHD_args = [PMMHD_exe]
    out_path = os.path.join(script_dir, "sim_out.txt")
    err_path = os.path.join(script_dir, "sim_err.txt")

    t1 = time.time()
    with open(out_path, "w") as log:
        with open(err_path, "w") as err_log:
            proc = subprocess.run(PMMHD_args, cwd=script_dir, stdout=log, stderr=err_log)
    t2 = time.time()

    if proc.returncode == 0:
        print("done!")
    else:
        print(f"simulation stopped with error code {proc.returncode}")

    with open(out_path, "a") as log:
        log.write(f"\ncompute time {t2 - t1:.03f}s\n")

    print(f"simulation compute time {t2 - t1:.03f}s")


if __name__ == "__main__":
    main()