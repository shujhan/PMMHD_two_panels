#!/usr/bin/env python3
import os
import json
import glob
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon


def load_deck_species(deck_path):
    with open(deck_path, "r") as f:
        d = json.load(f)
    return [sp["name"] for sp in d.get("initial_list", []) if isinstance(sp, dict) and "name" in sp]


def read_bin_doubles(path):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Missing file: {path}")
    return np.fromfile(path, dtype=np.float64)


def pick_latest_tag(species_dir):
    """
    Look under simulation_output/<species>/{xs,ys,w0s,j0s,q0s,u1s,u2s,b1s,b2s}/
    and pick a tag that exists in all needed folders.
    """
    need_dirs = ["xs", "ys", "w0s", "j0s", "q0s", "u1s", "u2s", "b1s", "b2s"]

    xs_dir = os.path.join(species_dir, "xs")
    xs_files = sorted(glob.glob(os.path.join(xs_dir, "xs_*")), key=os.path.getmtime)
    if not xs_files:
        raise FileNotFoundError(f"No xs_* files found in {xs_dir}")

    for candidate in reversed(xs_files):
        tag = os.path.basename(candidate)[len("xs_"):]
        ok = True
        for dname in need_dirs:
            prefix = dname[:-1] if dname.endswith("s") else dname
            # 这里实际文件名都是目录名同名，例如 w0s/w0s_0, b1s/b1s_0
            path = os.path.join(species_dir, dname, f"{dname}_{tag}")
            if not os.path.isfile(path):
                ok = False
                break
        if ok:
            return tag

    return os.path.basename(xs_files[-1])[len("xs_"):]


def build_paths(base_dir, species, time):
    """
    Read data from simulation_output/<species>/...
    """
    sp_dir = os.path.join(base_dir, "simulation_output", species)
    return {
        "xs":  os.path.join(sp_dir, "xs",  f"xs_{time}"),
        "ys":  os.path.join(sp_dir, "ys",  f"ys_{time}"),
        "w0s": os.path.join(sp_dir, "w0s", f"w0s_{time}"),
        "j0s": os.path.join(sp_dir, "j0s", f"j0s_{time}"),
        "q0s": os.path.join(sp_dir, "q0s", f"q0s_{time}"),
        "u1s": os.path.join(sp_dir, "u1s", f"u1s_{time}"),
        "u2s": os.path.join(sp_dir, "u2s", f"u2s_{time}"),
        "b1s": os.path.join(sp_dir, "b1s", f"b1s_{time}"),
        "b2s": os.path.join(sp_dir, "b2s", f"b2s_{time}"),
    }


def find_panels_path(base_dir, species, time):
    """
    Try a few common locations for panel file.
    """
    candidates = [
        os.path.join(base_dir, "simulation_output", species, "panels", f"leaf_point_inds_{time}"),
        os.path.join(base_dir, "simulation_output", "panels", f"leaf_point_inds_{time}"),
        os.path.join(base_dir, "simulation_output", f"leaf_point_inds_{time}"),
        os.path.join(base_dir, f"leaf_point_inds_{time}"),
    ]

    for path in candidates:
        if os.path.isfile(path):
            return path

    raise FileNotFoundError(
        f"Cannot find leaf_point_inds_{time}. Checked:\n" + "\n".join(candidates)
    )


def plot_species(ax, x, y, f, species, time, field_name=None):
    n = min(len(x), len(y), len(f))
    if n == 0:
        raise ValueError(f"No particles to plot for {species} (time {time}).")

    x = x[:n]
    y = y[:n]
    f = f[:n]

    sc = ax.scatter(x, y, c=f, s=3, alpha=0.8)
    ax.set_aspect("equal", adjustable="box")
    if field_name is None:
        ax.set_title(f"{species}   [{n} particles]   time: {time}")
    else:
        ax.set_title(f"{field_name}   [{n} particles]   time: {time}")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    plt.colorbar(sc, ax=ax)


def plot_phase_space(ax, xs, ys, fs, panels, species, step_ii, dt, symmetric=False, clim=None, field_name=None):
    xs = np.asarray(xs)
    ys = np.asarray(ys)
    fs = np.asarray(fs)
    panels = np.asarray(panels)

    simtime = int(step_ii) * dt

    num_panels = int(panels.size / 9)
    panels = np.reshape(panels, (num_panels, 9))
    panels_fs = np.zeros(4 * num_panels)

    patches = []
    for ii, panel in enumerate(panels):
        panel_xs = xs[panel]
        panel_ys = ys[panel]
        panel_fs = fs[panel]

        p0 = [0, 1, 4, 3]
        panels_fs[4 * ii] = 0.25 * sum(panel_fs[p0])
        rect_pts = np.vstack([panel_xs[p0], panel_ys[p0]]).T
        patches.append(Polygon(rect_pts))

        p1 = [1, 2, 5, 4]
        panels_fs[4 * ii + 1] = 0.25 * sum(panel_fs[p1])
        rect_pts = np.vstack([panel_xs[p1], panel_ys[p1]]).T
        patches.append(Polygon(rect_pts))

        p2 = [3, 4, 7, 6]
        panels_fs[4 * ii + 2] = 0.25 * sum(panel_fs[p2])
        rect_pts = np.vstack([panel_xs[p2], panel_ys[p2]]).T
        patches.append(Polygon(rect_pts))

        p3 = [4, 5, 8, 7]
        panels_fs[4 * ii + 3] = 0.25 * sum(panel_fs[p3])
        rect_pts = np.vstack([panel_xs[p3], panel_ys[p3]]).T
        patches.append(Polygon(rect_pts))

    if clim is None:
        if symmetric:
            m = np.percentile(np.abs(panels_fs), 95)
            if m == 0:
                m = np.max(np.abs(panels_fs)) if np.max(np.abs(panels_fs)) > 0 else 1.0
            clim = (-m, m)
        else:
            vmin = np.nanmin(fs)
            vmax = np.nanmax(fs)
            if vmin == vmax:
                vmin, vmax = float(np.min(panels_fs)), float(np.max(panels_fs))
                if vmin == vmax:
                    vmin, vmax = vmin - 1.0, vmax + 1.0
            clim = (vmin, vmax)

    p = PatchCollection(patches, cmap=matplotlib.cm.jet)
    p.set_array(panels_fs)
    p.set_clim(clim)
    ax.add_collection(p)

    ax.set_xlim(np.min(xs), np.max(xs))
    ax.set_ylim(np.min(ys), np.max(ys))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    if field_name is None:
        ax.set_title(f"{species}  time: {simtime}")
    else:
        ax.set_title(f"{field_name}  time: {simtime}")

    fig = ax.figure
    fig.colorbar(p, ax=ax)

    return p


def save_txt_outputs(base_dir, time, data_dict):
    for key, arr in data_dict.items():
        np.savetxt(os.path.join(base_dir, f"{key}_{time}.txt"), arr, fmt="%.6e")


def make_points_figure(species_dir, species, time, xs, ys, w0s, j0s, q0s, u1s, u2s, b1s, b2s):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.ravel()

    plot_species(axes[0], xs, ys, w0s, species, time, field_name="w0")
    plot_species(axes[1], xs, ys, j0s, species, time, field_name="j0")
    plot_species(axes[2], xs, ys, q0s, species, time, field_name="q0")
    plot_species(axes[3], xs, ys, u1s, species, time, field_name="u1")
    plot_species(axes[4], xs, ys, u2s, species, time, field_name="u2")
    plot_species(axes[5], xs, ys, b1s, species, time, field_name="b1")
    plot_species(axes[6], xs, ys, b2s, species, time, field_name="b2")
    axes[7].axis("off")

    fig.suptitle(f"{species} particle fields at time {time}", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    out_name = os.path.join(species_dir, f"particles_{time}.png")
    fig.savefig(out_name, dpi=200)
    plt.close(fig)
    print(f"[OK] Saved figure to {out_name}")


def make_phase_space_figure(species_dir, base_dir, species, time, dt, xs, ys, w0s, j0s, q0s, u1s, u2s, b1s, b2s, panels):
    fig, axes = plt.subplots(2, 4, figsize=(20, 9))
    axes = axes.ravel()

    plot_phase_space(axes[0], xs, ys, w0s, panels, species, time, dt, field_name="w0")
    plot_phase_space(axes[1], xs, ys, j0s, panels, species, time, dt, field_name="j0")
    plot_phase_space(axes[2], xs, ys, q0s, panels, species, time, dt, field_name="q0")
    plot_phase_space(axes[3], xs, ys, u1s, panels, species, time, dt, field_name="u1")
    plot_phase_space(axes[4], xs, ys, u2s, panels, species, time, dt, field_name="u2")
    plot_phase_space(axes[5], xs, ys, b1s, panels, species, time, dt, field_name="b1")
    plot_phase_space(axes[6], xs, ys, b2s, panels, species, time, dt, field_name="b2")
    axes[7].axis("off")

    fig.suptitle(f"{species} phase-space fields at time {time}", fontsize=16)
    fig.tight_layout(rect=[0, 0, 1, 0.96])

    if species == 'vorticity':
        out_name = os.path.join(base_dir, f"q_plus_phase_space_{time}.png")
    else:
        out_name = os.path.join(base_dir, f"q_minus_phase_space_{time}.png")
    fig.savefig(out_name, dpi=200)
    plt.close(fig)
    print(f"[OK] Saved figure to {out_name}")


def main():
    p = argparse.ArgumentParser(description="Read PMMHD particle dumps and plot fields for each species.")
    p.add_argument("--base-dir", default=".", help="Project root (default: current dir)")
    p.add_argument("--deck", default="deck.json", help="Deck path (default: deck.json in base dir)")
    p.add_argument("--time", default=None, help="File time suffix. If omitted, auto-pick latest per species.")
    p.add_argument("--plot_points", "-pts", action="store_true", help="plot point scatter")
    p.add_argument("--plot_phase_space", "-ps", action="store_true", help="plot panel phase space")
    args = p.parse_args()

    base_dir = os.path.abspath(args.base_dir)
    deck_path = os.path.join(base_dir, args.deck)
    if not os.path.isfile(deck_path):
        raise FileNotFoundError(f"Deck not found: {deck_path}")

    species_list = load_deck_species(deck_path)
    if not species_list:
        raise ValueError("No species found in deck.json: initial_list is empty.")

    with open(deck_path, "r") as f:
        d = json.load(f)
    dt = d.get("dt", 1.0)

    for species in species_list:
        species_dir = os.path.join(base_dir, "simulation_output", species)

        if args.time is None:
            time = pick_latest_tag(species_dir)
        else:
            time = str(args.time)

        paths = build_paths(base_dir, species, time)

        xs = read_bin_doubles(paths["xs"])
        ys = read_bin_doubles(paths["ys"])
        w0s = read_bin_doubles(paths["w0s"])
        j0s = read_bin_doubles(paths["j0s"])
        q0s = read_bin_doubles(paths["q0s"])
        u1s = read_bin_doubles(paths["u1s"])
        u2s = read_bin_doubles(paths["u2s"])
        b1s = read_bin_doubles(paths["b1s"])
        b2s = read_bin_doubles(paths["b2s"])

        # save_txt_outputs(base_dir, time, {
        #     "xs": xs,
        #     "ys": ys,
        #     "w0s": w0s,
        #     "j0s": j0s,
        #     "q0s": q0s,
        #     "u1s": u1s,
        #     "u2s": u2s,
        #     "b1s": b1s,
        #     "b2s": b2s,
        # })

        if args.plot_points:
            make_points_figure(species_dir, species, time, xs, ys, w0s, j0s, q0s, u1s, u2s, b1s, b2s)

        if args.plot_phase_space:
            panels_path = find_panels_path(base_dir, species, time)
            panels = np.fromfile(panels_path, dtype=np.int32)
            make_phase_space_figure(species_dir, base_dir, species, time, dt, xs, ys, w0s, j0s, q0s, u1s, u2s, b1s, b2s, panels)

    if (not args.plot_points) and (not args.plot_phase_space):
        print("[INFO] Nothing to do. Use -pts and/or -ps.")


if __name__ == "__main__":
    main()