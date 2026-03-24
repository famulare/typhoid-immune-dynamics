"""
Phase 2: Extract structured data from Figure S11 panels.

Panel boundaries and axis calibrations derived from Phase 1 inspection.
"""

import json
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_extraction")
RAW_DIR = OUT_DIR / "raw"

# =============================================================================
# Panel geometry (from Phase 1 grid lines and axis lines)
# =============================================================================

# Panel A: upper-left (individual Vi IgG trajectories)
PANEL_A = {
    "x_min": 69.9, "x_max": 266.8,
    "y_min": 150.3, "y_max": 236.1,
    # X-axis: LINEAR (verified: 5 ticks fit perfectly)
    # Grid lines: x=82.4→30d, x=103.4→90d, x=135.0→180d, x=198.2→360d, x=261.3→540d
    # Linear fit: pdf_x = 0.3508 * days + 71.9
    "x_ticks": [(82.4, 30), (103.4, 90), (135.0, 180), (198.2, 360), (261.3, 540)],
    # Y-axis: LOG10
    # Grid lines: y=175.1→1000 EU, y=219.0→100 EU
    # pdf_y = -43.9 * log10(EU) + 306.8
    "y_ticks": [(175.1, 1000), (219.0, 100)],
}

# Panel B: upper-right (box plots by age/catchment)
PANEL_B = {
    "x_min": 303.9, "x_max": 500.8,
    "y_min": 150.3, "y_max": 236.1,
    # X-axis: categorical
    "x_categories": {357.6: "Kathmandu", 447.1: "Kavre"},
    # Y-axis: LOG10 (same as Panel A)
    "y_ticks": [(175.1, 1000), (219.0, 100)],
}

# Panel C: bottom (boosting ratio box plots)
PANEL_C = {
    "x_min": 69.6, "x_max": 500.8,
    "y_min": 285.7, "y_max": 453.2,
    # X-axis: categorical
    "x_categories": {
        119.4: "HlyE IgA", 202.3: "HlyE IgG",
        285.2: "LPS IgA", 368.1: "LPS IgG", 451.0: "Vi IgG",
    },
    # Y-axis: LINEAR (0-50)
    # Grid lines: y=445.6→0, y=415.2→10, ..., y=293.3→50
    "y_ticks": [(445.6, 0), (415.2, 10), (384.7, 20),
                (354.2, 30), (323.8, 40), (293.3, 50)],
}


# =============================================================================
# Coordinate transforms
# =============================================================================

def fit_linear_transform(tick_pairs):
    """Fit pdf_coord = a * data_value + b from (pdf_coord, data_value) pairs."""
    pdf_coords = np.array([t[0] for t in tick_pairs])
    data_vals = np.array([t[1] for t in tick_pairs])
    # pdf = a * data + b  =>  solve by least squares
    A = np.column_stack([data_vals, np.ones_like(data_vals)])
    result = np.linalg.lstsq(A, pdf_coords, rcond=None)
    a, b = result[0]
    residuals = pdf_coords - (a * data_vals + b)
    return a, b, residuals


def fit_log_transform(tick_pairs):
    """Fit pdf_coord = a * log10(data_value) + b."""
    pdf_coords = np.array([t[0] for t in tick_pairs])
    log_vals = np.log10(np.array([t[1] for t in tick_pairs]))
    A = np.column_stack([log_vals, np.ones_like(log_vals)])
    result = np.linalg.lstsq(A, pdf_coords, rcond=None)
    a, b = result[0]
    residuals = pdf_coords - (a * log_vals + b)
    return a, b, residuals


def pdf_to_data_linear(pdf_coord, a, b):
    """Convert PDF coordinate to data value (linear axis)."""
    return (pdf_coord - b) / a


def pdf_to_data_log(pdf_coord, a, b):
    """Convert PDF coordinate to data value (log10 axis)."""
    log_val = (pdf_coord - b) / a
    return 10 ** log_val


def nearest_category(pdf_x, categories, tolerance=40):
    """Map a PDF x-coordinate to the nearest categorical label."""
    best_dist = float('inf')
    best_label = None
    for cat_x, label in categories.items():
        dist = abs(pdf_x - cat_x)
        if dist < best_dist:
            best_dist = dist
            best_label = label
    return best_label if best_dist < tolerance else None


def color_to_hex(c):
    if c is None:
        return None
    return "#{:02x}{:02x}{:02x}".format(
        int(c[0] * 255), int(c[1] * 255), int(c[2] * 255)
    )


def in_panel(d, panel, margin=5):
    """Check if a drawing's bounding box falls within a panel (with margin)."""
    r = d["rect"]
    if r is None:
        return False
    cx = (r[0] + r[2]) / 2
    cy = (r[1] + r[3]) / 2
    return (panel["x_min"] - margin <= cx <= panel["x_max"] + margin and
            panel["y_min"] - margin <= cy <= panel["y_max"] + margin)


# =============================================================================
# Extraction functions
# =============================================================================

def extract_line_points(d):
    """Extract ordered (x, y) points from a drawing with line segments."""
    points = []
    for item in d["items"]:
        if item["op"] == "line":
            if not points or (points[-1] != item["p1"]):
                points.append(item["p1"])
            points.append(item["p2"])
    return points


def extract_circle_center(d):
    """Extract center of a small ellipse/circle (2 or 4 bezier curves)."""
    if d["n_items"] in (2, 4) and all(item["op"] == "curve" for item in d["items"]):
        # Check it's small (data point sized, not a large shape)
        r = d["rect"]
        if r is None:
            return None
        w, h = r[2] - r[0], r[3] - r[1]
        if w > 10 or h > 10:
            return None
        # Center = midpoint of bounding rect
        return ((r[0] + r[2]) / 2, (r[1] + r[3]) / 2)
    return None


def extract_panel_a(drawings):
    """Extract individual trajectories from Panel A."""
    # Calibrate axes
    ax, bx, rx = fit_linear_transform(PANEL_A["x_ticks"])
    ay, by, ry = fit_log_transform(PANEL_A["y_ticks"])

    print("Panel A axis calibration:")
    print(f"  X (linear): pdf_x = {ax:.4f} * days + {bx:.2f}")
    print(f"    Residuals: {rx}")
    print(f"  Y (log10): pdf_y = {ay:.2f} * log10(EU) + {by:.2f}")
    print(f"    Residuals: {ry}")

    # Find all line-containing drawings in Panel A
    trajectories = []
    overall_line = None

    for d in drawings:
        if not in_panel(d, PANEL_A):
            continue

        # Skip non-line content (rects, circles)
        has_lines = any(item["op"] == "line" for item in d["items"])
        if not has_lines:
            continue

        # Skip grid lines and axes (gray #eaeaea or black #000000)
        hex_color = color_to_hex(d["stroke_color"])
        if hex_color in ("#eaeaea", "#000000", None):
            continue

        points = extract_line_points(d)
        if len(points) < 2:
            continue

        # Classify by line width
        w = d["line_width"]
        if w > 2.5:
            # Bold overall/median line
            traj_type = "overall"
        elif w > 1.5:
            traj_type = "age_stratum"
        else:
            traj_type = "individual"

        # Convert to data coordinates
        data_points = []
        for px, py in points:
            days = pdf_to_data_linear(px, ax, bx)
            eu = pdf_to_data_log(py, ay, by)
            data_points.append((days, eu))

        entry = {
            "primitive_id": d["primitive_id"],
            "trajectory_type": traj_type,
            "color_hex": hex_color,
            "line_width": w,
            "n_points": len(data_points),
            "points_pdf": points,
            "points_data": data_points,
        }

        if traj_type == "overall":
            overall_line = entry
        else:
            trajectories.append(entry)

    print(f"\n  Extracted {len(trajectories)} trajectories "
          f"({'+ overall line' if overall_line else 'no overall line'})")

    # Count types
    types = {}
    for t in trajectories:
        types.setdefault(t["trajectory_type"], 0)
        types[t["trajectory_type"]] += 1
    for k, v in types.items():
        print(f"    {k}: {v}")

    return trajectories, overall_line, {"x": (ax, bx, rx), "y": (ay, by, ry)}


def extract_panel_bc_points(drawings, panel, panel_label, y_scale="log"):
    """Extract individual data points (circles) from Panel B or C."""
    # Calibrate axes
    if y_scale == "log":
        ay, by, ry = fit_log_transform(panel["y_ticks"])
    else:
        ay, by, ry = fit_linear_transform(panel["y_ticks"])

    print(f"\nPanel {panel_label} axis calibration:")
    print(f"  Y ({'log10' if y_scale == 'log' else 'linear'}): pdf_y = {ay:.4f} * val + {by:.2f}")
    print(f"    Residuals (max): {max(abs(ry)):.4f}")

    points = []
    for d in drawings:
        if not in_panel(d, panel):
            continue

        center = extract_circle_center(d)
        if center is None:
            continue

        cx, cy = center
        hex_color = color_to_hex(d["fill_color"])

        # Take only filled versions (stroked duplicates have fill=None)
        if d["fill_color"] is None:
            continue

        # Skip axis-colored items
        if hex_color in ("#eaeaea", "#000000", "#ffffff", "#4d4d4d"):
            continue

        # Convert y to data value
        if y_scale == "log":
            data_y = pdf_to_data_log(cy, ay, by)
        else:
            data_y = pdf_to_data_linear(cy, ay, by)

        # Map x to category
        category = nearest_category(cx, panel["x_categories"])

        points.append({
            "primitive_id": d["primitive_id"],
            "pdf_x": cx,
            "pdf_y": cy,
            "category": category,
            "data_y": data_y,
            "color_hex": hex_color,
        })

    print(f"  Extracted {len(points)} data points")
    cats = {}
    for p in points:
        cats.setdefault(p["category"], 0)
        cats[p["category"]] += 1
    for k, v in sorted(cats.items(), key=lambda x: str(x)):
        print(f"    {k}: {v}")

    return points, {"y": (ay, by, ry)}


# =============================================================================
# Output
# =============================================================================

def write_panel_a_csv(trajectories, overall_line):
    """Write Panel A trajectories to CSV."""
    path = OUT_DIR / "panel_A_trajectories.csv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["trajectory_id", "trajectory_type", "color_hex",
                         "point_index", "days_since_fever_onset", "elisa_units"])
        tid = 0
        all_entries = trajectories + ([overall_line] if overall_line else [])
        for entry in all_entries:
            for i, (days, eu) in enumerate(entry["points_data"]):
                writer.writerow([tid, entry["trajectory_type"], entry["color_hex"],
                                 i, f"{days:.1f}", f"{eu:.1f}"])
            tid += 1
    print(f"\nWrote {path}")
    return path


# Panel B age group color mapping (from Phase 1 x-position clustering)
AGE_GROUP_COLORS = {
    "#f38efc": "<5",
    "#9133be": "5-15",
    "#000080": "16+",
}


def write_panel_b_csv(points):
    """Write Panel B data points to CSV with age group from color."""
    path = OUT_DIR / "panel_B_datapoints.csv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["point_id", "catchment", "age_group", "color_hex", "elisa_units"])
        for i, pt in enumerate(points):
            age_group = AGE_GROUP_COLORS.get(pt["color_hex"], "unknown")
            writer.writerow([i, pt["category"], age_group, pt["color_hex"],
                             f"{pt['data_y']:.1f}"])
    print(f"Wrote {path}")
    return path


def write_panel_c_csv(points):
    """Write Panel C data points to CSV."""
    path = OUT_DIR / "panel_C_datapoints.csv"
    with open(path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["point_id", "antigen_isotype", "color_hex", "boosting_ratio"])
        for i, pt in enumerate(points):
            writer.writerow([i, pt["category"], pt["color_hex"],
                             f"{pt['data_y']:.2f}"])
    print(f"Wrote {path}")
    return path


def generate_verification_plots(trajectories, overall_line, panel_b_points,
                                 panel_c_points, page_image_path):
    """Generate verification plots comparing extracted data to original."""
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))

    # Original figure image
    ax = axes[0, 0]
    img = plt.imread(str(page_image_path))
    ax.imshow(img)
    ax.set_title("Original Figure S11")
    ax.axis("off")

    # Panel A verification
    ax = axes[0, 1]
    for t in trajectories:
        pts = t["points_data"]
        days = [p[0] for p in pts]
        eus = [p[1] for p in pts]
        ax.plot(days, eus, color=t["color_hex"], linewidth=0.8, alpha=0.6)
    if overall_line:
        pts = overall_line["points_data"]
        days = [p[0] for p in pts]
        eus = [p[1] for p in pts]
        ax.plot(days, eus, color=overall_line["color_hex"], linewidth=2)
    ax.set_xlabel("Days since fever onset")
    ax.set_ylabel("ELISA units")
    ax.set_yscale("log")
    ax.set_title(f"Panel A: {len(trajectories)} trajectories extracted")
    ax.set_xlim(0, 600)
    ax.set_ylim(50, 2000)
    ax.grid(True, alpha=0.3)

    # Panel B verification — box/whisker by age group within each catchment
    ax = axes[1, 0]
    age_groups = ["<5", "5-15", "16+"]
    age_colors = {"<5": "#f38efc", "5-15": "#9133be", "16+": "#000080"}
    catchments = ["Kathmandu", "Kavre"]

    # Build grouped data
    b_data = {}  # (catchment, age_group) -> [values]
    for pt in panel_b_points:
        ag = AGE_GROUP_COLORS.get(pt["color_hex"], "unknown")
        key = (pt["category"], ag)
        b_data.setdefault(key, []).append(pt["data_y"])

    positions = []
    box_data = []
    box_colors = []
    labels = []
    pos = 0
    for ci, catch in enumerate(catchments):
        for ai, ag in enumerate(age_groups):
            key = (catch, ag)
            vals = b_data.get(key, [])
            if vals:
                positions.append(pos)
                box_data.append(vals)
                box_colors.append(age_colors[ag])
                labels.append(f"{ag}")
                # Jittered scatter underneath
                jitter = np.random.RandomState(int(pos * 10)).uniform(-0.15, 0.15, len(vals))
                ax.scatter(pos + jitter, vals, color=age_colors[ag],
                           s=12, alpha=0.4, zorder=1)
            pos += 1
        pos += 0.5  # gap between catchments

    bp = ax.boxplot(box_data, positions=positions, widths=0.6, patch_artist=True,
                    showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    # Label catchments
    catch_centers = [1, 4.5]  # midpoints of each group
    ax.set_xticks(catch_centers)
    ax.set_xticklabels(catchments)
    # Add age group sub-labels
    for p, lbl in zip(positions, labels):
        ax.text(p, ax.get_ylim()[0] * 0.85, lbl, ha="center", va="top",
                fontsize=6, rotation=45)
    ax.set_ylabel("ELISA units")
    ax.set_yscale("log")
    ax.set_title(f"Panel B: {len(panel_b_points)} points, by age group")
    ax.set_ylim(50, 2000)
    ax.grid(True, alpha=0.3, axis="y")
    # Legend
    from matplotlib.patches import Patch
    legend_patches = [Patch(facecolor=age_colors[ag], alpha=0.5, label=ag)
                      for ag in age_groups]
    ax.legend(handles=legend_patches, fontsize=7, title="Age", title_fontsize=7)

    # Panel C verification — box/whisker by antigen-isotype
    ax = axes[1, 1]
    cat_labels = ["HlyE IgA", "HlyE IgG", "LPS IgA", "LPS IgG", "Vi IgG"]
    # Group data by category
    c_data = {}
    c_colors = {}
    for pt in panel_c_points:
        cat = pt["category"]
        c_data.setdefault(cat, []).append(pt["data_y"])
        c_colors[cat] = pt["color_hex"]

    c_box_data = []
    c_box_colors = []
    for cat in cat_labels:
        vals = c_data.get(cat, [])
        c_box_data.append(vals)
        c_box_colors.append(c_colors.get(cat, "#999999"))
        # Jittered scatter
        if vals:
            idx = cat_labels.index(cat)
            jitter = np.random.RandomState(idx).uniform(-0.2, 0.2, len(vals))
            ax.scatter(idx + jitter, vals, color=c_colors.get(cat, "#999999"),
                       s=12, alpha=0.4, zorder=1)

    bp = ax.boxplot(c_box_data, positions=range(len(cat_labels)), widths=0.5,
                    patch_artist=True, showfliers=False, zorder=2,
                    medianprops=dict(color="black", linewidth=1.5))
    for patch, color in zip(bp["boxes"], c_box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    ax.set_xticks(range(len(cat_labels)))
    ax.set_xticklabels(cat_labels, rotation=30, ha="right")
    ax.set_ylabel("Ratio (1-month peak / population mean)")
    ax.set_title(f"Panel C: {len(panel_c_points)} points extracted")
    ax.set_ylim(-2, 55)
    ax.grid(True, alpha=0.3, axis="y")

    # Vi IgG annotation
    vi_vals = c_data.get("Vi IgG", [])
    if vi_vals:
        med = np.median(vi_vals)
        q25, q75 = np.percentile(vi_vals, [25, 75])
        ax.annotate(f"Vi IgG: median={med:.1f}\nIQR [{q25:.1f}, {q75:.1f}]",
                    xy=(4, med), xytext=(3.2, 35),
                    arrowprops=dict(arrowstyle="->", color="blue", alpha=0.5),
                    fontsize=8, color="blue")

    plt.tight_layout()
    out_path = OUT_DIR / "verification_plots.png"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Wrote {out_path}")


def write_calibration_report(cal_a, cal_b, cal_c):
    """Write calibration parameters to JSON."""
    report = {
        "panel_A": {
            "x_axis": {"type": "linear", "a": cal_a["x"][0], "b": cal_a["x"][1],
                        "residuals": cal_a["x"][2].tolist(),
                        "max_residual_pts": float(max(abs(cal_a["x"][2])))},
            "y_axis": {"type": "log10", "a": cal_a["y"][0], "b": cal_a["y"][1],
                        "residuals": cal_a["y"][2].tolist()},
        },
        "panel_B": {
            "y_axis": {"type": "log10", "a": cal_b["y"][0], "b": cal_b["y"][1],
                        "residuals": cal_b["y"][2].tolist()},
        },
        "panel_C": {
            "y_axis": {"type": "linear", "a": cal_c["y"][0], "b": cal_c["y"][1],
                        "residuals": cal_c["y"][2].tolist(),
                        "max_residual_pts": float(max(abs(cal_c["y"][2])))},
        },
    }
    path = OUT_DIR / "calibration_report.json"
    with open(path, "w") as f:
        json.dump(report, f, indent=2)
    print(f"Wrote {path}")


def main():
    # Load raw data
    with open(RAW_DIR / "all_drawings.json") as f:
        drawings = json.load(f)
    print(f"Loaded {len(drawings)} drawings")

    # Extract Panel A
    print("\n" + "=" * 60)
    print("PANEL A: Individual Vi IgG trajectories")
    print("=" * 60)
    trajectories, overall_line, cal_a = extract_panel_a(drawings)
    write_panel_a_csv(trajectories, overall_line)

    # Extract Panel B
    print("\n" + "=" * 60)
    print("PANEL B: Vi IgG by age group and catchment")
    print("=" * 60)
    panel_b_points, cal_b = extract_panel_bc_points(
        drawings, PANEL_B, "B", y_scale="log")
    write_panel_b_csv(panel_b_points)

    # Extract Panel C
    print("\n" + "=" * 60)
    print("PANEL C: Boosting ratios by antigen-isotype")
    print("=" * 60)
    panel_c_points, cal_c = extract_panel_bc_points(
        drawings, PANEL_C, "C", y_scale="linear")
    write_panel_c_csv(panel_c_points)

    # Validation summaries
    print("\n" + "=" * 60)
    print("VALIDATION CHECKS")
    print("=" * 60)
    vi_ratios = [pt["data_y"] for pt in panel_c_points if pt["category"] == "Vi IgG"]
    if vi_ratios:
        print(f"Vi IgG boosting ratios: n={len(vi_ratios)}, "
              f"median={np.median(vi_ratios):.2f}, "
              f"IQR=[{np.percentile(vi_ratios, 25):.2f}, {np.percentile(vi_ratios, 75):.2f}]")
        print(f"  Expected from paper: median ~1.5, IQR ~1.0-2.1")

    if trajectories:
        all_eus = [eu for t in trajectories for _, eu in t["points_data"]]
        print(f"Panel A ELISA range: {min(all_eus):.0f} - {max(all_eus):.0f} EU")
        print(f"  Expected: ~100 - 1000+ EU")

    # Calibration report
    write_calibration_report(cal_a, cal_b, cal_c)

    # Verification plots
    generate_verification_plots(
        trajectories, overall_line, panel_b_points, panel_c_points,
        RAW_DIR / "page_image.png")

    print("\nPhase 2 complete.")


if __name__ == "__main__":
    main()
