"""
Phase 1: Raw extraction of Figure S11 from Aiemjoy 2022 supplement.

Extracts all vector primitives and text from the page containing Figure S11,
saves raw inventory, and generates diagnostic overlay for inspection.
"""

import json
import os
from collections import Counter
from pathlib import Path

import fitz  # PyMuPDF
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Paths
PDF_PATH = Path("metastudies/waning_after_infection/input_papers/Aiemjoy et al._2022_Estimating typhoid incidence from community-based serosurveys a multicohort study -- supplementary appendix 4.pdf")
OUT_DIR = Path("metastudies/waning_after_infection/Aiemjoy2022_reanalysis/fig_s11_extraction/raw")
OUT_DIR.mkdir(parents=True, exist_ok=True)


def find_figure_page(doc):
    """Find the page containing Figure S11 by text search + drawings."""
    # Find pages with "S11" text AND significant drawing content
    candidates = []
    for i, page in enumerate(doc):
        hits = page.search_for("S11")
        n_drawings = len(page.get_drawings())
        if hits and n_drawings > 100:
            candidates.append((i, n_drawings, hits))
            print(f"Candidate page {i}: {n_drawings} drawings, {len(hits)} S11 hit(s)")
    if candidates:
        # Pick the candidate with the most drawings (the actual figure, not TOC)
        best = max(candidates, key=lambda x: x[1])
        print(f"Selected page {best[0]} (0-indexed) with {best[1]} drawings")
        return best[0]
    raise ValueError("Could not find Figure S11 in the PDF")


def serialize_drawing(idx, d):
    """Convert a PyMuPDF drawing dict to a JSON-serializable record."""
    items = []
    for item in d.get("items", []):
        op = item[0]
        if op == "l":  # line
            p1, p2 = item[1], item[2]
            items.append({"op": "line", "p1": [p1.x, p1.y], "p2": [p2.x, p2.y]})
        elif op == "c":  # cubic bezier
            pts = [{"x": item[j].x, "y": item[j].y} for j in range(1, 5)]
            items.append({"op": "curve", "points": pts})
        elif op == "re":  # rectangle
            r = item[1]
            items.append({"op": "rect", "rect": [r.x0, r.y0, r.x1, r.y1]})
        elif op == "qu":  # quad
            q = item[1]
            items.append({"op": "quad", "points": [[q[j].x, q[j].y] for j in range(4)]})
        else:
            items.append({"op": str(op), "raw": str(item[1:])})

    rect = d.get("rect")
    return {
        "primitive_id": idx,
        "items": items,
        "n_items": len(items),
        "stroke_color": list(d["color"]) if d.get("color") else None,
        "fill_color": list(d["fill"]) if d.get("fill") else None,
        "line_width": d.get("width"),
        "opacity": d.get("stroke_opacity", d.get("opacity")),
        "fill_opacity": d.get("fill_opacity"),
        "even_odd": d.get("even_odd"),
        "closePath": d.get("closePath"),
        "rect": [rect.x0, rect.y0, rect.x1, rect.y1] if rect else None,
        "layer": d.get("layer"),
    }


def serialize_text(text_dict):
    """Flatten PyMuPDF text dict into a list of span records."""
    spans = []
    for block in text_dict.get("blocks", []):
        if block["type"] != 0:  # skip image blocks
            continue
        for line in block.get("lines", []):
            for span in line.get("spans", []):
                spans.append({
                    "text": span["text"],
                    "bbox": list(span["bbox"]),
                    "font": span["font"],
                    "size": span["size"],
                    "color": span["color"],
                    "origin": list(span["origin"]),
                    "flags": span["flags"],
                })
    return spans


def color_to_hex(c):
    """Convert RGB float tuple to hex string."""
    if c is None:
        return "#000000"
    return "#{:02x}{:02x}{:02x}".format(
        int(c[0] * 255), int(c[1] * 255), int(c[2] * 255)
    )


def summarize_primitives(drawings):
    """Print summary statistics about extracted primitives."""
    print(f"\n{'='*60}")
    print(f"PRIMITIVE SUMMARY")
    print(f"{'='*60}")
    print(f"Total drawing entries: {len(drawings)}")

    # Count operations
    op_counts = Counter()
    for d in drawings:
        for item in d["items"]:
            op_counts[item["op"]] += 1
    print(f"\nOperation counts:")
    for op, count in op_counts.most_common():
        print(f"  {op}: {count}")

    # Item counts per drawing (path complexity)
    item_counts = [d["n_items"] for d in drawings]
    print(f"\nItems per drawing (path complexity):")
    print(f"  min: {min(item_counts)}, max: {max(item_counts)}, "
          f"median: {sorted(item_counts)[len(item_counts)//2]}")
    complexity_dist = Counter()
    for c in item_counts:
        if c == 1:
            complexity_dist["1 (simple)"] += 1
        elif c <= 4:
            complexity_dist["2-4 (small path)"] += 1
        elif c <= 10:
            complexity_dist["5-10 (medium path)"] += 1
        else:
            complexity_dist[f"11+ (complex path)"] += 1
    for k, v in sorted(complexity_dist.items()):
        print(f"  {k}: {v}")

    # Color distribution
    stroke_colors = Counter()
    for d in drawings:
        c = d["stroke_color"]
        if c:
            stroke_colors[color_to_hex(c)] += 1
    print(f"\nStroke color distribution (top 15):")
    for color, count in stroke_colors.most_common(15):
        print(f"  {color}: {count}")

    fill_colors = Counter()
    for d in drawings:
        c = d["fill_color"]
        if c:
            fill_colors[color_to_hex(c)] += 1
    print(f"\nFill color distribution (top 15):")
    for color, count in fill_colors.most_common(15):
        print(f"  {color}: {count}")

    # Line width distribution
    widths = Counter()
    for d in drawings:
        w = d["line_width"]
        if w is not None:
            widths[round(w, 3)] += 1
    print(f"\nLine width distribution:")
    for w, count in sorted(widths.items()):
        print(f"  {w}: {count}")

    # Multi-segment lines (potential trajectories)
    multi_line_paths = [d for d in drawings if d["n_items"] > 1
                        and all(item["op"] == "line" for item in d["items"])]
    print(f"\nMulti-segment line paths (potential trajectories): {len(multi_line_paths)}")
    if multi_line_paths:
        segs = [d["n_items"] for d in multi_line_paths]
        print(f"  Segment counts: {sorted(segs)}")


def generate_diagnostic_overlay(page, drawings, text_spans):
    """Generate a diagnostic overlay showing all extracted primitives."""
    # Render page as background
    pix = page.get_pixmap(dpi=150)
    img_data = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.h, pix.w, pix.n)

    fig, axes = plt.subplots(1, 2, figsize=(24, 14))

    # Left: page image with drawing overlays
    ax = axes[0]
    ax.imshow(img_data, extent=[0, pix.w, pix.h, 0])
    ax.set_title(f"Extracted drawings (n={len(drawings)})")

    scale = 150 / 72  # DPI / PDF points per inch

    for d in drawings:
        color = color_to_hex(d["stroke_color"]) if d["stroke_color"] else "#cccccc"
        lw = max(0.3, (d["line_width"] or 0.5) * scale * 0.5)
        for item in d["items"]:
            if item["op"] == "line":
                x = [item["p1"][0] * scale, item["p2"][0] * scale]
                y = [item["p1"][1] * scale, item["p2"][1] * scale]
                ax.plot(x, y, color=color, linewidth=lw, alpha=0.7)
            elif item["op"] == "rect":
                r = item["rect"]
                rect_patch = mpatches.Rectangle(
                    (r[0] * scale, r[1] * scale),
                    (r[2] - r[0]) * scale, (r[3] - r[1]) * scale,
                    linewidth=lw, edgecolor=color, facecolor='none', alpha=0.5
                )
                ax.add_patch(rect_patch)
            elif item["op"] == "curve":
                pts = item["points"]
                xs = [p["x"] * scale for p in pts]
                ys = [p["y"] * scale for p in pts]
                ax.plot(xs, ys, color=color, linewidth=lw, alpha=0.5, linestyle='--')

    ax.set_xlim(0, pix.w)
    ax.set_ylim(pix.h, 0)

    # Right: page image with text overlays
    ax2 = axes[1]
    ax2.imshow(img_data, extent=[0, pix.w, pix.h, 0])
    ax2.set_title(f"Extracted text spans (n={len(text_spans)})")

    for span in text_spans:
        bbox = span["bbox"]
        x0, y0, x1, y1 = [c * scale for c in bbox]
        rect_patch = mpatches.Rectangle(
            (x0, y0), x1 - x0, y1 - y0,
            linewidth=0.5, edgecolor='red', facecolor='yellow', alpha=0.3
        )
        ax2.add_patch(rect_patch)
        ax2.text(x0, y0 - 2, span["text"], fontsize=4, color='red', alpha=0.8)

    ax2.set_xlim(0, pix.w)
    ax2.set_ylim(pix.h, 0)

    plt.tight_layout()
    out_path = OUT_DIR / "diagnostic_overlay.png"
    fig.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved diagnostic overlay: {out_path}")


def main():
    doc = fitz.open(PDF_PATH)
    print(f"Opened PDF: {PDF_PATH}")
    print(f"Total pages: {len(doc)}")

    page_idx = find_figure_page(doc)
    page = doc[page_idx]
    print(f"\nPage dimensions: {page.rect.width:.1f} x {page.rect.height:.1f} points")

    # Extract drawings
    print("\nExtracting drawings...")
    raw_drawings = page.get_drawings()
    drawings = [serialize_drawing(i, d) for i, d in enumerate(raw_drawings)]
    print(f"Extracted {len(drawings)} drawing entries")

    # Extract text
    print("Extracting text...")
    text_dict = page.get_text("dict")
    text_spans = serialize_text(text_dict)
    print(f"Extracted {len(text_spans)} text spans")

    # Save raw JSON
    drawings_path = OUT_DIR / "all_drawings.json"
    with open(drawings_path, "w") as f:
        json.dump(drawings, f, indent=2)
    print(f"Saved: {drawings_path}")

    text_path = OUT_DIR / "all_text_blocks.json"
    with open(text_path, "w") as f:
        json.dump(text_spans, f, indent=2)
    print(f"Saved: {text_path}")

    # Save page image
    pix = page.get_pixmap(dpi=300)
    img_path = OUT_DIR / "page_image.png"
    pix.save(str(img_path))
    print(f"Saved: {img_path}")

    # Summary stats
    if not drawings:
        print("\nWARNING: No drawings found on this page!")
        print("The figure may be rasterized or on a different page.")
        doc.close()
        return
    summarize_primitives(drawings)

    # Print text content for inspection
    print(f"\n{'='*60}")
    print("TEXT CONTENT")
    print(f"{'='*60}")
    for span in text_spans:
        bbox = span["bbox"]
        print(f"  [{bbox[0]:.0f},{bbox[1]:.0f}]-[{bbox[2]:.0f},{bbox[3]:.0f}] "
              f"size={span['size']:.1f} '{span['text']}'")

    # Diagnostic overlay
    print("\nGenerating diagnostic overlay...")
    generate_diagnostic_overlay(page, drawings, text_spans)

    doc.close()
    print("\nPhase 1 complete.")


if __name__ == "__main__":
    main()
