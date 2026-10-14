#!/usr/bin/env python3
"""
figS5_realignment_accuracy.py - Figure S5: is lifting a BAM over as good as
re-aligning the reads?

Figure S2 shows liftover is ~100x faster and ~300x lighter than re-alignment.
That only matters if the coordinates it produces are as good. This figure is
the accuracy half of that claim: the same simulated reads go through BWA-MEM
against hg38 from scratch, and the two coordinate sets are compared read by
read, so every number here is a property of the method rather than of a
particular file. Two scales, 1M and 10M read pairs, because a concordance
number that only holds at the smaller scale is not a property of the method.

Two panels, answering the two questions a reader actually has:

  (a) the retention cascade -- of every read pair both tools mapped, what
      fraction survives each successively stricter test? Three rungs, read left
      to right, each drawn as a bar and each carrying its own value: all mapped
      -> same chromosome -> the exact same start base. The bars are steeped
      (both scales side by side, which is why the two bars of a rung look
      identical) and each rung is annotated with what it costs, in percentage
      points, so the shape of the loss is visible rather than just its endpoint.
      Drawn as irregular steps rather than a zero-based bar chart because all
      three values sit above 98.5%: from 0 they would be three identical-looking
      full-height bars, which shows nothing.

  (b) the discrepancy attribution -- where the reads that fail the last step
      actually went. This is the panel that answers the question a reviewer
      asks next: "is the 1.44% that does not match exactly evidence of drift in
      your coordinate arithmetic?" The answer is 99.5% no: the disagreements
      are almost entirely reads that were placed somewhere else *entirely* --
      57% land on a different chromosome, 43% land more than 500 bp away, both
      the signature of a read sitting in a region the two assemblies order
      differently or of a paralogue BWA-MEM resolved the other way. Only 0.5%
      of disagreements (0.008% of all reads) sit within 500 bp of the liftover
      coordinate, which is the band a buggy offset calculation would produce.

The offset bands are computed on same-chromosome pairs only (the benchmark
script defines an offset only when both coordinates are on one chromosome), so
"different chromosome" is counted separately rather than as an infinite offset.
The three categories are exhaustive and disjoint by construction: they are
read off the same count, and the script asserts the buckets re-sum to the
chromosome-concordant total before drawing anything.

WHY STRAND IS NOT A STAGE IN THE CASCADE. Strand concordance (99.58%) is
recorded alongside the other two, and the obvious thing to do is add it as a
fourth rung. It cannot be: it is a separate marginal over all both-mapped pairs,
not a subset of chromosome concordance, and it is numerically *higher* than
chromosome concordance (99.58% > 99.18%), so no nesting is possible. A joint
tabulation from the 1M BAMs shows why -- of the 1M pairs, 13,595 reads (0.68%)
land on a different chromosome while keeping their strand, while 5,757 (0.29%)
stay on the chromosome and flip. The two properties are near-independent. Any
"same chromosome AND same strand" rung would therefore have to be a *joint*
count, which the benchmark JSON does not record for the 10M scale, and drawing
the 1M joint beside 10M marginals would mix two different denominators in one
cascade. The strand number is reported as a side note instead, where it is
unambiguously a marginal. (That joint also shows every exact-position match is
a same-strand match: strand is a strict prerequisite of the last rung, which is
why the last rung needs no strand qualifier.)

Axis note for (a): the y-axis starts at 98%, not 0. Every value in the cascade
is above 98.5%, so a zero-based axis would draw three bars at indistinguishable
lengths. The truncation is labeled on the panel, and each stage's value is
printed at its marker, so nothing has to be read off the axis by eye.

Input:   paper/results/accuracy_realignment_1M.json
         paper/results/accuracy_realignment_10M.json
Output:  paper/figures/figS5_realignment_accuracy.{pdf,png}
         paper/results/realignment_accuracy_summary.json

Usage: python paper/plot/figS5_realignment_accuracy.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path("paper/results")
FIGURES_DIR = Path("paper/figures")
FIGURES_DIR.mkdir(parents=True, exist_ok=True)

SOURCES = [
    ("1M",  RESULTS_DIR / "accuracy_realignment_1M.json",  "#1f77b4"),
    ("10M", RESULTS_DIR / "accuracy_realignment_10M.json", "#ff7f0e"),
]

# the three rungs of the cascade, top to bottom. Each is a share of the SAME
# denominator (all both-mapped pairs), so the sequence is genuinely nested:
# an exact-start match is necessarily a same-chromosome match.
STAGES = [
    ("both mapped", 100.0),
    ("same chromosome", "__chrom__"),
    ("exact start (0 bp offset)", "__exact__"),
]

# everything the cascade loses between "both mapped" and "exact start", split
# into the three mutually exclusive places a read can end up. Order is by
# reading importance: the two "far" categories first, the benign band last.
CATEGORIES = [
    ("different chromosome", "#9467bd"),
    ("same chromosome, >500 bp", "#c44e52"),
    ("same chromosome, 1-500 bp", "#55a868"),
]

# the offset bands the benchmark records for same-chromosome pairs
NONZERO_BUCKETS = ["1-5", "6-10", "11-50", "51-100", "101-500", ">500"]

# the third category is a sliver and is marked rather than filled, so it needs
# an edge colour that reads against both the red segment beside it and the
# white page
SLIVER_EDGE = "#1b4d2b"

# the buckets are counts while the concordance rates are rounded to 3 decimals,
# so the two can differ by a handful of reads out of millions. Anything larger
# means the file's own numbers disagree with each other.
SUM_TOLERANCE = 0.001   # 0.1% of the pair count


def load():
    out = []
    for label, path, colour in SOURCES:
        if not path.exists():
            raise SystemExit(f"missing input: {path}")
        with open(path) as fh:
            doc = json.load(fh)
        out.append((label, colour, doc["result"]))
    return out


def cascade(res) -> list[float]:
    """The three retention shares, each a percentage of all both-mapped pairs."""
    exact = res["concordance_at_threshold"]["<=0bp"]
    return [100.0, res["chrom_concordance"], exact]


def attribution(res) -> dict:
    """
    The not-exact reads, split into three disjoint categories (counts).

    same_chrom_total is the sum of this format's own offset buckets, which the
    benchmark computes over same-chromosome pairs only; "different chromosome"
    is what is left of the both-mapped total, so the three parts sum back to
    (both_mapped - exact) exactly.
    """
    n = res["both_mapped"]
    dist = res["offset_distribution"]
    same_chrom = sum(dist.get(b, 0) for b in ["0"] + NONZERO_BUCKETS)
    far_within = dist.get(">500", 0)
    return {
        "both_mapped": n,
        "exact": dist.get("0", 0),
        "different chromosome": n - same_chrom,
        "same chromosome, >500 bp": far_within,
        "same chromosome, 1-500 bp": same_chrom - dist.get("0", 0) - far_within,
    }


def check_claim(datasets):
    """
    Assert the two things the figure asserts before drawing it.

    1. The file's own numbers are self-consistent: the offset buckets must
       re-sum to the chromosome-concordant count implied by chrom_concordance.
       They are computed from the same comparison, so a mismatch means the JSON
       was assembled from two different runs and nothing drawn from it is safe.
    2. The claim panel (b) makes -- that the 1-500 bp band, the one that would
       indicate bad coordinate arithmetic, is negligible -- actually holds.
    """
    problems = []
    for label, _, res in datasets:
        n = res["both_mapped"]
        attr = attribution(res)
        implied = n * res["chrom_concordance"] / 100.0
        observed = n - attr["different chromosome"]
        if abs(implied - observed) > SUM_TOLERANCE * n:
            problems.append(
                f"{label}: offset buckets sum to {observed:,} but chrom_concordance "
                f"{res['chrom_concordance']}% implies {implied:,.0f}"
            )
        not_exact = n - attr["exact"]
        band = attr["same chromosome, 1-500 bp"]
        if band / not_exact > 0.01:
            problems.append(
                f"{label}: the 1-500 bp band is {100 * band / not_exact:.2f}% of "
                "discrepancies -- too large to call the disagreement non-local"
            )
        # the two panels must be complements: panel (a) draws the exact-match
        # share, panel (b) the complement of it, both as a percentage of the
        # same denominator. They are read off different quantities (a stored
        # rate vs a count), so they can drift apart at the last digit and look
        # like a contradiction in print.
        pct_not_exact = 100.0 * not_exact / n
        if abs(res["concordance_at_threshold"]["<=0bp"] + pct_not_exact - 100.0) > 0.002:
            problems.append(
                f"{label}: panel (a) shows "
                f"{res['concordance_at_threshold']['<=0bp']:.3f}% exact while panel "
                f"(b) shows {pct_not_exact:.3f}% not exact -- they do not sum to 100"
            )
    if problems:
        raise SystemExit(
            "refusing to draw Figure S5:\n  " + "\n  ".join(problems)
        )
    print("claim holds: each file is self-consistent, and the 1-500 bp band is "
          "under 1% of discrepancies at both scales")


def panel_cascade(ax, datasets):
    """(a) the retention cascade, both scales, as grouped bars per rung."""
    x = np.arange(len(STAGES))
    n = len(datasets)
    width = 0.34
    base = 98.40   # the y-floor every bar is drawn from (the axis is truncated)

    for i, (label, colour, res) in enumerate(datasets):
        ys = cascade(res)
        offs = x + (i - (n - 1) / 2) * width
        ax.bar(offs, [y - base for y in ys], width, bottom=base,
               color=colour, edgecolor="white", lw=0.8,
               label=f"{label} read pairs", zorder=3)
        for xi, y in zip(offs, ys):
            # each bar carries its own value: the bars differ by well under 1 pp,
            # so the number is the result and the bar is only there to hold it.
            # THREE decimals, not two: at two the two scales round to the same
            # string (98.56 / 98.56) while panel (b) reports their complements
            # as 1.445 and 1.442 -- which reads as a contradiction between the
            # panels. At three the two scales differ where the data differs,
            # and each panel's value plus the other's is exactly 100.
            ax.text(xi, y + 0.03, f"{y:.3f}", ha="center", va="bottom",
                    fontsize=8, color="#333333", zorder=4)

    # what each rung costs, in percentage points, averaged over the two runs
    # (they agree to 0.004 pp, so the two bars of a rung differ by a hair). The
    # loss belongs in the tick label: anywhere inside the axes it would land on
    # a bar or on the note.
    chrom = np.mean([cascade(res)[1] for _, _, res in datasets])
    exact = np.mean([cascade(res)[2] for _, _, res in datasets])
    losses = [None, 100.0 - chrom, chrom - exact]

    labels = []
    for (name, _), loss in zip(STAGES, losses):
        labels.append(name if loss is None else f"{name}\n(−{loss:.2f} pp)")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_xlim(-0.6, len(STAGES) - 0.4)
    ax.set_ylim(base, 100.55)
    ax.set_ylabel("retained (% of all read pairs both tools mapped)")
    ax.set_title("(a) retention of both-mapped reads", fontsize=10, loc="left")
    ax.grid(alpha=0.25, axis="y", lw=0.5)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.9)


def panel_attribution(ax, datasets):
    """(b) where the reads that miss the last rung actually went."""
    short = {"different chromosome": "different chromosome",
             "same chromosome, >500 bp": "same chromosome,\n>500 bp away",
             "same chromosome, 1-500 bp": "1-500 bp"}
    for row, (label, _, res) in enumerate(datasets):
        attr = attribution(res)
        not_exact = attr["both_mapped"] - attr["exact"]
        left = 0.0
        for name, colour in CATEGORIES:
            pct = 100.0 * attr[name] / not_exact
            ax.barh(row, pct, left=left, height=0.72, color=colour,
                    edgecolor="white", lw=0.9, zorder=3)
            # the two wide segments carry their own label, so the panel needs no
            # legend; the third is a sliver and is boxed instead
            if name != "same chromosome, 1-500 bp":
                ax.text(left + pct / 2, row, f"{short[name]}\n{pct:.1f}%",
                        ha="center", va="center", fontsize=8.5, color="white",
                        fontweight="bold", zorder=4)
            left += pct

        # the 1-500 bp band is ~0.5% of the bar -- under a pixel of ink, and
        # indistinguishable from the bar's own edge. Box it so a reader can see
        # that there is a third category and that the two shares shown are not
        # the whole: 56.8 + 42.6 = 99.4, and this is the rest. The box is drawn
        # taller than the bar so its top and bottom edges land on white and the
        # rectangle reads as a deliberate mark rather than as part of the bar.
        band_pct = 100.0 * attr["same chromosome, 1-500 bp"] / not_exact
        ax.add_patch(plt.Rectangle(
            (100.0 - band_pct, row - 0.46), band_pct, 0.92,
            fill=False, edgecolor=SLIVER_EDGE, lw=1.4,
            linestyle=(0, (3, 2)), zorder=6))

        ax.text(0, row + 0.42, f"{label}: {not_exact:,} not-exact reads "
                               f"({100 * not_exact / attr['both_mapped']:.3f}% of all)",
                fontsize=8, color="#333333", va="bottom", ha="left")

    # name the boxed sliver in the free band under the bars; a box with no label
    # would just be an unexplained mark, and the segment is far too narrow to
    # hold its own text. The percentages are panel (b)'s own denominator -- the
    # share of the reads that are not an exact match, not of all reads -- and
    # the counts are given per scale because the box appears on both bars.
    band_txt = " and ".join(
        f"{attr['same chromosome, 1-500 bp']:,} at {lab}"
        for lab, _, res in datasets
        for attr in [attribution(res)]
    )
    a0 = attribution(datasets[0][2])
    band_pct0 = 100.0 * a0["same chromosome, 1-500 bp"] / (a0["both_mapped"] - a0["exact"])
    ax.text(100.0, -0.60,
            f"<500 bp band: {band_pct0:.1f}% of the "
            f"disagreements ({band_txt} reads)",
            ha="right", va="center", fontsize=7.5, color=SLIVER_EDGE)

    # what the bars *mean* -- that a long-range relocation is a difference in how
    # the two assemblies order the region, or a paralogue BWA-MEM resolved the
    # other way, rather than drift in the coordinate arithmetic -- is prose and
    # belongs in the caption; the panel carries data only.

    # a little room to the right of 100 so the boxed sliver is not read as part
    # of the right spine
    ax.set_xlim(0, 103)
    ax.set_ylim(-0.78, 1.75)
    ax.set_yticks(range(len(datasets)))
    ax.set_yticklabels([lab for lab, _, _ in datasets], fontsize=9)
    ax.set_xlabel("share of the reads that are not an exact match (%)")
    ax.set_title("(b) the disagreements are almost all long-range relocations",
                 fontsize=10, loc="left")
    ax.grid(alpha=0.25, axis="x", lw=0.5)
    ax.set_axisbelow(True)


def main():
    datasets = load()
    check_claim(datasets)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.6, 4.7),
                                   gridspec_kw=dict(width_ratios=[1.0, 1.05]))
    panel_cascade(ax1, datasets)
    panel_attribution(ax2, datasets)

    fig.suptitle("Figure S5. BAM liftover vs. BWA-MEM re-alignment on the same "
                 "simulated reads (truth carried in the read names)",
                 fontsize=11, y=1.00)
    fig.tight_layout()

    out_pdf = FIGURES_DIR / "figS5_realignment_accuracy.pdf"
    out_png = FIGURES_DIR / "figS5_realignment_accuracy.png"
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    print(f"wrote {out_pdf}")
    print(f"wrote {out_png}")

    # machine-readable copy of exactly what was drawn
    payload = {"cascade_denominator": "share of all read pairs both tools mapped",
               "attribution_denominator": "share of the pairs that are not an "
                                          "exact-start match",
               "strand_note": "strand_concordance is a separate marginal over "
                              "all both-mapped pairs, not a rung of the cascade; "
                              "it is not nested inside chromosome concordance",
               "scales": {}}
    for label, _, res in datasets:
        attr = attribution(res)
        not_exact = attr["both_mapped"] - attr["exact"]
        payload["scales"][label] = {
            "both_mapped": attr["both_mapped"],
            "cascade_pct": {name: round(v, 3)
                            for (name, _), v in zip(STAGES, cascade(res))},
            "strand_concordance_pct": res["strand_concordance"],
            "not_exact": not_exact,
            "attribution_count": {k: attr[k] for k, _ in CATEGORIES},
            "attribution_pct": {k: round(100 * attr[k] / not_exact, 3)
                                for k, _ in CATEGORIES},
        }
    path = RESULTS_DIR / "realignment_accuracy_summary.json"
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"wrote {path}")

    print()
    for label, _, res in datasets:
        attr = attribution(res)
        not_exact = attr["both_mapped"] - attr["exact"]
        c = cascade(res)
        print(f"{label:>4}: mapped {attr['both_mapped']:,}  "
              f"same chrom {c[1]:.3f}%  exact start {c[2]:.3f}%  "
              f"strand(separate) {res['strand_concordance']:.3f}%")
        for name, _ in CATEGORIES:
            print(f"        {name:<26} {attr[name]:>9,}  "
                  f"{100 * attr[name] / not_exact:6.2f}% of the {not_exact:,} "
                  f"not-exact reads")
        far = attr["different chromosome"] + attr["same chromosome, >500 bp"]
        print(f"        long-range total           {far:>9,}  "
              f"{100 * far / not_exact:6.2f}%")


if __name__ == "__main__":
    main()
