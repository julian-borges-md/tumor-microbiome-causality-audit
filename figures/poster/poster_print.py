#!/usr/bin/env python3
"""
RO-2026-008 | Poster mock builder v1.0
Composites the real Figures 1-7 into a 48 x 36 inch landscape poster.

Outputs POSTER_RO-2026-008_44x44.png (100 dpi preview) and POSTER_RO-2026-008_44x44.pdf.
Layout parameters mirror POSTER_SPEC.md v1.0. Figure widths respect the
170 ppi floor documented there.

    python3 poster_mock.py
"""
import textwrap
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import FancyBboxPatch, Rectangle

# --- asset resolution -------------------------------------------------------
# Figures live in ../ ; the de-titled crops and the QR are derived here on
# demand into _derived/ (gitignored) so `make poster` works from a clean clone.
import os as _os, subprocess as _sp
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_FIGS = _os.path.normpath(_os.path.join(_HERE, ".."))
_DER = _os.path.join(_HERE, "_derived")
_os.makedirs(_DER, exist_ok=True)

REPO_URL = "https://github.com/julian-borges-md/tumor-microbiome-causality-audit"


def _crop_title(name):
    """Remove the figure's internal title band. Returns path to the crop."""
    import numpy as _np
    from PIL import Image as _Im
    dst = _os.path.join(_DER, name.replace(".png", "_crop.png"))
    if _os.path.exists(dst):
        return dst
    im = _Im.open(_os.path.join(_FIGS, name))
    a = _np.array(im.convert("L"))
    ink = (a < 240).sum(axis=1)
    end = None
    for i, v in enumerate(ink):
        if v > 0:
            j = i
            while j < len(ink) and ink[j] > 0:
                j += 1
            end = j
            break
    im.crop((0, end + 22, im.size[0], im.size[1])).save(dst)
    return dst


def _hero_panel_b():
    """Split cropped Figure 4 at its widest internal whitespace, keep panel b."""
    import numpy as _np
    from PIL import Image as _Im
    dst = _os.path.join(_DER, "Fig4b_hpylori.png")
    if _os.path.exists(dst):
        return dst
    im = _Im.open(_crop_title("Figure4_real_data.png"))
    a = _np.array(im.convert("L"))
    ink = (a < 240).sum(axis=0)
    W = len(ink)
    lo, hi = int(W * 0.35), int(W * 0.65)
    bands, s = [], None
    for i in range(lo, hi):
        if ink[i] == 0 and s is None:
            s = i
        if ink[i] > 0 and s is not None:
            bands.append((s, i, i - s))
            s = None
    bands.sort(key=lambda b: -b[2])
    split = (bands[0][0] + bands[0][1]) // 2
    im.crop((split, 0, W, im.size[1])).save(dst)
    return dst


def _qr(url=REPO_URL):
    dst = _os.path.join(_DER, "repo_qr.png")
    if not _os.path.exists(dst):
        import qrcode
        q = qrcode.QRCode(box_size=20, border=1,
                          error_correction=qrcode.constants.ERROR_CORRECT_M)
        q.add_data(url)
        q.make(fit=True)
        q.make_image(fill_color="#12100E", back_color="white").save(dst)
    return dst


def _asset(name):
    """Resolve a figure name to a path, deriving crops when asked for."""
    if name.endswith("_crop.png"):
        return _crop_title(name.replace("_crop.png", ".png"))
    if name == "Fig4b_hpylori.png":
        return _hero_panel_b()
    if name == "repo_qr.png":
        return _qr()
    return _os.path.join(_FIGS, name)
# ---------------------------------------------------------------------------

W, H = 44.0, 44.0
MARGIN = 1.5
GUTTER = 1.0
NCOL = 3
COLW = (W - 2 * MARGIN - (NCOL - 1) * GUTTER) / NCOL   # 14.33
COLX = [MARGIN + i * (COLW + GUTTER) for i in range(NCOL)]

INK = "#1A1A1A"
MUTED = "#4A4A4A"
ACCENT = "#CC0000"          # BU scarlet
PANEL_BG = "#F4F2EF"
RULE = "#D8D4CE"

TITLE_PT = 76
SUB_PT = 34
HEAD_PT = 48
BODY_PT = 31
NUM_PT = 42
CAP_PT = 23
FOOT_PT = 21

fig = plt.figure(figsize=(W, H), dpi=100)
fig.patch.set_facecolor("white")


def inches(x, y):
    """Convert inches (origin bottom-left) to figure fraction."""
    return x / W, y / H


def text(x, y, s, pt, color=INK, weight="normal", va="top", ha="left",
         wrap_chars=None, leading=1.30, style="normal"):
    """Draw wrapped text at (x, y) in inches. Returns height consumed."""
    lines = []
    if wrap_chars:
        for para in s.split("\n"):
            lines.extend(textwrap.wrap(para, wrap_chars) or [""])
    else:
        lines = s.split("\n")
    line_h = pt * leading / 72.0
    for i, ln in enumerate(lines):
        fx, fy = inches(x, y - i * line_h)
        fig.text(fx, fy, ln, fontsize=pt, color=color, weight=weight,
                 va=va, ha=ha, style=style, family="DejaVu Sans")
    return len(lines) * line_h


def rule(x, y, w, lw=2.5, color=RULE):
    fx, fy = inches(x, y)
    fig.add_artist(Rectangle((fx, fy), w / W, lw / 72.0 / H,
                             facecolor=color, edgecolor="none",
                             transform=fig.transFigure))


def panel_bg(x, y_top, w, h, color=PANEL_BG):
    fx, fy = inches(x - 0.35, y_top - h)
    fig.add_artist(FancyBboxPatch((fx, fy), (w + 0.7) / W, h / H,
                                  boxstyle="round,pad=0,rounding_size=0.004",
                                  facecolor=color, edgecolor="none",
                                  transform=fig.transFigure, zorder=0))


def image(path, x, y_top, width):
    """Place an image with its top-left at (x, y_top). Returns height."""
    im = mpimg.imread(_asset(path))
    h = width * im.shape[0] / im.shape[1]
    fx, fy = inches(x, y_top - h)
    ax = fig.add_axes([fx, fy, width / W, h / H])
    ax.imshow(im)
    ax.axis("off")
    return h


def headline(x, y, n, s, w):
    """Numbered panel headline with accent bar. Returns height consumed."""
    fx, fy = inches(x - 0.35, y - 0.62)
    fig.add_artist(Rectangle((fx, fy), 0.16 / W, 0.62 / H,
                             facecolor=ACCENT, edgecolor="none",
                             transform=fig.transFigure))
    text(x, y, f"{n}   {s}", HEAD_PT, weight="bold", wrap_chars=30)
    n_lines = len(textwrap.wrap(f"{n}   {s}", 30))
    return n_lines * HEAD_PT * 1.22 / 72.0 + 0.22


# ----------------------------------------------------------------- header
text(MARGIN, H - MARGIN,
     "Cross-sectional intratumoral microbial abundance cannot establish causation",
     TITLE_PT, weight="bold", wrap_chars=64)
text(MARGIN, H - MARGIN - 3.05,
     "A calibration study using an established carcinogen as a natural control",
     44, color=MUTED, wrap_chars=88, style="italic")
text(MARGIN, H - MARGIN - 4.35,
     "Julian Borges, MD, MS   |   Department of Computer Science, Boston University Metropolitan College",
     SUB_PT, weight="bold")
text(MARGIN, H - MARGIN - 5.00,
     "ORCID 0009-0001-9929-3135   |   jyborges@bu.edu",
     28, color=MUTED)
rule(MARGIN, H - MARGIN - 5.70, W - 2 * MARGIN, lw=5, color=ACCENT)

TOP = H - MARGIN - 6.45
BOTTOM = MARGIN + 2.4

# --------------------------------------------------------------- column 1
y = TOP
y -= headline(COLX[0], y, "1", "The question the retraction did not answer", COLW)
y -= text(COLX[0], y,
          "The 2024 retraction of a pan-cancer tumor microbiome analysis "
          "established that reported signal can be human read misclassification "
          "and batch effect. The field responded with decontamination pipelines "
          "and curated resources. Those address one question: is a detected "
          "taxon really present. They do not address a second: if it is present, "
          "is it causal, and in which direction does the association run. The "
          "dominant design, cross-sectional tumor versus adjacent normal, detects "
          "a difference in abundance well and interprets it poorly, because it "
          "fixes neither the cause nor the temporal order.",
          BODY_PT, wrap_chars=55) + 0.50

y -= headline(COLX[0], y, "2", "Standard tests pass on data with no biology", COLW)
text(COLX[0], y, "10 of 10", NUM_PT, weight="bold", color=ACCENT)
text(COLX[0] + 3.9, y,
     "seeds in which a cohort built to contain no\nbiological signal passed both standard tests",
     27, color=MUTED)
y -= 1.55
y -= text(COLX[0], y,
          "Comparison against a no-information rate and label permutation both "
          "passed on zero-signal data. Only a confounder baseline and "
          "within-batch cross-validation discriminated. Batch-outcome "
          "confounding alone drove accuracy to 2.5 times chance on data "
          "with no biology, as confounding rose to 0.95.",
          BODY_PT, wrap_chars=55) + 0.25
y -= image("Figure1_audit_validation_crop.png", COLX[0] + 0.65, y, 11.7) + 0.14
y -= text(COLX[0], y,
          "Fig 1  Pass rate by test across ten seeds. T1 and T2 pass in both "
          "cohorts; only T3 and T5a separate them.",
          CAP_PT, color=MUTED, wrap_chars=74) + 0.30
y -= image("Figure2_confounding_sweep_crop.png", COLX[0] + 0.65, y, 11.7) + 0.14
y -= text(COLX[0], y,
          "Fig 2  Accuracy against confounding strength. Within-batch delta stays "
          "near zero for zero-signal data at every level.",
          CAP_PT, color=MUTED, wrap_chars=74)

# --------------------------------------------------------------- column 2
y = TOP
y -= headline(COLX[1], y, "3", "An established carcinogen runs backwards", COLW)
text(COLX[1], y, "8x", NUM_PT, weight="bold", color=ACCENT)
text(COLX[1] + 2.0, y,
     "depletion of H. pylori at the gastric tumor site\n"
     "39 matched pairs, diff -0.99 log units, p 4.5e-4",
     27, color=MUTED)
y -= 1.55
y -= text(COLX[1], y,
          "Helicobacter pylori is an IARC Group 1 carcinogen and the accepted "
          "cause of gastric adenocarcinoma. In a decontaminated reference cohort "
          "it is depleted, not enriched, at the tumor site. The mechanism is "
          "general: gastric carcinogenesis proceeds through atrophy and "
          "intestinal metaplasia, which eliminate the acid-adapted mucosa the "
          "organism requires. It causes the disease and is then displaced by it. "
          "For the one organism whose causal role is known, cross-sectional "
          "abundance gives the wrong direction.",
          BODY_PT, wrap_chars=55) + 0.25
y -= image("Figure4_real_data_crop.png", COLX[1] + 0.55, y, 11.9) + 0.14
y -= text(COLX[1], y,
          "Fig 4  (a) Within-tissue tumor versus adjacent-normal discrimination, "
          "weak in all four evaluable cancers. (b) Paired H. pylori depletion.",
          CAP_PT, color=MUTED, wrap_chars=74) + 0.50

y -= headline(COLX[1], y, "4", "Two ways clean data still misleads", COLW)
h1 = image("Figure3_detection_floor_crop.png", COLX[1] + 0.05, y, 6.2)
h2 = image("Figure5_nesting_crop.png", COLX[1] + 6.8, y, 6.2)
y -= max(h1, h2) + 0.14
text(COLX[1] + 0.1, y,
     "Fig 3  Detection floor is about 0.8 log units. A null must be reported as "
     "no signal above the floor.",
     CAP_PT, color=MUTED, wrap_chars=34)
text(COLX[1] + 7.35, y,
     "Fig 5  Discovery counts vary up to 5.7-fold with the taxonomic redundancy "
     "rule, on identical data.",
     CAP_PT, color=MUTED, wrap_chars=34)
y -= 1.55

m_h = 6.4
panel_bg(COLX[1], y, COLW, m_h, color="#F1F1F1")
y -= 0.48
text(COLX[1], y, "DATA AND METHODS", 27, weight="bold", color=MUTED)
y -= 0.92
text(COLX[1], y,
     "Simulation  Two synthetic cohorts, six cancer types, three batches, 150 "
     "taxa, 60 samples per type. One carries genuine signal, one carries none. "
     "Batch confounded with outcome at 0.75.\n"
     "Audit  T1 no-information rate, T2 label permutation, T3 confounder "
     "baseline, T5a within-batch cross-validation. Random forests.\n"
     "Real data  The Cancer Microbiome Atlas: 611 samples, 14,492 taxa, five "
     "TCGA projects, three centers.\n"
     "MR  MiBioGen, 211 taxa, 18,340 individuals. IVW primary, plus MR-Egger, "
     "weighted median, Cochran Q, BH-FDR.",
     27, wrap_chars=62)


# --------------------------------------------------------------- column 3
y = TOP
y -= headline(COLX[2], y, "5", "Germline anchoring: the causal screen", COLW)
text(COLX[2], y, "100,204", NUM_PT, weight="bold", color=ACCENT)
text(COLX[2] + 4.9, y,
     "case independent replication cohort\n211 gut taxa, two-sample MR",
     27, color=MUTED)
y -= 1.55
y -= text(COLX[2], y,
          "Every FinnGen lead failed replication. The one FDR-significant hit, "
          "phylum Cyanobacteria, reads as pleiotropy: Egger intercept p 0.05, "
          "an instrument at NOS2 (a colorectal carcinogenesis gene, violating "
          "the exclusion restriction), and a direction that flipped on "
          "replication. A Bifidobacterium signal collapsed on removal of "
          "rs182549 at LCT. Zero taxa survive FDR at scale.",
          BODY_PT, wrap_chars=55) + 0.25
y -= image("Figure6_MR_colorectal_crop.png", COLX[2], y, 13.0) + 0.14
y -= text(COLX[2], y,
          "Fig 6  Two-sample MR of 211 gut taxa against colorectal cancer, "
          "FinnGen R12.",
          CAP_PT, color=MUTED, wrap_chars=74) + 0.50

y -= headline(COLX[2], y, "6", "The signal is real but diffuse", COLW)
y -= image("Figure7_cross_cohort.png", COLX[2], y, 13.0) + 0.14
y -= text(COLX[2], y,
          "Fig 7  Cross-cohort concordance, 210 taxa, sign recovered. Agreement "
          "61.4 percent (129/210, binomial p 0.0011) against 50 percent under "
          "noise, robust to the tie rule.",
          CAP_PT, color=MUTED, wrap_chars=74) + 0.20
y -= text(COLX[2], y,
          "Alistipes points protective in three of three tests, nominal in two of "
          "three. Colorectal FinnGen is null at p 0.578. No test survives FDR. A "
          "lead requiring species-level instrumentation, not a finding.",
          BODY_PT, wrap_chars=55) + 0.40

# Conclusion box: height computed from the wrapped text, not hardcoded.
CONCLUSION = ("Decontamination is necessary and not sufficient. Cross-sectional "
              "abundance cannot establish causation and, on the one organism whose "
              "causal role is known, gets the direction backwards. Causal inference "
              "here requires germline anchoring or temporal ordering, and both "
              "require an audit first.")
_lines = textwrap.wrap(CONCLUSION, 55)
_lh = BODY_PT * 1.30 / 72.0
box_h = 0.36 + 0.76 + len(_lines) * _lh + 0.26
panel_bg(COLX[2], y, COLW, box_h, color="#EFE7E7")
y -= 0.36
text(COLX[2], y, "CONCLUSION", 25, weight="bold", color=ACCENT)
y -= 0.76
text(COLX[2], y, CONCLUSION, BODY_PT, wrap_chars=55)
y -= len(_lines) * _lh + 0.42

# QR_BLOCK_V1
# Anchored to the running column cursor, and shrunk if the space left above the
# footer is tight, so it can never collide with the conclusion box.
_qr_top = y - 0.18
_avail = _qr_top - (MARGIN + 1.78)
_qr_size = min(1.45, max(0.95, _avail))
_qr_x = COLX[2] + COLW - _qr_size
_qax = fig.add_axes([_qr_x / W, (_qr_top - _qr_size) / H, _qr_size / W, _qr_size / H])
_qax.imshow(mpimg.imread(_asset("repo_qr.png")))
_qax.axis("off")
text(COLX[2], _qr_top - 0.04,
     "Scan for the repository:\ncode, results, figures,\nand the corrections log.", 26,
     color=MUTED, leading=1.24)

fig.savefig("POSTER_RO-2026-008_44x44_preview.png", dpi=72, facecolor="white")
fig.savefig("POSTER_RO-2026-008_44x44.pdf", facecolor="white")
print("wrote POSTER_RO-2026-008_44x44.png and .pdf")
