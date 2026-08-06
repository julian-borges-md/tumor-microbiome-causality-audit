#!/usr/bin/env python3
"""
RO-2026-008 | Social card builder v1.1

Two adaptations of the poster, not reductions of it. A 44 x 44 poster shrunk
to a feed image is unreadable, so each card carries one hero figure, three
anchor numbers, and a single claim.

  SOCIAL_linkedin_1920x1005.png   1.91:1 landscape, LinkedIn feed
  SOCIAL_instagram_1080x1080.png  1:1 square, Instagram feed

Layout is a top-down cursor and the hero is sized from whatever space
remains, so no element can overrun the canvas.

    python3 social_cards.py
"""
import textwrap
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle

INK = "#12100E"
MUTED = "#5A544D"
ACCENT = "#CC0000"
BG = "#FBFAF8"
HERO = "Fig4b_hpylori.png"

EYEBROW = "TUMOR MICROBIOME   |   CAUSAL INFERENCE"
TITLE = "A known carcinogen\nruns backwards"
CLAIM = ("Cross-sectional tumor microbiome abundance recovers the wrong "
         "causal direction.")
NUMBERS = [
    ("10 / 10", "seeds where standard tests\npassed on zero-signal data"),
    ("8x", "depletion of H. pylori at\nthe gastric tumor site"),
    ("100,204", "case independent\nreplication cohort"),
]
EVENT = ("BU Health Data Science & AI Showcase   |   15 September 2026   |   "
         "Hiebert Lounge")
BYLINE = ("Julian Borges, MD, MS   |   Department of Computer Science, "
          "Boston University")


def build(path, px_w, px_h, dpi=150, show_claim=True, title_pt=34):
    W, H = px_w / dpi, px_h / dpi
    portrait_or_square = px_h >= px_w
    fig = plt.figure(figsize=(W, H), dpi=dpi)
    fig.patch.set_facecolor(BG)

    def t(x, y, s, pt, color=INK, weight="normal", wrap=None, leading=1.22):
        lines = []
        for para in s.split("\n"):
            lines.extend(textwrap.wrap(para, wrap) if wrap else [para])
        lh = pt * leading / 72.0
        for i, ln in enumerate(lines):
            fig.text(x / W, (y - i * lh) / H, ln, fontsize=pt, color=color,
                     weight=weight, ha="left", va="top", family="DejaVu Sans")
        return len(lines) * lh

    def bar(x, y, w, h, color=ACCENT):
        fig.add_artist(Rectangle((x / W, y / H), w / W, h / H,
                                 facecolor=color, edgecolor="none",
                                 transform=fig.transFigure))

    def hero(cx, y_top, avail_h, max_w):
        im = mpimg.imread(HERO)
        aspect = im.shape[1] / im.shape[0]
        h = avail_h
        w = h * aspect
        if w > max_w:
            w, h = max_w, max_w / aspect
        ax = fig.add_axes([(cx - w / 2) / W, (y_top - h) / H, w / W, h / H])
        ax.imshow(im)
        ax.axis("off")
        return h

    M = 0.50
    FOOT_H = 0.78

    if portrait_or_square:
        y = H - M
        bar(M, y - 0.09, 1.45, 0.09)
        y -= 0.26
        y -= t(M, y, EYEBROW, 11, color=ACCENT, weight="bold") + 0.14
        y -= t(M, y, TITLE, title_pt, weight="bold", leading=1.12) + 0.16
        if show_claim:
            y -= t(M, y, CLAIM, 14, color=MUTED, wrap=58) + 0.28
        else:
            y -= 0.14

        stats_h = 0.92
        floor = M + FOOT_H + 0.22 + stats_h + 0.26
        hero(W / 2, y, y - floor, W - 2 * M - 0.4)

        ys = M + FOOT_H + 0.22 + stats_h
        colw = (W - 2 * M) / 3
        for i, (big, small) in enumerate(NUMBERS):
            x = M + i * colw
            t(x, ys, big, 26, weight="bold", color=ACCENT)
            t(x, ys - 0.42, small, 10, color=MUTED, leading=1.25, wrap=26)
    else:
        LEFT_R = W * 0.47
        y = H - M
        bar(M, y - 0.09, 1.45, 0.09)
        y -= 0.26
        y -= t(M, y, EYEBROW, 11, color=ACCENT, weight="bold") + 0.14
        y -= t(M, y, TITLE, 34, weight="bold", leading=1.12) + 0.16
        y -= t(M, y, CLAIM, 14, color=MUTED, wrap=42) + 0.40

        for big, small in NUMBERS:
            t(M, y, big, 23, weight="bold", color=ACCENT)
            t(M + 1.75, y + 0.03, small, 10, color=MUTED, leading=1.25, wrap=34)
            y -= 0.64

        hx0 = LEFT_R + 0.25
        hero((hx0 + (W - M)) / 2, H - M - 0.05,
             H - 2 * M - FOOT_H - 0.05, (W - M) - hx0)

    bar(M, M + FOOT_H - 0.18, W - 2 * M, 0.016, color="#D8D2CA")
    t(M, M + FOOT_H - 0.30, EVENT, 10, color=MUTED)
    t(M, M + FOOT_H - 0.56, BYLINE, 10, color=MUTED)

    fig.savefig(path, dpi=dpi, facecolor=BG)
    plt.close(fig)
    print("wrote", path)


build("SOCIAL_linkedin_1920x1005.png", 1920, 1005)
build("SOCIAL_instagram_1080x1350.png", 1080, 1350)
build("SOCIAL_instagram_1080x1080.png", 1080, 1080,
      show_claim=False, title_pt=30)
