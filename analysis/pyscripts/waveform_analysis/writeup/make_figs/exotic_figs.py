#!/usr/bin/env python3
"""Crop page 2 (waveform panels) of the two exotic-event displays for the writeup.

Reads   ../../pileup_analysis/exotic_events/event_display_run24590_ev*_W*.pdf  (read-only)
Writes  ../figs/exotic_w121_baseline_burst.pdf, ../figs/exotic_w493_wide_pulse.pdf
Needs pypdf.  Boxes are (left, bottom, right, top) in points on the 1080 x 828 pt page.
"""
import os
from pypdf import PdfReader, PdfWriter, Transformation, PageObject

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "..", "..", "pileup_analysis", "exotic_events")
OUT = os.path.join(HERE, "..", "figs")
JOBS = [
    ("event_display_run24590_ev4084847_W121.pdf", "exotic_w121_baseline_burst.pdf", (345, 268, 975, 540)),
    ("event_display_run24590_ev2069998_W493.pdf", "exotic_w493_wide_pulse.pdf", (345, 265, 1080, 765)),
]
for src, out, (l, b, r, t) in JOBS:
    page = PdfReader(os.path.join(SRC, src)).pages[1]
    new = PageObject.create_blank_page(width=r - l, height=t - b)
    new.merge_transformed_page(page, Transformation().translate(-l, -b), expand=False)
    w = PdfWriter()
    w.add_page(new)
    w.write(os.path.join(OUT, out))
    print("wrote", out)
