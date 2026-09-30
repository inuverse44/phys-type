"""Digitize the equilateral f_NL curve of Fig. 7 of Dias, Frazer, Mulryne and
Seery, arXiv:1609.00379v2 (JCAP 12 (2016) 033), axion-quartic model.

Usage:
    python3 -m venv venv && venv/bin/pip install pypdf
    curl -L https://arxiv.org/pdf/1609.00379v2 -o dfms.pdf
    venv/bin/python docs/derivations/digitize_dfms_figure7.py dfms.pdf

The figure is a vector graphic, so no image processing is needed. The bottom
right panel of Fig. 7 is the Form XObject /Im19 on page 42 (index 41). Its
content stream contains the axis ticks and the curve as one polyline. The
axis mapping is read from the tick marks in the same stream:

- the vertical grid lines for N = 10 and 20 e-folds lie at x = 81.74286 and
  133.2286, so N = 10 + (x - 81.74286) / 5.148571;
- the horizontal ticks for f_NL = -50, ..., 20 lie at y = 46.26803 + 33.10455 n,
  and f_NL = 0 at y = 211.7909, so f_NL = (y - 211.7909) / 3.310455.

The printed points are the plotted samples themselves (no interpolation). The
repository test `axion_quartic_bispectrum_matches_dfms_and_numerical_delta_n`
uses nine of them after N = 59.
"""
import re
import sys

from pypdf import PdfReader


def main(path):
    page = PdfReader(path).pages[41]
    form = page["/Resources"]["/XObject"]["/Im19"].get_object()
    stream = form.get_data().decode("latin1")
    start = stream.index("SC 82.07237 203.4762 m")
    end = stream.index(" S", start)
    points = re.findall(r"([-\d.]+) ([-\d.]+) [ml]", stream[start + 3:end])
    for x, y in points:
        e_folds = 10.0 + (float(x) - 81.74286) / 5.148571
        f_nl = (float(y) - 211.7909) / 3.310455
        print(f"{e_folds:.4f} {f_nl:.4f}")


if __name__ == "__main__":
    main(sys.argv[1])
