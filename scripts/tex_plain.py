# /// script
# requires-python = ">=3.11"
# ///
"""Plain text of a TeX fragment, for the Japanese prose checks of yomiyasu.

yomiyasu's lint and diff expect Markdown or plain text. This drops
displayed equations, tables and figures, replaces inline math with ``X``
and references with ``R``, keeps section titles as Markdown headings, and
joins the one-sentence-per-line source into paragraphs. Run it on the
original and the rewritten fragment, then give both to
``yomiyasu_check.py run --source ... --final ...``:

    uv run scripts/tex_plain.py original.tex > original.md
"""

import re
import sys

DISPLAYED = r"equation|align|gather|center|figure|tabular"


def plain(source: str) -> str:
    s = re.sub(r"(?m)^%.*\n", "", source)
    s = re.sub(rf"\\begin\{{({DISPLAYED})\*?\}}.*?\\end\{{\1\*?\}}", "", s, flags=re.S)
    s = re.sub(r"\\\[.*?\\\]", "", s, flags=re.S)
    s = re.sub(r"\\texorpdfstring\{([^{}]*(?:\{[^{}]*\}[^{}]*)*)\}\{[^{}]*\}", r"\1", s)
    s = re.sub(
        r"\\((?:sub)*)section\{(.*?)\}\n",
        lambda m: "#" * (2 + len(m.group(1)) // 3) + " " + m.group(2) + "\n\n",
        s,
    )
    s = re.sub(r"\$[^$]*\$", "X", s)
    s = re.sub(r"\\footnote\{%?", "（", s)
    s = re.sub(r"\\(label|ref|eqref|cite)\{[^}]*\}", "R", s)
    s = re.sub(r"\\href\{[^}]*\}\{([^}]*)\}", r"\1", s)
    s = re.sub(r"\\(texttt|code|repo|emph)\{([^}]*)\}", r"\2", s)
    s = re.sub(r"\\item\s*", "- ", s)
    s = re.sub(r"\\(begin|end)\{[^}]*\}", "", s)
    s = s.replace("\\S", "§").replace("~", "").replace("\\%", "%")
    s = s.replace("\\_", "_").replace("\\#", "#").replace("--", "–")
    s = re.sub(r"\\[a-zA-Z]+", "", s)
    s = s.replace("{", "").replace("}", "")
    paragraphs = (p.replace("\n", "").strip() for p in re.split(r"\n\s*\n", s))
    return "\n\n".join(p for p in paragraphs if p)


if __name__ == "__main__":
    with open(sys.argv[1], encoding="utf-8") as file:
        print(plain(file.read()))
