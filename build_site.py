"""Render the course report into the static Project 4 viewer page."""

from html import escape
from pathlib import Path
import re

import markdown2

ROOT = Path(__file__).resolve().parent


def render_report(source):
    heading, separator, report = source.partition("\n")
    if not separator or not heading.startswith("# "):
        raise ValueError("index.md must begin with a level-one heading")
    math = []

    def protect(match):
        display = match.group(1) is not None
        tex = escape(match.group(1) if display else match.group(2), quote=False)
        delimiters = (r"\[", r"\]") if display else (r"\(", r"\)")
        math.append(delimiters[0] + tex.strip() + delimiters[1])
        return f"PROJECTFOURMATH{len(math) - 1:04d}TOKEN"

    report = re.sub(r"\$\$(.+?)\$\$|(?<!\$)\$([^$\n]+)\$(?!\$)", protect, report, flags=re.S)
    body = markdown2.markdown(report, extras=["tables", "fenced-code-blocks"])
    for index, equation in enumerate(math):
        body = body.replace(f"PROJECTFOURMATH{index:04d}TOKEN", equation)
    return escape(heading[2:]), body


def build():
    title, body = render_report((ROOT / "index.md").read_text(encoding="utf-8"))
    template = (ROOT / "site_template.html").read_text(encoding="utf-8")
    page = template.replace("<!-- TITLE -->", title).replace("<!-- REPORT -->", body)
    (ROOT / "index.html").write_text(page, encoding="utf-8", newline="\n")


if __name__ == "__main__":
    build()
