"""Render report.md -> report.html -> ../report.pdf (WeasyPrint)."""
from pathlib import Path

import markdown
import weasyprint

HERE = Path(__file__).resolve().parent
md = (HERE / "report.md").read_text()

# Python-Markdown needs a blank line before a list that follows a paragraph.
import re
LIST = re.compile(r"^\s*(?:[-*]|\d+\.)\s")
lines, prev = [], ""
for line in md.splitlines():
    if LIST.match(line) and prev.strip() and not LIST.match(prev) and not prev.startswith((" ", "|", ">")):
        lines.append("")
    lines.append(line)
    prev = line
md = "\n".join(lines)
body = markdown.markdown(md, extensions=["tables", "attr_list", "md_in_html", "sane_lists"])

CSS = """
@page { size: A4; margin: 16mm 14mm 18mm 14mm;
        @bottom-center { content: "Page " counter(page) " of " counter(pages); font-size: 8pt; color: #666; } }
body { font-family: "DejaVu Sans", Arial, sans-serif; font-size: 9pt; line-height: 1.38; color: #111; }
h1 { font-size: 15pt; border-bottom: 2px solid #234; padding-bottom: 3px; margin-top: 10px; color: #123; }
h2 { font-size: 11.5pt; margin: 14px 0 4px; color: #234; }
table { border-collapse: collapse; width: 100%; margin: 6px 0 4px; font-size: 7.6pt; page-break-inside: avoid; }
th, td { border: 1px solid #bbb; padding: 2px 4px; vertical-align: top; }
th { background: #e9eef4; }
tr:nth-child(even) td { background: #fafafa; }
img { max-width: 100%; max-height: 120mm; display: block; margin: 6px auto 0; }
p em:first-child { color: #333; }
blockquote { border-left: 4px solid #d08a00; background: #fff6e5; margin: 6px 0; padding: 4px 8px; }
code { font-size: 8pt; background: #f2f2f2; padding: 0 2px; }
.pagebreak { page-break-after: always; }
"""
html = f"<!doctype html><html><head><meta charset='utf-8'><title>Assignment 3 report</title><style>{CSS}</style></head><body>{body}</body></html>"
(HERE / "report.html").write_text(html)
weasyprint.HTML(string=html, base_url=str(HERE)).write_pdf(HERE.parent / "report.pdf")
print("wrote", HERE.parent / "report.pdf")
