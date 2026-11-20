#!/usr/bin/env python3
"""Minimal Markdown -> HTML for the report (headings, paragraphs, nested lists, tables, images, inline styles)."""
import html, re, sys

src, dst, title = sys.argv[1], sys.argv[2], sys.argv[3]
lines = open(src, encoding="utf-8").read().split("\n")


def inline(t):
    t = html.escape(t, quote=False)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<![\w*])\*(?!\s)(.+?)(?<!\s)\*(?![\w*])", r"<em>\1</em>", t)
    return t


out, i = [], 0
while i < len(lines):
    ln = lines[i]
    if not ln.strip():
        i += 1; continue
    m = re.match(r"^(#{1,6})\s+(.*)", ln)
    if m:
        n = len(m.group(1)); out.append(f"<h{n}>{inline(m.group(2))}</h{n}>"); i += 1; continue
    m = re.match(r"^!\[(.*?)\]\((.*?)\)\s*$", ln)
    if m:
        cap = f"<figcaption>{inline(m.group(1))}</figcaption>" if m.group(1) else ""
        out.append(f'<figure><img src="{m.group(2)}">{cap}</figure>'); i += 1; continue
    if ln.lstrip().startswith("|"):
        rows = []
        while i < len(lines) and lines[i].lstrip().startswith("|"):
            rows.append([c.strip() for c in lines[i].strip().strip("|").split("|")]); i += 1
        align = ["right" if a.endswith(":") and not a.startswith(":") else "left" for a in rows[1]]
        h = "".join(f"<th>{inline(c)}</th>" for c in rows[0])
        body = "".join("<tr>" + "".join(f'<td style="text-align:{align[k] if k < len(align) else "left"}">{inline(c)}</td>'
                                        for k, c in enumerate(r)) + "</tr>" for r in rows[2:])
        out.append(f"<table><thead><tr>{h}</tr></thead><tbody>{body}</tbody></table>"); continue
    if re.match(r"^\s*(-|\d+\.)\s+", ln):
        stack = []  # (indent, tag)
        while i < len(lines) and (re.match(r"^\s*(-|\d+\.)\s+", lines[i]) or
                                  (lines[i].startswith("  ") and lines[i].strip())):
            l2 = lines[i]; m = re.match(r"^(\s*)(-|\d+\.)\s+(.*)", l2)
            if not m:  # continuation line
                out[-1] = out[-1][:-5] + " " + inline(l2.strip()) + "</li>"; i += 1; continue
            ind, tag = len(m.group(1)), ("ul" if m.group(2) == "-" else "ol")
            while stack and ind < stack[-1][0]:
                out.append(f"</{stack.pop()[1]}>")
                if stack: out.append("</li>")
            if not stack or ind > stack[-1][0]:
                if stack and out[-1].endswith("</li>"): out[-1] = out[-1][:-5]
                stack.append((ind, tag)); out.append(f"<{tag}>")
            out.append(f"<li>{inline(m.group(3))}</li>"); i += 1
        while stack:
            out.append(f"</{stack.pop()[1]}>")
            if stack: out.append("</li>")
        continue
    para = []
    while i < len(lines) and lines[i].strip() and not re.match(r"^(#|\||!\[|\s*-\s|\s*\d+\.\s)", lines[i]):
        para.append(lines[i].strip()); i += 1
    out.append(f"<p>{inline(' '.join(para))}</p>")

css = """body{font-family:'DejaVu Sans',Arial,sans-serif;font-size:10pt;line-height:1.45;color:#0b0b0b;max-width:185mm;margin:auto}
h1{font-size:17pt}h2{font-size:13.5pt;border-bottom:1px solid #d0cfca;padding-bottom:3px;margin-top:22px}h3{font-size:11.5pt;margin-top:18px}
table{border-collapse:collapse;margin:8px 0;font-size:8.5pt;page-break-inside:avoid}th,td{border:1px solid #d0cfca;padding:3px 6px;vertical-align:top}
th{background:#f2f1ed}figure{margin:10px 0;page-break-inside:avoid}img{max-width:100%;max-height:228mm;width:auto;height:auto;display:block;margin:auto}
figcaption{font-size:9pt;color:#52514e;font-style:italic;margin-top:4px}code{font-size:8.5pt;background:#f2f1ed;padding:0 2px}
li{margin:2px 0}@page{size:A4;margin:16mm 14mm}"""
open(dst, "w", encoding="utf-8").write(f'<!doctype html><html><head><meta charset="utf-8"><title>{title}</title>'
                                       f"<style>{css}</style></head><body>{''.join(out)}</body></html>")
