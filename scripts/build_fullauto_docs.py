#!/usr/bin/env python3
"""Render the maintained fullauto Markdown guides to standalone HTML and TeX."""
from pathlib import Path
import subprocess
import tempfile
import re

ROOT = Path(__file__).resolve().parents[1]
STYLE = '''<style>
:root{color-scheme:light;--ink:#163047;--accent:#086f82;--line:#d6e1e8}
body{font:17px/1.65 system-ui,sans-serif;color:var(--ink);max-width:1100px;margin:36px auto;padding:0 28px;background:#fbfdff}
h1,h2,h3{line-height:1.25;color:#13394b;scroll-margin-top:20px}h1{font-size:2.3em}h2{margin-top:2.2em;padding-top:.6em;border-top:1px solid var(--line)}
a{color:var(--accent)}pre{background:#edf3f7;border:1px solid var(--line);border-radius:8px;padding:18px;overflow:auto;line-height:1.45}code{font-size:.88em;overflow-wrap:anywhere}
table{border-collapse:collapse;width:100%;display:block;overflow:auto;margin:22px 0}th,td{padding:11px 14px;border:1px solid var(--line);text-align:left;vertical-align:top}th{background:#e8f1f5}
nav{background:#f0f6f8;padding:18px 24px;border-radius:10px}nav ul{margin:.3em 0}math[display=block]{overflow:auto;padding:12px 0}p{margin:.8em 0}
@media print{body{font-size:10pt;max-width:none;margin:0}pre,table{break-inside:avoid}nav{display:none}}
</style>'''

def main():
    with tempfile.TemporaryDirectory(prefix='mbs-manual-') as d:
        header=Path(d)/'style.html';header.write_text(STYLE)
        for source,target,title in [
            ('FULLAUTO_RU.md','FULLAUTO_ADAPTIVE_RU.html','MBS-fast: подробный fullauto'),
            ('FULLAUTO_QUICK_RU.md','MINI_FULLAUTO_RU.html','MBS-fast: краткий fullauto')]:
            subprocess.run(['pandoc',str(ROOT/'docs'/source),'--standalone','--toc','--toc-depth=2',
                '--mathml','--no-highlight','--include-in-header',str(header),
                '--metadata','title='+title,'--metadata','lang=ru','-o',str(ROOT/'docs'/target)],check=True)
    for source,target in [('FULLAUTO.md','FULLAUTO.inc.tex'),('FULLAUTO_RU.md','FULLAUTO_RU.inc.tex')]:
        subprocess.run(['pandoc',str(ROOT/'docs'/source),'-t','latex','--no-highlight',
            '-o',str(ROOT/'docs'/target)],check=True)
        p=ROOT/'docs'/target
        def wrap(match):
            value=match.group(1)
            return '\\nolinkurl{'+value.replace('\\_', '_')+'}' if ('\\ ' not in value and ('_' in value or '/' in value or value.startswith('--'))) else match.group(0)
        text=re.sub(r'\\texttt\{([^{}]*)\}',wrap,p.read_text())
        text=text.replace('\\begin{verbatim}', '\\noindent\\begin{minipage}{\\linewidth}\n\\begin{verbatim}').replace('\\end{verbatim}', '\\end{verbatim}\n\\end{minipage}')
        p.write_text(text)

if __name__=='__main__':
    main()
