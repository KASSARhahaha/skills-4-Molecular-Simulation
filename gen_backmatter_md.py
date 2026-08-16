#!/usr/bin/env python3
"""生成在线版的缩略语表 docs/acronyms.md 与符号表 docs/symbols.md。

数据源和纸质版同一份：译稿仓 tools/backmatter.py 里的 ACRONYMS / SYMBOLS，
免得两个版本对不上。符号栏的 LaTeX 直接交给页面上的 MathJax 渲染。

用法：python3 gen_backmatter_md.py
"""
import os
import sys

TEX_DIR = ("/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/"
           "Chinese-Understanding Molecular Simulation")
DOCS_DIR = "/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/molsim-online/docs"

sys.path.insert(0, os.path.join(TEX_DIR, "tools"))
from backmatter import ACRONYMS, SYMBOLS, gloss   # noqa: E402


def write(path, title, intro, header, rows):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("# %s\n\n%s\n\n%s\n" % (title, intro, header))
        fh.write("".join(rows))
    print("Wrote %s（%d 条）" % (os.path.basename(path), len(rows)))


write(os.path.join(DOCS_DIR, "acronyms.md"), "缩略语表",
      "按原书 Acronyms 逐条译出，收全书缩略语 %d 条，顺序与原书一致。"
      "本页由译稿仓的 `tools/backmatter.py` 生成。" % len(ACRONYMS),
      "| 缩略语 | 含义 |\n|---|---|",
      ["| **%s** | %s |\n" % (a, gloss(en, zh)) for a, en, zh in ACRONYMS])

write(os.path.join(DOCS_DIR, "symbols.md"), "符号表",
      "按原书 Glossary 逐条译出，收全书主要符号 %d 条，顺序与原书一致。"
      "符号的字体照原书还原：粗体为矢量或矩阵，花体为算符或泛函。"
      "同一字母出现多次的，表示原书在不同场合用它表示不同的量，以释义区分。"
      "本页由译稿仓的 `tools/backmatter.py` 生成。" % len(SYMBOLS),
      "| 符号 | 含义 |\n|---|---|",
      ["| $%s$ | %s |\n" % (tex, gloss(en, zh)) for tex, en, zh in SYMBOLS])
