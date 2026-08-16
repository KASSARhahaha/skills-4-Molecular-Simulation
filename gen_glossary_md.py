#!/usr/bin/env python3
"""生成在线版的英汉术语对照表 docs/glossary.md。

数据源和纸质版同一份：译稿仓 tools/terms.py + tools/names.py，
筛选规则直接复用 tools/glossary.py 里的 is_term/norm/cap，
免得两个版本的词条对不上。

用法：python3 gen_glossary_md.py
"""
import os
import sys

TEX_DIR = ("/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/"
           "Chinese-Understanding Molecular Simulation")
DOCS_DIR = "/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/molsim-online/docs"

sys.path.insert(0, os.path.join(TEX_DIR, "tools"))
import glossary as G          # noqa: E402  复用同一套筛选与大小写还原
from names import KEEP, PERSON  # noqa: E402

items = G.collect()

L = ["# 英汉术语对照表", "",
     "词条与主题索引同源，译名均取自正文实际用词，可与英文原著逐条对照查阅。",
     "本页由译稿仓的 `tools/terms.py` 与 `tools/names.py` 生成。", "",
     "## 人名类术语的写法", "",
     "人名充当普通科技术语的构词成分时用规范汉译；人名充当算法、模型、力场或软件的"
     "专名时保留拉丁原名。同一人名全书只用一种写法。参考文献与人名索引一律用拉丁"
     "原名，因为那里索引的是人，不是术语。", "",
     "### 用汉译的人名", "", "| 原名 | 汉译 |", "|---|---|"]
for en in sorted(PERSON, key=str.lower):
    L.append(f"| {en} | {PERSON[en]} |")

L += ["", "### 保留拉丁原名的专名", "", "| 原名 | 说明 |", "|---|---|"]
for en in sorted(KEEP, key=str.lower):
    L.append("| %s | %s |" % (en.replace("\\'{e}", "é"), KEEP[en]))

L += ["", "## 英汉术语对照", ""]
cur = None
for en, cn in items:
    g = G.group_letter(en)
    if g != cur:
        cur = g
        L += ["", f"### {g}", "", "| 英文 | 中文 |", "|---|---|"]
    L.append(f"| {G.cap(en)} | {cn} |")

out = os.path.join(DOCS_DIR, "glossary.md")
open(out, "w", encoding="utf-8").write("\n".join(L) + "\n")
print(f"Wrote glossary.md（{len(items)} 条，人名 {len(PERSON)}+{len(KEEP)}）")
