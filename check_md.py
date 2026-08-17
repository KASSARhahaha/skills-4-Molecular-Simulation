# -*- coding: utf-8 -*-
r"""在线版闸门。两趟：mkdocs 之前量 Markdown，之后量成品 HTML。

    在线版保真 30 章 + index；裸宏 0 处；逐章 图 101/101、编号公式 1076/1076
    在线版成品 未渲染的行间公式 0 块；code/pre 里的数学 0 处；裸反斜杠 0 处

**为什么非两趟不可。** 纸质版有 37 道闸门，在线版一道也没有——而它同样是
已发布的交付物。第一版闸门只量 .md，当场就漏掉了本轮最大的一处：全书 108 块
行间公式在 .md 里写法**完全正常**（`$$` 独占一行、内容无误、`\tag` 齐全），
只是缩进比该有的少了一格，arithmatex 的块处理器不认，`$$…$$` 原样落进
Markdown 的行内处理器，`_{\mathrm{sum}}` 的下划线被当成强调，页面上印出

    \mathbf{v}<em _mathrm_sum="\mathrm{sum">{\mathrm{sum}}(1, j) = \mathbf{v}</em>(t)

——既是裸 LaTeX、又被撕坏，裹着 76 条编号公式。`.md` 干净、`mkdocs
build --strict` 通过、内容量两侧相符，三样全绿，东西还是坏的。**只有到成品
HTML 上才量得出来。**

── 第一趟：Markdown（mkdocs 之前）────────────────────────────────────

判据一：生成物里不得有裸的反斜杠命令。先把代码块、`$$…$$`、行内 `$…$`、
        HTML 注释剥掉，剩下的 `\字母` 一律算漏网。**范围只框在 tex2md 的
        生成物上**（`CHAPTERS` 的 30 章 + index.md）——`docs/skills/` 下是
        手写文档，不是这条产线的产物，混进来只会制造噪声。

判据二：内容量两侧相符。逐章比 `\includegraphics` 与 `![`、比 `\tag{`。
        这两样是内容承载单元，转换器吃掉一张图或一条编号公式，页面上就是
        实打实的缺失，而字面比对做不了（一边 LaTeX 一边 Markdown）。

判据三：数学定界符配平。`$` 落单会让整段以源码示人。

判据四：行内代码里不得有 `$`。**MathJax 默认跳过 `code`/`pre`**，
        `` `for $1 \le x \le y$ do` `` 会连着美元号原样印出来，而纸质版上
        那里是 `for 1 ≤ x ≤ y do`。

── 第二趟：成品 HTML（mkdocs 之后，`--html`）─────────────────────────

判据五：不得有字面 `$$`。arithmatex 接走了就变成 `<div class="arithmatex">`，
        HTML 里还留着 `$$` 就是没接走。

判据六：`code`/`pre` 里不得有 `$…$`。同判据四，从成品那侧再量一遍。

判据七：剥掉全部数学之后，不得有裸的反斜杠命令。

判据八：成品上的图张数、编号公式条数，与 Markdown 侧相符——渲染这一段
        不许把内容弄丢。

**章节清单跟转换器同一份真相源**（`from tex2md import CHAPTERS`），不自己
glob——闸门用 glob 猜「该查哪些文件」迟早会多读或少读，这个坑记过一次。
"""
import io
import os
import re
import sys

SITE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SITE)
from tex2md import CHAPTERS, TEX_DIR, DOCS_DIR              # noqa: E402

TEXCOMMENT = re.compile(r"(?<!\\)%.*$", re.M)
FENCE = re.compile(r"^```.*?^```", re.S | re.M)
DISPMATH = re.compile(r"\$\$.*?\$\$", re.S)
INLMATH = re.compile(r"\$[^$\n]*\$")
HTMLCOM = re.compile(r"<!--.*?-->", re.S)
CMD = re.compile(r"\\([A-Za-z@]+)")
MDCODE = re.compile(r"(?<!`)`([^`\n]+)`(?!`)")

SCRIPT = re.compile(r"<script.*?</script>|<style.*?</style>", re.S)
HMATH = re.compile(r"\\\(.*?\\\)|\\\[.*?\\\]|\$\$.*?\$\$|\$[^$\n]*\$", re.S)
HCODE = re.compile(r"<(code|pre)[^>]*>(.*?)</\1>", re.S)
HDOLLAR = re.compile(r"\$[^$]{1,80}\$", re.S)
HIMG = re.compile(r"<img\s[^>]*src=\"[^\"]*images/")   # \s 不能省：<img[^>]+
# 会把 <imgX …> 也算进来，点火时正是这条让「渲染把图弄丢了」那一支没红
TAGCMD = re.compile(r"\\tag\{")

MIN_FIG, MIN_TAG = 90, 1000     # 两侧都少于这个数，就是清单或路径出了问题

bad = []


def die(m):
    bad.append(m)


def strip(md):
    """剥掉代码块与数学，剩下的反斜杠命令才是漏网的。"""
    t = FENCE.sub(" ", md)
    t = HTMLCOM.sub(" ", t)
    t = DISPMATH.sub(" ", t)
    return INLMATH.sub(" ", t)


def unbalanced(md):
    """行内 `$` 落单的行号。$$ 块先剥掉。"""
    t = DISPMATH.sub(" ", FENCE.sub(" ", md))
    return [i for i, ln in enumerate(t.split("\n"), 1) if ln.count("$") % 2]


def md_files():
    f = [(s, t, os.path.join(DOCS_DIR, "chapters", s + ".md"))
         for s, t, _ in CHAPTERS]
    f.append(("index", None, os.path.join(DOCS_DIR, "index.md")))
    return f


def check_md():
    fig = tag = mfig = mtag = leak = 0
    for slug, tex_file, mp in md_files():
        if not os.path.exists(mp) or not os.path.getsize(mp):
            die("%s.md 不存在或是空的——这一章没生成出来" % slug)
            continue
        md = io.open(mp, encoding="utf-8").read()
        clean = strip(md)
        for m in CMD.finditer(clean):
            leak += 1
            die("%s.md 留下裸宏 \\%s——Markdown 不认反斜杠命令，会当普通文字"
                "印在页面上：…%s…"
                % (slug, m.group(1),
                   clean[max(0, m.start() - 30):m.end() + 24].replace("\n", " ")))
        for i in unbalanced(md):
            die("%s.md 第 %d 行的 $ 落单——整段会以源码示人" % (slug, i))
        for m in MDCODE.finditer(md):
            if "$" in m.group(1):
                die("%s.md 行内代码里有数学 `%s`——MathJax 跳过 code/pre，"
                    "美元号会原样印出来" % (slug, m.group(1)[:46]))
        if tex_file is None:
            continue
        tp = os.path.join(TEX_DIR, tex_file)
        if not os.path.exists(tp):
            die("%s 找不到书稿 %s——闸门量不了内容量" % (slug, tex_file))
            continue
        t = TEXCOMMENT.sub("", io.open(tp, encoding="utf-8").read())
        a, b = len(re.findall(r"\\includegraphics", t)), len(re.findall(r"!\[", md))
        c, e = len(TAGCMD.findall(t)), len(TAGCMD.findall(md))
        fig += a
        mfig += b
        tag += c
        mtag += e
        if a != b:
            die("%s 图数不符：书稿 %d 张，网页 %d 张" % (slug, a, b))
        if c != e:
            die("%s 编号公式不符：书稿 %d 条，网页 %d 条" % (slug, c, e))
    if fig < MIN_FIG or tag < MIN_TAG:
        die("书稿那侧只数出图 %d 张、编号公式 %d 条（该有 %d+/%d+）——"
            "清单或路径出了问题，本项不作数" % (fig, tag, MIN_FIG, MIN_TAG))
    print("在线版保真 %d 章 + index（只量 tex2md 的生成物，docs/skills/ 下是手写"
          "文档不在此列）；裸宏 %d 处；逐章 图 %d/%d、编号公式 %d/%d"
          % (len(CHAPTERS), leak, fig, mfig, tag, mtag))
    return mfig, mtag


def check_html(want_fig, want_tag):
    root = os.path.join(SITE, "site", "chapters")
    if not os.path.isdir(root):
        die("找不到 %s——还没跑过 mkdocs build，成品这侧量不了，本项不作数" % root)
        print("在线版成品 未查（没有 site/）")
        return
    pages = raw = incode = leak = fig = tag = 0
    for d in sorted(os.listdir(root)):
        p = os.path.join(root, d, "index.html")
        if not os.path.exists(p):
            continue
        pages += 1
        t = SCRIPT.sub(" ", io.open(p, encoding="utf-8").read())
        fig += len(HIMG.findall(t))
        tag += len(TAGCMD.findall(t))
        n = t.count("$$")
        if n:
            raw += n // 2
            die("%s 有 %d 块行间公式没渲染——HTML 里还留着字面 $$，"
                "arithmatex 没接走，页面上是一框源码" % (d, n // 2))
        for m in HCODE.finditer(t):
            for k in HDOLLAR.finditer(m.group(2)):
                incode += 1
                die("%s 的 code/pre 里有数学 %s——MathJax 跳过这两种元素，"
                    "会原样印出来" % (d, k.group(0)[:46].replace("\n", " ")))
        for m in CMD.finditer(HMATH.sub(" ", t)):
            leak += 1
            die("%s 页面上留着裸宏 \\%s" % (d, m.group(1)))
    if pages < len(CHAPTERS):
        die("成品只有 %d 页，章节清单有 %d 章——site/ 是旧的或没建全，本项不作数"
            % (pages, len(CHAPTERS)))
    if fig != want_fig:
        die("成品图 %d 张，Markdown 侧 %d 张——渲染这一段把图弄丢了"
            % (fig, want_fig))
    if tag != want_tag:
        die("成品编号公式 %d 条，Markdown 侧 %d 条——渲染这一段把公式弄丢了"
            % (tag, want_tag))
    print("在线版成品 %d 页；未渲染的行间公式 %d 块；code/pre 里的数学 %d 处；"
          "裸反斜杠 %d 处；图 %d/%d、编号公式 %d/%d"
          % (pages, raw, incode, leak, fig, want_fig, tag, want_tag))


def main():
    want = check_md()
    if "--html" in sys.argv:
        check_html(*want)
    for m in bad:
        print("  ✗ %s" % m)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
