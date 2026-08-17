#!/usr/bin/env python3
"""Convert Chinese TeX source files to Markdown for MkDocs Material site.

Handles: chapters, sections, equations (MathJax), figures, tables,
algorithms, citations, cross-references, lists, quotes, footnotes, etc.
"""

import os
import re
import sys

TEX_DIR = "/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/Chinese-Understanding Molecular Simulation"
DOCS_DIR = "/home/yangkai/1-AI-CCUS/8-AI-CCUS/molsim2026/molsim-online/docs"
FIGURES_DIR = os.path.join(TEX_DIR, "figures")
AUX_FILE = os.path.join(TEX_DIR, "main_cn.aux")

# label -> 编号，由 main_cn.aux 的 \newlabel 填充（见 load_labels）
LABELS = {}


def load_labels():
    """从 main_cn.aux 读取 \\newlabel，建立 label -> 编号 的权威映射。

    式号例外：正文用 \\tag{} 手工钉死三级号，aux 里记录的是 LaTeX 自动计数器
    （不可用），而 eq 标签名本身就是式号（如 eq:11.2.15），故直接取标签名。
    """
    if not os.path.exists(AUX_FILE):
        print(f"警告：找不到 {AUX_FILE}，交叉引用将退化为无编号文本")
        return
    with open(AUX_FILE, encoding="utf-8", errors="replace") as f:
        aux = f.read()
    for m in re.finditer(r"\\newlabel\{([^}]+)\}\{\{([^}]*)\}\{", aux):
        name, num = m.group(1), m.group(2).strip()
        if name.startswith("eq:") or not num:
            continue
        LABELS[name] = num
    for m in re.finditer(r"\\newlabel\{(eq:[^}]+)\}", aux):
        LABELS[m.group(1)] = m.group(1)[3:]
    print(f"载入 {len(LABELS)} 条标签映射（来自 main_cn.aux）")

# Mapping from tex slug -> (input_file, title)
CHAPTERS = [
    ("preface-berend", "preface_berend_cn.tex", "译者序"),
    ("preface1", "preface1_cn.tex", "第一版前言"),
    ("preface2", "preface2_cn.tex", "第二版前言"),
    ("preface3", "preface3_cn.tex", "第三版前言"),
    ("ch01-introduction", "ch01_introduction_cn.tex", "第1章 绪论"),
    ("ch02-thermo-statmech", "ch02_thermo_statmech_cn.tex", "第2章 统计力学基础"),
    ("ch03-monte-carlo", "ch03_monte_carlo_cn.tex", "第3章 蒙特卡洛模拟"),
    ("ch04-md", "ch04_md_cn.tex", "第4章 分子动力学模拟"),
    ("ch05-computer-exp", "ch05_computer_exp_cn.tex", "第5章 计算机实验"),
    ("ch06-mc-ensembles", "ch06_mc_ensembles_cn.tex", "第6章 蒙特卡洛系综"),
    ("ch07-md-ensembles", "ch07_md_ensembles_cn.tex", "第7章 分子动力学系综"),
    ("ch08-free-energy", "ch08_free_energy_cn.tex", "第8章 自由能计算"),
    ("ch09-free-energy-solids", "ch09_free_energy_solids_cn.tex", "第9章 固体自由能"),
    ("ch10-free-energy-chains", "ch10_free_energy_chains_cn.tex", "第10章 链状分子自由能"),
    ("ch11-long-ranged", "ch11_long_ranged_cn.tex", "第11章 长程相互作用"),
    ("ch12-cbmc", "ch12_cbmc_cn.tex", "第12章 构型偏倚蒙特卡洛"),
    ("ch13-accel-mc", "ch13_accel_mc_cn.tex", "第13章 加速蒙特卡洛方法"),
    ("ch14-timescale", "ch14_timescale_cn.tex", "第14章 时间尺度加速方法"),
    ("ch15-rare-events", "ch15_rare_events_cn.tex", "第15章 稀有事件"),
    ("ch16-mesoscopic", "ch16_mesoscopic_cn.tex", "第16章 介观模拟方法"),
    ("appA-lagrangian", "appA_lagrangian_cn.tex", "附录A 拉格朗日力学"),
    ("appB-non-hamiltonian", "appB_non_hamiltonian_cn.tex", "附录B 非哈密顿动力学"),
    ("appC-kirkwood-buff", "appC_kirkwood_buff_cn.tex", "附录C Kirkwood-Buff理论"),
    ("appD-non-eq-thermo", "appD_non_eq_thermo_cn.tex", "附录D 非平衡热力学"),
    ("appE-non-eq-work", "appE_non_eq_work_cn.tex", "附录E 非平衡功"),
    ("appF-linear-response", "appF_linear_response_cn.tex", "附录F 线性响应理论"),
    ("appG-committor", "appG_committor_cn.tex", "附录G Committor分析"),
    ("appH-sdpd", "appH_sdpd_cn.tex", "附录H 确定性颗粒动力学"),
    ("appI-saving-cpu", "appI_saving_cpu_cn.tex", "附录I 节省CPU时间"),
    ("appJ-algorithms", "appJ_algorithms_cn.tex", "附录J 算法汇总"),
]


def read_tex(filename):
    path = os.path.join(TEX_DIR, filename)
    with open(path, encoding="utf-8") as f:
        return f.read()


def write_md(slug, content):
    os.makedirs(os.path.join(DOCS_DIR, "chapters"), exist_ok=True)
    path = os.path.join(DOCS_DIR, "chapters", f"{slug}.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write(content)
    return path


def process_labels(text):
    """Remove \label{...} commands.

    索引标记 \sidx{排序键}{显示形式}、\aidx{姓名} 只服务纸质版的索引页，
    在线版靠站内搜索，直接丢掉（\sidx 是两个参数，别只吃掉一个）。
    """
    text = re.sub(r"\\sidx\{[^{}]*\}\{[^{}]*\}", "", text)
    text = re.sub(r"\\aidx\{[^{}]*\}", "", text)
    # \markboth{左}{右} 与 \markright{右} 只管纸质版的书眉，
    # 在线版没有书眉，直接丢掉，否则裸宏会漏进 Markdown。
    text = re.sub(r"\\markboth\{[^{}]*\}\{[^{}]*\}", "", text)
    text = re.sub(r"\\markright\{[^{}]*\}", "", text)
    # \figkey{...} 是题注末尾的「图内英文标注」中译，纸质版排成
    # 「（图内标注：…）」，在线版同样要出，不能留下裸宏。
    text = re.sub(r"\\figkey\{([^{}]*)\}", r"（图内标注：\1）", text)
    return re.sub(r"\\label\{[^}]*\}", "", text)


# 当前章的图号前缀（"3"、"A" …），由 main() 每章设置，供 \captionof{figure} 用
CURRENT_CHAP = ""


def brace_arg(text, start):
    """从 text[start] 处的 '{' 取出配对花括号内的内容，返回 (内容, 右括号后位置)。"""
    if start >= len(text) or text[start] != "{":
        return None, start
    depth, i = 0, start
    while i < len(text):
        if text[i] == "{":
            depth += 1
        elif text[i] == "}":
            depth -= 1
            if depth == 0:
                return text[start + 1:i], i + 1
        i += 1
    return None, start



def web_image(path):
    """把正文里的图名换成网页能显示的文件名。

    .tex 引的是 tools/refigure.py 裁出的矢量 PDF（印刷用），浏览器
    <img> 显示不了；同一个脚本按 1400 px 渲了同名 .png，这里换过去。
    """
    name = path.replace("figures/", "")
    return name[:-4] + ".png" if name.endswith(".pdf") else name

def process_boxed_figures(text):
    """转换 \\begin{center} 里的插图。

    例/例证框是 shaded 环境，浮动体不能嵌套其中，故书中框内插图写作
    center + \\includegraphics + \\setcounter{figure} + \\captionof{figure}{...}。
    """
    def center_replace(m):
        inner = m.group(1)
        img = re.search(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", inner)
        if not img:
            # 框内表格：center + \captionof{table}{...} + tabular
            if "\\begin{tabular}" in inner:
                cap = ""
                cm = re.search(r"\\captionof\{table\}", inner)
                if cm:
                    cap, _ = brace_arg(inner, cm.end())
                    cap = (cap or "").strip()
                body = process_tables(
                    "\\begin{table}" + re.sub(r"\\captionof\{table\}", "", inner)
                    + "\\end{table}"
                )
                return (f"\n*{cap}*\n{body}\n" if cap else f"\n{body}\n")
            return m.group(0)
        img_path = web_image(img.group(1))

        caption = ""
        cm = re.search(r"\\captionof\{figure\}", inner)
        if cm:
            caption, _ = brace_arg(inner, cm.end())
            caption = (caption or "").strip()
        num = re.search(r"\\setcounter\{figure\}\{(\d+)\}", inner)
        if num and CURRENT_CHAP:
            caption = f"图 {CURRENT_CHAP}.{int(num.group(1)) + 1}　{caption}".rstrip("　")

        # alt 只放短图号：题注含括号/引号/数学，塞进 alt 会破坏 Markdown 链接
        alt = f"图 {CURRENT_CHAP}.{int(num.group(1)) + 1}" if (num and CURRENT_CHAP) else ""
        out = f"\n![{alt}](../images/{img_path})\n"
        if caption:
            out += f"\n*{caption}*\n"
        return out

    return re.sub(
        r"\\begin\{center\}(.*?)\\end\{center\}",
        center_replace,
        text,
        flags=re.DOTALL,
    )


def resolve_label(name):
    """label -> 编号字符串；查不到时退回标签名里的数字部分。"""
    if name in LABELS:
        return LABELS[name]
    tail = name.split(":", 1)[-1]
    # alg:mc_verlet 这类别名查不到就返回 None，由调用方决定如何降级
    return tail.replace("_", ".") if re.fullmatch(r"[0-9A-Za-z._]+", tail) else None


def process_refs(text):
    """把 \\ref{...} 解析为原书编号（而非丢弃）。"""
    # \pageref 在网页上没有意义——线上版没有页码。整个括注换成按章节定位的说法，
    # 指向第 1 章的「算法」小节，也就是印本里「第 7 页」那段伪代码约定。
    text = re.sub(
        r"（一般说明参见第~?\\pageref\{[^}]*\}~?页）",
        "（一般说明见第 1 章「算法」）",
        text,
    )
    text = re.sub(r"\\pageref\{[^}]*\}", "", text)   # 兜底：其余 pageref 一律去掉

    def sub_ref(m):
        num = resolve_label(m.group(1))
        return num if num else ""

    # \eqref{eq:...} -> (式号)
    text = re.sub(
        r"\\eqref\{(eq:[^}]*)\}",
        lambda m: f"({resolve_label(m.group(1)) or ''})",
        text,
    )
    # 其余一律按标签解析为编号，保留正文里原有的「图/表/式/算法/第…节」字样
    text = re.sub(r"\\ref\{([^}]*)\}", sub_ref, text)
    return text


def _compress(nums):
    """[3,4,5,9] -> '3–5,9'，与 LaTeX 端 cite 宏包的压缩规则一致。"""
    nums = sorted(set(nums))
    out, i = [], 0
    while i < len(nums):
        j = i
        while j + 1 < len(nums) and nums[j + 1] == nums[j] + 1:
            j += 1
        if j - i >= 2:
            out.append(f"{nums[i]}–{nums[j]}")
        else:
            out.extend(str(n) for n in nums[i:j + 1])
        i = j + 1
    return ",".join(out)


def process_citations(text):
    r"""\cite{3,4,5} -> 链接到参考文献页的 [3–5]。

    正文里的裸 [N] 早已在 LaTeX 源中统一成 \cite{}，此处只兜底处理残留。
    """
    def sub_cite(m):
        nums = [int(x) for x in re.findall(r"\d+", m.group(1))]
        if not nums:
            return m.group(0)
        label = _compress(nums)
        return f"[[{label}]](references.md#ref-{nums[0]})"

    text = re.sub(r"\\cite\{([\d,\s\-]+)\}", sub_cite, text)
    return text


def process_sections(text):
    """Convert LaTeX section commands to Markdown headers."""
    # \texorpdfstring{显示形式}{书签形式} -> 显示形式（须先于标题正则展开，
    # 否则标题里的嵌套花括号会把 [^}]* 截断）
    out, i = [], 0
    while True:
        j = text.find("\\texorpdfstring{", i)
        if j == -1:
            out.append(text[i:])
            break
        out.append(text[i:j])
        shown, k = brace_arg(text, j + len("\\texorpdfstring"))
        _, k = brace_arg(text, k)
        out.append(shown or "")
        i = k
    text = "".join(out)
    # \chapter*{...} and \chapter{...}
    text = re.sub(r"\\chapter\*\{([^}]*)\}", r"# \1", text)
    text = re.sub(r"\\chapter\{([^}]*)\}", r"# \1", text)
    # \section*{...} and \section{...}
    text = re.sub(r"\\section\*\{([^}]*)\}", r"## \1", text)
    text = re.sub(r"\\section\{([^}]*)\}", r"## \1", text)
    # \subsection*{...} and \subsection{...}
    text = re.sub(r"\\subsection\*\{([^}]*)\}", r"### \1", text)
    text = re.sub(r"\\subsection\{([^}]*)\}", r"### \1", text)
    # \subsubsection*{...} and \subsubsection{...}
    text = re.sub(r"\\subsubsection\*\{([^}]*)\}", r"#### \1", text)
    text = re.sub(r"\\subsubsection\{([^}]*)\}", r"#### \1", text)
    # \paragraph*{...} and \paragraph{...}
    text = re.sub(r"\\paragraph\*?\{([^}]*)\}", r"**\1**", text)
    # \addcontentsline
    text = re.sub(r"\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", text)
    return text


def tighten_math(inner):
    r"""收紧数学块内部的空行。

    \label{...} 被剥掉后原地留下一个空行。Markdown 在空行处断段，$$ 块被
    拆成两个 <p>，arithmatex 匹配不到，公式就以 `$$ ... \tag{} $$` 的字面量
    印在页面上。全站 142 个显示式栽在这上面。
    """
    return re.sub(r"\n[ \t]*\n+", "\n", inner).strip()


def split_to_aligned(inner):
    r"""裸 split 换成 aligned。

    equation / \[ \] 的外壳在转换时被剥掉，split 就没有外层环境了。LaTeX 与
    MathJax 都要求 split 嵌在 equation/align 之类里面，独立使用会报
    "Erroneous nesting of equation structures"。aligned 是它的独立版，
    对齐点与换行行为一致。align 分支不需要这个：那里 split 仍是嵌套的。
    """
    return inner.replace("\\begin{split}", "\\begin{aligned}") \
                .replace("\\end{split}", "\\end{aligned}")


def process_display_equations(text):
    """Convert display equations, preserving LaTeX math for MathJax."""
    # \begin{equation}...\end{equation} -> $$...$$
    def equation_replace(m):
        inner = m.group(1)
        inner = process_labels(inner)
        inner = strip_resizebox(inner)
        inner = split_to_aligned(inner)
        return f"\n$$\n{tighten_math(inner)}\n$$\n"

    text = re.sub(
        r"\\begin\{equation\*?\}(.*?)\\end\{equation\*?\}",
        equation_replace,
        text,
        flags=re.DOTALL,
    )

    # \begin{align}...\end{align} -> $$\begin{aligned}...\end{aligned}$$
    def align_replace(m):
        inner = m.group(1)
        inner = process_labels(inner)
        # 保留 align（不降级为 aligned）：全书 87 个 align 中 71 个带 \tag{}，
        # aligned 里 MathJax 不产生行内编号，会丢掉原书式号
        env = "align" if "\\tag{" in inner else "aligned"
        if env == "aligned":
            inner = re.sub(r"\\nonumber|\\notag", "", inner)
        return f"\n$$\n\\begin{{{env}}}\n{tighten_math(inner)}\n\\end{{{env}}}\n$$\n"

    text = re.sub(
        r"\\begin\{align\*?\}(.*?)\\end\{align\*?\}",
        align_replace,
        text,
        flags=re.DOTALL,
    )

    # \begin{multline}...\end{multline} -> $$...$$（\tag 原样保留）
    text = re.sub(
        r"\\begin\{multline\*?\}(.*?)\\end\{multline\*?\}",
        lambda m: "\n$$\n" + tighten_math(process_labels(m.group(1))) + "\n$$\n",
        text,
        flags=re.DOTALL,
    )

    # \begin{gather}...\end{gather}
    def gather_replace(m):
        inner = m.group(1)
        inner = process_labels(inner)
        inner = re.sub(r"\\nonumber", "", inner)
        return f"\n$$\n\\begin{{gathered}}\n{tighten_math(inner)}\n\\end{{gathered}}\n$$\n"

    text = re.sub(
        r"\\begin\{gather\*?\}(.*?)\\end\{gather\*?\}",
        gather_replace,
        text,
        flags=re.DOTALL,
    )

    # \[...\] -> $$...$$
    def display_replace(m):
        inner = m.group(1)
        inner = strip_resizebox(inner)
        inner = split_to_aligned(inner)
        return f"\n$$\n{tighten_math(inner)}\n$$\n"

    text = re.sub(r"\\\[(.*?)\\\]", display_replace, text, flags=re.DOTALL)

    return text


def strip_resizebox(text):
    """Remove \resizebox{\linewidth}{!}{...} wrapper."""
    # \resizebox{\linewidth}{!}{$\displaystyle ... $}
    text = re.sub(
        r"\\resizebox\{[^}]*\}\{[^}]*\}\{\s*\$\s*\\displaystyle\s*(.*?)\s*\$\s*\}",
        r"\1",
        text,
        flags=re.DOTALL,
    )
    # \resizebox{...}{...}{...}
    text = re.sub(r"\\resizebox\{[^}]*\}\{[^}]*\}\{", "", text)
    return text


def process_figures(text):
    """Convert figure environments to Markdown image syntax."""
    def figure_replace(m):
        inner = m.group(1)
        # Extract \includegraphics path
        img_match = re.search(
            r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", inner
        )
        if not img_match:
            # 无图的 figure 浮动体是算法续排框（如算法 10 的下半部分），
            # 保留正文，交给后续 process_code_blocks 转成代码块，切勿丢弃
            return "\n" + re.sub(r"\\(?:centering|label\{[^}]*\})", "", inner) + "\n"
        img_path = web_image(img_match.group(1))
        # 题注常含 $...$ 与嵌套花括号，必须配对提取，不能用 [^}]*
        caption = ""
        cm = re.search(r"\\caption\{", inner)
        if cm:
            caption, _ = brace_arg(inner, cm.end() - 1)
            caption = (caption or "").strip()
        # 图号取自 main_cn.aux（与原书一致）
        lm = re.search(r"\\label\{(fig:[^}]*)\}", inner)
        num = LABELS.get(lm.group(1)) if lm else None
        if num:
            caption = f"图 {num}　{caption}".rstrip("　")
        # 题注里有括号/引号/数学，放进 alt 或 title 都会破坏 Markdown 链接，
        # 故图下另起一行作斜体题注
        alt = f"图 {num}" if num else ""
        out = f"\n![{alt}](../images/{img_path})\n"
        if caption:
            out += f"\n*{caption}*\n"
        return out

    text = re.sub(
        r"\\begin\{figure\}(.*?)\\end\{figure\}",
        figure_replace,
        text,
        flags=re.DOTALL,
    )
    return text


def process_tables(text):
    """Convert simple tabular environments to Markdown tables."""
    def table_replace(m):
        inner = m.group(0)
        # Extract caption
        caption_match = re.search(r"\\caption\{([^}]*)\}", inner)
        caption = caption_match.group(1) if caption_match else ""
        # Extract tabular content
        tab_match = re.search(
            r"\\begin\{tabular\}\{[^}]*\}(.*?)\\end\{tabular\}", inner, re.DOTALL
        )
        if not tab_match:
            return ""
        tab_content = tab_match.group(1).strip()
        rows = tab_content.split("\\\\")
        md_rows = []
        for row in rows:
            row = row.strip()
            if not row:
                continue
            # Remove \hline
            row = row.replace("\\hline", "").strip()
            if not row:
                continue
            cells = [c.strip() for c in row.split("&")]
            md_rows.append("| " + " | ".join(cells) + " |")

        if len(md_rows) == 0:
            return ""

        # Add header separator after first row
        if len(md_rows) >= 1:
            ncols = len(md_rows[0].split("|")) - 2
            md_rows.insert(1, "|" + "|".join([" --- "] * ncols) + "|")

        result = "\n"
        if caption:
            result += f"**{caption}**\n\n"
        result += "\n".join(md_rows) + "\n"
        return result

    text = re.sub(
        r"\\begin\{table\}.*?\\end\{table\}",
        table_replace,
        text,
        flags=re.DOTALL,
    )
    return text


def process_algorithms(text):
    """Convert custom algorithm environments to code blocks."""
    def algo_replace(m):
        inner = m.group(1)
        # Extract caption
        caption_match = re.search(r"\\caption\{([^}]*)\}", inner)
        caption = caption_match.group(1) if caption_match else ""
        # 算法号由 \setcounter{algorithm}{N-1} 钉死（全书连续编号 1--41）
        num_match = re.search(r"\\setcounter\{algorithm\}\{(\d+)\}", inner)
        if num_match:
            caption = f"算法 {int(num_match.group(1)) + 1}　{caption}".rstrip("　")

        # lstlisting / verbatim 两种代码体
        code_match = re.search(
            r"\\begin\{(?:lstlisting|verbatim)\}(?:\[[^\]]*\])?(.*?)"
            r"\\end\{(?:lstlisting|verbatim)\}",
            inner,
            re.DOTALL,
        )
        if code_match:
            code = code_match.group(1).strip()
        else:
            # Extract algorithmic content
            algo_match = re.search(
                r"\\begin\{algorithmic\*?\}(.*?)\\end\{algorithmic\*?\}",
                inner,
                re.DOTALL,
            )
            if algo_match:
                code = algo_match.group(1).strip()
            else:
                code = inner.strip()

            # 去掉浮动体样板（位置参数、计数器、题注、标签、居中）
            code = re.sub(r"^\s*\[[a-zA-Z!]+\]", "", code)
            code = re.sub(r"\\setcounter\{[^}]*\}\{[^}]*\}", "", code)
            code = re.sub(r"\\caption\{[^}]*\}", "", code)
            code = re.sub(r"\\label\{[^}]*\}", "", code)
            code = re.sub(r"\\centering", "", code)

            # Convert algorithm commands to pseudocode
            code = re.sub(r"\\program\{([^}]*)\}", r"program \1", code)
            code = re.sub(r"\\endprogram", "end program", code)
            code = re.sub(r"\\function\{([^}]*)\}\{([^}]*)\}", r"function \1(\2)", code)
            code = re.sub(r"\\endfunction", "end function", code)
            code = re.sub(r"\\FOR\{([^}]*)\}", r"for \1", code)
            code = re.sub(r"\\For\{([^}]*)\}", r"for \1", code)
            code = re.sub(r"\\ENDFOR|\\EndFor", "end for", code)
            code = re.sub(r"\\WHILE\{([^}]*)\}", r"while \1", code)
            code = re.sub(r"\\While\{([^}]*)\}", r"while \1", code)
            code = re.sub(r"\\ENDWHILE|\\EndWhile", "end while", code)
            code = re.sub(r"\\IF\{([^}]*)\}", r"if \1", code)
            code = re.sub(r"\\If\{([^}]*)\}", r"if \1", code)
            code = re.sub(r"\\ELSE|\\Else", "else", code)
            code = re.sub(r"\\ENDIF|\\EndIf", "end if", code)
            code = re.sub(r"\\RETURN\b|\\Return\b", "return", code)
            code = re.sub(r"\\REQUIRE\b", "Input:", code)
            code = re.sub(r"\\ENSURE\b", "Output:", code)
            code = re.sub(r"\\STATE\b|\\State\b", "", code)
            code = re.sub(r"\\item\b", "", code)
            code = re.sub(r"\\comment\{([^}]*)\}", r"// \1", code)
            code = re.sub(r"\\Comment\{([^}]*)\}", r"// \1", code)
            # Clean up extra whitespace
            code = re.sub(r"\n\s*\n", "\n", code)

        result = "\n"
        if caption:
            result += f"**{caption}**\n\n"
        result += f"```\n{code}\n```\n"
        return result

    text = re.sub(
        r"\\begin\{algorithm\}(.*?)\\end\{algorithm\}",
        algo_replace,
        text,
        flags=re.DOTALL,
    )
    return text


def process_lists(text):
    """Convert itemize and enumerate to Markdown lists."""
    # \begin{enumerate}...\end{enumerate}
    def enum_replace(m):
        inner = m.group(1)
        lines = inner.split("\n")
        md_lines = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            # \item[(a)] -> 1.
            item_match = re.match(r"\\item\[([^\]]+)\]\s*(.*)", line)
            if item_match:
                md_lines.append(f"{item_match.group(1)} {item_match.group(2)}")
                continue
            # \item -> numbered
            item_match2 = re.match(r"\\item\s+(.*)", line)
            if item_match2:
                md_lines.append(f"1. {item_match2.group(1)}")
                continue
            # continuation line
            if md_lines:
                md_lines.append(f"   {line}")
        return "\n" + "\n".join(md_lines) + "\n"

    # 由内向外逐层替换：非贪婪匹配会把嵌套 enumerate 的内层 \begin 吞掉，
    # 留下孤立的外层 \end。这里先匹配不含嵌套的最内层，循环到不再变化。
    innermost_enum = re.compile(
        r"\\begin\{enumerate\}"
        r"((?:(?!\\begin\{enumerate\}|\\end\{enumerate\}).)*?)"
        r"\\end\{enumerate\}",
        re.DOTALL,
    )
    while True:
        new = innermost_enum.sub(enum_replace, text)
        if new == text:
            break
        text = new

    # \begin{itemize}...\end{itemize}
    def itemize_replace(m):
        inner = m.group(1)
        lines = inner.split("\n")
        md_lines = []
        for line in lines:
            line = line.strip()
            if not line:
                continue
            item_match = re.match(r"\\item\s+(.*)", line)
            if item_match:
                md_lines.append(f"- {item_match.group(1)}")
                continue
            if md_lines:
                md_lines.append(f"  {line}")
        return "\n" + "\n".join(md_lines) + "\n"

    innermost_item = re.compile(
        r"\\begin\{itemize\}"
        r"((?:(?!\\begin\{itemize\}|\\end\{itemize\}).)*?)"
        r"\\end\{itemize\}",
        re.DOTALL,
    )
    while True:
        new = innermost_item.sub(itemize_replace, text)
        if new == text:
            break
        text = new

    return text


def process_quotes(text):
    """Convert quote environments to Markdown blockquotes."""
    def quote_replace(m):
        inner = m.group(1).strip()
        lines = inner.split("\n")
        return "\n" + "\n".join(f"> {l}" for l in lines) + "\n"

    text = re.sub(
        r"\\begin\{quote\}(.*?)\\end\{quote\}",
        quote_replace,
        text,
        flags=re.DOTALL,
    )
    return text


def _indent_block(body):
    """把整块内容缩进 4 空格，供 Material admonition 使用（空行不加尾随空格）。"""
    return "\n".join("    " + l if l.strip() else "" for l in body.split("\n"))


def process_custom_envs(text):
    """Convert custom environments (example, boxtext, etc.)."""
    # \begin{exbox}{N}{标题} ... \end{exbox}  -> 例 N（标题）
    # \begin{illbox}{N}{标题} ... \end{illbox} -> 例证 N（标题）
    def box_replace(kind):
        def _sub(m):
            num, title, inner = m.group(1), m.group(2), m.group(3).strip()
            head = f"{kind} {num}（{title}）" if title else f"{kind} {num}"
            return f'\n\n???+ example "{head}"\n\n{_indent_block(inner)}\n\n'
        return _sub

    text = re.sub(
        r"\\begin\{exbox\}\{([^}]*)\}\{([^}]*)\}(.*?)\\end\{exbox\}",
        box_replace("例"),
        text,
        flags=re.DOTALL,
    )
    text = re.sub(
        r"\\begin\{illbox\}\{([^}]*)\}\{([^}]*)\}(.*?)\\end\{illbox\}",
        box_replace("例证"),
        text,
        flags=re.DOTALL,
    )

    # \begin{example}[title]...\end{example}
    def example_replace(m):
        title = m.group(1) or ""
        inner = m.group(2).strip()
        result = '\n???+ example "例'
        if title:
            result += f"：{title}"
        result += '"\n'
        result += f"    {inner}\n"
        return result

    text = re.sub(
        r"\\begin\{example\}(?:\[([^\]]*)\])?(.*?)\\end\{example\}",
        example_replace,
        text,
        flags=re.DOTALL,
    )

    # \begin{boxtext}[title]...\end{boxtext}
    text = re.sub(
        r"\\begin\{boxtext\}(?:\[([^\]]*)\])?(.*?)\\end\{boxtext\}",
        lambda m: f'\n> **{m.group(1) or ""}**\n> {m.group(2).strip()}\n',
        text,
        flags=re.DOTALL,
    )

    # \begin{illustration}[title]...\end{illustration}
    text = re.sub(
        r"\\begin\{illustration\}(?:\[([^\]]*)\])?(.*?)\\end\{illustration\}",
        lambda m: f'\n> **说明：{m.group(1) or ""}**\n> {m.group(2).strip()}\n',
        text,
        flags=re.DOTALL,
    )

    # \begin{question}[title]...\end{question}
    text = re.sub(
        r"\\begin\{question\}(?:\[([^\]]*)\])?(.*?)\\end\{question\}",
        lambda m: f'\n???+ question "问题：{m.group(1) or ""}"\n    {m.group(2).strip()}\n',
        text,
        flags=re.DOTALL,
    )

    # \begin{exercise}[title]...\end{exercise}
    text = re.sub(
        r"\\begin\{exercise\}(?:\[([^\]]*)\])?(.*?)\\end\{exercise\}",
        lambda m: f'\n???+ exercise "练习：{m.group(1) or ""}"\n    {m.group(2).strip()}\n',
        text,
        flags=re.DOTALL,
    )

    return text


def process_code_blocks(text):
    """Convert verbatim and lstlisting to code blocks."""
    # \begin{lstlisting}[options]...\end{lstlisting}
    text = re.sub(
        r"\\begin\{lstlisting\}(?:\[[^\]]*\])?(.*?)\\end\{lstlisting\}",
        lambda m: f"\n```\n{m.group(1).strip()}\n```\n",
        text,
        flags=re.DOTALL,
    )

    # \begin{verbatim}...\end{verbatim} (only standalone, not inside algorithm)
    text = re.sub(
        r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}",
        lambda m: f"\n```\n{m.group(1).strip()}\n```\n",
        text,
        flags=re.DOTALL,
    )
    return text


def process_inline(text):
    """Process inline LaTeX formatting commands."""
    # \textbf{...} -> **...**
    text = re.sub(r"\\textbf\{([^}]*)\}", r"**\1**", text)
    # \textit{...} -> *...*
    text = re.sub(r"\\textit\{([^}]*)\}", r"*\1*", text)
    # \emph{...} -> *...*
    text = re.sub(r"\\emph\{([^}]*)\}", r"*\1*", text)
    # \texttt{...} -> `...`
    # 反引号内是逐字文本，LaTeX 的 \_ \& \% \# \$ 转义要脱掉，
    # 否则 new\_vlist 会照着反斜杠一起显示出来。
    text = re.sub(
        r"\\texttt\{([^}]*)\}",
        lambda m: "`%s`" % re.sub(r"\\([_&%#$])", r"\1", m.group(1)),
        text,
    )
    # \textsuperscript{...} -> ^(...)
    text = re.sub(r"\\textsuperscript\{([^}]*)\}", r"^\1", text)
    # \ldots -> ...
    text = re.sub(r"\\ldots", "...", text)
    # \noindent
    text = re.sub(r"\\noindent", "", text)
    # \centering
    text = re.sub(r"\\centering", "", text)
    # Non-breaking space ~ -> space
    text = text.replace("~", " ")
    # \url{...} -> link
    text = re.sub(r"\\url\{([^}]*)\}", r"[\1](\1)", text)
    # \footnote{...} -> [^N] (collect footnotes, support nested braces)
    footnotes = []

    def extract_brace_content(text, start):
        """Extract content inside matching braces starting at position start.
        start should point to the opening { after \footnote.
        Returns (content, end_position) or (None, original_pos)."""
        depth = 0
        i = start
        while i < len(text):
            if text[i] == '{':
                depth += 1
            elif text[i] == '}':
                depth -= 1
                if depth == 0:
                    return text[start+1:i], i
            i += 1
        return None, start

    footnote_positions = []
    pos = 0
    while True:
        idx = text.find("\\footnote{", pos)
        if idx == -1:
            break
        content, end = extract_brace_content(text, idx + len("\\footnote"))
        if content is not None:
            footnote_positions.append((idx, end + 1, content))
            pos = end + 1
        else:
            pos = idx + 1

    # Replace from end to start to preserve positions
    fn_idx = len(footnote_positions)
    for start, end, content in reversed(footnote_positions):
        fn_idx_current = fn_idx
        fn_idx -= 1
        footnotes.insert(0, content)
        text = text[:start] + f"[^{fn_idx_current}]" + text[end:]

    # Append footnotes at end
    if footnotes:
        text += "\n\n---\n\n"
        for i, fn in enumerate(footnotes, 1):
            text += f"[^{i}]: {fn}\n"
    return text


def process_math_commands(text):
    """Convert custom math shortcut commands to standard LaTeX."""
    # \ket{x} -> |x\rangle
    text = re.sub(r"\\ket\{([^}]*)\}", r"|\\rangle", text)
    # \bra{x} -> \langle x|
    text = re.sub(r"\\bra\{([^}]*)\}", r"\\langle |", text)
    # \braket{x}{y} -> \langle x|y\rangle
    text = re.sub(
        r"\\braket\{([^}]*)\}\{([^}]*)\}", r"\\langle \1|\2\\rangle", text
    )
    # \expect{x} -> \langle x\rangle
    text = re.sub(r"\\expect\{([^}]*)\}", r"\\langle \1\\rangle", text)
    # \muVT -> \mu VT
    text = text.replace("\\muVT", "\\mu VT")
    # \dd -> \mathrm{d} (only standalone \dd, not \ddot)
    text = re.sub(r"\\dd(?![a-zA-Z])", r"\\mathrm{d}", text)
    # \vv -> \mathbf{v} (only standalone)
    text = re.sub(r"\\vv(?![a-zA-Z])", r"\\mathbf{v}", text)
    # \rr -> \mathbf{r}
    text = re.sub(r"\\rr(?![a-zA-Z])", r"\\mathbf{r}", text)
    # \ff -> \mathbf{f}
    text = re.sub(r"\\ff(?![a-zA-Z])", r"\\mathbf{f}", text)
    # \OO -> \mathcal{O}
    text = re.sub(r"\\OO(?![a-zA-Z])", r"\\mathcal{O}", text)
    # \fitdisplay{...} -> just the math
    text = re.sub(r"\\fitdisplay\{([^}]*)\}", r"\1", text)
    # \bm{...} -> \boldsymbol{...}
    text = re.sub(r"\\bm\{([^}]*)\}", r"\\boldsymbol{\1}", text)
    # \sideset{}{'}\sum -> \sum'
    text = text.replace("\\sideset{}{'}\\sum", "\\sum'")
    # \tag{...} 保留：全书 1076 个式号靠它钉死三级号，MathJax 原生支持
    return text


def cleanup(text):
    """Final cleanup of remaining LaTeX artifacts."""
    # Remove remaining \hline, \toprule, \midrule, \bottomrule
    text = re.sub(r"\\[th]?line|\\toprule|\\midrule|\\bottomrule", "", text)
    # Remove \vspace{...}, \hspace{...}
    text = re.sub(r"\\[vh]space\{[^}]*\}", "", text)
    # Remove \small, \footnotesize, \normalsize, \large etc
    text = re.sub(r"\\(?:tiny|scriptsize|footnotesize|small|normalsize|large|Large|LARGE|huge|Huge)", "", text)
    # Remove \rmfamily, \sffamily, \ttfamily, \cjkfont
    text = re.sub(r"\\(?:rm|sf|tt)family", "", text)
    text = re.sub(r"\\cjkfont\{[^}]*\}", "", text)
    # Remove \newpage, \clearpage, \pagebreak
    text = re.sub(r"\\(?:newpage|clearpage|pagebreak)", "", text)
    # \footnotemark 的脚注号已由 \footnote 处理流程给出，这里去掉命令残留
    text = re.sub(r"\\footnotemark\b", "", text)
    # 重音/特殊字母
    text = re.sub(r"\\AA\{\}|\\AA\b", "Å", text)
    text = re.sub(r"\\o\{\}|\\o(?![a-zA-Z])", "ø", text)
    # Remove remaining braces from simple commands
    # Clean multiple blank lines
    text = re.sub(r"\n{3,}", "\n\n", text)
    # Remove trailing whitespace
    text = re.sub(r" +\n", "\n", text)
    return text.strip()


# 数学、代码、链接目标、HTML 标签、attr_list 花括号里不能动
PANGU_SKIP = re.compile(
    r"```.*?```|~~~.*?~~~|\$\$.*?\$\$|\$[^$\n]*\$|`[^`\n]*`"
    r"|\]\([^)]*\)|<[^>\n]*>|\{[^}\n]*\}",
    re.S)
PANGU_CJK_FIRST = re.compile(r"(?<=[\u4e00-\u9fff])(?=[A-Za-z0-9])")
PANGU_LAT_FIRST = re.compile(r"(?<=[A-Za-z0-9\)\]])(?=[\u4e00-\u9fff])")


def pangu(text):
    """中英/中数之间补一个空格（俗称「盘古之白」）。

    LaTeX 端交给 xeCJK，源文件里空格可有可无；Markdown 端没有自动字距，
    这里统一补齐，顺带把译稿里本来就不一致的那 4000 处也一并抹平。
    """
    out, last = [], 0
    for m in PANGU_SKIP.finditer(text):
        out.append(_pangu_plain(text[last:m.start()]))
        out.append(m.group(0))
        last = m.end()
    out.append(_pangu_plain(text[last:]))
    return "".join(out)


# LaTeX 的 en dash；三个及以上是 markdown 的分隔线/表头，不能碰
EN_DASH = re.compile(r"(?<!-)--(?!-)")


def _pangu_plain(s):
    s = EN_DASH.sub("\u2013", s)
    return PANGU_LAT_FIRST.sub(" ", PANGU_CJK_FIRST.sub(" ", s))


def convert_tex_to_md(tex_content):
    """Full conversion pipeline."""
    text = tex_content

    # Remove comment lines
    text = re.sub(r"^\s*%.*$", "", text, flags=re.MULTILINE)

    # 索引标记只服务纸质版索引页，在线版靠站内搜索。必须在最前面剥掉：
    # 图注、例证框等处理会先吃掉一层大括号，到后面就只剩半截匹配不上了。
    text = re.sub(r"\\sidx\{[^{}]*\}\{[^{}]*\}", "", text)
    text = re.sub(r"\\aidx\{[^{}]*\}", "", text)

    # Process in order: blocks first, then inline
    text = process_display_equations(text)
    text = process_figures(text)
    text = process_boxed_figures(text)
    text = process_tables(text)
    text = process_algorithms(text)
    text = process_lists(text)
    text = process_quotes(text)
    text = process_custom_envs(text)
    text = process_code_blocks(text)
    text = process_sections(text)
    text = process_labels(text)
    text = process_refs(text)
    text = process_math_commands(text)
    text = process_inline(text)
    text = process_citations(text)
    text = cleanup(text)
    text = pangu(text)

    return text


def sync_images():
    """把各章实际引用到的图拷进 docs/images/（CI 上没有源目录，必须入库）。"""
    import shutil
    images_dir = os.path.join(DOCS_DIR, "images")
    os.makedirs(images_dir, exist_ok=True)
    used, copied, missing = set(), 0, []
    for md in sorted(os.listdir(os.path.join(DOCS_DIR, "chapters"))):
        if not md.endswith(".md"):
            continue
        with open(os.path.join(DOCS_DIR, "chapters", md), encoding="utf-8") as f:
            used |= set(re.findall(r"\]\(\.\./images/([^\s)]+)", f.read()))
    for name in sorted(used):
        src = os.path.join(FIGURES_DIR, name)
        dst = os.path.join(images_dir, name)
        if not os.path.exists(src):
            missing.append(name)
            continue
        if not os.path.exists(dst) or os.path.getsize(src) != os.path.getsize(dst):
            shutil.copy2(src, dst)
            copied += 1
    print(f"图片：正文引用 {len(used)} 张，新增/更新 {copied} 张，缺失 {len(missing)} 张")
    for name in missing:
        print(f"  缺失: {name}")
    return missing


def generate_references():
    r"""把 references_cn.tex 的 \bibitem 转成带锚点的参考文献页。

    正文里的 \cite{N} 会链到本页的 #ref-N，锚点用 attr_list 语法挂在条目段落上。
    """
    path = os.path.join(TEX_DIR, "references_cn.tex")
    if not os.path.exists(path):
        print(f"警告：找不到 {path}，跳过参考文献页")
        return 0
    with open(path, encoding="utf-8") as f:
        tex = f.read()
    body = tex.split(r"\begin{thebibliography}", 1)[-1]
    body = body.split(r"\end{thebibliography}", 1)[0]

    entries = []
    for m in re.finditer(r"\\bibitem\{(\d+)\}(.*?)(?=\\bibitem\{|\Z)", body, re.S):
        num, txt = int(m.group(1)), m.group(2)
        txt = re.sub(r"%.*", "", txt)                       # 行末注释
        txt = re.sub(r"\\(textit|textbf|emph|texttt)\{([^{}]*)\}", r"\2", txt)
        txt = txt.replace("~", " ").replace(r"\&", "&").replace("\\%", "%")
        txt = re.sub(r"(?<!-)--(?!-)", "\u2013", txt)
        txt = " ".join(txt.split())
        # 译稿仓把 DOI 包进了 \url{}（LaTeX 靠它断行），这边先脱掉，
        # 否则下面那步会把大括号一起卷进链接地址
        txt = re.sub(r"\\url\{([^}]*)\}", r"\1", txt)
        # 裸 DOI/URL 转成链接
        txt = re.sub(r"(https?://[^\s,]+?)(\.?)(?=$|[\s,])", r"[\1](\1)\2", txt)
        if txt:
            entries.append((num, txt))

    entries.sort()
    lines = [
        "# 参考文献",
        "",
        f"原书参考文献共 {len(entries)} 条，编号与英文原版一致；",
        "正文中的 \\[N] 均链接到此处对应条目。",
        "",
    ]
    for num, txt in entries:
        lines.append(f"**[{num}]** {txt}")
        lines.append("{ #ref-%d }" % num)
        lines.append("")
    out = os.path.join(DOCS_DIR, "chapters", "references.md")
    with open(out, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"Wrote references.md（{len(entries)} 条）")
    return len(entries)


def main():
    os.makedirs(os.path.join(DOCS_DIR, "chapters"), exist_ok=True)
    load_labels()

    # Write index page
    index_content = """# 理解分子模拟：从算法到应用（第三版）

!!! info "关于本书"

    本书是 Frenkel & Smit《Understanding Molecular Simulation: From Algorithms to Applications》第三版的中文翻译版。

    **原作者：** Daan Frenkel, Berend Smit

    **译者：** 杨凯

---

## 目录

### 前言

- [译者序](chapters/preface-berend.md)
- [第一版前言](chapters/preface1.md)
- [第二版前言](chapters/preface2.md)
- [第三版前言](chapters/preface3.md)

### 正文

| 章节 | 标题 |
|------|------|
| [第1章](chapters/ch01-introduction.md) | 绪论 |
| [第2章](chapters/ch02-thermo-statmech.md) | 统计力学基础 |
| [第3章](chapters/ch03-monte-carlo.md) | 蒙特卡洛模拟 |
| [第4章](chapters/ch04-md.md) | 分子动力学模拟 |
| [第5章](chapters/ch05-computer-exp.md) | 计算机实验 |
| [第6章](chapters/ch06-mc-ensembles.md) | 蒙特卡洛系综 |
| [第7章](chapters/ch07-md-ensembles.md) | 分子动力学系综 |
| [第8章](chapters/ch08-free-energy.md) | 自由能计算 |
| [第9章](chapters/ch09-free-energy-solids.md) | 固体自由能 |
| [第10章](chapters/ch10-free-energy-chains.md) | 链状分子自由能 |
| [第11章](chapters/ch11-long-ranged.md) | 长程相互作用 |
| [第12章](chapters/ch12-cbmc.md) | 构型偏倚蒙特卡洛 |
| [第13章](chapters/ch13-accel-mc.md) | 加速蒙特卡洛方法 |
| [第14章](chapters/ch14-timescale.md) | 时间尺度加速方法 |
| [第15章](chapters/ch15-rare-events.md) | 稀有事件 |
| [第16章](chapters/ch16-mesoscopic.md) | 介观模拟方法 |

### 附录

| 附录 | 标题 |
|------|------|
| [附录A](chapters/appA-lagrangian.md) | 拉格朗日力学 |
| [附录B](chapters/appB-non-hamiltonian.md) | 非哈密顿动力学 |
| [附录C](chapters/appC-kirkwood-buff.md) | Kirkwood-Buff理论 |
| [附录D](chapters/appD-non-eq-thermo.md) | 非平衡热力学 |
| [附录E](chapters/appE-non-eq-work.md) | 非平衡功 |
| [附录F](chapters/appF-linear-response.md) | 线性响应理论 |
| [附录G](chapters/appG-committor.md) | Committor分析 |
| [附录H](chapters/appH-sdpd.md) | 确定性颗粒动力学 |
| [附录I](chapters/appI-saving-cpu.md) | 节省CPU时间 |
| [附录J](chapters/appJ-algorithms.md) | 算法汇总 |

## 参考文献

[参考文献总表](chapters/references.md)（786 条，编号与英文原版一致，正文引用可直接跳转）

## 术语

[英汉术语对照表](glossary.md)（456 条，与主题索引同源；含人名类术语的写法规则）

---

## 版权与授权状态

本书译自 Daan Frenkel and Berend Smit, *Understanding Molecular Simulation:
From Algorithms to Applications*, Third Edition, Academic Press（Elsevier
旗下品牌）, 2023, ISBN 978-0-323-90292-2。

Original English language edition copyright © 2023 Elsevier Inc.
All rights reserved. 原著版权 © 2023 Elsevier Inc.，保留所有权利。
中文译文版权 © 2026 杨凯。

**本译本尚未取得 Elsevier 的中文翻译出版授权。** 本站为译稿的在线阅读版本，
仅供学习与审校参考，不作商业用途；正式出版须由中文出版社完成版权引进。
若权利人认为本站内容不当，请联系译者，将立即下线。

**翻译范围**：对应原著纸质版第 1–16 章与附录 A–J。原著附录 K–P 为出版商
在线补充材料（页码带 e 前缀），不在本译本范围内。
"""
    with open(os.path.join(DOCS_DIR, "index.md"), "w", encoding="utf-8") as f:
        f.write(index_content)
    print(f"Wrote index.md")

    # Convert each chapter
    global CURRENT_CHAP
    for slug, tex_file, title in CHAPTERS:
        print(f"Converting {tex_file} -> {slug}.md ...")
        m = re.match(r"ch0?(\d+)-", slug) or re.match(r"app([A-J])-", slug)
        CURRENT_CHAP = m.group(1) if m else ""
        tex = read_tex(tex_file)
        md = convert_tex_to_md(tex)
        path = write_md(slug, md)
        lines = md.count("\n") + 1
        print(f"  -> {path} ({lines} lines)")

    generate_references()

    print(f"\nDone! {len(CHAPTERS)} chapters converted.")
    sync_images()


if __name__ == "__main__":
    main()
