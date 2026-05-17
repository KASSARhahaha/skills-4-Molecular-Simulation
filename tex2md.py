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
    """Remove \label{...} commands."""
    return re.sub(r"\\label\{[^}]*\}", "", text)


def process_refs(text):
    """Convert \ref{...} to a placeholder readable form."""
    # 图~\ref{fig:...} -> 图 X
    text = re.sub(r"[~ ]*\\ref\{fig:([^}]*)\}", "", text)
    # 表~\ref{tab:...}
    text = re.sub(r"[~ ]*\\ref\{tab:([^}]*)\}", "", text)
    # 第~X~节 -> section ref
    text = re.sub(r"第~?\\ref\{sec:([^}]*)\}~?节", "", text)
    text = re.sub(r"\\ref\{sec:([^}]*)\}", "", text)
    # 式~(\ref{eq:...})
    text = re.sub(r"[~ ]*\(?~?\\ref\{eq:([^}]*)\}~?\)?", "", text)
    # 算法~\ref{alg:...}
    text = re.sub(r"[~ ]*\\ref\{alg:([^}]*)\}", "", text)
    # chapter ref
    text = re.sub(r"\\ref\{ch:([^}]*)\}", "", text)
    # app ref
    text = re.sub(r"\\ref\{app:([^}]*)\}", "", text)
    # generic \ref
    text = re.sub(r"\\ref\{([^}]*)\}", "", text)
    return text


def process_citations(text):
    """Convert [N] style citations to superscript."""
    # [1], [2,3], [3--5], [21--24,26,105]
    text = re.sub(r"\[(\d+(?:[,\-\d]*\d)*)\]", r"^[\\ref{\1}]", text)
    # \cite{...}
    text = re.sub(r"\\cite\{([^}]*)\}", "", text)
    return text


def process_sections(text):
    """Convert LaTeX section commands to Markdown headers."""
    # \chapter*{...} and \chapter{...}
    text = re.sub(r"\\chapter\*\{([^}]*)\}", r"# \1", text)
    text = re.sub(r"\\chapter\{([^}]*)\}", r"# \1", text)
    # \section*{...} and \section{...}
    text = re.sub(r"\\section\*\{([^}]*)\}", r"## \1", text)
    text = re.sub(r"\\section\{([^}]*)\}", r"## \1", text)
    # \subsection*{...} and \subsection{...}
    text = re.sub(r"\\subsection\*\{([^}]*)\}", r"### \1", text)
    text = re.sub(r"\\subsection\{([^}]*)\}", r"### \1", text)
    # \subsubsection{...}
    text = re.sub(r"\\subsubsection\{([^}]*)\}", r"#### \1", text)
    # \paragraph{...}
    text = re.sub(r"\\paragraph\{([^}]*)\}", r"**\1**", text)
    # \addcontentsline
    text = re.sub(r"\\addcontentsline\{[^}]*\}\{[^}]*\}\{[^}]*\}", "", text)
    return text


def process_display_equations(text):
    """Convert display equations, preserving LaTeX math for MathJax."""
    # \begin{equation}...\end{equation} -> $$...$$
    def equation_replace(m):
        inner = m.group(1)
        inner = process_labels(inner)
        inner = strip_resizebox(inner)
        return f"\n$$\n{inner.strip()}\n$$\n"

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
        inner = re.sub(r"\\nonumber", "", inner)
        inner = re.sub(r"\\notag", "", inner)
        return f"\n$$\n\\begin{{aligned}}\n{inner.strip()}\n\\end{{aligned}}\n$$\n"

    text = re.sub(
        r"\\begin\{align\*?\}(.*?)\\end\{align\*?\}",
        align_replace,
        text,
        flags=re.DOTALL,
    )

    # \begin{gather}...\end{gather}
    def gather_replace(m):
        inner = m.group(1)
        inner = process_labels(inner)
        inner = re.sub(r"\\nonumber", "", inner)
        return f"\n$$\n\\begin{{gathered}}\n{inner.strip()}\n\\end{{gathered}}\n$$\n"

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
        return f"\n$$\n{inner.strip()}\n$$\n"

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
            return ""
        img_path = img_match.group(1)
        # Normalize path: strip figures/ prefix, make relative
        img_path = img_path.replace("figures/", "")
        # Extract caption
        caption_match = re.search(r"\\caption\{([^}]*)\}", inner)
        caption = caption_match.group(1) if caption_match else ""
        # Use relative path from docs/chapters/ to ../../images/
        return f'\n![{caption}](../../images/{img_path} "{caption}")\n'

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

        # Check for verbatim content
        verbatim_match = re.search(
            r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}", inner, re.DOTALL
        )
        if verbatim_match:
            code = verbatim_match.group(1).strip()
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

    text = re.sub(
        r"\\begin\{enumerate\}(.*?)\\end\{enumerate\}",
        enum_replace,
        text,
        flags=re.DOTALL,
    )

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

    text = re.sub(
        r"\\begin\{itemize\}(.*?)\\end\{itemize\}",
        itemize_replace,
        text,
        flags=re.DOTALL,
    )

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


def process_custom_envs(text):
    """Convert custom environments (example, boxtext, etc.)."""
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
    text = re.sub(r"\\texttt\{([^}]*)\}", r"`\1`", text)
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
    # \footnote{...} -> [^N] (collect footnotes)
    footnotes = []

    def footnote_replace(m):
        idx = len(footnotes) + 1
        footnotes.append(m.group(1))
        return f"[^{idx}]"

    text = re.sub(r"\\footnote\{([^}]*)\}", footnote_replace, text)
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
    # \dd -> \mathrm{d}
    text = text.replace("\\dd", "\\mathrm{d}")
    # \vv -> \mathbf{v}
    text = text.replace("\\vv", "\\mathbf{v}")
    # \rr -> \mathbf{r}
    text = text.replace("\\rr", "\\mathbf{r}")
    # \ff -> \mathbf{f}
    text = text.replace("\\ff", "\\mathbf{f}")
    # \OO -> \mathcal{O}
    text = text.replace("\\OO", "\\mathcal{O}")
    # \fitdisplay{...} -> just the math
    text = re.sub(r"\\fitdisplay\{([^}]*)\}", r"\1", text)
    # \bm{...} -> \boldsymbol{...}
    text = re.sub(r"\\bm\{([^}]*)\}", r"\\boldsymbol{\1}", text)
    # \sideset{}{'}\sum -> \sum'
    text = text.replace("\\sideset{}{'}\\sum", "\\sum'")
    # \tag{...} - remove tags
    text = re.sub(r"\\tag\{[^}]*\}", "", text)
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
    # Remove remaining braces from simple commands
    # Clean multiple blank lines
    text = re.sub(r"\n{3,}", "\n\n", text)
    # Remove trailing whitespace
    text = re.sub(r" +\n", "\n", text)
    return text.strip()


def convert_tex_to_md(tex_content):
    """Full conversion pipeline."""
    text = tex_content

    # Remove comment lines
    text = re.sub(r"^\s*%.*$", "", text, flags=re.MULTILINE)

    # Process in order: blocks first, then inline
    text = process_display_equations(text)
    text = process_figures(text)
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

    return text


def main():
    os.makedirs(os.path.join(DOCS_DIR, "chapters"), exist_ok=True)

    # Copy figures to docs/images/
    images_dir = os.path.join(DOCS_DIR, "images")
    if os.path.exists(FIGURES_DIR):
        os.makedirs(images_dir, exist_ok=True)
        # Create symlink instead of copying
        link_path = images_dir
        if not os.path.exists(link_path):
            os.symlink(FIGURES_DIR, link_path)
            print(f"Linked figures: {FIGURES_DIR} -> {link_path}")

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
"""
    with open(os.path.join(DOCS_DIR, "index.md"), "w", encoding="utf-8") as f:
        f.write(index_content)
    print(f"Wrote index.md")

    # Convert each chapter
    for slug, tex_file, title in CHAPTERS:
        print(f"Converting {tex_file} -> {slug}.md ...")
        tex = read_tex(tex_file)
        md = convert_tex_to_md(tex)
        path = write_md(slug, md)
        lines = md.count("\n") + 1
        print(f"  -> {path} ({lines} lines)")

    print(f"\nDone! {len(CHAPTERS)} chapters converted.")


if __name__ == "__main__":
    main()
