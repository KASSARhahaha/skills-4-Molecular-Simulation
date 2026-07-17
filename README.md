# Skills for Molecular Simulation · 分子模拟技能集

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Deploy](https://github.com/KASSARhahaha/skills-4-Molecular-Simulation/actions/workflows/ci.yml/badge.svg)](https://github.com/KASSARhahaha/skills-4-Molecular-Simulation/actions)
[![Site](https://img.shields.io/badge/site-live-teal)](https://kassarhahaha.github.io/skills-4-Molecular-Simulation/)
[![Language: zh-CN](https://img.shields.io/badge/lang-zh--CN-red)](docs/)
[![MkDocs Material](https://img.shields.io/badge/MkDocs-Material-02aff6)](https://squidfunk.github.io/mkdocs-material/)

> 开放的中文分子模拟 / CCUS / 吸附教学知识库。把 6 本领域标准教材转成可搜索、可交互、带 LaTeX 公式的网页，并提供 Claude Code skill 包供研究者本地查阅。

🌐 **公开站**: <https://kassarhahaha.github.io/skills-4-Molecular-Simulation/>

---

## 项目背景

分子模拟、碳捕集与封存（CCUS）、吸附过程模拟是计算化学、化学工程与能源研究的基础工具。但中文领域缺乏系统、开放、高质量的教学参考资源 — 大多数经典教材只有英文原版，研究生入门门槛高。

本仓库把以下 **6 本书** 的核心内容转化为结构化衍生品（速查 / 模式 / 公式 / 推导 / 术语），并配合 Frenkel & Smit《理解分子模拟》第三版的完整中文译本（译者：杨凯，482 页 XeLaTeX 排版），形成完整的学习参考体系。

## 收录书目

| Slug | 原书 | 语言 | 章节 | 公式数 |
|---|---|---|---|---|
| `molsim-cn` | Frenkel & Smit《理解分子模拟：从算法到应用》（第三版中译本） | 中文 | 16 章 + 10 附录 | 7421 OMML |
| `molsim-tutorial-exercises` | *Lab Course MolSim 2026* (EPFL/UC Berkeley 配套练习) | English | 12 章 | — |
| `ccs-intro` | *Introduction to Carbon Capture and Sequestration* (Berkeley Lectures Vol 1) | English | 11 章 | — |
| `adsorption-sim` | *Simulation of Adsorption Process* (+ DOCX) | English | 10 章 | — |
| `adsorption-sim-notes` | *Simulation of Adsorption Process — Notes* | English | 7 章 | — |
| `aspen-adsorption` | 《Aspen Adsorption 在气体吸附过程模拟方面的应用》 | 中文 | 5 章 | — |

公开站目前只挂高密度衍生品（速查 / 模式 / 公式 / 推导 / 术语），**不包含书的原文**，避免版权风险。

## `molsim-cn` 三层知识结构（试点）

`molsim-cn` 已完成深度升级，作为整套仓库的质量标杆：

| 层级 | 文件 | 用途 |
|---|---|---|
| 📖 **完整译本** | `docs/chapters/` (26 章) | Frenkel & Smit 第三版中文译本，XeLaTeX 排版 |
| ⚡ **速查** | `skills/molsim-cn/{cheatsheet,formulas,glossary}.md` | 公式速查 / 决策表 / 中英术语 |
| 🎓 **教学** | `skills/molsim-cn/derivations.md` | 从第一性原理推导 10 个核心算法（Metropolis、Verlet、Nosé-Hoover、LJ 尾部校正、Ewald、FEP、TI、BAR、Crooks/Jarzynski、Green-Kubo） |
| 🧩 **模式** | `skills/molsim-cn/patterns.md` | 设计模式 / 反模式 |

其余 5 本待按相同标准升级。

## 本地使用

### 浏览在线内容

```bash
git clone https://github.com/KASSARhahaha/skills-4-Molecular-Simulation.git
cd skills-4-Molecular-Simulation
pip install mkdocs-material
mkdocs serve
# 浏览器打开 http://localhost:8000
```

### 安装为 Claude Code skill（含逐章摘要）

完整 skill（71 个 chapter files + 24 个 supporting files）在作者本地维护。如需自行生成：

```bash
# 使用 book-to-skill 工具
git clone https://github.com/virgiliojr94/book-to-skill.git ~/.claude/skills/book-to-skill
# 在 Claude Code 会话里
/book-to-skill ~/path/to/book.pdf <skill-slug>
```

### 把中文 .tex 转 Word（保留 OMML 公式）

使用 [tex2word](https://github.com/yfyang86/tex2word)（独立项目，99.96% 公式转 OMML）：

```bash
uv tool install --python 3.12 "tex2word[pdf,mathml,csl]"
tex2word convert ch03_monte_carlo_cn.tex -o ch03.docx
```

## 项目结构

```
skills-4-Molecular-Simulation/
├── docs/
│   ├── chapters/        # 中文译本 26 章（molsim-cn）
│   ├── images/          # 图文件
│   ├── skills/          # 6 本书的衍生品
│   │   ├── index.md     # 总入口
│   │   ├── molsim-cn/   # 含 formulas.md / derivations.md
│   │   ├── ccs-intro/
│   │   ├── adsorption-sim/
│   │   ├── adsorption-sim-notes/
│   │   ├── aspen-adsorption/
│   │   └── molsim-tutorial-exercises/
│   └── index.md         # 站点首页
├── .github/workflows/ci.yml   # GitHub Pages 自动部署
├── mkdocs.yml          # MathJax + arithmatex 已配置
└── LICENSE             # MIT
```

## 路线图

- [x] **molsim-cn 深度升级**：保留 LaTeX 公式 + 手工策展 formulas/derivations
- [x] **公开站 CI**：MkDocs strict build + GitHub Actions 自动部署
- [x] **6 本书的 cheatsheet/patterns/glossary** 衍生品
- [ ] **扩展深度升级到其余 5 本**（CCUS、吸附、Aspen、MolSim tutorial）
- [ ] **整书合并 docx**：跨章 `\ref` 解析 + 单一 .docx 输出
- [ ] **双语搜索**：中英术语对齐 + 跨书交叉引用
- [ ] **集成 EPFL/UC Berkeley 课程练习**（Day 0–9）
- [ ] **docx 下载**：把章节 docx 挂到公开站供读者下载

## 贡献

欢迎 issue / PR：
- 中文译本校对（typo、术语统一、公式排版）
- 新书扩展（提供 PDF/源文件 + 主题）
- 公式 / 推导的补充与勘误
- 翻译其他语种（如需合作出英中对照版）

请先开 issue 讨论方向，再提交 PR。所有改动需通过 `mkdocs build --strict` 零警告。

## 致谢

- **Daan Frenkel & Berend Smit** — 原书作者
- **EPFL / UC Berkeley MolSim 2026 课程团队** — 配套练习
- [virgiliojr94/book-to-skill](https://github.com/virgiliojr94/book-to-skill) — skill 自动生成工具
- [tex2word](https://github.com/yfyang86/tex2word) — LaTeX → Word 转换
- [MkDocs Material](https://squidfunk.github.io/mkdocs-material/) — 文档主题

## 许可证

[MIT](LICENSE) — 代码、衍生品、中文译本（杨凯译作）。

原书及原图版权归 Academic Press / Elsevier 所有；本仓库仅发布衍生教学摘要，不包含原书正文与原图。

---

⭐ 如果项目对你有帮助，欢迎 star — 这对我申请 Codex for Open Source 维护者计划有帮助。
