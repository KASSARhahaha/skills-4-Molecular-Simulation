#!/bin/bash
# 构建在线版：书稿 .tex -> Markdown -> 静态站点。退出码 0 表示可发布。
#
# 闸门跑两趟：mkdocs 之前量 Markdown，之后量成品 HTML。**两趟都要**——.md
# 干净不等于渲染出来是对的。全书 108 块行间公式曾经在 .md 里写法完全正常，
# 只是缩进差一格，arithmatex 接不走，成品上是一框源码；那一轮 .md 侧全绿、
# mkdocs --strict 也通过。

set -uo pipefail
cd "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)" || exit 1

echo "=== [1/4] 书稿 -> Markdown ==="
python3 tex2md.py || exit 1

echo
echo "=== [2/4] 闸门（Markdown 侧）==="
md_bad=0
md_out=$(python3 check_md.py 2>&1) || md_bad=1
printf '  %s\n' "$(printf '%s' "$md_out" | grep -o '在线版保真 .*' | tail -1)"
if [ "$md_bad" -ne 0 ]; then
    printf '%s\n' "$md_out" | grep '✗' | head -20 | sed 's/^/    /'
    echo "✗ 在线版未达标"; exit 1
fi

echo
echo "=== [3/4] 静态站点 ==="
.venv/bin/mkdocs build --strict > /tmp/mkdocs_build.log 2>&1 || {
    tail -20 /tmp/mkdocs_build.log; echo "✗ mkdocs 构建失败"; exit 1; }
printf '  %s\n' "$(grep -o 'Documentation built in .*' /tmp/mkdocs_build.log | tail -1)"

echo
echo "=== [4/4] 闸门（成品 HTML 侧）==="
h_bad=0
h_out=$(python3 check_md.py --html 2>&1) || h_bad=1
printf '  %s\n' "$(printf '%s' "$h_out" | grep -o '在线版成品 .*' | tail -1)"
if [ "$h_bad" -ne 0 ]; then
    printf '%s\n' "$h_out" | grep '✗' | head -20 | sed 's/^/    /'
    echo "✗ 成品未达标"; exit 1
fi

echo
echo "✓ 在线版构建通过（可发布）"
