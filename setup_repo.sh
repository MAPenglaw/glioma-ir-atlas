#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# 把 reproducibility/ 推成公开 GitHub 仓库，并引导拿到 Zenodo DOI
#
# 用法（在本文件所在目录执行）：
#   GH_USER=你的GitHub用户名 REPO_NAME=glioma-intron-retention-atlas ./setup_repo.sh
#
# 只做本地提交不推送：
#   DRY_RUN=1 GH_USER=xxx ./setup_repo.sh
#
# 前置：brew install gh  &&  gh auth login
# ---------------------------------------------------------------------------
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GH_USER="${GH_USER:-}"
REPO_NAME="${REPO_NAME:-glioma-intron-retention-atlas}"
DRY_RUN="${DRY_RUN:-0}"
TAG="v1.0"

cd "$REPO_DIR"

echo "== 步骤 0 · 前置检查 =="
command -v git >/dev/null || { echo "缺 git"; exit 1; }
if [ "$DRY_RUN" = "0" ]; then
  command -v gh  >/dev/null || { echo "缺 gh，请先 brew install gh"; exit 1; }
  gh auth status >/dev/null 2>&1 || { echo "gh 未登录，请先 gh auth login"; exit 1; }
  [ -z "$GH_USER" ] && { echo "请设置 GH_USER=你的GitHub用户名"; exit 1; }
fi

echo
echo "== 步骤 1 · 隐私复查（必须逐条看过）=="
echo "-- 本机绝对路径（可接受，README 已声明需替换 ROOT/DATA）:"
grep -rn "/Users/" code/ README.md | wc -l | sed 's/^/   命中行数: /'
echo "-- 邮箱 / 密钥（通讯邮箱白名单外必须为空）:"
# 通讯作者邮箱是有意公开的，加入白名单；其余邮箱与任何密钥都视为风险。
CONTACT_EMAIL="${CONTACT_EMAIL:-mapenglaw@163.com}"
HITS=$(grep -rniE "[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Za-z]{2,}|api[_-]?key|password|secret|token|sk-[A-Za-z0-9]{10,}" \
         code/ README.md .zenodo.json 2>/dev/null | grep -v -F "$CONTACT_EMAIL" || true)
if [ -n "$HITS" ]; then
  echo "$HITS"
  echo "   !! 发现疑似敏感信息（非白名单邮箱或密钥），先清掉再继续"; exit 1
else
  echo "    无（已放行通讯邮箱 ${CONTACT_EMAIL}）"
fi
echo "-- 残留 [[TO FILL]] 标记（推上去前建议先填好作者信息）:"
grep -rn "\[\[TO FILL" README.md .zenodo.json || echo "    无"

echo
echo "== 步骤 2 · 本地提交 =="
if [ ! -d .git ]; then
  git init -q
  git config user.name  "${GIT_NAME:-${GH_USER:-author}}"
  git config user.email "${GIT_EMAIL:-author@example.org}"
fi
git add -A
git commit -q -m "Analysis code for: global intron retention burden in glioma (v1.0)" || echo "    (无新改动可提交)"
git tag -f "$TAG" >/dev/null 2>&1 || true
echo "    已提交并打 tag $TAG"

if [ "$DRY_RUN" = "1" ]; then
  echo; echo "DRY_RUN=1，到此为止，未推送。"; exit 0
fi

echo
echo "== 步骤 3 · 建仓并推送 =="
if ! git remote get-url origin >/dev/null 2>&1; then
  gh repo create "$REPO_NAME" --public --source=. \
     --description "Analysis code for: global intron retention burden in glioma (TCGA GBM/LGG, pan-cancer validation, independent-pipeline replication)" \
     --push
else
  git push -u origin main 2>/dev/null || git push -u origin master
fi
git push -f origin "$TAG"
REPO_URL="https://github.com/${GH_USER}/${REPO_NAME}"
echo "    仓库: $REPO_URL"

cat <<MSG

== 步骤 4 · 拿 Zenodo DOI（网页操作，脚本做不了）==

1. 打开 https://zenodo.org/account/settings/github/ ，登录（可用 GitHub 账号）
2. 在仓库列表里把 ${GH_USER}/${REPO_NAME} 的开关打开（Enable）
3. 回到 GitHub 仓库页 → Releases → 用刚推的 tag ${TAG} 发一个 Release
4. Zenodo 会在几分钟内抓取并分配 DOI，形如 10.5281/zenodo.XXXXXXX
   （若超过 30 分钟没动静：进 Zenodo → Upload → GitHub 页签 → 手动点仓库名）

== 步骤 5 · 回填 DOI 与仓库地址（共 5 处）==

  sed -i '' 's|\\[\\[TO FILL: GitHub repository URL\\]\\]|$REPO_URL|' \\
      ../manuscript/MANUSCRIPT_DRAFT.md ../manuscript/COVER_LETTER.md
  sed -i '' 's|\\[\\[TO FILL: Zenodo DOI, e.g. 10.5281/zenodo.XXXXXXX\\]\\]|10.5281/zenodo.XXXXXXX|' \\
      ../manuscript/MANUSCRIPT_DRAFT.md ../manuscript/COVER_LETTER.md

  然后重跑：
      cd ../manuscript && python3 build_jtm.py && python3 md2docx_jtm.py

  另外 README.md 的 Contact 与 .zenodo.json 的 creators 也要填真名（步骤 1 会提示）。
MSG
