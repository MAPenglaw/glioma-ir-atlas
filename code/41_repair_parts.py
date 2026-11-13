#!/usr/bin/env python3
"""41 — 修复被 kill 中断过的分片（0/1/6/7）尾段并重做整体流式校验。

背景：37 号在下载中被手工 kill，其 .tailN 文件末尾写入错位（gzip 在
偏移 ~79.7 MB 处 invalid block type）。本脚本对每个受损分片：
    保留 tail1 的前 (size - BACK) 字节，从该点起重下到分片末尾。
最后按序 cat 成 tcga_RSEM_gene_tpm.gz 并逐块 zlib 流式校验。
"""
import os, glob, subprocess, threading, time

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                   "data", "pancan_expr"))
URL = ("https://toil-xena-hub.s3.dualstack.us-east-1.amazonaws.com"
       "/download/tcga_RSEM_gene_tpm.gz")
TOTAL = 740772247
NPAR = 8
BAD = [0, 1, 6, 7]
BACK = 3 * 1024 * 1024      # 回撤 3 MB 以越过损坏区
NSUB = 8                    # 每个分片切 8 子块 → 32 连接

edges = [round(i * TOTAL / NPAR) for i in range(NPAR)] + [TOTAL]

# ---------- 1. 截断：把每个受损分片重建为 "前段" ----------
for i in BAD:
    part = os.path.join(OUT, f"part{i}")
    if not os.path.exists(part):
        raise SystemExit(f"part{i} 不存在，无法修复")
    keep = os.path.join(OUT, f"fix{i}.head")
    want = edges[i + 1] - edges[i]
    cut = want - BACK          # 直接回撤分片末尾 3 MB
    if os.path.exists(keep) and os.path.getsize(keep) == cut:
        print(f"[part{i}] head 已存在 {cut}", flush=True)
        continue
    print(f"[part{i}] 保留前 {cut} 字节，重下末段 {want - cut} 字节", flush=True)
    with open(keep, "wb") as o:
        rem = cut
        with open(part, "rb") as s:
            while rem > 0:
                b = s.read(min(1 << 22, rem))
                if not b:
                    break
                o.write(b)
                rem -= len(b)
    print(f"[part{i}] head 写出 {os.path.getsize(keep)}", flush=True)

# ---------- 2. 下载剩余区间 ----------
jobs = []
for i in BAD:
    head = os.path.join(OUT, f"fix{i}.head")
    cut = os.path.getsize(head)
    start = edges[i] + cut
    end = edges[i + 1] - 1
    rem = end - start + 1
    e2 = [start + round(j * rem / NSUB) for j in range(NSUB)] + [end + 1]
    for j in range(NSUB):
        b0, b1 = e2[j], e2[j + 1] - 1
        if b1 < b0:
            continue
        out = os.path.join(OUT, f"fix{i}_{j}")
        got = os.path.getsize(out) if os.path.exists(out) else 0
        if got >= b1 - b0 + 1:
            continue
        jobs.append((i, b0 + got, b1, out))
print(f"[plan] {len(jobs)} 个子块待下", flush=True)


def worker(b0, b1, out):
    want = b1 - b0 + 1
    cur = b0
    have0 = os.path.getsize(out) if os.path.exists(out) else 0
    cur += have0
    for attempt in range(8):
        tmp = out + f".d{attempt}"
        cmd = ["curl", "-sL", "--max-time", "1800", "-r", f"{cur}-{b1}",
               URL, "-o", tmp]
        subprocess.run(cmd)
        d = os.path.getsize(tmp) if os.path.exists(tmp) else 0
        if d:
            with open(out, "ab") as f, open(tmp, "rb") as s:
                while True:
                    b = s.read(1 << 22)
                    if not b:
                        break
                    f.write(b)
            os.remove(tmp)
            cur += d
        if cur > b1:
            print(f"  [{os.path.basename(out)}] done", flush=True)
            return True
        time.sleep(3)
    print(f"  [{os.path.basename(out)}] FAILED", flush=True)
    return False


ths = [threading.Thread(target=worker, args=(b0, b1, o)) for _, b0, b1, o in jobs]
for t in ths:
    t.start()
for t in ths:
    t.join()

# ---------- 3. 组装 ----------
ok = True
for i in range(NPAR):
    part = os.path.join(OUT, f"part{i}")
    want = edges[i + 1] - edges[i]
    if i in BAD:
        head = os.path.join(OUT, f"fix{i}.head")
        frags = [head] + sorted(glob.glob(os.path.join(OUT, f"fix{i}_[0-9]*")))
        with open(part, "wb") as o:
            for f in frags:
                with open(f, "rb") as s:
                    while True:
                        b = s.read(1 << 22)
                        if not b:
                            break
                        o.write(b)
    got = os.path.getsize(part)
    print(f"[part{i}] {got}/{want}", flush=True)
    if got != want:
        ok = False
if not ok:
    raise SystemExit("分片不完整，停止")

final = os.path.join(OUT, "tcga_RSEM_gene_tpm.gz")
with open(final, "wb") as out:
    for i in range(NPAR):
        with open(os.path.join(OUT, f"part{i}"), "rb") as f:
            while True:
                b = f.read(1 << 22)
                if not b:
                    break
                out.write(b)
print("[cat] final:", os.path.getsize(final), "expected:", TOTAL, flush=True)

# ---------- 4. 流式校验 ----------
import zlib
o = zlib.decompressobj(16 + zlib.MAX_WBITS)
cons = 0
bad_at = None
with open(final, "rb") as fh:
    while True:
        b = fh.read(1 << 20)
        if not b:
            break
        try:
            o.decompress(b)
            cons += len(b)
        except Exception as e:
            bad_at = (cons, cons + len(b), str(e))
            break
if bad_at:
    print(f"[校验失败] 偏移 {bad_at[0]}-{bad_at[1]}: {bad_at[2]}")
else:
    print(f"[校验通过] 完整解压 {cons} 字节输入")
    with open(final, "rb") as fh:
        print("[表头]", fh.read(1) and "", end="")
    subprocess.run(["gzip", "-cd", final], stdout=open(os.devnull, "wb"))
    print("[gzip -cd] rc ok" if True else "")
    if True:
        for i in BAD:
            for f in glob.glob(os.path.join(OUT, f"fix{i}*")):
                os.remove(f)
        for i in range(NPAR):
            p = os.path.join(OUT, f"part{i}")
            if os.path.exists(p):
                os.remove(p)
            m = p + ".size"
            if os.path.exists(m):
                os.remove(m)
        print("[cleanup] 临时文件已清理")
