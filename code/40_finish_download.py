#!/usr/bin/env python3
"""40 — 收尾补丁：对 37 号未完成的分片，把剩余区间再细分并行抢完。

前提: part0..7 中部分已完成（partN >= 目标大小），未完成的以 partN 当前大小
      为起点，把剩余字节切成 SUB 块并行下载到 resN_M，全部完成后按序追加。
"""
import os, subprocess, threading, time, glob

OUTDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                      "data", "pancan_expr"))
URL = ("https://toil-xena-hub.s3.dualstack.us-east-1.amazonaws.com"
       "/download/tcga_RSEM_gene_tpm.gz")
TOTAL = 740772247
NPAR, SUB = 8, 4

edges = [round(i * TOTAL / NPAR) for i in range(NPAR)] + [TOTAL]
jobs = []           # (part_idx, b0, b1, outfile)
for i in range(NPAR):
    part = os.path.join(OUTDIR, f"part{i}")
    have = os.path.getsize(part) if os.path.exists(part) else 0
    want = edges[i + 1] - edges[i]
    if have >= want:
        print(f"[part{i}] complete", flush=True)
        continue
    # 清理旧 tail，重算实际起点
    for t in glob.glob(os.path.join(OUTDIR, f"part{i}.tail*")):
        have += os.path.getsize(t)
    start = edges[i] + have
    end = edges[i + 1] - 1
    if start > end:
        print(f"[part{i}] effective complete", flush=True)
        continue
    rem = end - start + 1
    e2 = [start + round(j * rem / SUB) for j in range(SUB)] + [end + 1]
    for j in range(SUB):
        b0, b1 = e2[j], e2[j + 1] - 1
        if b1 < b0:
            continue
        out = os.path.join(OUTDIR, f"res{i}_{j}")
        have_j = os.path.getsize(out) if os.path.exists(out) else 0
        if have_j >= b1 - b0 + 1:
            continue
        jobs.append((i, b0 + have_j, b1, out))
print(f"[plan] 剩余 {len(jobs)} 个子块", flush=True)


def worker(idx, b0, b1, out):
    want = b1 - b0 + 1
    for attempt in range(6):
        cmd = ["curl", "-sL", "--max-time", "3600", "-r", f"{b0}-{b1}",
               URL, "-o", out + f".d{attempt}"]
        r = subprocess.run(cmd)
        dsize = os.path.getsize(out + f".d{attempt}") if os.path.exists(out + f".d{attempt}") else 0
        if dsize >= want:
            os.rename(out + f".d{attempt}", out)
            print(f"[{out.split('/')[-1]}] done ({dsize})", flush=True)
            return True
        # 保存部分数据
        if dsize:
            with open(out, "ab") as f, open(out + f".d{attempt}", "rb") as s:
                f.write(s.read())
            os.remove(out + f".d{attempt}")
            b0 += dsize
        print(f"[{out.split('/')[-1]}] attempt {attempt+1}: got {dsize}/{want}, retry", flush=True)
        time.sleep(4)
    return False


ths = [threading.Thread(target=worker, args=j) for j in jobs]
for t in ths: t.start()
for t in ths: t.join()

# 合并：part i = part i + tail* + res i_*
ok = True
for i in range(NPAR):
    part = os.path.join(OUTDIR, f"part{i}")
    want = edges[i + 1] - edges[i]
    have = os.path.getsize(part)
    frags = sorted(glob.glob(os.path.join(OUTDIR, f"part{i}.tail*"))) + \
            sorted(glob.glob(os.path.join(OUTDIR, f"res{i}_*")))
    for f in frags:
        have += os.path.getsize(f)
    print(f"[part{i}] have {have}/{want}", flush=True)
    if have >= want:
        with open(part, "ab") as o:
            for f in frags:
                with open(f, "rb") as s:
                    while True:
                        b = s.read(1 << 22)
                        if not b: break
                        o.write(b)
        for f in frags:
            os.remove(f)
        print(f"[part{i}] merged -> {os.path.getsize(part)}", flush=True)
    else:
        ok = False

if ok:
    final = os.path.join(OUTDIR, "tcga_RSEM_gene_tpm.gz")
    with open(final, "wb") as out:
        for i in range(NPAR):
            with open(os.path.join(OUTDIR, f"part{i}"), "rb") as f:
                while True:
                    b = f.read(1 << 22)
                    if not b: break
                    out.write(b)
    print("[cat] final size:", os.path.getsize(final), flush=True)
    r = subprocess.run(["gzip", "-t", final])
    print("[gzip -t] rc =", r.returncode, flush=True)
    if r.returncode == 0:
        for i in range(NPAR):
            p = os.path.join(OUTDIR, f"part{i}")
            if os.path.exists(p): os.remove(p)
            m = p + ".size"
            if os.path.exists(m): os.remove(m)
        print("[cleanup] 分片已清理", flush=True)
