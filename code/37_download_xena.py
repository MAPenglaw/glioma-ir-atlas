#!/usr/bin/env python3
"""37 — 分块并行下载 Xena 泛癌 TCGA 表达矩阵 (log2 TPM)。

数据源: https://toil.xenahubs.net/download/tcga_RSEM_gene_tpm.gz
        10,535 TCGA 样本 × ~60k 基因, 约 740 MB gzip

策略: 8 个 Range 分片并行下载 -> 按序拼接 -> gzip -t 校验 -> 流式解压抽取目标基因行。
带块级续传: 每个 part 单独记录已完成字节, 失败只需重下该 part。
"""
import os, sys, subprocess, time, threading, queue

URL = ("https://toil-xena-hub.s3.dualstack.us-east-1.amazonaws.com"
       "/download/tcga_RSEM_gene_tpm.gz")
OUTDIR = os.path.join(os.path.dirname(__file__), "..", "data", "pancan_expr")
OUTDIR = os.path.abspath(OUTDIR)
TOTAL = 740772247
NPAR = 8
FINAL = os.path.join(OUTDIR, "tcga_RSEM_gene_tpm.gz")

os.makedirs(OUTDIR, exist_ok=True)

# ---- 分片区间 ----
edges = [round(i * TOTAL / NPAR) for i in range(NPAR)] + [TOTAL]
parts = [(edges[i], edges[i + 1] - 1, i) for i in range(NPAR)]


def size_of(p):
    return os.path.getsize(p) if os.path.exists(p) else 0


def part_bytes(i):
    """分片 i 的已下载字节数。

    必须包含正在写入的 .tailN 文件 —— 只统计 part{i} 会低估约 100%
    （curl 期间字节都在 .tail 里），这是 CGGA 下载时就踩过的坑。
    """
    import glob
    tot = size_of(os.path.join(OUTDIR, f"part{i}"))
    for f in glob.glob(os.path.join(OUTDIR, f"part{i}.tail*")):
        tot += size_of(f)
    return tot


def worker(idx, start, end, q):
    part = os.path.join(OUTDIR, f"part{idx}")
    meta = os.path.join(OUTDIR, f"part{idx}.size")
    want = end - start + 1
    for attempt in range(1, 6):
        have = size_of(part)
        if have >= want:
            q.put(f"[part{idx}] already complete ({have}/{want})")
            return True
        b0 = start + have
        if b0 > end:
            break
        q.put(f"[part{idx}] attempt {attempt}: bytes {b0}-{end} "
              f"(have {have}/{want})")
        cmd = ["curl", "-sL", "--max-time", "7200", "-r", f"{b0}-{end}",
               URL, "-o", part + f".tail{attempt}"]
        r = subprocess.run(cmd)
        if r.returncode != 0:
            q.put(f"[part{idx}] curl rc={r.returncode}, retry in 8s")
            time.sleep(8)
            continue
        # 追加到 part
        with open(part + f".tail{attempt}", "rb") as f, open(part, "ab") as o:
            while True:
                b = f.read(1 << 20)
                if not b:
                    break
                o.write(b)
        os.remove(part + f".tail{attempt}")
        have = size_of(part)
        open(meta, "w").write(str(have))
        if have >= want:
            q.put(f"[part{idx}] DONE ({have}/{want})")
            return True
        time.sleep(5)
    q.put(f"[part{idx}] FAILED after retries, have={size_of(part)}/{want}")
    return False


q = queue.Queue()


def logger():
    t0 = time.time()
    while True:
        try:
            m = q.get(timeout=2)
        except queue.Empty:
            done = sum(1 for i in range(NPAR)
                       if size_of(os.path.join(OUTDIR, f"part{i}"))
                       >= (edges[i + 1] - edges[i]))
            if done == NPAR:
                break
            tot = sum(part_bytes(i) for i in range(NPAR))
            el = time.time() - t0
            sp = tot / el / 1e6 if el > 0 else 0
            eta = (TOTAL - tot) / (tot / el) / 60 if tot > 0 else 0
            bar = "#" * int(30 * tot / TOTAL) + "." * (30 - int(30 * tot / TOTAL))
            print(f"[{el/60:6.1f}m] [{bar}] {tot/1e6:7.1f}/{TOTAL/1e6:.1f} MB "
                  f"({100*tot/TOTAL:5.1f}%) {sp:.2f} MB/s  ETA {eta:.0f}m "
                  f"done={done}/{NPAR}", flush=True)
            continue
        print(m, flush=True)


th = threading.Thread(target=logger, daemon=True)
th.start()

threads = [threading.Thread(target=worker, args=(i, s, e, q))
           for s, e, i in parts]
for t in threads:
    t.start()
for t in threads:
    t.join()
time.sleep(3)

ok = all(size_of(os.path.join(OUTDIR, f"part{i}")) >= (edges[i + 1] - edges[i])
         for i in range(NPAR))
print("[status] parts ok =", ok)
if ok:
    with open(FINAL, "wb") as out:
        for i in range(NPAR):
            with open(os.path.join(OUTDIR, f"part{i}"), "rb") as f:
                while True:
                    b = f.read(1 << 22)
                    if not b:
                        break
                    out.write(b)
    print("[cat] final size:", os.path.getsize(FINAL))
    r = subprocess.run(["gzip", "-t", FINAL])
    print("[gzip -t] rc =", r.returncode)
