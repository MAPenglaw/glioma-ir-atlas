#!/usr/bin/env python3
"""42 — gzip 流式自动迭代修复。

循环：组装 → 流式 zlib 校验 → 定位第一个坏点 → 判断所属分片 →
回退 1 MB 重下该分片尾段 → 重新组装校验，最多 12 轮。
"""
import os, glob, subprocess, threading, time, zlib

OUT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..",
                                   "data", "pancan_expr"))
URL = ("https://toil-xena-hub.s3.dualstack.us-east-1.amazonaws.com"
       "/download/tcga_RSEM_gene_tpm.gz")
TOTAL = 740772247
NPAR = 8
NSUB = 16
BACK = 1 * 1024 * 1024
edges = [round(i * TOTAL / NPAR) for i in range(NPAR)] + [TOTAL]
FINAL = os.path.join(OUT, "tcga_RSEM_gene_tpm.gz")


def assemble():
    with open(FINAL, "wb") as out:
        for i in range(NPAR):
            with open(os.path.join(OUT, f"part{i}"), "rb") as f:
                while True:
                    b = f.read(1 << 22)
                    if not b:
                        break
                    out.write(b)


def verify():
    """返回 (ok, bad_offset)"""
    o = zlib.decompressobj(16 + zlib.MAX_WBITS)
    cons = 0
    with open(FINAL, "rb") as fh:
        while True:
            b = fh.read(1 << 20)
            if not b:
                break
            try:
                o.decompress(b)
                cons += len(b)
            except Exception as e:
                return False, cons, str(e)
    return True, cons, ""


def dl(b0, b1, out):
    cur = b0
    if os.path.exists(out):
        cur += os.path.getsize(out)
    for a in range(8):
        if cur > b1:
            return True
        tmp = f"{out}.d{a}"
        subprocess.run(["curl", "-sL", "--max-time", "1800",
                        "-r", f"{cur}-{b1}", URL, "-o", tmp])
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
        time.sleep(2)
    return cur > b1


def repair_part(i, keep):
    """保留 part{i} 前 keep 字节，重下 [edges[i]+keep, edges[i+1]-1]"""
    part = os.path.join(OUT, f"part{i}")
    start = edges[i] + keep
    end = edges[i + 1] - 1
    rem = end - start + 1
    tag = f"rp{i}_{keep}"
    outs = []
    e2 = [start + round(j * rem / NSUB) for j in range(NSUB)] + [end + 1]
    jobs = []
    for j in range(NSUB):
        b0, b1 = e2[j], e2[j + 1] - 1
        if b1 < b0:
            continue
        out = os.path.join(OUT, f"{tag}_{j}")
        outs.append(out)
        if os.path.exists(out) and os.path.getsize(out) >= b1 - b0 + 1:
            continue
        jobs.append((b0, b1, out))
    print(f"  重下 part{i} 尾段 {rem/1e6:.1f} MB（{len(jobs)} 块）", flush=True)
    ths = [threading.Thread(target=dl, args=j) for j in jobs]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    return outs


tried = []
for rnd in range(12):
    assemble()
    ok, cons, err = verify()
    if ok:
        print(f"[round {rnd}] 校验通过，解压输入 {cons} 字节", flush=True)
        break
    print(f"[round {rnd}] 坏点 offset={cons} ({err})", flush=True)
    # 定位分片
    i = min(range(NPAR), key=lambda k: abs(edges[k] - cons) if edges[k] <= cons else 1e18)
    if edges[i] > cons:
        i -= 1
    local = cons - edges[i]
    keep = max(0, local - BACK)
    # 重复坏点 → 加大回撤
    while (i, keep) in tried:
        keep = max(0, keep - 5 * 1024 * 1024)
        print(f"  重复坏点，加大回撤至 {keep}", flush=True)
        if keep == 0:
            break
    print(f"  -> 属于 part{i}，分片内偏移 {local}，保留前 {keep}", flush=True)
    # 备份旧 part 的前 keep 字节
    part = os.path.join(OUT, f"part{i}")
    head = os.path.join(OUT, f"part{i}.keep")
    with open(head, "wb") as o:
        rem = keep
        with open(part, "rb") as s:
            while rem > 0:
                b = s.read(min(1 << 22, rem))
                if not b:
                    break
                o.write(b)
                rem -= len(b)
    outs = repair_part(i, keep)
    # 合并 head + 新块
    with open(part, "wb") as o:
        with open(head, "rb") as s:
            while True:
                b = s.read(1 << 22)
                if not b:
                    break
                o.write(b)
        for f in sorted(outs, key=lambda x: int(x.rsplit("_", 1)[1])):
            if not os.path.exists(f):
                print(f"  [warn] 缺块 {os.path.basename(f)}", flush=True)
                continue
            with open(f, "rb") as s:
                while True:
                    b = s.read(1 << 22)
                    if not b:
                        break
                    o.write(b)
    got = os.path.getsize(part)
    want = edges[i + 1] - edges[i]
    print(f"  part{i} 重建 -> {got}/{want}", flush=True)
    if got != want:
        print(f"  [error] part{i} 长度不符，本轮中止重下", flush=True)
        raise SystemExit(1)
    tried.append((i, keep))

# 最终
assemble()
ok, cons, err = verify()
print("[final]", "OK" if ok else f"FAILED at {cons}: {err}")
if ok:
    r = subprocess.run(["gzip", "-t", FINAL])
    print("gzip -t rc =", r.returncode)
    if r.returncode == 0:
        for i in range(NPAR):
            p = os.path.join(OUT, f"part{i}")
            if os.path.exists(p):
                os.remove(p)
            m = p + ".size"
            if os.path.exists(m):
                os.remove(m)
        for f in glob.glob(os.path.join(OUT, "rp*")) + glob.glob(os.path.join(OUT, "fix*")):
            os.remove(f)
        print("[cleanup] done")
        print("[final size]", os.path.getsize(FINAL))
