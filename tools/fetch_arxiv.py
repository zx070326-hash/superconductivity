"""
从 library/index/Dn.json 读取文献元数据，下载 arXiv PDF 到 library/pdf/Dn/，
并生成 library/index.md 总索引。

用法（在 superconductivity/ 目录下）：
    python tools/fetch_arxiv.py               # 下载 core + important
    python tools/fetch_arxiv.py --tier all    # 包括 reference
    python tools/fetch_arxiv.py --only D3     # 只处理 D3
    python tools/fetch_arxiv.py --dry-run     # 只列出不下载

JSON 条目格式（agent 产出）：
    {"title", "first_author", "authors_short", "year", "venue",
     "arxiv_id" | null, "doi" | null, "tier": "core|important|reference",
     "one_line", "why"}
"""
import argparse
import json
import re
import sys
import time
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
INDEX_DIR = ROOT / "library" / "index"
PDF_DIR = ROOT / "library" / "pdf"
INDEX_MD = ROOT / "library" / "index.md"

TIER_ORDER = {"core": 0, "important": 1, "reference": 2}
UA = "Mozilla/5.0 (research library sync; contact: local)"
SLEEP_S = 3.0  # arXiv 礼貌抓取间隔


def slug(s: str, n: int = 40) -> str:
    # 只保留 ASCII 字母数字，避免中文占位符（"作者未核实"）进文件名
    s = re.sub(r"[^A-Za-z0-9\-]+", "_", s).strip("_")
    return s[:n] or "unknown"


def pdf_name(e: dict) -> str:
    aid = (e.get("arxiv_id") or "noarxiv").replace("/", "-")
    return f"{e.get('year', 'xxxx')}_{slug(e.get('first_author', 'unknown'), 20)}_{aid}.pdf"


def arxiv_url(aid: str) -> str:
    aid = aid.strip()
    aid = re.sub(r"v\d+$", "", aid)  # 去版本号，取最新
    return f"https://arxiv.org/pdf/{aid}"


def download(url: str, dest: Path) -> bool:
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r:
            data = r.read()
    except Exception as ex:  # noqa: BLE001
        print(f"    FAIL {url}: {ex}")
        return False
    if not data.startswith(b"%PDF"):
        print(f"    FAIL {url}: not a PDF ({len(data)} bytes)")
        return False
    dest.write_bytes(data)
    return True


def load_indexes(only: str | None) -> dict[str, list[dict]]:
    out = {}
    for p in sorted(INDEX_DIR.glob("D*.json")):
        d = p.stem.split("_")[0]  # D3 或 D3_xxx → D3
        if only and d != only:
            continue
        with p.open(encoding="utf-8") as f:
            entries = json.load(f)
        for e in entries:
            e["_dir"] = d
            e["_file"] = pdf_name(e)
        out[d] = entries
    return out


def write_index_md(indexes: dict[str, list[dict]]):
    lines = ["# 文献总索引", "", f"自动生成于 {time.strftime('%Y-%m-%d')}，来源 `library/index/*.json`。", ""]
    for d, entries in indexes.items():
        lines += [f"## {d}", "", "| 级别 | 年份 | 第一作者 | 标题 | 出处 | arXiv | 本地 | 一句话 |", "|---|---|---|---|---|---|---|---|"]
        for e in sorted(entries, key=lambda x: (TIER_ORDER.get(x.get("tier"), 9), -int(x.get("year") or 0))):
            local = PDF_DIR / d / e["_file"]
            has = "✅" if local.exists() else ""
            aid = e.get("arxiv_id") or ""
            alink = f"[{aid}](https://arxiv.org/abs/{aid})" if aid else ""
            lines.append(
                f"| {e.get('tier','')} | {e.get('year','')} | {e.get('first_author','')} | {e.get('title','')} "
                f"| {e.get('venue','')} | {alink} | {has} | {e.get('one_line','')} |"
            )
        lines.append("")
    INDEX_MD.write_text("\n".join(lines), encoding="utf-8")
    print(f"wrote {INDEX_MD}")


def main():
    # Windows 控制台默认 GBK，打印中文占位符会崩
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    ap = argparse.ArgumentParser()
    ap.add_argument("--tier", default="important", choices=["core", "important", "all"],
                    help="下载到哪一级（含更高级）")
    ap.add_argument("--only", default=None, help="只处理某个方向，如 D3")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    max_tier = {"core": 0, "important": 1, "all": 2}[args.tier]
    indexes = load_indexes(args.only)
    if not indexes:
        print(f"no index json under {INDEX_DIR}")
        sys.exit(1)

    n_ok = n_skip = n_fail = n_noid = 0
    for d, entries in indexes.items():
        (PDF_DIR / d).mkdir(parents=True, exist_ok=True)
        print(f"== {d}: {len(entries)} entries")
        for e in entries:
            if TIER_ORDER.get(e.get("tier"), 9) > max_tier:
                continue
            dest = PDF_DIR / d / e["_file"]
            aid = e.get("arxiv_id")
            if not aid:
                n_noid += 1
                print(f"  [no arxiv] {e.get('year')} {e.get('first_author')}: {e.get('title')}")
                continue
            if dest.exists():
                n_skip += 1
                continue
            print(f"  GET {aid} -> {dest.name}")
            if args.dry_run:
                continue
            if download(arxiv_url(aid), dest):
                n_ok += 1
            else:
                n_fail += 1
            time.sleep(SLEEP_S)

    write_index_md(indexes)
    print(f"done: downloaded={n_ok} skipped(existing)={n_skip} failed={n_fail} no_arxiv={n_noid}")


if __name__ == "__main__":
    main()
