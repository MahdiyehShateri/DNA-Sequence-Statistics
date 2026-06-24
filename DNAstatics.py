#!/usr/bin/env python3
"""
dna_stats.py — DNA Sequence Statistics Analyzer
=================================================
A modular, production-quality bioinformatics script for analyzing DNA sequences.

Features
--------
* Basic statistics    : length, nucleotide counts, GC/AT content, GC skew
* Codon usage         : frequency table for all 64 codons with amino acid mapping
* k-mer frequencies   : dinucleotide (k=2) and trinucleotide (k=3) tables
* Motif search        : all occurrences (overlapping) with 0-based and 1-based positions
* Sequence complexity : Shannon entropy, linguistic complexity, sliding-window GC plot

Installation
------------
    pip install biopython matplotlib rich

Usage examples
--------------
# Analyze sequences from a FASTA file
    python dna_stats.py --fasta sequences.fasta

# Analyze a single inline sequence
    python dna_stats.py --sequence ATGCGTACGATCGATCG

# Search for motifs while analyzing a FASTA file
    python dna_stats.py --fasta sequences.fasta --motifs ATG TATAAT GAATTC

# Save all plots to a custom directory
    python dna_stats.py --fasta sequences.fasta --outdir ./results

Author  : Bioinformatics Learning Project
License : MIT
"""

import argparse
import math
import re
import sys
from collections import Counter
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# ── Third-party imports ────────────────────────────────────────────────────────
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqUtils import gc_fraction
except ImportError:
    sys.exit(
        "Biopython is required.  Install it with:  pip install biopython"
    )

try:
    import matplotlib
    matplotlib.use("Agg")          # non-interactive backend; works without a display
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("[WARNING] matplotlib not found — plots will be skipped.  "
          "Install with:  pip install matplotlib")

try:
    from rich.console import Console
    from rich.table import Table
    from rich import box
    HAS_RICH = True
    console = Console()
except ImportError:
    HAS_RICH = False
    console = None
    print("[WARNING] rich not found — plain-text output will be used.  "
          "Install with:  pip install rich")

# ── Standard genetic code (codon → single-letter amino acid) ──────────────────
GENETIC_CODE: Dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

# ─────────────────────────────────────────────────────────────────────────────
# 1. BASIC STATISTICS
# ─────────────────────────────────────────────────────────────────────────────

def compute_basic_stats(sequence: str) -> Dict:
    """
    Compute fundamental statistics for a DNA sequence.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence string.

    Returns
    -------
    dict
        Keys: length, counts (A/T/C/G/N/other), gc_content, at_content, gc_skew.
    """
    length = len(sequence)
    if length == 0:
        raise ValueError("Sequence is empty.")

    # Raw nucleotide counts
    counter = Counter(sequence)
    a = counter.get("A", 0)
    t = counter.get("T", 0)
    c = counter.get("C", 0)
    g = counter.get("G", 0)
    n = counter.get("N", 0)
    other = length - a - t - c - g - n

    gc = g + c
    at = a + t
    valid = gc + at           # exclude N and ambiguous bases for percentages

    gc_content = (gc / valid * 100) if valid else 0.0
    at_content = (at / valid * 100) if valid else 0.0

    # GC skew: (G − C) / (G + C); positive = more G than C
    gc_skew = ((g - c) / gc) if gc else 0.0

    return {
        "length": length,
        "counts": {"A": a, "T": t, "C": c, "G": g, "N": n, "other": other},
        "gc_content": gc_content,
        "at_content": at_content,
        "gc_skew": gc_skew,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 2. CODON USAGE
# ─────────────────────────────────────────────────────────────────────────────

def compute_codon_usage(sequence: str) -> Dict[str, Dict]:
    """
    Count all non-overlapping codons starting at position 0 (reading frame +1).

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence.

    Returns
    -------
    dict
        Mapping codon → {"amino_acid": str, "count": int, "frequency": float}.
    """
    # Initialise all 64 codons with zero counts
    all_codons = ["".join(c) for c in product("ACGT", repeat=3)]
    codon_counts: Dict[str, int] = {cod: 0 for cod in all_codons}

    # Walk the sequence in triplets (frame +1)
    triplets = [sequence[i: i + 3] for i in range(0, len(sequence) - 2, 3)]
    valid_triplets = [t for t in triplets if set(t).issubset(set("ACGT"))]

    for triplet in valid_triplets:
        codon_counts[triplet] += 1

    total = sum(codon_counts.values())

    usage = {}
    for codon, count in codon_counts.items():
        usage[codon] = {
            "amino_acid": GENETIC_CODE.get(codon, "?"),
            "count": count,
            "frequency": (count / total * 100) if total else 0.0,
        }
    return usage


# ─────────────────────────────────────────────────────────────────────────────
# 3. k-MER FREQUENCIES
# ─────────────────────────────────────────────────────────────────────────────

def compute_kmer_frequencies(sequence: str, k: int) -> Dict[str, float]:
    """
    Calculate the relative frequency (%) of every possible k-mer.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence.
    k : int
        k-mer length (2 for dinucleotides, 3 for trinucleotides, …).

    Returns
    -------
    dict
        kmer → frequency (%).  Every possible k-mer is present (zero if absent).
    """
    # All possible k-mers over {A, C, G, T}
    all_kmers = ["".join(c) for c in product("ACGT", repeat=k)]
    counts: Dict[str, int] = {km: 0 for km in all_kmers}

    for i in range(len(sequence) - k + 1):
        kmer = sequence[i: i + k]
        if kmer in counts:          # skip k-mers containing N or other chars
            counts[kmer] += 1

    total = sum(counts.values())
    return {km: (cnt / total * 100) if total else 0.0 for km, cnt in counts.items()}


# ─────────────────────────────────────────────────────────────────────────────
# 4. MOTIF SEARCH
# ─────────────────────────────────────────────────────────────────────────────

def search_motifs(sequence: str, motifs: List[str]) -> Dict[str, List[int]]:
    """
    Find all (overlapping) occurrences of each motif in the sequence.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence.
    motifs : list of str
        List of motif strings to search for (case-insensitive).

    Returns
    -------
    dict
        motif → list of 0-based start positions.
    """
    results: Dict[str, List[int]] = {}
    for motif in motifs:
        pattern = motif.upper()
        # re.finditer with a lookahead captures overlapping matches
        positions = [m.start() for m in re.finditer(f"(?={re.escape(pattern)})", sequence)]
        results[pattern] = positions
    return results


# ─────────────────────────────────────────────────────────────────────────────
# 5. SEQUENCE COMPLEXITY
# ─────────────────────────────────────────────────────────────────────────────

def compute_shannon_entropy(sequence: str) -> float:
    """
    Compute the Shannon entropy (bits) of a DNA sequence.

    H = -Σ p(x) · log₂(p(x))   for each unique character x

    Maximum entropy for a 4-symbol alphabet is 2.0 bits (perfectly random).

    Parameters
    ----------
    sequence : str
        DNA sequence string.

    Returns
    -------
    float
        Entropy value in bits (0–2).
    """
    if not sequence:
        return 0.0
    counter = Counter(sequence)
    length = len(sequence)
    entropy = -sum(
        (freq / length) * math.log2(freq / length)
        for freq in counter.values()
        if freq > 0
    )
    return entropy


def compute_linguistic_complexity(sequence: str) -> float:
    """
    Estimate linguistic complexity as the ratio of observed unique k-mers
    to the maximum possible unique k-mers, averaged over k = 1 … min(|seq|, 4).

    A value close to 1.0 indicates a highly complex (low-repeat) sequence.
    A value close to 0.0 suggests a highly repetitive sequence.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence.

    Returns
    -------
    float
        Linguistic complexity score in [0, 1].
    """
    if len(sequence) == 0:
        return 0.0

    scores = []
    for k in range(1, min(len(sequence), 5)):           # k = 1, 2, 3, 4
        observed = len({sequence[i: i + k] for i in range(len(sequence) - k + 1)})
        max_possible = min(4**k, len(sequence) - k + 1)
        if max_possible > 0:
            scores.append(observed / max_possible)

    return sum(scores) / len(scores) if scores else 0.0


def compute_sliding_gc(sequence: str, window: int = 100, step: int = 10) -> Tuple[List[int], List[float]]:
    """
    Calculate GC content in a sliding window along the sequence.

    Parameters
    ----------
    sequence : str
        Uppercase DNA sequence.
    window : int
        Window size in base pairs (default 100).
    step : int
        Step size between consecutive windows (default 10).

    Returns
    -------
    tuple
        (positions, gc_values) — mid-point positions and GC% for each window.
    """
    positions, gc_values = [], []
    for start in range(0, len(sequence) - window + 1, step):
        chunk = sequence[start: start + window]
        gc = (chunk.count("G") + chunk.count("C"))
        valid = sum(1 for b in chunk if b in "ACGT")
        gc_pct = (gc / valid * 100) if valid else 0.0
        positions.append(start + window // 2)           # mid-point of window
        gc_values.append(gc_pct)
    return positions, gc_values


# ─────────────────────────────────────────────────────────────────────────────
# PLOTTING FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def plot_nucleotide_composition(counts: Dict[str, int], seq_id: str, outdir: Path) -> None:
    """Bar chart of nucleotide composition (A, T, C, G, N)."""
    if not HAS_MATPLOTLIB:
        return

    bases = ["A", "T", "C", "G", "N"]
    values = [counts.get(b, 0) for b in bases]
    colors = ["#4CAF50", "#F44336", "#2196F3", "#FF9800", "#9E9E9E"]

    fig, ax = plt.subplots(figsize=(7, 4))
    bars = ax.bar(bases, values, color=colors, edgecolor="black", linewidth=0.6)

    # Annotate each bar with its count
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + max(values) * 0.01,
                str(val), ha="center", va="bottom", fontsize=9)

    ax.set_title(f"Nucleotide Composition — {seq_id}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Nucleotide")
    ax.set_ylabel("Count")
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda x, _: f"{int(x):,}"))
    plt.tight_layout()

    out = outdir / f"{_safe_filename(seq_id)}_nucleotide_composition.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    _print(f"  [Plot saved] {out}")


def plot_gc_window(positions: List[int], gc_values: List[float],
                   seq_id: str, outdir: Path, window: int) -> None:
    """Line plot of GC content in sliding windows."""
    if not HAS_MATPLOTLIB or not positions:
        return

    fig, ax = plt.subplots(figsize=(10, 3.5))
    ax.plot(positions, gc_values, color="#1565C0", linewidth=0.8, alpha=0.85)
    ax.axhline(sum(gc_values) / len(gc_values), color="#E53935",
               linestyle="--", linewidth=1, label="Mean GC%")
    ax.fill_between(positions, gc_values, alpha=0.15, color="#1565C0")

    ax.set_title(f"Sliding-Window GC Content ({window} bp) — {seq_id}",
                 fontsize=12, fontweight="bold")
    ax.set_xlabel("Position (bp)")
    ax.set_ylabel("GC Content (%)")
    ax.set_ylim(0, 100)
    ax.legend(fontsize=9)
    plt.tight_layout()

    out = outdir / f"{_safe_filename(seq_id)}_gc_window.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    _print(f"  [Plot saved] {out}")


def plot_codon_usage(codon_usage: Dict[str, Dict], seq_id: str, outdir: Path, top_n: int = 20) -> None:
    """Horizontal bar chart of the top-N most-used codons."""
    if not HAS_MATPLOTLIB:
        return

    # Sort by count descending, take top N
    sorted_codons = sorted(codon_usage.items(), key=lambda x: x[1]["count"], reverse=True)[:top_n]
    labels = [f"{cod}\n({info['amino_acid']})" for cod, info in sorted_codons]
    values = [info["count"] for _, info in sorted_codons]

    fig, ax = plt.subplots(figsize=(8, max(4, top_n * 0.35)))
    bars = ax.barh(labels[::-1], values[::-1], color="#5C6BC0", edgecolor="black", linewidth=0.4)

    for bar, val in zip(bars, values[::-1]):
        ax.text(bar.get_width() + max(values) * 0.005, bar.get_y() + bar.get_height() / 2,
                str(val), va="center", fontsize=7)

    ax.set_title(f"Top {top_n} Codon Usage — {seq_id}", fontsize=12, fontweight="bold")
    ax.set_xlabel("Count")
    plt.tight_layout()

    out = outdir / f"{_safe_filename(seq_id)}_codon_usage.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    _print(f"  [Plot saved] {out}")


# ─────────────────────────────────────────────────────────────────────────────
# REPORT PRINTING HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def _print(msg: str, style: str = "") -> None:
    """Unified print that uses Rich if available, otherwise plain print."""
    if HAS_RICH:
        console.print(msg, style=style)
    else:
        # Strip any Rich markup before printing
        clean = re.sub(r"\[/?[a-zA-Z0-9_ ]+\]", "", msg)
        print(clean)


def _safe_filename(name: str) -> str:
    """Convert a sequence ID to a filesystem-safe string."""
    return re.sub(r"[^\w\-]", "_", name)[:60]


def _rich_table(title: str, columns: List[str], rows: List[List[str]]) -> None:
    """Print a formatted table using Rich or plain-text fallback."""
    if HAS_RICH:
        tbl = Table(title=title, box=box.SIMPLE_HEAD, show_lines=False,
                    title_style="bold cyan")
        for col in columns:
            tbl.add_column(col, style="white", no_wrap=True)
        for row in rows:
            tbl.add_row(*row)
        console.print(tbl)
    else:
        # Plain-text table
        col_widths = [max(len(col), max((len(r[i]) for r in rows), default=0))
                      for i, col in enumerate(columns)]
        sep = "  ".join("-" * w for w in col_widths)
        header = "  ".join(col.ljust(col_widths[i]) for i, col in enumerate(columns))
        print(f"\n{title}")
        print(header)
        print(sep)
        for row in rows:
            print("  ".join(cell.ljust(col_widths[i]) for i, cell in enumerate(row)))


# ─────────────────────────────────────────────────────────────────────────────
# REPORT FUNCTIONS
# ─────────────────────────────────────────────────────────────────────────────

def report_basic_stats(stats: Dict, seq_id: str) -> None:
    """Pretty-print basic sequence statistics."""
    _print(f"\n[bold green]━━━ Basic Statistics: {seq_id} ━━━[/bold green]")
    counts = stats["counts"]
    rows = [
        ["Sequence Length", f"{stats['length']:,} bp"],
        ["A", f"{counts['A']:,}  ({counts['A']/stats['length']*100:.2f}%)"],
        ["T", f"{counts['T']:,}  ({counts['T']/stats['length']*100:.2f}%)"],
        ["C", f"{counts['C']:,}  ({counts['C']/stats['length']*100:.2f}%)"],
        ["G", f"{counts['G']:,}  ({counts['G']/stats['length']*100:.2f}%)"],
        ["N / ambiguous", f"{counts['N']:,}"],
        ["Other", f"{counts['other']:,}"],
        ["GC Content", f"{stats['gc_content']:.2f}%"],
        ["AT Content", f"{stats['at_content']:.2f}%"],
        ["GC Skew (G−C)/(G+C)", f"{stats['gc_skew']:.4f}"],
    ]
    _rich_table("", ["Metric", "Value"], [[r[0], r[1]] for r in rows])


def report_codon_usage(codon_usage: Dict[str, Dict], seq_id: str, top_n: int = 10) -> None:
    """Print a sorted codon usage table (top-N and bottom-N by count)."""
    _print(f"\n[bold green]━━━ Codon Usage: {seq_id} ━━━[/bold green]")

    sorted_by_count = sorted(codon_usage.items(), key=lambda x: x[1]["count"], reverse=True)

    # Top N
    rows = [
        [cod, info["amino_acid"], str(info["count"]), f"{info['frequency']:.3f}%"]
        for cod, info in sorted_by_count[:top_n]
    ]
    _rich_table(f"Top {top_n} Most-Used Codons", ["Codon", "AA", "Count", "Frequency"], rows)

    # Bottom N (exclude zero-count codons to keep output meaningful)
    nonzero = [(c, i) for c, i in sorted_by_count if i["count"] > 0]
    if nonzero:
        least = nonzero[-top_n:][::-1]
        rows_bot = [
            [cod, info["amino_acid"], str(info["count"]), f"{info['frequency']:.3f}%"]
            for cod, info in least
        ]
        _rich_table(f"Bottom {top_n} Least-Used Codons (non-zero)", ["Codon", "AA", "Count", "Frequency"], rows_bot)


def report_kmer_frequencies(frequencies: Dict[str, float], k: int, seq_id: str, top_n: int = 10) -> None:
    """Print top-N and bottom-N k-mer frequencies."""
    label = {2: "Dinucleotide", 3: "Trinucleotide"}.get(k, f"{k}-mer")
    _print(f"\n[bold green]━━━ {label} Frequencies: {seq_id} ━━━[/bold green]")

    sorted_freq = sorted(frequencies.items(), key=lambda x: x[1], reverse=True)

    top_rows = [[km, f"{freq:.4f}%"] for km, freq in sorted_freq[:top_n]]
    _rich_table(f"Top {top_n} {label}s", ["k-mer", "Frequency"], top_rows)

    bot_rows = [[km, f"{freq:.4f}%"] for km, freq in sorted_freq[-top_n:][::-1]]
    _rich_table(f"Bottom {top_n} {label}s", ["k-mer", "Frequency"], bot_rows)


def report_motif_search(results: Dict[str, List[int]], seq_id: str) -> None:
    """Print motif search results with 0-based and 1-based positions."""
    _print(f"\n[bold green]━━━ Motif Search: {seq_id} ━━━[/bold green]")
    if not results:
        _print("  No motifs specified.")
        return

    for motif, positions in results.items():
        count = len(positions)
        _print(f"\n  [cyan]Motif:[/cyan] {motif}  — {count} occurrence(s)")
        if count == 0:
            _print("    Not found in sequence.")
        else:
            rows = []
            for pos in positions[:50]:          # cap display at first 50 hits
                rows.append([str(pos), str(pos + 1), str(pos + len(motif) - 1)])
            if len(positions) > 50:
                rows.append(["…", "…", f"(+{len(positions)-50} more)"])
            _rich_table("", ["0-based start", "1-based start", "0-based end"], rows)


def report_complexity(entropy: float, lc: float, seq_id: str) -> None:
    """Print sequence complexity metrics."""
    _print(f"\n[bold green]━━━ Sequence Complexity: {seq_id} ━━━[/bold green]")
    rows = [
        ["Shannon Entropy (bits)", f"{entropy:.4f}", "Max = 2.0 (uniform ACGT)"],
        ["Linguistic Complexity", f"{lc:.4f}", "1.0 = maximally complex"],
    ]
    _rich_table("", ["Metric", "Value", "Note"], rows)


# ─────────────────────────────────────────────────────────────────────────────
# SEQUENCE LOADING
# ─────────────────────────────────────────────────────────────────────────────

def load_sequences(fasta: Optional[str], sequence: Optional[str]) -> List[Tuple[str, str]]:
    """
    Load sequences from a FASTA file or an inline string.

    Parameters
    ----------
    fasta : str or None
        Path to a FASTA file.
    sequence : str or None
        Raw DNA sequence string (used when --sequence is given).

    Returns
    -------
    list of (seq_id, uppercase_sequence)
    """
    seqs: List[Tuple[str, str]] = []

    if fasta:
        path = Path(fasta)
        if not path.exists():
            sys.exit(f"[ERROR] FASTA file not found: {fasta}")
        try:
            for record in SeqIO.parse(str(path), "fasta"):
                seqs.append((record.id, str(record.seq).upper()))
        except Exception as exc:
            sys.exit(f"[ERROR] Failed to parse FASTA file: {exc}")

    if sequence:
        clean = re.sub(r"\s+", "", sequence).upper()
        seqs.append(("inline_sequence", clean))

    if not seqs:
        sys.exit("[ERROR] No sequences provided.  Use --fasta or --sequence.")

    return seqs


# ─────────────────────────────────────────────────────────────────────────────
# MAIN PIPELINE
# ─────────────────────────────────────────────────────────────────────────────

def analyse_sequence(
    seq_id: str,
    sequence: str,
    motifs: Optional[List[str]],
    outdir: Path,
    window: int = 100,
) -> None:
    """
    Run the full analysis pipeline for a single sequence.

    Parameters
    ----------
    seq_id : str
        Identifier for the sequence (used in titles and filenames).
    sequence : str
        Uppercase DNA sequence.
    motifs : list of str or None
        Motifs to search for.
    outdir : Path
        Directory where plot images are saved.
    window : int
        Sliding-window size for GC content plot.
    """
    _print(f"\n[bold magenta]{'═'*60}[/bold magenta]")
    _print(f"[bold magenta]  Analysing: {seq_id}  ({len(sequence):,} bp)[/bold magenta]")
    _print(f"[bold magenta]{'═'*60}[/bold magenta]")

    # ── 1. Basic statistics ──────────────────────────────────────────────────
    stats = compute_basic_stats(sequence)
    report_basic_stats(stats, seq_id)
    plot_nucleotide_composition(stats["counts"], seq_id, outdir)

    # ── 2. Codon usage ───────────────────────────────────────────────────────
    codon_usage = compute_codon_usage(sequence)
    report_codon_usage(codon_usage, seq_id)
    plot_codon_usage(codon_usage, seq_id, outdir)

    # ── 3. k-mer frequencies ─────────────────────────────────────────────────
    di_freq = compute_kmer_frequencies(sequence, k=2)
    report_kmer_frequencies(di_freq, k=2, seq_id=seq_id)

    tri_freq = compute_kmer_frequencies(sequence, k=3)
    report_kmer_frequencies(tri_freq, k=3, seq_id=seq_id)

    # ── 4. Motif search ───────────────────────────────────────────────────────
    if motifs:
        motif_results = search_motifs(sequence, motifs)
        report_motif_search(motif_results, seq_id)
    else:
        _print("\n[dim]Tip: use --motifs ATG TATAAT GAATTC to search for motifs.[/dim]")

    # ── 5. Sequence complexity ────────────────────────────────────────────────
    entropy = compute_shannon_entropy(sequence)
    lc = compute_linguistic_complexity(sequence)
    report_complexity(entropy, lc, seq_id)

    # ── Sliding-window GC plot ────────────────────────────────────────────────
    effective_window = min(window, len(sequence))
    positions, gc_values = compute_sliding_gc(sequence, window=effective_window)
    plot_gc_window(positions, gc_values, seq_id, outdir, effective_window)


# ─────────────────────────────────────────────────────────────────────────────
# ARGUMENT PARSING
# ─────────────────────────────────────────────────────────────────────────────

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="dna_stats.py",
        description="DNA Sequence Statistics Analyzer using Biopython",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples
--------
  python dna_stats.py --sequence ATGCGTACGATCGATCG
  python dna_stats.py --fasta genome.fasta --motifs ATG GAATTC
  python dna_stats.py --fasta sequences.fasta --outdir ./results --window 200
        """,
    )
    src = parser.add_mutually_exclusive_group()
    src.add_argument("--fasta", metavar="FILE",
                     help="Path to a FASTA file (one or more sequences).")
    src.add_argument("--sequence", metavar="SEQ",
                     help="Single DNA sequence string (uppercase or lowercase).")
    parser.add_argument("--motifs", nargs="+", metavar="MOTIF",
                        help="One or more motifs to search for (e.g. ATG TATAAT).")
    parser.add_argument("--outdir", default=".", metavar="DIR",
                        help="Output directory for plots (default: current directory).")
    parser.add_argument("--window", type=int, default=100, metavar="INT",
                        help="Sliding-window size for GC content plot (default: 100).")
    # Allow running with no args to show help
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)
    return parser.parse_args()


# ─────────────────────────────────────────────────────────────────────────────
# ENTRY POINT
# ─────────────────────────────────────────────────────────────────────────────

def main() -> None:
    args = parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    sequences = load_sequences(args.fasta, args.sequence)
    _print(f"\n[bold]Loaded {len(sequences)} sequence(s).[/bold]")

    for seq_id, seq in sequences:
        try:
            analyse_sequence(
                seq_id=seq_id,
                sequence=seq,
                motifs=args.motifs,
                outdir=outdir,
                window=args.window,
            )
        except ValueError as exc:
            _print(f"[red][SKIP] {seq_id}: {exc}[/red]")

    _print(f"\n[bold green]✔  Analysis complete.  Plots saved to: {outdir.resolve()}[/bold green]\n")


if __name__ == "__main__":
    main()


# ─────────────────────────────────────────────────────────────────────────────
# EXAMPLE USAGE (run directly for a quick demo)
# ─────────────────────────────────────────────────────────────────────────────
# To run a quick built-in demo without any external files, execute:
#
#   python - <<'EOF'
#   from dna_stats import *
#   from pathlib import Path
#
#   demo_seq = (
#       "ATGAAACCCGGGTTTTAA"           # Met-Lys-Pro-Gly-Phe-Stop
#       "GCTAGCTAGCTAGCTAGC"           # repeat region
#       "GAATTCAAGCTTGGATCC"           # EcoRI + HindIII + BamHI sites
#       "TATAATAGCGATCGATCG" * 5       # TATA-box region × 5
#   )
#
#   stats   = compute_basic_stats(demo_seq)
#   report_basic_stats(stats, "demo")
#
#   usage   = compute_codon_usage(demo_seq)
#   report_codon_usage(usage, "demo", top_n=5)
#
#   di      = compute_kmer_frequencies(demo_seq, k=2)
#   report_kmer_frequencies(di, k=2, seq_id="demo", top_n=5)
#
#   motifs  = search_motifs(demo_seq, ["ATG", "GAATTC", "TATAAT"])
#   report_motif_search(motifs, "demo")
#
#   entropy = compute_shannon_entropy(demo_seq)
#   lc      = compute_linguistic_complexity(demo_seq)
#   report_complexity(entropy, lc, "demo")
#   EOF
