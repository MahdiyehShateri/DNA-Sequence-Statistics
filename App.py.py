import streamlit as st
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from collections import Counter
import plotly.express as px
import plotly.graph_objects as go
import math
import re
from itertools import product
import base64
from io import BytesIO
import zipfile

st.set_page_config(
    page_title="DNA Stats Pro",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Dark mode by default
st.markdown("""
<style>
    .stApp {
        background-color: #0E1117;
        color: #FAFAFA;
    }
</style>
""", unsafe_allow_html=True)

# ====================== HELPER FUNCTIONS ======================

def compute_basic_stats(seq: str):
    seq = seq.upper()
    length = len(seq)
    counter = Counter(seq)
    a = counter.get("A", 0)
    t = counter.get("T", 0)
    c = counter.get("C", 0)
    g = counter.get("G", 0)
    n = counter.get("N", 0)
    other = length - (a + t + c + g + n)
    
    gc = g + c
    at = a + t
    valid = gc + at
    gc_content = (gc / valid * 100) if valid else 0
    at_content = (at / valid * 100) if valid else 0
    gc_skew = ((g - c) / gc) if gc else 0
    
    return {
        "length": length,
        "gc_content": round(gc_content, 2),
        "at_content": round(at_content, 2),
        "gc_skew": round(gc_skew, 4),
        "counts": {"A": a, "T": t, "C": c, "G": g, "N": n, "Other": other}
    }

def compute_codon_usage(seq: str):
    seq = seq.upper()
    codons = [seq[i:i+3] for i in range(0, len(seq)-2, 3)]
    counter = Counter([c for c in codons if len(c) == 3 and set(c).issubset("ACGT")])
    
    genetic_code = {
        "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
        "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
        "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
        "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
        "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
        "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
        "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
        "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G",
    }
    
    total = sum(counter.values())
    data = []
    for codon, count in sorted(counter.items(), key=lambda x: x[1], reverse=True):
        aa = genetic_code.get(codon, "?")
        freq = (count / total * 100) if total else 0
        data.append([codon, aa, count, round(freq, 3)])
    return pd.DataFrame(data, columns=["Codon", "Amino Acid", "Count", "Frequency (%)"])

def compute_kmer_freq(seq: str, k=3):
    seq = seq.upper()
    kmers = [seq[i:i+k] for i in range(len(seq)-k+1)]
    counter = Counter([kmer for kmer in kmers if set(kmer).issubset("ACGT")])
    total = sum(counter.values())
    df = pd.DataFrame(counter.items(), columns=["k-mer", "Count"])
    df["Frequency (%)"] = round(df["Count"] / total * 100, 4)
    return df.sort_values("Frequency (%)", ascending=False)

def search_motifs(seq: str, motifs):
    seq = seq.upper()
    results = {}
    for motif in motifs:
        motif = motif.strip().upper()
        if not motif:
            continue
        positions = [m.start() for m in re.finditer(f"(?={re.escape(motif)})", seq)]
        results[motif] = positions
    return results

def compute_shannon_entropy(seq):
    seq = seq.upper()
    counter = Counter(seq)
    length = len(seq)
    entropy = -sum((count/length) * math.log2(count/length) for count in counter.values() if count > 0)
    return round(entropy, 4)

# ====================== STREAMLIT APP ======================

st.title("🧬 DNA Stats Pro")
st.markdown("**Professional DNA Sequence Statistics Analyzer** | Built for learning & research")

# Sidebar
with st.sidebar:
    st.header("Input Sequence")
    
    input_method = st.radio("Choose input method:", ["Upload FASTA", "Paste Sequence"])
    
    if input_method == "Upload FASTA":
        uploaded_file = st.file_uploader("Upload FASTA file", type=["fasta", "fa", "txt"])
        sequence = ""
        seq_id = ""
        if uploaded_file:
            records = list(SeqIO.parse(uploaded_file, "fasta"))
            if records:
                sequence = str(records[0].seq)
                seq_id = records[0].id
    else:
        sequence = st.text_area("Paste DNA sequence here:", height=200)
        seq_id = "Pasted_Sequence"
    
    st.divider()
    
    # Example sequences
    st.subheader("Example Sequences")
    col1, col2 = st.columns(2)
    with col1:
        if st.button("Lambda Phage"):
            sequence = "GGGCGGCGACCTCGCGGGTTTTCGCTATTTATGAAAATTTTCCGGTTTAAGGCGTTTCCGTTCTTCTTCGTCATAACTTAATGTTTTTATTTAAAATACCCTCTGAAAAGAAAGGAAACGACAGGTGCTGAAAGCGAGGCTTTTTGGCCTCTGTCGTTTCCTTTCTCTGTTTTTGTCCGTGGAATGAACAATGGAAGTCAACAAAAAGCAGCTGGCTGACATTTTCGGTGCGAGTATCCGTACCATTCAGAACTGGCAGGAACAGGGAATGCCCGTTCTGCGAGGCGGTGGCAAGG"
            seq_id = "Lambda_Phage"
    with col2:
        if st.button("E. coli example"):
            sequence = "ATGAAACCCGGGTTTTAA" * 20
            seq_id = "Ecoli_Example"
    
    motifs_input = st.text_area("Motifs to search (one per line)", 
                               "ATG\nTATAAT\nGAATTC\nGGATCC", 
                               height=100)
    motifs = [m.strip() for m in motifs_input.split("\n") if m.strip()]
    
    window_size = st.selectbox("Sliding Window Size (bp)", [50, 100, 200, 500], index=1)
    
    run_button = st.button("🚀 Run Analysis", type="primary", use_container_width=True)

# Main Area
if run_button and sequence:
    with st.spinner("Analyzing DNA sequence..."):
        seq = sequence.upper().replace(" ", "").replace("\n", "")
        
        if len(seq) < 10:
            st.error("Sequence is too short!")
            st.stop()
        
        # Basic Stats
        stats = compute_basic_stats(seq)
        
        # Tabs
        tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
            "📊 Overview", 
            "🔠 Nucleotide Composition", 
            "🧬 Codon Usage", 
            "📏 k-mer Analysis",
            "🔍 Motif Search", 
            "📈 Complexity & GC Profile"
        ])
        
        with tab1:
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Length", f"{stats['length']:,} bp")
            with col2:
                st.metric("GC Content", f"{stats['gc_content']}%")
            with col3:
                st.metric("AT Content", f"{stats['at_content']}%")
            with col4:
                st.metric("GC Skew", stats['gc_skew'])
            
            st.success("Analysis completed successfully!")
        
        with tab2:
            counts = stats["counts"]
            df_counts = pd.DataFrame.from_dict(counts, orient='index', columns=['Count'])
            fig = px.bar(x=df_counts.index, y=df_counts['Count'], 
                        labels={'x':'Nucleotide', 'y':'Count'},
                        color=df_counts.index,
                        title="Nucleotide Composition")
            st.plotly_chart(fig, use_container_width=True)
            st.dataframe(df_counts, use_container_width=True)
        
        with tab3:
            codon_df = compute_codon_usage(seq)
            st.dataframe(codon_df.head(30), use_container_width=True)
            
            fig = px.bar(codon_df.head(20), x="Codon", y="Count", 
                        color="Amino Acid", title="Top 20 Codon Usage")
            st.plotly_chart(fig, use_container_width=True)
        
        with tab4:
            col1, col2 = st.columns(2)
            with col1:
                st.subheader("Dinucleotides (k=2)")
                di_df = compute_kmer_freq(seq, k=2)
                st.dataframe(di_df.head(15), use_container_width=True)
            with col2:
                st.subheader("Trinucleotides (k=3)")
                tri_df = compute_kmer_freq(seq, k=3)
                st.dataframe(tri_df.head(15), use_container_width=True)
        
        with tab5:
            if motifs:
                results = search_motifs(seq, motifs)
                for motif, positions in results.items():
                    if positions:
                        st.success(f"**{motif}** found {len(positions)} times")
                        pos_df = pd.DataFrame({
                            "0-based": positions,
                            "1-based": [p+1 for p in positions],
                            "End": [p + len(motif) for p in positions]
                        })
                        st.dataframe(pos_df, use_container_width=True)
                    else:
                        st.warning(f"Motif **{motif}** not found.")
            else:
                st.info("No motifs provided.")
        
        with tab6:
            entropy = compute_shannon_entropy(seq)
            st.metric("Shannon Entropy", f"{entropy} bits")
            st.caption("Higher value = more random/complex sequence (max ≈ 2.0)")
            
            # Simple GC window plot
            positions = list(range(0, len(seq)-window_size+1, window_size//2))
            gc_values = [gc_fraction(seq[i:i+window_size])*100 for i in range(0, len(seq)-window_size+1, window_size//2)]
            
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=positions, y=gc_values, mode='lines', name='GC%'))
            fig.update_layout(title=f"Sliding Window GC Content ({window_size} bp)", 
                            xaxis_title="Position (bp)", yaxis_title="GC Content (%)")
            st.plotly_chart(fig, use_container_width=True)

else:
    st.info("👈 Upload a FASTA file or paste a DNA sequence in the sidebar and click **Run Analysis**")

# Footer
st.caption("DNA Stats Pro — Educational Bioinformatics Tool")