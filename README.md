# DNA Stats Pro 🧬

**Professional Web Application for DNA Sequence Statistical Analysis**

My first bioinformatics project — built to learn and demonstrate skills in Computational Biology and Medical Genetics.

## Features

- Basic DNA statistics (Sequence length, GC content, AT content, GC skew)
- Codon usage bias analysis with amino acid translation
- Dinucleotide (k=2) and Trinucleotide (k=3) frequency analysis
- Motif search (e.g. ATG, TATA-box, restriction sites such as GAATTC)
- Sequence complexity measurement (Shannon Entropy)
- Interactive sliding window GC content profile
- Modern, clean and responsive dark-themed dashboard

## Technologies Used

- **Python**
- **Streamlit** (Web Framework)
- **Biopython**
- **Plotly** (Interactive visualizations)
- **Pandas**

## How to Run Locally

```bash
git clone https://github.com/YOUR-USERNAME/DNA-Sequence-Statistics.git
cd DNA-Sequence-Statistics
pip install -r requirements.txt
streamlit run app.py
