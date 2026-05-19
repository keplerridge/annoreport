# annoreport

A command-line tool for summarizing and visualizing gene annotation results from metagenome-assembled genome (MAG) workflows. Supports output from [Prokka](https://github.com/tseemann/prokka) and [Bakta](https://github.com/oschwengers/bakta), enriches results via the [UniProt REST API](https://www.uniprot.org/help/api), and produces a polished interactive HTML report alongside a TSV summary table.

[![Bioconda](https://img.shields.io/conda/vn/bioconda/annoreport.svg)](https://anaconda.org/bioconda/annoreport)

---

## Features

- **Auto-detects** Prokka or Bakta output from directory contents
- **Counts and ranks** the top N most common annotated gene products across all bins
- **Separates hypothetical proteins** from annotated CDS and reports them independently
- **UniProt enrichment** — looks up gene names and one-line function descriptions for each top gene product
- **Functional clustering** — groups genes into biological categories (DNA Metabolism, Translation, Energy & Metabolism, Stress & Chaperones, etc.) based on UniProt keywords
- **Interactive HTML report** with:
  - Summary stat cards (bins, contigs, assembly size, CDS counts, annotation rate, RNA features)
  - Feature type and RNA feature tables
  - Functional cluster cards with per-gene CDS counts
  - Functional category distribution bar chart
  - Searchable top-N gene product table
- **TSV output** for downstream analysis in R, Python, or Excel
- **`--no_uniprot` flag** for offline/fast runs — skips UniProt lookup and omits clustering sections

---

## Requirements

- Python 3.9+
- No external dependencies — uses Python standard library only
- Internet access required for UniProt enrichment (unless `--no_uniprot` is used)

---

## Installation

### Bioconda (recommended)
```bash
conda install -c bioconda annoreport
```

### From source
```bash
git clone https://github.com/keplerridge/annoreport.git
cd annoreport
```

Or copy `annotation_report.py` directly into your project's `scripts/` directory.

---

## Usage

### Basic (auto-detect tool)

```bash
python3 annotation_report.py \
    --annotation_dir results/prokka \
    --outdir results/annotation_summary
```

### Bakta output

```bash
python3 annotation_report.py \
    --annotation_dir results/bakta \
    --outdir results/annotation_summary
```

### Force tool type

```bash
python3 annotation_report.py \
    --annotation_dir results/bakta \
    --outdir results/annotation_summary \
    --tool bakta
```

### Skip UniProt lookup (offline / fast mode)

```bash
python3 annotation_report.py \
    --annotation_dir results/prokka \
    --outdir results/annotation_summary \
    --no_uniprot
```

### Change number of top genes reported

```bash
python3 annotation_report.py \
    --annotation_dir results/bakta \
    --outdir results/annotation_summary \
    --top_n 50
```

---

## Arguments

| Argument | Default | Description |
|---|---|---|
| `--annotation_dir` | *(required)* | Path to Prokka or Bakta output directory |
| `--outdir` | `annotation_summary` | Directory for output files |
| `--top_n` | `100` | Number of top gene products to report |
| `--tool` | auto-detect | Force tool type: `prokka` or `bakta` |
| `--no_uniprot` | off | Skip UniProt lookup; omits clustering and gene/function columns |

---

## Output

Two files are written to `--outdir`:

### `annotation_gene_summary.html`
An interactive HTML report containing:
- **Summary cards** — bins/MAGs, contigs, assembly size, total CDS, annotation rate, RNA features
- **Hypothetical protein callout** — count and percentage of CDS with no known function
- **Feature type summary** — counts of CDS, tRNA, rRNA, tmRNA, and other features
- **RNA features table** — breakdown of non-coding RNA annotations
- **Functional clusters** *(UniProt mode only)* — top genes grouped by biological function with CDS counts per gene
- **Functional category distribution** *(UniProt mode only)* — bar chart of CDS counts per category
- **Top N gene products table** — searchable table with gene product name, gene name, UniProt function description, CDS count, and percentage

### `annotation_gene_summary.tsv`
Tab-separated summary with columns:
```
rank  count  percent_of_total_cds  product  gene_name  function  keywords
```
Plus a feature type summary appended at the bottom.

---

## Supported Annotation Tools

| Tool | File Types Used | Notes |
|---|---|---|
| **Prokka** | `.tsv`, `.gff` | Reads EC numbers and COG categories if present |
| **Bakta** | `.tsv`, `.gff3`, `.json` | Reads database cross-references; auto-skips `hypotheticals.tsv` and `inference.tsv`|

Tool auto-detection checks for `.gff3` or `.json` files (Bakta) versus `.gff` or `.tsv` only (Prokka).

---

## Functional Categories

When UniProt lookup is enabled, genes are assigned to one of the following categories based on UniProt keyword matching (first match wins):

| Category | Example keywords |
|---|---|
| DNA Metabolism | DNA replication, DNA repair, DNA-binding |
| Transcription | Transcription, RNA-binding, Sigma factor |
| Translation & Ribosomes | Protein biosynthesis, Ribosomal protein, Elongation factor |
| Energy & Metabolism | ATP synthesis, Oxidoreductase, TCA cycle, Glycolysis |
| Transport & Membrane | Transport, Membrane, ABC transporter, Porin |
| Stress & Chaperones | Chaperone, Heat shock, Oxidative stress, Protease |
| Cell Division & Structure | Cell division, Peptidoglycan, Cell wall |
| Signaling & Regulation | Kinase, Two-component regulatory system, Signal transduction |
| Nucleotide Binding | ATP-binding, GTP-binding, Isomerase, Hydrolase |
| Other / Unclassified | No matching keywords found |

Each gene is assigned to exactly one category. Priority follows the order above.

---

## Example Workflow

This tool is designed to run after a MAG annotation step in a Snakemake workflow:

```python (snakemake)
rule annotation_report:
    input:
        annotation_dir = 'results/bakta'
    output:
        html = 'results/annotation_summary/annotation_gene_summary.html',
        tsv  = 'results/annotation_summary/annotation_gene_summary.tsv'
    params:
        outdir = 'results/annotation_summary'
    conda:
        'envs/annotation_report.yaml'
    threads: 1
    resources:
        mem_mb=4000,
        runtime=30
    shell:
        """
        python3 scripts/annotation_report.py \
            --annotation_dir {input.annotation_dir} \
            --outdir {params.outdir} \
            --top_n 100
        """
```

---

## Notes

- UniProt queries use the reviewed (Swiss-Prot) database only for high-quality annotations
- A 0.2 second delay is applied between UniProt API calls to respect rate limits
- For 100 gene products, the UniProt lookup phase takes approximately 30–40 seconds
- The `--no_uniprot` flag is recommended for quick runs or environments without internet access

---

## License

MIT License — see [LICENSE](LICENSE) for details.

---

## Citation

If you use this tool in your research, please cite it as:

> annoreport: a summary and visualization tool for MAG annotation workflows. 
> https://github.com/keplerridge/annoreport

> annoreport: a summary and visualization tool for MAG annotation workflows. https://github.com/keplerridge/annoreport
