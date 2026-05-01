#!/usr/bin/env python3
"""
annotation_report.py

Summarizes annotated gene products across all bins from Prokka or Bakta output.
Auto-detects the annotation tool from directory contents. Reports hypothetical
proteins separately, enriches top hits via UniProt, clusters genes by functional
category, and outputs both a TSV and a polished HTML report.

Supported tools:
    - Prokka  (reads .tsv, .gff)
    - Bakta   (reads .tsv, .gff3, .json for richer cross-references)

Usage:
    python3 annotation_report.py \
        --annotation_dir results/prokka \
        --outdir results/annotation_summary \
        --top_n 100

    python3 annotation_report.py \
        --annotation_dir results/bakta \
        --outdir results/annotation_summary \
        --tool bakta          # optional: override auto-detection

    python3 annotation_report.py \
        --annotation_dir results/bakta \
        --outdir results/annotation_summary \
        --no_uniprot          # skip UniProt lookup, omit clustering and gene/function columns
"""

import argparse
import os
import sys
import glob
import time
import json
import re
from collections import Counter, defaultdict
import html as html_module
import urllib.request
import urllib.parse
import urllib.error
from datetime import datetime

__version__ = '0.1.0'

# ─────────────────────────────────────────────
# Tool detection
# ─────────────────────────────────────────────

def detect_tool(annotation_dir):
    gff3_files = glob.glob(os.path.join(annotation_dir, '**', '*.gff3'), recursive=True)
    json_files = glob.glob(os.path.join(annotation_dir, '**', '*.json'), recursive=True)
    gff_files  = glob.glob(os.path.join(annotation_dir, '**', '*.gff'),  recursive=True)
    tsv_files  = glob.glob(os.path.join(annotation_dir, '**', '*.tsv'),  recursive=True)

    if gff3_files or json_files:
        print("Auto-detected: Bakta output (found .gff3 / .json files)")
        return 'bakta'
    elif gff_files or tsv_files:
        print("Auto-detected: Prokka output (found .gff / .tsv files)")
        return 'prokka'
    else:
        sys.exit(
            f"ERROR: No annotation files found under '{annotation_dir}'.\n"
            "Expected .tsv/.gff (Prokka) or .tsv/.gff3/.json (Bakta)."
        )


# ─────────────────────────────────────────────
# GFF / GFF3 parsing
# ─────────────────────────────────────────────

def parse_gff_files(annotation_dir, extension='gff'):
    pattern   = os.path.join(annotation_dir, '**', f'*.{extension}')
    gff_files = glob.glob(pattern, recursive=True)

    contig_lengths = {}
    rna_counter    = Counter()

    for gff_path in gff_files:
        with open(gff_path, errors='replace') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('##sequence-region'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            contig_id  = parts[1]
                            seq_length = int(parts[3]) - int(parts[2]) + 1
                            contig_lengths[contig_id] = seq_length
                        except (ValueError, IndexError):
                            pass
                    continue
                if line.startswith('#') or not line:
                    continue
                cols = line.split('\t')
                if len(cols) < 3:
                    continue
                feature = cols[2]
                if feature in ('tRNA', 'rRNA', 'tmRNA', 'ncRNA', 'repeat_region',
                               'CRISPR', 'misc_RNA', 'regulatory'):
                    rna_counter[feature] += 1

    return contig_lengths, rna_counter


# ─────────────────────────────────────────────
# TSV parsing
# ─────────────────────────────────────────────

def parse_tsvs(annotation_dir, tool):
    tsv_files = glob.glob(os.path.join(annotation_dir, '**', '*.tsv'), recursive=True)
    tsv_files = [f for f in tsv_files if 'hypothetical' not in os.path.basename(f).lower()]

    if not tsv_files:
        sys.exit(f"ERROR: No .tsv files found under '{annotation_dir}'.")

    print(f"Found {len(tsv_files)} TSV file(s) for tool={tool}.")

    cds_counter        = Counter()
    hypothetical_count = 0
    feature_counter    = Counter()
    total_cds          = 0
    bin_count          = 0
    ec_counter         = Counter()
    cog_counter        = Counter()
    dbxref_counter     = Counter()

    for tsv_path in tsv_files:
        bin_count += 1
        with open(tsv_path, errors='replace') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                parts = line.split('\t')

                if tool == 'prokka':
                    if parts[0] == 'locus_tag':
                        continue
                    if len(parts) < 4:
                        continue
                    feature = parts[1].strip()
                    product = parts[6].strip() if len(parts) > 6 else (parts[3].strip() if len(parts) > 3 else 'unknown')
                    ec      = parts[4].strip() if len(parts) > 4 else ''
                    cog     = parts[5].strip() if len(parts) > 5 else ''

                    feature_counter[feature] += 1

                    if feature == 'CDS':
                        total_cds += 1
                        if product.lower() in ('hypothetical protein', 'hypothetical_protein', ''):
                            hypothetical_count += 1
                        else:
                            cds_counter[product] += 1
                            if ec:
                                ec_counter[ec] += 1
                            if cog:
                                cog_counter[cog] += 1

                elif tool == 'bakta':
                    if parts[0].startswith('#') or parts[0] == 'Sequence Id':
                        continue
                    if len(parts) < 4:
                        continue
                    feature = parts[1].strip()
                    product = parts[7].strip() if len(parts) > 7 else 'unknown'
                    dbxrefs = parts[8].strip() if len(parts) > 8 else ''

                    feature_counter[feature.lower()] += 1

                    if feature.lower() == 'cds':
                        total_cds += 1
                        if product.lower() in ('hypothetical protein', 'hypothetical_protein', ''):
                            hypothetical_count += 1
                        else:
                            cds_counter[product] += 1
                            if dbxrefs:
                                for ref in dbxrefs.split(','):
                                    ref = ref.strip()
                                    db  = ref.split(':')[0] if ':' in ref else ref
                                    if db:
                                        dbxref_counter[db] += 1

    print(f"Processed {bin_count} bin(s).")
    print(f"Total CDS: {total_cds} ({hypothetical_count} hypothetical, "
          f"{total_cds - hypothetical_count} annotated)")

    return (cds_counter, hypothetical_count, feature_counter,
            total_cds, bin_count, ec_counter, cog_counter, dbxref_counter)


# ─────────────────────────────────────────────
# Bakta JSON enrichment
# ─────────────────────────────────────────────

def parse_bakta_json(annotation_dir):
    json_files = glob.glob(os.path.join(annotation_dir, '**', '*.json'), recursive=True)
    meta = {}

    for jf in json_files[:1]:
        try:
            with open(jf) as f:
                data = json.load(f)
            meta['bakta_version'] = data.get('version', '')
            meta['db_version']    = data.get('db', {}).get('version', '')
            stats = data.get('stats', {})
            meta['genome_size']   = stats.get('size', 0)
            meta['gc_content']    = stats.get('gc', 0.0)
            meta['no_sequences']  = stats.get('no_sequences', 0)
        except Exception:
            pass

    return meta


# ─────────────────────────────────────────────
# UniProt lookup
# ─────────────────────────────────────────────

def lookup_uniprot(product_name):
    try:
        params = urllib.parse.urlencode({
            'query': f'protein_name:{product_name} AND reviewed:true',
            'fields': 'gene_names,keyword,cc_function',
            'format': 'json',
            'size': '1'
        })
        url = f"https://rest.uniprot.org/uniprotkb/search?{params}"
        req = urllib.request.Request(url, headers={'User-Agent': 'annotation_report/2.0'})
        with urllib.request.urlopen(req, timeout=10) as response:
            data = json.loads(response.read().decode())

        results = data.get('results', [])
        if not results:
            return '—', '—', []

        entry = results[0]

        # Gene name
        gene_names = entry.get('genes', [])
        gene_name  = '—'
        if gene_names:
            primary   = gene_names[0].get('geneName', {})
            gene_name = primary.get('value', '—')

        # Function
        function_text = '—'
        for comment in entry.get('comments', []):
            if comment.get('commentType') == 'FUNCTION':
                texts = comment.get('texts', [])
                if texts:
                    raw      = texts[0].get('value', '')
                    sentence = raw.split('.')[0].strip()
                    if len(sentence) > 200:
                        sentence = sentence[:197] + '...'
                    function_text = sentence + '.'
                    break

        # Keywords
        keywords = [k['name'] for k in entry.get('keywords', [])]

        return gene_name, function_text, keywords

    except Exception:
        return '—', '—', []


def enrich_with_uniprot(top_genes):
    enriched = []
    total    = len(top_genes)
    print(f"\nQuerying UniProt for {total} gene products...")

    for i, (product, count) in enumerate(top_genes, 1):
        print(f"  [{i}/{total}] {product[:70]}", end='\r', flush=True)
        gene_name, function_text, keywords = lookup_uniprot(product)
        enriched.append((product, count, gene_name, function_text, keywords))
        time.sleep(0.2)

    print(f"\nUniProt lookup complete.{' ' * 60}")
    return enriched


# ─────────────────────────────────────────────
# Functional clustering
# ─────────────────────────────────────────────

CATEGORY_MAP = [
    ('DNA Metabolism',           ['DNA replication', 'DNA repair', 'DNA-binding',
                                  'DNA-directed DNA polymerase', 'Nucleotidyltransferase',
                                  'DNA recombination', 'DNA packaging', 'Chromosome']),
    ('Transcription',            ['Transcription', 'Transcription regulation',
                                  'RNA-binding', 'Sigma factor', 'mRNA']),
    ('Translation & Ribosomes',  ['Protein biosynthesis', 'Ribosomal protein',
                                  'Ribosome', 'Elongation factor', 'Initiation factor',
                                  'tRNA', 'Aminoacyl-tRNA']),
    ('Energy & Metabolism',      ['ATP synthesis', 'Oxidoreductase', 'Electron transport',
                                  'Hydrogen ion transport', 'NAD', 'FAD', 'TCA cycle',
                                  'Glycolysis', 'Photosynthesis', 'Carbon fixation',
                                  'Lipid metabolism', 'Fatty acid']),
    ('Transport & Membrane',     ['Transport', 'Ion transport', 'Membrane', 'Transmembrane',
                                  'Transporter', 'ABC transporter', 'Porin', 'Secretion']),
    ('Stress & Chaperones',      ['Chaperone', 'Stress response', 'Heat shock',
                                  'Oxidative stress', 'SOS response', 'Protease',
                                  'Unfolded protein response']),
    ('Cell Division & Structure',['Cell division', 'Cell cycle', 'Cytoskeleton',
                                  'Cell shape', 'Peptidoglycan', 'Cell wall',
                                  'Cell membrane']),
    ('Signaling & Regulation',   ['Kinase', 'Phosphoprotein', 'Two-component regulatory system',
                                  'Signal transduction', 'Allosteric enzyme']),
    ('Nucleotide Binding',       ['ATP-binding', 'GTP-binding', 'Nucleotide-binding',
                                  'Isomerase', 'Ligase', 'Transferase', 'Hydrolase']),
    ('Other / Unclassified',     []),
]

CLUSTER_COLORS = [
    '#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6',
    '#06b6d4', '#f97316', '#ec4899', '#84cc16', '#6b7280'
]


def cluster_genes(enriched_genes):
    clusters = {cat: [] for cat, _ in CATEGORY_MAP}

    for entry in enriched_genes:
        product, count, gene_name, function_text, keywords = entry
        assigned = False
        kw_set   = set(keywords)

        for category, kw_list in CATEGORY_MAP[:-1]:
            if kw_set & set(kw_list):
                clusters[category].append(entry)
                assigned = True
                break

        if not assigned:
            clusters['Other / Unclassified'].append(entry)

    return {k: v for k, v in clusters.items() if v}


# ─────────────────────────────────────────────
# TSV output
# ─────────────────────────────────────────────

def write_tsv(outpath, enriched_genes, hypothetical_count, feature_counter,
              total_cds, bin_count, top_n, tool):
    annotated = total_cds - hypothetical_count
    with open(outpath, 'w') as f:
        f.write(f"# Annotation Gene Product Summary\n")
        f.write(f"# Tool: {tool.capitalize()}\n")
        f.write(f"# Bins analyzed: {bin_count}\n")
        f.write(f"# Total CDS: {total_cds}\n")
        f.write(f"# Hypothetical proteins: {hypothetical_count} "
                f"({100*hypothetical_count/total_cds:.1f}% of CDS)\n")
        f.write(f"# Annotated CDS: {annotated}\n")
        f.write(f"# Top {top_n} annotated gene products shown\n#\n")
        f.write("rank\tcount\tpercent_of_total_cds\tproduct\tgene_name\tfunction\tkeywords\n")

        for rank, (product, count, gene_name, func, keywords) in enumerate(enriched_genes, 1):
            pct    = 100 * count / total_cds
            kw_str = '|'.join(keywords) if keywords else '—'
            f.write(f"{rank}\t{count}\t{pct:.2f}%\t{product}\t{gene_name}\t{func}\t{kw_str}\n")

        f.write(f"\n# Feature type summary\n")
        f.write("feature_type\tcount\n")
        for feature, count in feature_counter.most_common():
            f.write(f"{feature}\t{count}\n")

    print(f"TSV written to: {outpath}")


# ─────────────────────────────────────────────
# HTML output
# ─────────────────────────────────────────────

def write_html(outpath, enriched_genes, hypothetical_count, feature_counter,
               total_cds, bin_count, top_n, tool, rna_counter,
               contig_lengths, ec_counter, cog_counter, dbxref_counter,
               bakta_meta, annotation_dir, no_uniprot=False):

    annotated = total_cds - hypothetical_count
    hyp_pct   = 100 * hypothetical_count / total_cds if total_cds else 0
    ann_pct   = 100 * annotated / total_cds if total_cds else 0

    total_contigs = len(contig_lengths)
    total_bp      = sum(contig_lengths.values())
    total_mbp     = total_bp / 1_000_000

    # N50
    n50 = 0
    if contig_lengths:
        sorted_lens = sorted(contig_lengths.values(), reverse=True)
        cumsum = 0
        for ln in sorted_lens:
            cumsum += ln
            if cumsum >= total_bp / 2:
                n50 = ln
                break

    total_rna  = sum(rna_counter.values())
    tool_label = tool.capitalize()
    tool_color = '#10b981' if tool == 'bakta' else '#3b82f6'

    # ── Bakta metadata pills ──
    bakta_meta_html = ''
    if tool == 'bakta' and bakta_meta:
        items = []
        if bakta_meta.get('bakta_version'):
            items.append(f"Bakta v{html_module.escape(str(bakta_meta['bakta_version']))}")
        if bakta_meta.get('db_version'):
            items.append(f"DB v{html_module.escape(str(bakta_meta['db_version']))}")
        if bakta_meta.get('gc_content'):
            items.append(f"GC {bakta_meta['gc_content']:.1f}%")
        if items:
            bakta_meta_html = (
                '<div class="meta-pill-row">'
                + ''.join(f'<span class="meta-pill">{i}</span>' for i in items)
                + '</div>'
            )

    # ── Stat cards ──
    def stat_card(number, label, sub='', color=''):
        color_attr = f'style="color:{color}"' if color else ''
        sub_html   = f'<div class="stat-sub">{sub}</div>' if sub else ''
        return f'''
        <div class="stat-card">
            <div class="stat-num" {color_attr}>{number}</div>
            <div class="stat-label">{label}</div>
            {sub_html}
        </div>'''

    mbp_display = f'{total_mbp:.1f} Mbp' if total_mbp >= 1 else f'{total_bp:,} bp'
    stats_html = (
        stat_card(f'{bin_count:,}',      'Bins / MAGs')
        + stat_card(f'{total_contigs:,}','Contigs', sub=f'N50: {n50:,} bp' if n50 else '')
        + stat_card(mbp_display,          'Assembly Size')
        + stat_card(f'{total_cds:,}',    'Total CDS')
        + stat_card(f'{ann_pct:.1f}%',   'Annotation Rate',
                    color=('#10b981' if ann_pct >= 60 else '#f59e0b'))
        + stat_card(f'{total_rna:,}',    'RNA Features')
    )

    # ── RNA rows ──
    rna_rows = ''
    for feat, cnt in rna_counter.most_common():
        rna_rows += f'<tr><td>{html_module.escape(feat)}</td><td>{cnt:,}</td></tr>'
    if not rna_rows:
        rna_rows = '<tr><td colspan="2" class="empty">No RNA features found in GFF</td></tr>'

    # ── Feature type rows ──
    feat_total   = sum(feature_counter.values())
    feature_rows = ''
    for feat, cnt in feature_counter.most_common():
        pct = 100 * cnt / feat_total if feat_total else 0
        bar = int(pct * 1.5)
        feature_rows += f'''
        <tr>
            <td class="feat-name">{html_module.escape(feat)}</td>
            <td class="feat-count">{cnt:,}</td>
            <td class="feat-bar">
                <div class="mini-bar-wrap">
                    <div class="mini-bar" style="width:{bar}px"></div>
                    <span>{pct:.1f}%</span>
                </div>
            </td>
        </tr>'''

    # ── Sidebar tags ──
    def top_tags(counter, n=8, label=''):
        if not counter:
            return ''
        rows    = ''
        total_c = sum(counter.values())
        for item, cnt in counter.most_common(n):
            pct = 100 * cnt / total_c
            rows += f'''
            <div class="tag-row">
                <span class="tag-name">{html_module.escape(str(item))}</span>
                <span class="tag-count">{cnt:,}</span>
                <div class="tag-bar-bg"><div class="tag-bar" style="width:{min(pct*2,100):.0f}%"></div></div>
            </div>'''
        return f'<div class="tag-section"><div class="tag-title">{label}</div>{rows}</div>'

    sidebar_html = ''
    if tool == 'prokka':
        sidebar_html += top_tags(ec_counter,     label='Top EC Numbers')
        sidebar_html += top_tags(cog_counter,    label='Top COG Categories')
    else:
        sidebar_html += top_tags(dbxref_counter, label='Database Cross-References')

    # ── Functional clusters (only when UniProt was used) ──
    cluster_section_html = ''
    if not no_uniprot:
        clusters = cluster_genes(enriched_genes)

        cluster_cards_html = ''
        for i, (category, genes) in enumerate(clusters.items()):
            color     = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
            total_cat = sum(g[1] for g in genes)
            gene_tags = ''.join(
                f'<span class="cluster-tag" title="{html_module.escape(g[2] if g[2] != chr(8212) else g[0])}">'
                f'{html_module.escape(g[0][:45])}'
                f'<span class="cluster-tag-count">{g[1]:,}</span></span>'
                for g in sorted(genes, key=lambda x: x[1], reverse=True)
            )
            cluster_cards_html += f'''
            <div class="cluster-card" style="--cluster-color:{color}">
                <div class="cluster-header">
                    <span class="cluster-name">{html_module.escape(category)}</span>
                    <span class="cluster-stats">{len(genes)} gene{"s" if len(genes) != 1 else ""} &nbsp;·&nbsp; {total_cat:,} CDS</span>
                </div>
                <div class="cluster-tags">{gene_tags}</div>
            </div>'''

        max_cat_count    = max(sum(g[1] for g in genes) for genes in clusters.values())
        cluster_bar_html = ''
        for i, (category, genes) in enumerate(clusters.items()):
            color     = CLUSTER_COLORS[i % len(CLUSTER_COLORS)]
            total_cat = sum(g[1] for g in genes)
            pct       = 100 * total_cat / total_cds if total_cds else 0
            bar_w     = int(300 * total_cat / max_cat_count)
            cluster_bar_html += f'''
            <div class="cbar-row">
                <div class="cbar-label">{html_module.escape(category)}</div>
                <div class="cbar-bar" style="width:{bar_w}px; background:{color}"></div>
                <div class="cbar-count">{total_cat:,} <span class="cbar-pct">({pct:.1f}%)</span></div>
            </div>'''

        cluster_section_html = f'''
    <!-- Functional Clusters -->
    <div class="section">
        <div class="section-title">Functional Clusters</div>
        <p style="font-size:0.8em; color:var(--muted); margin-bottom:18px;">
            Top {top_n} genes grouped by biological function based on UniProt keywords.
            Each tag shows the gene product name and its total CDS count.
        </p>
        <div class="cluster-grid">
            {cluster_cards_html}
        </div>
    </div>

    <!-- Cluster distribution bar chart -->
    <div class="section">
        <div class="section-title">Functional Category Distribution</div>
        <p style="font-size:0.8em; color:var(--muted); margin-bottom:16px;">
            Total CDS counts per functional category across all bins.
        </p>
        <div style="padding: 4px 0;">
            {cluster_bar_html}
        </div>
    </div>'''

    # ── Gene table rows ──
    gene_rows = ''
    for rank, (product, count, gene_name, func, keywords) in enumerate(enriched_genes, 1):
        pct   = 100 * count / total_cds
        bar_w = min(int(pct * 25), 120)
        extra_cols = '' if no_uniprot else f'''
            <td class="g-gene">{html_module.escape(gene_name)}</td>
            <td class="g-func">{html_module.escape(func)}</td>'''
        gene_rows += f'''
        <tr>
            <td class="g-rank">{rank}</td>
            <td class="g-product">{html_module.escape(product)}</td>
            {extra_cols}
            <td class="g-count">{count:,}</td>
            <td class="g-pct">
                <div class="bar-wrap">
                    <div class="bar" style="width:{bar_w}px"></div>
                    <span>{pct:.2f}%</span>
                </div>
            </td>
        </tr>'''

    gene_table_extra_headers = '' if no_uniprot else '<th>Gene</th><th>Function (UniProt)</th>'

    generated  = datetime.now().strftime('%Y-%m-%d %H:%M')
    source_dir = html_module.escape(os.path.abspath(annotation_dir))

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Annotation Report — {tool_label}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link href="https://fonts.googleapis.com/css2?family=DM+Mono:ital,wght@0,300;0,400;0,500;1,400&family=Syne:wght@400;600;700;800&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg:       #0f1117;
            --surface:  #181c25;
            --surface2: #1e2330;
            --border:   #2a3040;
            --text:     #e2e8f0;
            --muted:    #64748b;
            --accent:   {tool_color};
            --accent2:  #f59e0b;
            --danger:   #ef4444;
            --success:  #10b981;
            --mono:     'DM Mono', monospace;
            --sans:     'Syne', sans-serif;
        }}

        * {{ box-sizing: border-box; margin: 0; padding: 0; }}

        body {{
            font-family: var(--sans);
            background: var(--bg);
            color: var(--text);
            min-height: 100vh;
        }}

        .header {{
            background: var(--surface);
            border-bottom: 1px solid var(--border);
            padding: 28px 40px 24px;
            display: flex;
            align-items: flex-start;
            justify-content: space-between;
            gap: 20px;
        }}

        .header-left h1 {{
            font-size: 1.6em;
            font-weight: 800;
            letter-spacing: -0.02em;
            color: #fff;
            line-height: 1.2;
        }}

        .header-left .source-path {{
            font-family: var(--mono);
            font-size: 0.72em;
            color: var(--muted);
            margin-top: 6px;
            word-break: break-all;
        }}

        .tool-badge {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: var(--accent);
            color: #fff;
            font-size: 0.78em;
            font-weight: 700;
            letter-spacing: 0.08em;
            text-transform: uppercase;
            padding: 5px 12px;
            border-radius: 20px;
            white-space: nowrap;
            margin-top: 4px;
        }}

        .tool-badge::before {{
            content: '';
            width: 7px; height: 7px;
            border-radius: 50%;
            background: rgba(255,255,255,0.7);
            display: inline-block;
        }}

        .meta-pill-row {{ display: flex; flex-wrap: wrap; gap: 6px; margin-top: 10px; }}
        .meta-pill {{
            background: var(--surface2);
            border: 1px solid var(--border);
            font-family: var(--mono);
            font-size: 0.72em;
            color: var(--muted);
            padding: 3px 9px;
            border-radius: 4px;
        }}

        .page-body {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 32px 40px;
        }}

        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(150px, 1fr));
            gap: 12px;
            margin-bottom: 28px;
        }}

        .stat-card {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 18px 20px;
            text-align: center;
        }}

        .stat-num {{
            font-size: clamp(1.1em, 2vw, 1.8em);
            font-weight: 800;
            color: var(--accent);
            letter-spacing: -0.03em;
            line-height: 1;
            word-break: break-word;
        }}

        .stat-label {{
            font-size: 0.72em;
            color: var(--muted);
            text-transform: uppercase;
            letter-spacing: 0.07em;
            margin-top: 5px;
        }}

        .stat-sub {{
            font-family: var(--mono);
            font-size: 0.68em;
            color: var(--muted);
            margin-top: 3px;
        }}

        .section {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 24px 28px;
            margin-bottom: 20px;
        }}

        .section-title {{
            font-size: 0.75em;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.1em;
            color: var(--muted);
            margin-bottom: 18px;
            padding-bottom: 12px;
            border-bottom: 1px solid var(--border);
        }}

        .hyp-callout {{
            background: #1c1a0e;
            border: 1px solid #3d3510;
            border-left: 4px solid var(--accent2);
            border-radius: 8px;
            padding: 18px 22px;
            display: flex;
            align-items: center;
            gap: 20px;
        }}

        .hyp-num  {{ font-size: 2.2em; font-weight: 800; color: var(--accent2); white-space: nowrap; letter-spacing: -0.03em; }}
        .hyp-text {{ font-size: 0.85em; color: #a8915a; line-height: 1.5; }}
        .hyp-pct  {{ font-family: var(--mono); font-size: 0.8em; margin-top: 3px; color: var(--muted); }}

        .two-col {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin-bottom: 20px;
        }}

        @media (max-width: 900px) {{ .two-col {{ grid-template-columns: 1fr; }} }}

        table {{ width: 100%; border-collapse: collapse; font-size: 0.84em; }}

        thead th {{
            background: var(--surface2);
            padding: 9px 12px;
            text-align: left;
            font-size: 0.72em;
            font-weight: 700;
            text-transform: uppercase;
            letter-spacing: 0.07em;
            color: var(--muted);
            border-bottom: 1px solid var(--border);
        }}

        tbody td {{
            padding: 9px 12px;
            border-bottom: 1px solid var(--border);
            vertical-align: middle;
        }}

        tbody tr:hover {{ background: var(--surface2); }}
        tbody tr:last-child td {{ border-bottom: none; }}
        .empty {{ color: var(--muted); font-style: italic; padding: 16px 12px; }}

        .feat-name  {{ font-family: var(--mono); color: var(--accent); font-size: 0.9em; }}
        .feat-count {{ text-align: right; font-variant-numeric: tabular-nums; color: var(--text); font-weight: 600; width: 70px; }}
        .feat-bar   {{ width: 200px; }}
        .mini-bar-wrap {{ display: flex; align-items: center; gap: 8px; }}
        .mini-bar {{ height: 6px; background: var(--accent); border-radius: 2px; min-width: 2px; opacity: 0.7; }}
        .mini-bar-wrap span {{ font-size: 0.78em; color: var(--muted); white-space: nowrap; font-family: var(--mono); }}

        .rna-table td:first-child {{ font-family: var(--mono); color: var(--success); }}
        .rna-table td:last-child  {{ text-align: right; color: var(--text); font-weight: 600; }}

        .tag-section {{ margin-bottom: 20px; }}
        .tag-title {{ font-size: 0.68em; font-weight: 700; text-transform: uppercase;
                      letter-spacing: 0.1em; color: var(--muted); margin-bottom: 10px; }}
        .tag-row {{ display: flex; align-items: center; gap: 8px; padding: 5px 0;
                    border-bottom: 1px solid var(--border); }}
        .tag-row:last-child {{ border-bottom: none; }}
        .tag-name {{ font-family: var(--mono); font-size: 0.8em; color: var(--accent);
                     flex: 0 0 auto; min-width: 80px; }}
        .tag-count {{ font-size: 0.78em; color: var(--muted); flex: 0 0 50px;
                      text-align: right; font-variant-numeric: tabular-nums; }}
        .tag-bar-bg {{ flex: 1; height: 4px; background: var(--border); border-radius: 2px; }}
        .tag-bar {{ height: 100%; background: var(--accent); border-radius: 2px; opacity: 0.6; }}

        /* ── Functional cluster cards ── */
        .cluster-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(360px, 1fr));
            gap: 16px;
        }}

        .cluster-card {{
            background: var(--surface2);
            border: 1px solid var(--border);
            border-top: 3px solid var(--cluster-color);
            border-radius: 8px;
            padding: 16px 18px;
        }}

        .cluster-header {{
            display: flex;
            justify-content: space-between;
            align-items: baseline;
            margin-bottom: 12px;
            gap: 8px;
        }}

        .cluster-name  {{ font-weight: 700; font-size: 0.88em; color: var(--cluster-color); }}
        .cluster-stats {{ font-family: var(--mono); font-size: 0.7em; color: var(--muted); white-space: nowrap; }}
        .cluster-tags  {{ display: flex; flex-wrap: wrap; gap: 6px; }}

        .cluster-tag {{
            background: var(--surface);
            border: 1px solid var(--border);
            border-radius: 4px;
            padding: 3px 8px;
            font-size: 0.74em;
            color: var(--text);
            display: inline-flex;
            align-items: center;
            gap: 5px;
            cursor: default;
            max-width: 100%;
        }}

        .cluster-tag-count {{
            background: var(--cluster-color);
            color: #fff;
            border-radius: 3px;
            padding: 1px 5px;
            font-family: var(--mono);
            font-size: 0.82em;
            opacity: 0.85;
            flex-shrink: 0;
        }}

        /* ── Cluster bar chart ── */
        .cbar-row {{
            display: flex;
            align-items: center;
            gap: 14px;
            padding: 7px 0;
            border-bottom: 1px solid var(--border);
        }}

        .cbar-row:last-child {{ border-bottom: none; }}
        .cbar-label {{ font-size: 0.8em; color: var(--text); width: 220px; flex-shrink: 0; }}
        .cbar-bar   {{ height: 12px; border-radius: 3px; opacity: 0.8; min-width: 2px; }}
        .cbar-count {{ font-family: var(--mono); font-size: 0.78em; color: var(--text); font-weight: 600; white-space: nowrap; }}
        .cbar-pct   {{ color: var(--muted); font-weight: 400; }}

        /* ── Gene table ── */
        .g-rank    {{ color: var(--muted); font-family: var(--mono); font-size: 0.8em; width: 36px; }}
        .g-product {{ font-weight: 600; color: var(--text); min-width: 160px; }}
        .g-gene    {{ font-family: var(--mono); color: #a78bfa; width: 80px; font-size: 0.88em; }}
        .g-func    {{ color: var(--muted); font-size: 0.82em; line-height: 1.5; min-width: 240px; }}
        .g-count   {{ text-align: right; color: var(--accent); font-weight: 700;
                      font-variant-numeric: tabular-nums; width: 70px; font-family: var(--mono); }}
        .g-pct     {{ width: 160px; }}

        .bar-wrap  {{ display: flex; align-items: center; gap: 8px; }}
        .bar       {{ height: 8px; background: var(--accent); border-radius: 2px; min-width: 2px; opacity: 0.75; }}
        .bar-wrap span {{ font-family: var(--mono); font-size: 0.78em; color: var(--muted); white-space: nowrap; }}

        .search-bar {{
            width: 100%;
            background: var(--surface2);
            border: 1px solid var(--border);
            border-radius: 6px;
            padding: 9px 14px;
            color: var(--text);
            font-family: var(--mono);
            font-size: 0.85em;
            margin-bottom: 14px;
            outline: none;
            transition: border-color 0.2s;
        }}
        .search-bar:focus {{ border-color: var(--accent); }}
        .search-bar::placeholder {{ color: var(--muted); }}

        .footer {{
            text-align: center;
            color: var(--muted);
            font-size: 0.75em;
            font-family: var(--mono);
            margin-top: 32px;
            padding-top: 20px;
            border-top: 1px solid var(--border);
        }}
    </style>
</head>
<body>

<div class="header">
    <div class="header-left">
        <h1>Genome Annotation Report</h1>
        <div class="source-path">{source_dir}</div>
        {bakta_meta_html}
    </div>
    <div>
        <div class="tool-badge">{tool_label}</div>
    </div>
</div>

<div class="page-body">

    <div class="stats-grid">
        {stats_html}
    </div>

    <div class="section" style="margin-bottom:20px;">
        <div class="section-title">Hypothetical Proteins</div>
        <div class="hyp-callout">
            <div class="hyp-num">{hypothetical_count:,}</div>
            <div>
                <div class="hyp-text">
                    CDS with no known function — excluded from the top-{top_n} gene product table.
                </div>
                <div class="hyp-pct">
                    {hyp_pct:.1f}% of all CDS &nbsp;·&nbsp;
                    {annotated:,} annotated CDS ({ann_pct:.1f}%)
                </div>
            </div>
        </div>
    </div>

    <div class="two-col">
        <div class="section">
            <div class="section-title">Feature Type Summary</div>
            <table>
                <thead>
                    <tr><th>Type</th><th style="text-align:right">Count</th><th>Distribution</th></tr>
                </thead>
                <tbody>{feature_rows}</tbody>
            </table>
        </div>
        <div class="section">
            <div class="section-title">RNA Features</div>
            <table class="rna-table">
                <thead><tr><th>Feature</th><th style="text-align:right">Count</th></tr></thead>
                <tbody>{rna_rows}</tbody>
            </table>
            {'<div class="tag-section" style="margin-top:20px;">' + sidebar_html + '</div>' if sidebar_html else ''}
        </div>
    </div>

    {cluster_section_html}

    <!-- Top N gene products -->
    <div class="section">
        <div class="section-title">Top {top_n} Annotated Gene Products</div>
        <p style="font-size:0.8em; color:var(--muted); margin-bottom:14px;">
            Percentages are of total CDS (including hypotheticals).
            {'Gene names and function descriptions from UniProt/Swiss-Prot reviewed entries.' if not no_uniprot else ''}
        </p>
        <input
            class="search-bar"
            type="text"
            id="geneSearch"
            placeholder="Filter by product name{',' if not no_uniprot else ''}{'gene, or function' if not no_uniprot else ''}..."
            oninput="filterTable()"
        >
        <table id="geneTable">
            <thead>
                <tr>
                    <th>#</th>
                    <th>Gene Product</th>
                    {gene_table_extra_headers}
                    <th style="text-align:right">Count</th>
                    <th>% of CDS</th>
                </tr>
            </thead>
            <tbody>{gene_rows}</tbody>
        </table>
    </div>

    <div class="footer">
        annotation_report.py &nbsp;·&nbsp; {tool_label} &nbsp;·&nbsp;
        Generated {generated}
        {'&nbsp;·&nbsp; UniProt/Swiss-Prot reviewed entries only' if not no_uniprot else '&nbsp;·&nbsp; UniProt lookup skipped'}
    </div>

</div>

<script>
function filterTable() {{
    const q    = document.getElementById('geneSearch').value.toLowerCase();
    const rows = document.querySelectorAll('#geneTable tbody tr');
    rows.forEach(row => {{
        row.style.display = row.textContent.toLowerCase().includes(q) ? '' : 'none';
    }});
}}
</script>

</body>
</html>"""

    with open(outpath, 'w') as f:
        f.write(html_content)

    print(f"HTML report written to: {outpath}")


# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            'Generate a gene annotation summary report from Prokka or Bakta output.\n'
            'Auto-detects the tool from directory contents; override with --tool.'
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--annotation_dir', required=True,
                        help='Path to annotation output directory (Prokka or Bakta)')
    parser.add_argument('--outdir', default='annotation_summary',
                        help='Output directory (default: annotation_summary)')
    parser.add_argument('--top_n', type=int, default=100,
                        help='Number of top gene products to report (default: 100)')
    parser.add_argument('--tool', choices=['prokka', 'bakta'], default=None,
                        help='Force tool type (default: auto-detect)')
    parser.add_argument('--no_uniprot', action='store_true',
                        help='Skip UniProt lookup — omits clustering sections and gene/function columns')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    tool = args.tool or detect_tool(args.annotation_dir)

    (cds_counter, hypothetical_count, feature_counter,
     total_cds, bin_count, ec_counter, cog_counter, dbxref_counter) = parse_tsvs(
        args.annotation_dir, tool
    )

    if total_cds == 0:
        sys.exit("ERROR: No CDS features found. Check your annotation directory.")

    ext = 'gff3' if tool == 'bakta' else 'gff'
    print(f"Parsing .{ext} files for contig and RNA data...")
    contig_lengths, rna_counter = parse_gff_files(args.annotation_dir, extension=ext)
    print(f"  Contigs: {len(contig_lengths):,}   RNA features: {sum(rna_counter.values()):,}")

    bakta_meta = {}
    if tool == 'bakta':
        bakta_meta = parse_bakta_json(args.annotation_dir)

    top_genes = cds_counter.most_common(args.top_n)
    if args.no_uniprot:
        print("Skipping UniProt lookup (--no_uniprot).")
        enriched_genes = [(p, c, '—', '—', []) for p, c in top_genes]
    else:
        enriched_genes = enrich_with_uniprot(top_genes)

    tsv_out  = os.path.join(args.outdir, 'annotation_gene_summary.tsv')
    html_out = os.path.join(args.outdir, 'annotation_gene_summary.html')

    write_tsv(tsv_out, enriched_genes, hypothetical_count, feature_counter,
              total_cds, bin_count, args.top_n, tool)

    write_html(html_out, enriched_genes, hypothetical_count, feature_counter,
               total_cds, bin_count, args.top_n, tool, rna_counter,
               contig_lengths, ec_counter, cog_counter, dbxref_counter,
               bakta_meta, args.annotation_dir, no_uniprot=args.no_uniprot)

    print(f"\nDone. Open: {html_out}")


if __name__ == '__main__':
    main()
