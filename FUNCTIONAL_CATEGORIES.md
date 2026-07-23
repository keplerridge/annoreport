# Functional Categories in annoreport

annoreport assigns annotated gene products to one of ten biological functional 
categories using a priority-ordered keyword matching scheme applied to UniProt 
annotations. Each gene is assigned to the first matching category based on the 
priority order below. If no keywords match, the gene is assigned to 
Other/Unclassified.

## Assignment Scheme

Keywords are matched against UniProt keyword fields returned by the Swiss-Prot 
reviewed database. Matching is exact and case-sensitive. Priority follows the 
order of the table below — a gene matching keywords in multiple categories is 
assigned to the highest-priority category.

## Category Definitions

| Priority | Category | Keywords |
|---|---|---|
| 1 | DNA Metabolism | DNA replication, DNA repair, DNA-binding, DNA-directed DNA polymerase, Nucleotidyltransferase, DNA recombination, DNA packaging, Chromosome |
| 2 | Transcription | Transcription, Transcription regulation, RNA-binding, Sigma factor, mRNA |
| 3 | Translation & Ribosomes | Protein biosynthesis, Ribosomal protein, Ribosome, Elongation factor, Initiation factor, tRNA, Aminoacyl-tRNA |
| 4 | Energy & Metabolism | ATP synthesis, Oxidoreductase, Electron transport, Hydrogen ion transport, NAD, FAD, TCA cycle, Glycolysis, Photosynthesis, Carbon fixation, Lipid metabolism, Fatty acid |
| 5 | Transport & Membrane | Transport, Ion transport, Membrane, Transmembrane, Transporter, ABC transporter, Porin, Secretion |
| 6 | Stress & Chaperones | Chaperone, Stress response, Heat shock, Oxidative stress, SOS response, Protease, Unfolded protein response |
| 7 | Cell Division & Structure | Cell division, Cell cycle, Cytoskeleton, Cell shape, Peptidoglycan, Cell wall, Cell membrane |
| 8 | Signaling & Regulation | Kinase, Phosphoprotein, Two-component regulatory system, Signal transduction, Allosteric enzyme |
| 9 | Nucleotide Binding | ATP-binding, GTP-binding, Nucleotide-binding, Isomerase, Ligase, Transferase, Hydrolase |
| 10 | Other / Unclassified | No keyword match — assigned by default |

## Design Rationale

These categories were selected to capture broad functional themes commonly used 
in microbial ecological and environmental genomics, broadly mirroring clusters 
of orthologous groups (COG) produced by Prokka and KEGG metabolic pathway 
classifications. The framework is designed to provide an interpretable ecological 
summary rather than a formal ontology-based annotation system, enabling rapid 
identification of dominant functional patterns across large MAG collections 
without replacing pathway reconstruction or curated metabolic inference pipelines.

## Notes

- UniProt queries use the reviewed Swiss-Prot database only
- Keywords are matched against the UniProt keyword field, not free text
- A gene is assigned to exactly one category — the highest priority match
- Assugnment is based on priority order, not keyword frequency. A gene matching one keyword in a higher-priority category will be assigned there even if it match multiple keywords in a lower-priority category
- Genes with no UniProt match or no keyword overlap are assigned to 
  Other/Unclassified
- The --no_uniprot flag disables UniProt annotation and functional clustering entirely
