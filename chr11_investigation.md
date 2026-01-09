# Chr11 Pangenome Investigation

## Overview

Analysis of the mouse pangenome reveals unusual characteristics on chromosome 11:
- **+122.2 Mb** excess sequence compared to mm10 reference (largest of all chromosomes)
- **~3.5 million nodes** (highest node count)
- **Visible duplication** signal in the pangenome graph

This document outlines the investigation strategy to understand these anomalies.

---

## Key Questions

### 1. **Where are these nodes localized on chr11?**
   - Are duplicates randomly distributed or clustered in specific regions?
   - Do they concentrate in particular genomic features?
   - Is there a gradient or hotspot pattern?

### 2. **Do they overlap with HOX clusters?**
   - **Hoxa cluster** (chr11:96,161,468-96,273,464 in mm10)
   - Check for paralogous sequences from other HOX clusters (Hoxb, Hoxc, Hoxd on other chromosomes)
   - Investigate potential mis-mappings or structural variants in HOX regions

### 3. **Are they near centromeric regions?**
   - Centromeric/pericentromeric regions (chr11:0-3,000,000)
   - These regions are repeat-rich and challenging for assembly
   - May cause graph fragmentation and apparent "duplication"

### 4. **What else should we investigate?**

---

## Additional Investigation Points

### A. Known Duplicated/Repetitive Regions on Chr11

#### **1. Immunoglobulin Lambda (IgL) Locus**
- **Location:** chr11:113,500,000-115,000,000 (approximate)
- **Characteristics:**
  - Contains multiple V, J, and C gene segments
  - High sequence similarity between segments
  - Expected to generate many nodes in pangenome
  - Structural variation across mouse strains

#### **2. Olfactory Receptor (OR) Gene Clusters**
- **Cluster 1:** chr11:49,000,000-52,000,000
- **Cluster 2:** chr11:58,000,000-60,000,000
- **Characteristics:**
  - Largest gene family in mammals
  - High copy number variation
  - Frequent duplication/deletion events
  - Poor assembly quality in many strains

#### **3. Segmental Duplications**
- Check UCSC Genome Browser for annotated segmental duplications
- These are >1kb sequences with >90% identity appearing in multiple locations
- Common source of pangenome complexity

### B. Assembly Quality Issues

#### **1. Supernova Assembly Artifacts**
- Long 10x contigs may contain unresolved repeats
- Pseudohaplotype collapse/expansion
- Check N50 statistics for chr11 across strains
- Investigate break points in assemblies

#### **2. Graph Construction Parameters**
- Review pggb parameters used:
  - `-s` (segment length): shorter = more nodes
  - `-p` (identity threshold): lower = more variation captured
  - `-k` (k-mer size)
  - `-B` (block size)
- May need chr11-specific parameter tuning

### C. Biological Variation

#### **1. Structural Variants (SVs)**
- Large duplications (>50kb) segregating in population
- Copy number variants (CNVs)
- Inversions creating apparent duplications
- Use `odgi stats` to check graph topology

#### **2. Strain-Specific Expansions**
- Some strains may have real chr11 expansions
- Check individual strain assemblies against reference
- Look for strain-specific branches in graph

### D. Technical Artifacts

#### **1. Read Mapping Errors**
- Sequences from other chromosomes mis-mapping to chr11
- Check wfmash/edyeet mapping statistics
- Investigate unmapped sequences

#### **2. Reference Gaps**
- mm10 reference has gaps on chr11
- Pangenome may correctly resolve these
- Check gap locations in reference

---

## Analysis Workflow

### Step 1: Run Initial Analysis
```bash
# Python analysis
python script/analyze_chr11_nodes.py DBmm10.nodes.tag.txt.gz

# R statistical analysis
Rscript script/chr11_statistical_analysis.R DBmm10.nodes.tag.txt.gz
```

**Outputs:**
- `chr11_analysis_overview.png` - Summary statistics
- `chr11_analysis_position.png` - Positional distribution
- `chr11_detailed_analysis.png` - R statistical plots
- `chr11_analysis_report.txt` - Text summary

### Step 2: Graph Topology Analysis
```bash
# If you have the GFA file
odgi stats -i chr11.og -S -m

# Check node degree distribution
odgi degree -i chr11.og > chr11_degree.txt

# Extract duplicate node sequences
odgi paths -i chr11.og -L | grep -A1 "<duplicate_node_id>"
```

### Step 3: Sequence-Level Investigation

#### For specific duplicate regions:
```bash
# Extract sequences from GFA
grep "^S" chr11.gfa | awk '$2 in duplicate_nodes' > chr11_duplicate_seqs.fa

# BLAST against mm10 to find origin
blastn -query chr11_duplicate_seqs.fa -db mm10 -outfmt 6 -max_target_seqs 10

# Check for repetitive elements
RepeatMasker chr11_duplicate_seqs.fa
```

### Step 4: Strain-Level Comparison
```bash
# Extract each strain's path through graph
for strain in $(cat strain_list.txt); do
  odgi paths -i chr11.og -f -P $strain > ${strain}_chr11_path.txt
done

# Compare path lengths
awk '{print FILENAME, NF}' *_chr11_path.txt | sort -k2 -rn
```

### Step 5: Compare to Other Pangenomes
- Check human HPRC pangenome - does chr2 (syntenic to mouse chr11) show similar pattern?
- Review literature on mouse structural variation
- Check other mouse pangenome projects (e.g., Sanger Mouse Genomes)

---

## Expected Findings

### Likely Explanations (in order of probability):

1. **Immunoglobulin Lambda Locus** - Most likely major contributor
   - High structural variation
   - Multiple paralogous segments
   - Known to be complex in all mammals

2. **Olfactory Receptor Clusters** - Secondary contributor
   - Copy number variation
   - Assembly difficulties
   - Strain-specific expansions

3. **Segmental Duplications** - Moderate contributor
   - Real biological duplications
   - May be collapsed in reference
   - Resolved in pangenome

4. **Centromeric Regions** - Minor contributor (if present)
   - Should be mostly filtered in preprocessing
   - Repeat-rich sequences

5. **Assembly Artifacts** - Variable contribution
   - Strain-dependent
   - Check assembly quality metrics

### Unexpected Findings to Watch For:

- **Novel structural variants** not present in reference
- **Contamination** from non-mouse sequences
- **Graph algorithm artifacts** from parameter settings
- **Mis-annotations** in reference genome

---

## Validation Steps

### 1. Cross-reference with Known Databases
- **UCSC Genome Browser:** Check for annotated features
- **Ensembl:** Verify gene annotations in duplicate regions
- **MGI (Mouse Genome Informatics):** Check for known variants
- **dbVar/DGVa:** Look for reported structural variants

### 2. Orthogonal Methods
- Call structural variants using traditional methods (Manta, Delly, Lumpy)
- Compare with short-read SV calls if available
- Use optical mapping data if available (Bionano)

### 3. Experimental Validation (if needed)
- qPCR for copy number
- FISH for large duplications
- Long-read sequencing of specific strains

---

## Reporting Findings

### Key Metrics to Report:
1. Total duplicated sequence (Mb) and percentage
2. Number of duplicate nodes
3. Genomic distribution (hotspots vs. distributed)
4. Overlap with known features (genes, repeats, SVs)
5. Strain-specific vs. shared duplications
6. Impact on downstream analyses (variant calling, etc.)

### Visualizations to Create:
1. ✅ Node count and length distributions
2. ✅ Positional heatmap along chr11
3. Circos plot showing duplications
4. Graph topology metrics (degree, component size)
5. Strain-specific presence/absence matrix
6. Comparison to other chromosomes

---

## Next Steps After Initial Analysis

Based on findings, consider:

1. **If duplicates are real biology:**
   - Document in supplementary materials
   - Consider separate analysis of complex regions
   - May need region-specific graph parameters

2. **If duplicates are technical artifacts:**
   - Adjust pggb parameters
   - Filter problematic strains/contigs
   - Re-run graph construction for chr11

3. **If duplicates are assembly errors:**
   - Identify problematic assemblies
   - Consider excluding or correcting them
   - Document limitations

---

## References

### Mouse Genome Resources:
- **mm10 reference:** [UCSC Genome Browser](https://genome.ucsc.edu/cgi-bin/hgGateway?db=mm10)
- **MGI:** [Mouse Genome Informatics](http://www.informatics.jax.org/)
- **Sanger Mouse Genomes:** [https://www.sanger.ac.uk/data/mouse-genomes-project/](https://www.sanger.ac.uk/data/mouse-genomes-project/)

### Pangenome Tools:
- **pggb:** [https://github.com/pangenome/pggb](https://github.com/pangenome/pggb)
- **odgi:** [https://github.com/pangenome/odgi](https://github.com/pangenome/odgi)
- **Bandage:** For graph visualization [https://rrwick.github.io/Bandage/](https://rrwick.github.io/Bandage/)

### Relevant Literature:
- Garrison et al. (2023) "Building pangenome graphs" - pggb methods
- Eizenga et al. (2020) "Pangenome graphs" - Graph genome concepts
- Mouse Genomes Project papers on structural variation

---

## Contact & Collaboration

If unusual findings are discovered:
- Consult with pangenome graph experts
- Share with pggb community (GitHub discussions)
- Consider collaboration with mouse genetics groups

---

## Summary Checklist

- [ ] Run Python analysis script
- [ ] Run R statistical analysis
- [ ] Review generated plots and reports
- [ ] Check overlap with HOX clusters
- [ ] Check overlap with centromeric regions
- [ ] Investigate Ig lambda locus
- [ ] Investigate olfactory receptor clusters
- [ ] Analyze graph topology
- [ ] Extract and BLAST duplicate sequences
- [ ] Compare strain-level paths
- [ ] Cross-reference with databases
- [ ] Document findings
- [ ] Decide on next steps
