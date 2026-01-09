# Chr11 Duplicate Analysis Scripts

## Quick Start

### Prerequisites
- Python 3.6+ with pandas, numpy, matplotlib, seaborn
- R with tidyverse, ggplot2, viridis, gridExtra

### Install Python dependencies:
```bash
pip install pandas numpy matplotlib seaborn
```

### Install R dependencies:
```R
install.packages(c("tidyverse", "ggplot2", "scales", "viridis", "gridExtra"))
```

## Usage

### 1. Python Analysis (Primary)
```bash
python script/analyze_chr11_nodes.py DBmm10.nodes.tag.txt.gz
```

**Outputs:**
- `chr11_analysis_overview.png` - Overview statistics and distributions
- `chr11_analysis_position.png` - Positional distribution along chromosome
- `chr11_analysis_report.txt` - Detailed text report

**Features:**
- Filters and analyzes chr11 nodes
- Identifies duplicate node patterns
- Checks overlap with:
  - HOX gene clusters
  - Centromeric regions
  - Immunoglobulin lambda locus
  - Olfactory receptor clusters
- Generates comprehensive visualizations

### 2. R Statistical Analysis (Complementary)
```bash
Rscript script/chr11_statistical_analysis.R DBmm10.nodes.tag.txt.gz
```

**Outputs:**
- `chr11_detailed_analysis.png` - Statistical plots and analyses

**Features:**
- Statistical tests (Kruskal-Wallis)
- Duplicate clustering analysis
- Duplication hotspot identification
- Publication-quality figures

## Expected Input Format

The input file should be a CSV (can be gzipped) with these columns:
```
nodeID,type,length,chr,pos,tag
1,mm10,16299,chrM,0,UNIQUE
2,B,16299,NA,NA,UNPLACED
...
```

**Column descriptions:**
- `nodeID`: Unique node identifier
- `type`: Node type (mm10, B, D, etc.)
- `length`: Node length in base pairs
- `chr`: Chromosome (chr11, chr1, etc., or NA)
- `pos`: Position on chromosome (or NA)
- `tag`: Classification (UNIQUE, DUPLICATED, UNPLACED)

## Interpreting Results

### High Duplication Regions to Investigate:

1. **~113.5-115 Mb**: Immunoglobulin lambda locus
   - Expected: High duplication due to V/J/C segments
   - Action: Verify sequences are Ig-related

2. **~49-52 Mb, ~58-60 Mb**: Olfactory receptor clusters
   - Expected: Copy number variation
   - Action: Check for assembly artifacts

3. **~96.2 Mb**: HOX-A cluster
   - Expected: Some complexity, but not extreme
   - Action: Look for paralogous mappings from other HOX clusters

4. **0-3 Mb**: Centromeric region
   - Expected: Repeats, potential artifacts
   - Action: Consider filtering if problematic

### Red Flags:

- Uniform distribution of duplicates (suggests technical issue)
- Extremely high node counts in non-repeat regions
- Duplicates only in specific strains (assembly problem)
- Sequences mapping to other chromosomes (mis-mapping)

## Troubleshooting

### "No chr11 nodes found"
- Check that your input file has `chr11` in the chr column
- Verify file format matches expected structure
- Try: `zcat DBmm10.nodes.tag.txt.gz | grep chr11 | head`

### Memory issues
- Large datasets may require more RAM
- Consider analyzing subsets or increasing swap space

### Missing visualizations
- Ensure matplotlib/ggplot2 are installed correctly
- Check that output directory is writable
- Look for error messages in console output

## Further Analysis

See `chr11_investigation.md` for:
- Detailed investigation strategy
- Graph topology analysis
- Sequence-level validation
- Cross-referencing with databases
- Next steps based on findings

## Contact

For issues or questions about these scripts, please open an issue in the repository.
