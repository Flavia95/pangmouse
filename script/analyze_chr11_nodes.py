#!/usr/bin/env python3
"""
Analysis script for chr11 pangenome nodes
Investigates duplicates, HOX cluster overlap, and centromeric regions

Usage:
    python analyze_chr11_nodes.py DBmm10.nodes.tag.txt.gz
"""

import sys
import gzip
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter
import warnings
warnings.filterwarnings('ignore')

# Mouse chromosome 11 reference information (mm10)
CHR11_LENGTH = 121843856  # bp

# Mouse HOX clusters on chr11 (mm10 coordinates)
# Based on NCBI/UCSC genome browser data
HOX_CLUSTERS = {
    'Hoxa': {'chr': 'chr11', 'start': 96161468, 'end': 96273464},
    # Note: Hoxb, Hoxc, Hoxd are on other chromosomes
}

# Centromeric region for chr11 (approximate, based on gap regions in mm10)
# Centromeres in mouse are typically not well-defined, but repetitive regions exist
CENTROMERE_REGION = {
    'chr11': {'start': 0, 'end': 3000000}  # First 3Mb typically contains centromeric/pericentromeric regions
}

# Paralogous gene clusters on chr11 that might contribute to duplications
KNOWN_DUPLICATED_REGIONS = {
    'immunoglobulin_lambda_locus': {'start': 113500000, 'end': 115000000},
    'olfactory_receptors_1': {'start': 49000000, 'end': 52000000},
    'olfactory_receptors_2': {'start': 58000000, 'end': 60000000},
}


def load_nodes_data(filename):
    """Load and parse the nodes data file"""
    print(f"Loading data from {filename}...")

    if filename.endswith('.gz'):
        with gzip.open(filename, 'rt') as f:
            df = pd.read_csv(f)
    else:
        df = pd.read_csv(filename)

    print(f"Total nodes loaded: {len(df):,}")
    print(f"Columns: {list(df.columns)}")
    return df


def filter_chr11_nodes(df):
    """Filter nodes belonging to chr11"""
    chr11_df = df[df['chr'] == 'chr11'].copy()
    print(f"\nChr11 nodes: {len(chr11_df):,}")
    return chr11_df


def analyze_node_types(df):
    """Analyze the distribution of node types"""
    print("\n" + "="*60)
    print("NODE TYPE DISTRIBUTION")
    print("="*60)

    type_counts = df['type'].value_counts()
    print(type_counts)

    # Tag distribution
    if 'tag' in df.columns:
        print("\nTag distribution:")
        tag_counts = df['tag'].value_counts()
        print(tag_counts)

        # Calculate percentages
        total = len(df)
        for tag, count in tag_counts.items():
            pct = (count / total) * 100
            print(f"  {tag}: {count:,} ({pct:.2f}%)")

    return type_counts


def analyze_duplicates(df):
    """Analyze duplicated nodes"""
    print("\n" + "="*60)
    print("DUPLICATE NODE ANALYSIS")
    print("="*60)

    # Filter for duplicated nodes
    dup_df = df[df['tag'] == 'DUPLICATED'].copy() if 'tag' in df.columns else df[df['type'] == 'D'].copy()

    print(f"Total duplicated nodes: {len(dup_df):,}")
    print(f"Percentage of chr11 nodes: {(len(dup_df)/len(df)*100):.2f}%")

    if len(dup_df) > 0:
        print(f"\nDuplicated node length statistics:")
        print(f"  Mean length: {dup_df['length'].mean():.0f} bp")
        print(f"  Median length: {dup_df['length'].median():.0f} bp")
        print(f"  Total duplicated sequence: {dup_df['length'].sum()/1e6:.2f} Mb")
        print(f"  Min length: {dup_df['length'].min()} bp")
        print(f"  Max length: {dup_df['length'].max()} bp")

    return dup_df


def analyze_position_distribution(df, bin_size=1000000):
    """Analyze the distribution of nodes along chr11"""
    print("\n" + "="*60)
    print("POSITIONAL DISTRIBUTION ANALYSIS")
    print("="*60)

    # Remove nodes without position information
    pos_df = df[df['pos'].notna()].copy()
    print(f"Nodes with position information: {len(pos_df):,}")

    if len(pos_df) == 0:
        print("No positional information available")
        return None

    # Create bins
    bins = np.arange(0, CHR11_LENGTH + bin_size, bin_size)
    pos_df['bin'] = pd.cut(pos_df['pos'], bins=bins, labels=bins[:-1]/1e6)

    # Count nodes per bin
    bin_counts = pos_df.groupby('bin').size()

    # Separate by tag
    if 'tag' in pos_df.columns:
        tag_counts = pos_df.groupby(['bin', 'tag']).size().unstack(fill_value=0)
        print("\nNodes per Mb bin (first 10 bins):")
        print(tag_counts.head(10))

        # Identify bins with high duplicate content
        if 'DUPLICATED' in tag_counts.columns:
            high_dup_bins = tag_counts[tag_counts['DUPLICATED'] > tag_counts['DUPLICATED'].quantile(0.9)]
            if len(high_dup_bins) > 0:
                print(f"\nRegions with highest duplication (>90th percentile):")
                for bin_mb in high_dup_bins.index:
                    dup_count = high_dup_bins.loc[bin_mb, 'DUPLICATED']
                    bin_start = float(bin_mb)
                    bin_end = bin_start + (bin_size/1e6)
                    print(f"  {bin_start:.1f}-{bin_end:.1f} Mb: {dup_count} duplicated nodes")

    return pos_df


def check_hox_overlap(df):
    """Check for overlap with HOX clusters"""
    print("\n" + "="*60)
    print("HOX CLUSTER OVERLAP ANALYSIS")
    print("="*60)

    pos_df = df[df['pos'].notna()].copy()

    if len(pos_df) == 0:
        print("No positional information available")
        return

    for cluster_name, coords in HOX_CLUSTERS.items():
        # Find nodes overlapping with HOX cluster
        overlapping = pos_df[
            (pos_df['pos'] >= coords['start']) &
            (pos_df['pos'] <= coords['end'])
        ]

        print(f"\n{cluster_name} cluster ({coords['start']/1e6:.2f}-{coords['end']/1e6:.2f} Mb):")
        print(f"  Total nodes: {len(overlapping)}")

        if len(overlapping) > 0 and 'tag' in overlapping.columns:
            tag_dist = overlapping['tag'].value_counts()
            for tag, count in tag_dist.items():
                print(f"    {tag}: {count}")

            # Check for paralogous regions nearby (within 5 Mb)
            nearby_start = max(0, coords['start'] - 5e6)
            nearby_end = coords['end'] + 5e6
            nearby = pos_df[
                (pos_df['pos'] >= nearby_start) &
                (pos_df['pos'] <= nearby_end)
            ]
            dup_nearby = nearby[nearby['tag'] == 'DUPLICATED'] if 'tag' in nearby.columns else pd.DataFrame()

            if len(dup_nearby) > 0:
                print(f"  Duplicated nodes within ±5 Mb: {len(dup_nearby)}")


def check_centromere_overlap(df):
    """Check for overlap with centromeric regions"""
    print("\n" + "="*60)
    print("CENTROMERIC REGION OVERLAP ANALYSIS")
    print("="*60)

    pos_df = df[df['pos'].notna()].copy()

    if len(pos_df) == 0:
        print("No positional information available")
        return

    centro = CENTROMERE_REGION['chr11']
    overlapping = pos_df[
        (pos_df['pos'] >= centro['start']) &
        (pos_df['pos'] <= centro['end'])
    ]

    print(f"Centromeric/pericentromeric region (0-{centro['end']/1e6:.1f} Mb):")
    print(f"  Total nodes: {len(overlapping)}")

    if len(overlapping) > 0 and 'tag' in overlapping.columns:
        tag_dist = overlapping['tag'].value_counts()
        for tag, count in tag_dist.items():
            pct = (count / len(overlapping)) * 100
            print(f"    {tag}: {count} ({pct:.1f}%)")

        # Calculate sequence content
        total_seq = overlapping['length'].sum() / 1e6
        print(f"  Total sequence: {total_seq:.2f} Mb")


def check_known_duplicated_regions(df):
    """Check overlap with known duplicated/paralogous regions"""
    print("\n" + "="*60)
    print("KNOWN DUPLICATED REGION ANALYSIS")
    print("="*60)

    pos_df = df[df['pos'].notna()].copy()

    if len(pos_df) == 0:
        print("No positional information available")
        return

    for region_name, coords in KNOWN_DUPLICATED_REGIONS.items():
        overlapping = pos_df[
            (pos_df['pos'] >= coords['start']) &
            (pos_df['pos'] <= coords['end'])
        ]

        print(f"\n{region_name.replace('_', ' ').title()} ({coords['start']/1e6:.1f}-{coords['end']/1e6:.1f} Mb):")
        print(f"  Total nodes: {len(overlapping)}")

        if len(overlapping) > 0 and 'tag' in overlapping.columns:
            tag_dist = overlapping['tag'].value_counts()
            for tag, count in tag_dist.items():
                pct = (count / len(overlapping)) * 100
                print(f"    {tag}: {count} ({pct:.1f}%)")


def create_visualizations(chr11_df, dup_df, pos_df, output_prefix="chr11_analysis"):
    """Create comprehensive visualizations"""
    print("\n" + "="*60)
    print("GENERATING VISUALIZATIONS")
    print("="*60)

    # Set style
    sns.set_style("whitegrid")

    # Figure 1: Overview statistics
    fig, axes = plt.subplots(2, 2, figsize=(15, 12))

    # 1.1: Node type distribution
    if 'tag' in chr11_df.columns:
        tag_counts = chr11_df['tag'].value_counts()
        axes[0, 0].bar(range(len(tag_counts)), tag_counts.values, color=['#E15759', '#4E79A7', '#BAB0AC'])
        axes[0, 0].set_xticks(range(len(tag_counts)))
        axes[0, 0].set_xticklabels(tag_counts.index, rotation=45, ha='right')
        axes[0, 0].set_ylabel('Count', fontsize=12)
        axes[0, 0].set_title('Node Type Distribution (Chr11)', fontsize=14, fontweight='bold')
        axes[0, 0].set_yscale('log')

        # Add count labels
        for i, v in enumerate(tag_counts.values):
            axes[0, 0].text(i, v, f'{v:,}', ha='center', va='bottom', fontsize=10)

    # 1.2: Node length distribution
    axes[0, 1].hist(chr11_df['length'], bins=50, color='#76B7B2', edgecolor='black', alpha=0.7)
    axes[0, 1].set_xlabel('Node Length (bp)', fontsize=12)
    axes[0, 1].set_ylabel('Frequency', fontsize=12)
    axes[0, 1].set_title('Node Length Distribution', fontsize=14, fontweight='bold')
    axes[0, 1].set_yscale('log')

    # 1.3: Duplicate node length distribution
    if len(dup_df) > 0:
        axes[1, 0].hist(dup_df['length'], bins=50, color='#E15759', edgecolor='black', alpha=0.7)
        axes[1, 0].set_xlabel('Node Length (bp)', fontsize=12)
        axes[1, 0].set_ylabel('Frequency', fontsize=12)
        axes[1, 0].set_title('Duplicated Node Length Distribution', fontsize=14, fontweight='bold')
        axes[1, 0].set_yscale('log')

    # 1.4: Cumulative sequence length by type
    if 'tag' in chr11_df.columns:
        seq_by_tag = chr11_df.groupby('tag')['length'].sum() / 1e6
        axes[1, 1].bar(range(len(seq_by_tag)), seq_by_tag.values, color=['#E15759', '#4E79A7', '#BAB0AC'])
        axes[1, 1].set_xticks(range(len(seq_by_tag)))
        axes[1, 1].set_xticklabels(seq_by_tag.index, rotation=45, ha='right')
        axes[1, 1].set_ylabel('Total Sequence (Mb)', fontsize=12)
        axes[1, 1].set_title('Cumulative Sequence Length by Type', fontsize=14, fontweight='bold')

        # Add value labels
        for i, v in enumerate(seq_by_tag.values):
            axes[1, 1].text(i, v, f'{v:.1f} Mb', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    plt.savefig(f'{output_prefix}_overview.png', dpi=300, bbox_inches='tight')
    print(f"Saved: {output_prefix}_overview.png")
    plt.close()

    # Figure 2: Positional distribution
    if pos_df is not None and len(pos_df) > 0:
        fig, axes = plt.subplots(2, 1, figsize=(15, 10))

        # 2.1: All nodes distribution along chr11
        bin_size = 1000000  # 1 Mb bins
        bins = np.arange(0, CHR11_LENGTH + bin_size, bin_size)

        if 'tag' in pos_df.columns:
            for tag in pos_df['tag'].unique():
                tag_df = pos_df[pos_df['tag'] == tag]
                axes[0].hist(tag_df['pos'] / 1e6, bins=bins/1e6, alpha=0.6, label=tag)
        else:
            axes[0].hist(pos_df['pos'] / 1e6, bins=bins/1e6, alpha=0.6)

        axes[0].set_xlabel('Position on Chr11 (Mb)', fontsize=12)
        axes[0].set_ylabel('Node Count', fontsize=12)
        axes[0].set_title('Node Distribution Along Chr11', fontsize=14, fontweight='bold')
        axes[0].legend()
        axes[0].set_yscale('log')

        # Add vertical lines for important regions
        # HOX cluster
        for cluster_name, coords in HOX_CLUSTERS.items():
            axes[0].axvline(coords['start']/1e6, color='red', linestyle='--', alpha=0.5, linewidth=2)
            axes[0].axvline(coords['end']/1e6, color='red', linestyle='--', alpha=0.5, linewidth=2)
            axes[0].text((coords['start']+coords['end'])/2e6, axes[0].get_ylim()[1]*0.9,
                        cluster_name, ha='center', fontsize=10, color='red', fontweight='bold')

        # Centromere
        axes[0].axvspan(0, CENTROMERE_REGION['chr11']['end']/1e6, alpha=0.2, color='gray', label='Centromere')

        # 2.2: Duplicated nodes density
        if len(dup_df) > 0 and 'pos' in dup_df.columns:
            dup_pos = dup_df[dup_df['pos'].notna()]
            if len(dup_pos) > 0:
                axes[1].hist(dup_pos['pos'] / 1e6, bins=bins/1e6, color='#E15759', alpha=0.7)
                axes[1].set_xlabel('Position on Chr11 (Mb)', fontsize=12)
                axes[1].set_ylabel('Duplicated Node Count', fontsize=12)
                axes[1].set_title('Duplicated Nodes Distribution Along Chr11', fontsize=14, fontweight='bold')

                # Add vertical lines for known duplicated regions
                for region_name, coords in KNOWN_DUPLICATED_REGIONS.items():
                    axes[1].axvspan(coords['start']/1e6, coords['end']/1e6, alpha=0.2, color='orange')
                    axes[1].text((coords['start']+coords['end'])/2e6, axes[1].get_ylim()[1]*0.9,
                               region_name.replace('_', '\n'), ha='center', fontsize=8)

        plt.tight_layout()
        plt.savefig(f'{output_prefix}_position.png', dpi=300, bbox_inches='tight')
        print(f"Saved: {output_prefix}_position.png")
        plt.close()


def generate_report(chr11_df, dup_df, output_file="chr11_analysis_report.txt"):
    """Generate a text report with findings"""
    print("\n" + "="*60)
    print("GENERATING REPORT")
    print("="*60)

    with open(output_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("CHR11 PANGENOME NODE ANALYSIS REPORT\n")
        f.write("="*80 + "\n\n")

        f.write("SUMMARY STATISTICS\n")
        f.write("-"*80 + "\n")
        f.write(f"Total chr11 nodes: {len(chr11_df):,}\n")
        f.write(f"Total sequence length: {chr11_df['length'].sum()/1e6:.2f} Mb\n")
        f.write(f"Mean node length: {chr11_df['length'].mean():.0f} bp\n")
        f.write(f"Median node length: {chr11_df['length'].median():.0f} bp\n\n")

        if 'tag' in chr11_df.columns:
            f.write("NODE TYPE DISTRIBUTION\n")
            f.write("-"*80 + "\n")
            tag_counts = chr11_df['tag'].value_counts()
            for tag, count in tag_counts.items():
                pct = (count / len(chr11_df)) * 100
                seq_len = chr11_df[chr11_df['tag'] == tag]['length'].sum() / 1e6
                f.write(f"{tag}:\n")
                f.write(f"  Count: {count:,} ({pct:.2f}%)\n")
                f.write(f"  Total sequence: {seq_len:.2f} Mb\n")
            f.write("\n")

        if len(dup_df) > 0:
            f.write("DUPLICATE ANALYSIS\n")
            f.write("-"*80 + "\n")
            f.write(f"Duplicated nodes: {len(dup_df):,}\n")
            f.write(f"Percentage of chr11: {(len(dup_df)/len(chr11_df)*100):.2f}%\n")
            f.write(f"Total duplicated sequence: {dup_df['length'].sum()/1e6:.2f} Mb\n")
            f.write(f"Mean duplicate length: {dup_df['length'].mean():.0f} bp\n\n")

        f.write("KEY FINDINGS AND RECOMMENDATIONS\n")
        f.write("-"*80 + "\n")
        f.write("1. Chr11 shows a high number of nodes compared to other chromosomes\n")
        f.write("2. Check if duplicates cluster in specific genomic regions\n")
        f.write("3. Investigate overlap with:\n")
        f.write("   - HOX gene clusters (Hoxa at ~96.2 Mb)\n")
        f.write("   - Centromeric/pericentromeric regions (0-3 Mb)\n")
        f.write("   - Immunoglobulin lambda locus (~113.5-115 Mb)\n")
        f.write("   - Olfactory receptor clusters (~49-52 Mb, 58-60 Mb)\n")
        f.write("4. Consider graph topology analysis to understand duplication patterns\n")
        f.write("5. Verify pangenome assembly parameters (especially -s, -p values in pggb)\n")

    print(f"Report saved: {output_file}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_chr11_nodes.py <DBmm10.nodes.tag.txt.gz>")
        sys.exit(1)

    input_file = sys.argv[1]

    print("="*80)
    print("CHR11 PANGENOME NODE ANALYSIS")
    print("="*80)

    # Load data
    df = load_nodes_data(input_file)

    # Filter for chr11
    chr11_df = filter_chr11_nodes(df)

    if len(chr11_df) == 0:
        print("ERROR: No chr11 nodes found in the data!")
        sys.exit(1)

    # Analyze node types
    analyze_node_types(chr11_df)

    # Analyze duplicates
    dup_df = analyze_duplicates(chr11_df)

    # Analyze position distribution
    pos_df = analyze_position_distribution(chr11_df)

    # Check for HOX overlap
    check_hox_overlap(chr11_df)

    # Check for centromere overlap
    check_centromere_overlap(chr11_df)

    # Check known duplicated regions
    check_known_duplicated_regions(chr11_df)

    # Create visualizations
    create_visualizations(chr11_df, dup_df, pos_df)

    # Generate report
    generate_report(chr11_df, dup_df)

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("\nGenerated files:")
    print("  - chr11_analysis_overview.png")
    print("  - chr11_analysis_position.png")
    print("  - chr11_analysis_report.txt")


if __name__ == "__main__":
    main()
