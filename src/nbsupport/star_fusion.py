import pandas as pd


def read_fusions(file_path):
    """Reads fusions from STAR-Fusion final candidate file."""

    col_names = [
        'fusion_name', 'support_junction', 'support_spanning', 'splice_type',
        'left_gene', 'left_breakpoint', 'right_gene', 'right_breakpoint',
        'junction_reads', 'spanning_reads', 'extra'
    ]

    return (pd.read_csv(
        str(file_path),
        sep='\t',
        header=None,
        skiprows=1,
        index_col=None,
        names=col_names).drop(
            'extra', axis=1))


def expand_fusion_name(fusions):
    """Splits fusion name into two additional gene columns."""

    expanded = fusions['fusion_name'].str.split('--', expand=True)
    expanded.columns = ['gene_a', 'gene_b']

    return pd.concat([fusions, expanded], axis=1)


def filter_sb_transposon(fusions):
    """Filters fusions with genes related to transposon sequences."""
    mask = fusions['fusion_name'].str.contains('En2|Foxf2')
    return fusions.loc[~mask]
