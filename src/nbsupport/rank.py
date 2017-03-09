"""Functions for comparing and plotting CTG/hit ranks."""

from collections import OrderedDict

from matplotlib import pyplot as plt
from matplotlib import patches as mpatches
import pandas as pd
import numpy as np
import seaborn as sns
import toolz

from imfusion.expression import test_de
from imfusion.model import Insertion


RANK_TYPES = ['Shared', 'Shared, no CTG', 'Shared, not DE',
              'Shared, no CTG, not DE', 'Not shared']  # yapf: disable


def compare_ranks(ranks_a, ranks_b, suffixes):
    """Compares RNA and DNA hit ranks."""

    # Add ranks if missing.
    ranks_a = _add_rank(ranks_a)
    ranks_b = _add_rank(ranks_b)

    # Merge rankings.
    merged = pd.merge(
        ranks_a[['gene_name', 'rank']],
        ranks_b[['gene_name', 'rank']],
        on='gene_name',
        how='outer',
        suffixes=suffixes)

    return merged


def _add_rank(hits):
    if 'rank' not in hits:
        ranks = np.arange(1, len(hits) + 1)
        hits = hits.assign(rank=ranks)
    return hits


# def annotate_de(ranks,
#                 rna_ctgs,
#                 threshold=0.05,
#                 gene_col='gene_name',
#                 col_name='is_de'):
#     """Annotate ranks for DE status."""

#     annot = rna_ctgs[[gene_col]].copy()
#     annot[col_name] = rna_ctgs['de_pvalue'] < threshold
#     merged = pd.merge(ranks, annot, on=gene_col, how='left')

#     return merged


def annotate_de(ranks,
                insertions,
                exon_counts,
                gene_id_map,
                threshold=0.05,
                col_name='is_de',
                fallback_to_gene=False):
    """Annotate ranks for DE status."""

    # Translate names to ids.
    count_ids = set(exon_counts.index.get_level_values(0))
    gene_ids = [gene_id_map[name] for name in ranks['gene_name']
                if name in gene_id_map]  # yapf: disable

    # Check if we are missing counts for any genes.
    missing_counts = {gene_id for gene_id in gene_ids
                      if gene_id not in count_ids}  # yapf: disable

    if len(missing_counts) > 0:
        print('Warning: missing counts for some genes ({})'
              .format(', '.join(missing_counts)))  #  yapf: disable
        gene_ids = list(set(gene_ids) - missing_counts)

    # Perform DE test.
    de_result = test_de(
        list(Insertion.from_frame(insertions)),
        exon_counts,
        gene_ids=gene_ids,
        fallback_to_gene=fallback_to_gene)

    # Add binary is_de column.
    de_result[col_name] = de_result['p_value'] <= threshold

    # Map gene_ids in result back to names.
    name_map = {v: k for k, v in gene_id_map.items()}
    de_result['gene_name'] = de_result['gene_id'].map(name_map)

    # Merge with ranks.
    merged = pd.merge(ranks, de_result[['gene_name', col_name]],
                      on='gene_name', how='left')  # yapf: disable

    return merged


def annotate_ctg(ranks, ctgs, col_name='is_ctg'):
    """Annotates ranks for CTG status."""
    return ranks.assign(
        **{col_name: ranks['gene_name'].isin(ctgs['gene_name'])})


def annotate_ins(ranks,
                 insertions,
                 gene_col='gene_name',
                 col_name='in_insertions'):
    """Annotates ranks for presence in insertions."""

    annot = ranks[[gene_col]].copy()
    annot[col_name] = annot[gene_col].isin(set(insertions[gene_col]))
    merged = pd.merge(ranks, annot, on=gene_col, how='left')

    return merged


def annotate_type(ranks):
    """Annotates ranks with their match type.

    Expects ranks to be annotated with the following columns: ctg_rna,
    ctg_dna, ins_rna, ins_dna and is_de. These annotations are added using
    the annotate_ctg, annotate_ins and annotate_de functions respectively.

    """

    ranks = ranks.copy()
    ranks['type'] = ranks.apply(_rank_type, axis=1)

    return ranks


def _rank_type(row):
    if row.ctg_rna and row.ctg_dna:
        if row.is_de:
            return 'Shared'
        else:
            return 'Shared, not DE'
    else:
        if (row.ctg_rna and row.ins_dna) or (row.ctg_dna and row.ins_rna):
            if row.is_de:
                return 'Shared, no CTG'
            else:
                return 'Shared, no CTG, not DE'
        else:
            return 'Not shared'


def plot_ranks(data,
               rank_a,
               rank_b,
               label,
               hue=None,
               palette=None,
               order=None,
               ax=None,
               offset=5,
               xlim=(-0.5, 1.5),
               label_kws=None,
               legend_kws=None):
    """Plots ranks using matplotlib."""

    if ax is None:
        _, ax = plt.subplots()

    color_map = _build_colormap(data, hue, palette, order=order)

    default_label_kws = {'va': 'center', 'textcoords': 'offset points'}
    label_kws = toolz.merge(default_label_kws, label_kws or {})

    # Plot lines.
    y_values = data[[rank_a, rank_b]].dropna().T
    ax.plot([0, 1], y_values, color='lightgrey')

    # Plot labels.
    for row in data.itertuples():
        if hue is not None:
            color = color_map.get(getattr(row, hue), 'black')
        else:
            color = 'black'

        ax.annotate(
            getattr(row, label),
            xy=(0, getattr(row, rank_a)),
            xytext=(-offset, 0),
            ha='right',
            color=color,
            **label_kws)

        ax.annotate(
            getattr(row, label),
            xy=(1, getattr(row, rank_b)),
            xytext=(offset, 0),
            ha='left',
            color=color,
            **label_kws)

    # Style axes.
    ax.set_xticks([0, 1])
    ax.set_xticklabels([rank_a, rank_b])

    ax.set_ylabel('Rank')
    ax.set_yticks([1] + list(np.arange(5, len(data), 5)))

    # Set lims.
    max_rank = max(data[rank_a].max(), data[rank_b].max())
    ax.set_ylim(max_rank + 1, 0)
    ax.set_xlim(*xlim)

    # Draw legend.
    if hue is not None:
        _draw_legend(color_map, ax=ax, **(legend_kws or {}))

    sns.despine(ax=ax)


def _build_colormap(data, hue, palette, order):
    """Builds a colormap."""

    if hue is None:
        color_map = {}
    else:
        if palette is None:
            palette = sns.color_palette()

        if order is None:
            order = data[hue].unique()

        color_map = OrderedDict(zip(order, palette))

    return color_map


def _draw_legend(color_map, ax, **kwargs):
    """Draws a (patch) legend for the given color map."""
    patches = [mpatches.Patch(color=color, label=label)
               for label, color in color_map.items()]  # yapf: disable
    return ax.legend(handles=patches, **kwargs)
