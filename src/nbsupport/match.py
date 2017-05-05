"""Functions for matching RNA insertions to DNA insertions."""

import itertools

import matplotlib_venn
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import toolz

from ngstk.tabix import GtfIterator

METHOD_PALETTE = dict(
    zip(['IM-Fusion-specific', 'Shared', 'ShearSplink-specific'],
        sns.color_palette()))


def match_rnaseq_insertions(rna_insertions,
                            dna_insertions,
                            keep_unmatched=False,
                            with_gene_name=False):
    """Matches rnaseq insertions to compatible dnaseq insertions."""

    # Generate matches.
    match_gen = _match_rnaseq_gen(rna_insertions, dna_insertions)

    # Convert to DataFrame.
    match_rows = itertools.chain.from_iterable(
        zip(
            itertools.cycle([rna_id]), dna_ids,
            itertools.cycle([sample]), itertools.cycle([gene_id]))
        for rna_id, sample, gene_id, dna_ids in match_gen)

    matches = pd.DataFrame.from_records(
        match_rows, columns=['rna_id', 'dna_id', 'sample', 'gene_id'])

    # Add unmatched insertions to frame, if requested.
    if keep_unmatched:

        # Identify unmatched insertions.
        unmatched_rna = (set(zip(rna_insertions.id, rna_insertions.gene_id)) -
                         set(zip(matches.rna_id, matches.gene_id)))

        unmatched_dna = (set(zip(dna_insertions.id, dna_insertions.gene_id)) -
                         set(zip(matches.dna_id, matches.gene_id)))

        # Convert to frames.
        unmatched_rna_df = _extract_unmatched(
            rna_insertions, unmatched_rna, is_rna=True)

        unmatched_dna_df = _extract_unmatched(
            dna_insertions, unmatched_dna, is_rna=False)

        # Combine frames.
        matches = pd.concat(
            [matches, unmatched_rna_df, unmatched_dna_df],
            axis=0,
            ignore_index=True)
        matches = matches[['rna_id', 'dna_id', 'sample', 'gene_id']]

    if with_gene_name:
        rna_mapping = dict(zip(rna_insertions['gene_id'],
                               rna_insertions['gene_name']))
        dna_mapping = dict(zip(dna_insertions['gene_id'],
                               dna_insertions['gene_name']))
        mapping = toolz.merge(rna_mapping, dna_mapping)

        matches['gene_name'] = matches['gene_id'].map(mapping)

    return matches


def _extract_unmatched(insertions, id_gene_indices, is_rna):
    if is_rna:
        set_id, null_id = 'rna_id', 'dna_id'
    else:
        set_id, null_id = 'dna_id', 'rna_id'

    insertions = insertions.set_index(['id', 'gene_id'])

    unmatched = (
        insertions.loc[list(id_gene_indices)]
        .reset_index()
        .rename(columns={
            'id': set_id
        }).assign(**{null_id: np.nan})
        [['rna_id', 'dna_id', 'sample', 'gene_id']]) # yapf: disable

    return unmatched


def _match_rnaseq_gen(rna_insertions, dna_insertions):
    """Generator matching RNA-seq insertions with compatible DNA insertions."""

    # Group DNA insertions by strand, sample and gene.
    grouped_dna = {
        key: grp
        for key, grp in dna_insertions.groupby(
            ['sample', 'strand', 'gene_id'])
    }

    for ins in rna_insertions.itertuples():
        # Select candidates with same sample, strand and gene.
        key = (ins.sample, ins.strand, ins.gene_id)
        dna_candidates = grouped_dna.get(key, None)

        if dna_candidates is None:
            # Return empty set.
            candidate_ids = set()
        else:
            # Filter for candidates with compatible location.
            dna_candidates = _filter_compatible_location(ins, dna_candidates)
            candidate_ids = set(dna_candidates.id)

        yield ins.id, ins.sample, ins.gene_id, candidate_ids


def _filter_compatible_location(rna_insertion, dna_candidates):
    """Filters DNA insertion candidates with incompatible location."""

    # First we determine which direction to look, correcting
    # for the strand of the involved gene.
    if rna_insertion.feature_type == 'SA':
        # Splice acceptor insertions must be downstream.
        direction = 1
    elif rna_insertion.feature_type == 'SD':
        # SD insertions must be upstream.
        direction = -1
    else:
        raise ValueError('Unknown feature type {}'
                         .format(rna_insertion.feature_type))

    direction *= rna_insertion.gene_strand

    # Then we select candidates that are in the correct relative position.
    dists = dna_candidates.position - rna_insertion.position
    mask = (np.sign(dists) * np.sign(direction)) >= 0

    return dna_candidates.ix[mask]


def annotate_dna_info(matches, dna_insertions):
    """Annotate matches with metadata from the DNA insertion."""

    cols = ['id', 'position', 'support']
    if 'clonality' in dna_insertions.columns:
        cols.append('clonality')

    subset = dna_insertions[cols]
    subset = subset.rename(columns=lambda c: 'dna_' + c)

    return pd.merge(matches, subset.drop_duplicates(), how='left', on='dna_id')


def annotate_rna_info(matches, rna_insertions):
    """Annotate matches with metadata from the RNA insertion."""

    subset = rna_insertions[['id', 'position', 'support']]
    subset = subset.rename(columns=lambda c: 'rna_' + c)

    return pd.merge(matches, subset.drop_duplicates(), how='left', on='rna_id')


def annotate_dna_gene_distance(matches, genes):
    """Annotate matches with distance of DNA insertion to gene."""

    if 'dna_position' not in matches.columns:
        raise ValueError('Matches must be annotated for DNA position first')

    # Calculate minimum distance to gene start/end.
    positions = matches.dna_position.values
    gene_starts = genes.ix[matches['gene_id'], 'start'].values
    gene_ends = genes.ix[matches['gene_id'], 'end'].values

    start_dist = np.abs(gene_starts - positions)
    end_dist = np.abs(gene_ends - positions)
    dist = np.minimum(start_dist, end_dist)

    # Set positions within genes to zero.
    with np.errstate(invalid='ignore'):
        within_mask = (gene_starts <= positions) & (positions < gene_ends)
        dist[within_mask] = 0.0

    return matches.assign(gene_distance=dist)


def annotate_dna_gene_position(matches, genes):
    """Annotates matches with relative position of DNA insertion within gene."""

    if 'dna_position' not in matches.columns:
        raise ValueError('Matches must be annotated for DNA position first')

    # Calculate relative position.
    gene_starts = genes.ix[matches['gene_id'], 'start'].values
    gene_ends = genes.ix[matches['gene_id'], 'end'].values

    gene_sizes = gene_ends - gene_starts

    positions = matches.dna_position.values
    rel_positions = (positions - gene_starts) / gene_sizes

    # Correct for strand of gene.
    is_reverse = (genes.ix[matches['gene_id'], 'strand'] == '-').values
    rel_positions[is_reverse] = 1. - rel_positions[is_reverse]

    # Replace values outside gene with nan.
    with np.errstate(invalid='ignore'):
        rel_positions[rel_positions < 0.] = np.nan
        rel_positions[rel_positions > 1.] = np.nan

    return matches.assign(dna_gene_position=rel_positions)


def annotate_type(matches):
    """Annotate matches with their match type."""

    # Determine masks.
    rna_mask = matches.rna_id.notnull()
    dna_mask = matches.dna_id.notnull()

    # Assign type.
    matches = matches.copy()

    matches['type'] = np.nan
    matches.ix[rna_mask & dna_mask, 'type'] = 'shared'
    matches.ix[rna_mask & ~dna_mask, 'type'] = 'rna-only'
    matches.ix[~rna_mask & dna_mask, 'type'] = 'dna-only'

    return matches


def annotate_expression(matches, expr):
    """Annotates matches with gene expression values."""

    expr_long = pd.melt(
        expr.ix[matches['gene_id'].unique()].reset_index(),
        id_vars='gene_id',
        var_name='sample',
        value_name='expr')

    return pd.merge(matches, expr_long, on=['gene_id', 'sample'], how='left')


def plot_insertion_overlap(matches,
                           ax=None,
                           fontsize=12,
                           label_offset=(50, 80),
                           line_width=1):
    """Plots overlap between different insertion match types."""

    if ax is None:
        _, ax = plt.subplots()

    # Determine set sizes.
    match_counts = matches.type.value_counts()

    set_sizes = {
        '11': match_counts.ix['Shared'],
        '01': match_counts.ix['IM-Fusion-specific'],
        '10': match_counts.ix['ShearSplink-specific']
    }

    # Draw venn diagram.
    venn = matplotlib_venn.venn2(
        set_sizes,
        set_labels=('ShearSplink', 'IM-Fusion'),
        set_colors=(sns.color_palette()[2], sns.color_palette()[0]),
        alpha=0.5)
    venn.get_patch_by_id('11').set_color(sns.color_palette()[1])

    for label in venn.subset_labels:
        label.set_fontsize(fontsize)

    matplotlib_venn.venn2_circles(set_sizes, linewidth=line_width)

    ax.set_title('Insertion overlap')

    # Adjust font sizes.
    for label in ('A', 'B'):
        venn.get_label_by_id(label).set_fontsize(fontsize)

    # Add arrow.
    ax.annotate(
        'Shared',
        xy=venn.get_label_by_id('11').get_position() - np.array([0, -0.05]),
        xytext=label_offset,
        ha='center',
        textcoords='offset points',
        fontsize=fontsize,
        arrowprops=dict(
            arrowstyle='->',
            connectionstyle='arc3,rad=0.2',
            color='black',
            linewidth=1))

    return ax


def plot_dna_position_bias(dna_matches, ax=None):
    """Plots distribution of relative insertion position per match type."""

    ax = _plot_distributions(dna_matches, var='dna_gene_position', ax=ax)

    ax.set_xticks(np.arange(0, 1.1, 0.2))
    ax.set_xlim(-0.2, 1.2)

    ax.set_title('Position bias')
    ax.set_xlabel('Relative position within gene')
    ax.set_ylabel('Density')

    sns.despine(ax=ax)

    return ax


def _plot_distributions(matches,
                        var,
                        ax=None,
                        log=False,
                        distplot_kws=None,
                        legend_kws=None):
    if ax is None:
        _, ax = plt.subplots()

    default_kws = dict(hist=False, kde_kws={'shade': True})
    distplot_kws = toolz.merge(default_kws, distplot_kws or {})

    for ins_type, grp in matches.groupby('type'):
        values = grp[var].dropna()

        if log:
            values = np.log2(values + 1)

        sns.distplot(
            values,
            color=METHOD_PALETTE[ins_type],
            label=ins_type,
            ax=ax,
            **distplot_kws)

        if legend_kws is not None:
            ax.legend(**legend_kws)

    return ax


def plot_dna_depth_bias(dna_matches, ax=None, var='dna_support'):
    """Plots distribution of DNA insertion depths per match type."""
    # TODO: Use clonality.

    ax = _plot_distributions(dna_matches, var=var, ax=ax)

    ax.set_title('Support bias (IM-Fusion)')
    ax.set_xlabel('Insertion support')
    ax.set_ylabel('Density')
    ax.set_xlim(0, None)

    sns.despine(ax=ax)

    return ax


def plot_expression_bias(dna_matches, ax=None, legend_kws=None):
    """Plots distribution of gene expression per match type."""

    ax = _plot_distributions(
        dna_matches, var='expr', ax=ax, legend_kws=legend_kws)

    ax.set_title('Expression bias')
    ax.set_xlabel('Expression')
    ax.set_ylabel('Density')
    ax.set_xlim(0, None)

    sns.despine(ax=ax)

    return ax


def plot_rna_depth_bias(rna_matches, ax=None):
    """Plots distribution of RNA insertion depths per match type."""

    distplot_kws = {'kde_kws': {'shade': True, 'bw': 0.4}}
    ax = _plot_distributions(
        rna_matches,
        var='rna_support',
        ax=ax,
        distplot_kws=distplot_kws,
        log=True)

    ax.set_title('Support bias (ShearSplink)')
    ax.set_xlabel('Insertion support (log2)')
    ax.set_ylabel('Density')
    ax.set_xlim(1, None)

    sns.despine(ax=ax)

    return ax
