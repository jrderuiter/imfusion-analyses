import itertools
import operator

from matplotlib import pyplot as plt
import numpy as np
import toolz
import seaborn as sns

from geneviz.tracks import (plot_tracks, FeatureTrack, BiomartTrack,
                            SpliceTrack, CoverageTrack)

from imfusion.expression import test_de_exon, test_de_exon_single
from imfusion.model import Insertion

from .util import read_star_junctions, tidy_shared_axes


def annotate_insertions(insertions, window, genes, select_closest=False):
    """Annotates insertions for genes within given window."""

    annotated = _annotate_insertions(insertions, window, genes)

    if select_closest:
        annotated = _select_closest_gene(annotated)

    yield from annotated


def _annotate_insertions(insertions, window, genes):
    for insertion in insertions:
        start = insertion.position - window
        end = insertion.position + window

        hits = genes.search(insertion.seqname, start, end)

        if len(hits) > 0:
            for hit in hits.itertuples():
                gene_meta = {
                    'gene_name': hit.gene_name,
                    'gene_id': hit.gene_id,
                    'gene_distance': _gene_distance(insertion, hit)
                }
                new_meta = toolz.merge(insertion.metadata, gene_meta)
                yield insertion._replace(metadata=new_meta)
        else:
            yield insertion


def _gene_distance(insertion, gene):
    """Calculates distance between insertion and gene."""

    position = insertion.position
    gene_start = gene.start
    gene_end = gene.end

    if gene_start <= position < gene_end:
        return 0
    elif position < gene_start:
        return gene_start - position
    else:
        return position - gene_end


def _select_closest_gene(insertions):
    """Selects closest annotated gene (returns multiple in case of tie)."""

    for _, grp in itertools.groupby(insertions, key=operator.attrgetter('id')):
        grp = list(grp)

        if len(grp) > 1:
            dists = np.array([ins.metadata['gene_distance'] for ins in grp])
            min_indices = np.where(dists == dists.min())[0]

            for ind in min_indices:
                yield grp[ind]
        else:
            yield grp[0]


def plot_insertion(insertion,
                   junction_path,
                   bam_path,
                   region,
                   reverse=False,
                   figsize=(12, None),
                   fontsize=12,
                   cov_resample_interval=100):
    """Visualizes insertion together with coverage/junction information."""

    # Coverage track.
    cov_track = CoverageTrack(
        bam_path, height=1.5, resample_interval=cov_resample_interval)

    # Splice track.
    junctions = read_star_junctions(junction_path)
    junctions = junctions.rename(columns={'unique_reads': 'score'})

    splice_track = SpliceTrack(
        data=junctions, patch_kws={'facecolor': 'white'}, height=1.5)

    # Insertion track.
    palette = [sns.color_palette()[0], sns.color_palette()[2]]

    ins_data = Insertion.to_frame([insertion])
    ins_data = ins_data.rename(columns={'seqname': 'chromosome'})

    width = (region[2] - region[1]) / 25

    feature_track = FeatureTrack.from_position(
        ins_data,
        width=width,
        height=0.5,
        hue='orientation',
        palette=palette,
        hue_order=['sense', 'antisense'],
        patch_kws={'linewidth': 0.5, 'edgecolor': 'black'})

    # Gene track.
    gene_track = BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        height=0.5,
        collapse='transcript',
        gene_id='gene_name',
        filter='gene_name == {!r}'.format(insertion.metadata['gene_name']),
        label_kws={'fontsize': fontsize,
                   'fontstyle': 'italic'})

    # Draw tracks.
    g = plot_tracks(
        [cov_track, splice_track, feature_track, gene_track],
        region=region,
        figsize=figsize,
        reverse=reverse)

    # Set ylabels for coverage/splice tracks.
    g.axes[0].set_ylabel('Coverage')
    g.axes[1].set_ylabel('Junctions')

    # Set yticks for coverage/splice tracks.
    ylim = g.axes[0].get_ylim()[1]
    g.axes[0].set_yticks(np.arange(ylim, step=100))

    ylim = g.axes[1].get_ylim()[1]
    g.axes[1].set_yticks(np.arange(50, ylim, step=50))

    # Flip splice track.
    g.axes[1].set_ylim(g.axes[1].get_ylim()[::-1])

    # Despine axes.
    sns.despine(ax=g.axes[0], left=False, bottom=False)
    sns.despine(ax=g.axes[1], left=False, bottom=True)
    sns.despine(ax=g.axes[2], left=True, bottom=True)
    sns.despine(ax=g.axes[3], left=True, bottom=False)

    # Reduce spacing between tracks.
    plt.subplots_adjust(hspace=0.01)

    return g


def plot_single_de(insertions, exon_counts, insertion, log=False, ax=None):
    """Plots DE for a single insertion."""

    # TODO: Add legend for colors?

    # Calculate DE.
    result = test_de_exon_single(insertions, exon_counts, insertion.id,
                                 insertion.metadata['gene_id'])

    # Plot result.
    ax = result.plot_sums(log=log, ax=ax)
    ax.set_title(insertion.metadata['gene_name'])

    sns.despine(ax=ax)

    return ax


def plot_feature_bias(insertions, genes=None, ax=None):
    """Plots frequency of transposon features per gene as barplot."""

    if ax is None:
        _, ax = plt.subplots()

    # Count feature occurrences per gene.
    feature_counts = (
        insertions.groupby(['gene_name', 'feature_name'])['sample'].nunique()
        .reset_index(name='count'))

    # Create barplot.
    sns.barplot(
        data=feature_counts,
        x='gene_name',
        y='count',
        order=genes,
        hue='feature_name',
        hue_order=['SD', 'SA', 'En2SA'],
        palette=['#86b0df', '#ef858f', '#7cba92'],
        ax=ax)

    ax.legend(title='Transposon feature')

    ax.set_xlabel('Gene')
    ax.set_ylabel('Count')

    plt.setp(ax.get_xticklabels(), fontstyle='italic')

    sns.despine(ax=ax)

    return ax


def plot_de_multiple(insertions,
                     exon_counts,
                     gene_ids,
                     gene_names=None,
                     figsize=(4, 3),
                     ncols=3,
                     **kwargs):
    """Plots groupwise DE for multiple genes."""

    nrows = len(gene_ids) // ncols
    if len(gene_ids) % ncols != 0:
        nrows += 1
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows, figsize=figsize)

    if gene_names is None:
        gene_names = gene_ids

    for gene_id, gene_name, ax in zip(gene_ids, gene_names, axes.flatten()):
        plot_de(insertions, exon_counts, gene_id,
                gene_name=gene_name, ax=ax, **kwargs)

    tidy_shared_axes(fig, axes, ylabel_x=-0.01)

    fig.tight_layout()
    sns.despine(fig)

    return fig


def plot_de(insertions,
            exon_counts,
            gene_id,
            gene_name=None,
            ax=None,
            box_kws=None,
            **kwargs):
    """Plots groupwise DE for a given gene."""

    if gene_name is None:
        gene_name = gene_id

    result = test_de_exon(insertions, exon_counts, gene_id=gene_id, **kwargs)

    ax = result.plot_boxplot(
        ax=ax, log=True, show_points=False,
        box_kws=toolz.merge({'color': 'lightgrey'}, box_kws or {}))

    pval_label = _format_pval(result.p_value)
    ax.set_title('{} ({})'.format(gene_name, pval_label), fontstyle='italic')

    return ax


def _format_pval(p_value):
    """Helper function for formatting p-value."""
    if p_value < 0.0001:
        return 'p < 0.0001'
    else:
        return 'p = {:5.4f}'.format(p_value)


def plot_insertion_track(insertions,
                         region,
                         gene=None,
                         transcript_id=None,
                         ins_ratio=1 / 25,
                         linewidth=0.5,
                         bm_host='http://www.ensembl.org',
                         bm_gene_name='external_gene_name',
                         **kwargs):
    ori_order = ['sense', 'antisense']
    ori_palette = [sns.color_palette()[0], sns.color_palette()[2]]

    width = (region[2] - region[1]) * ins_ratio

    if 'gene_orientation' in insertions:
        hue = 'gene_orientation'
    else:
        hue = None

    ins_track = FeatureTrack.from_position(
        data=insertions,
        width=width,
        height=0.25,
        hue=hue,
        hue_order=ori_order,
        palette=ori_palette,
        patch_kws={'edgecolor': 'black',
                   'linewidth': linewidth})

    if gene is not None:
        filter = 'gene_name == {!r}'.format(gene)
    elif transcript_id is not None:
        filter = 'transcript_id == {!r}'.format(transcript_id)
    else:
        filter = None

    gene_track = BiomartTrack(
        dataset='mmusculus_gene_ensembl',
        host=bm_host,
        height=0.4,
        collapse='transcript',
        filter=filter,
        gene_id='gene_name',
        patch_kws={'linewidth': linewidth},
        line_kws={'linewidth': 1},
        label_kws={'fontstyle': 'italic'},
        bm_gene_name=bm_gene_name)

    fig = plot_tracks(
        [ins_track, gene_track],
        region=region,
        despine=True,
        **kwargs)

    fig.axes[0].set_title('Insertions')
    fig.axes[-1].set_xlabel('Chromosome {}'.format(region[0]))

    return fig
