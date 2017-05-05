from io import StringIO

import numpy as np
import pandas as pd
import pybiomart
import pysam

from ngstk.tabix import GtfIterator
from ngstk.util.genomic import GenomicDataFrame


def read_bed(bed_path):
    """Reads bed file into a DataFrame."""

    bed_cols = [
        'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand',
        'thickStart', 'thickEnd', 'itemRgb', 'blockCount', 'blockSizes',
        'blockStarts'
    ]

    data = pd.read_csv(bed_path, sep='\t', header=None)
    data.columns = bed_cols[:data.shape[1]]

    return data


def read_gtf_genes(gtf_path):
    """Reads genes from a gtf file."""

    # Read gene records.
    gtf_iter = GtfIterator(gtf_path)

    records = gtf_iter.fetch(filters=[lambda rec: rec.feature == 'gene'])

    # Convert to frame.
    def _record_to_row(rec):
        return {
            'id': rec['gene_id'],
            'name': rec['gene_name'],
            'contig': rec.contig,
            'start': rec.start,
            'end': rec.end,
            'strand': rec.strand
        }

    columns = ['id', 'name', 'contig', 'start', 'end', 'strand']

    genes = pd.DataFrame.from_records(
        (_record_to_row(rec) for rec in records), columns=columns)

    return genes


def fetch_biomart_genes_mm9():
    """Fetches mm9 genes from Ensembl via biomart."""

    return _fetch_genes_biomart(
        host='http://may2012.archive.ensembl.org',
        gene_name_attr='external_gene_id')


def fetch_biomart_genes_mm10():
    """Fetches mm9 genes from Ensembl via biomart."""

    return _fetch_genes_biomart(
        host='http://www.ensembl.org', gene_name_attr='external_gene_name')


def _fetch_genes_biomart(host, gene_name_attr):
    # Fetch genes from biomart.
    dataset = pybiomart.Dataset(name='mmusculus_gene_ensembl', host=host)

    gene_data = dataset.query(
        attributes=[
            'ensembl_gene_id', 'chromosome_name', 'start_position',
            'end_position', 'strand', gene_name_attr, 'gene_biotype'
        ],
        use_attr_names=True)

    # Rename columns for convenience.
    gene_data = gene_data.rename(columns={
        'chromosome_name': 'chromosome',
        'start_position': 'start',
        'end_position': 'end',
        'ensembl_gene_id': 'gene_id',
        gene_name_attr: 'gene_name'
    })

    # Assemble GenomicDataFrame (which has fast lookups).
    genes = GenomicDataFrame(gene_data)

    return genes


def read_star_junctions(file_path):
    """Reads junctions from STAR tab file."""

    data = pd.read_csv(
        str(file_path),
        sep='\t',
        names=[
            'chromosome', 'start', 'end', 'strand', 'intron_motif',
            'annotated', 'unique_reads', 'multi_reads', 'overhang'
        ],
        dtype={'chromosome': str})

    data['strand'] = data['strand'].map({0: np.nan, 1: 1, 2: -1})

    return data


def tidy_shared_axes(fig, axes, ylabel_x=0):
    """Removes duplicate labels from shared axes."""

    if len(axes.shape) > 1:
        for row in axes[:-1]:
            for ax in row:
                ax.set_xticks([])

        ylabel = axes[0, 0].get_ylabel()

        for ax in axes.flatten():
            ax.set_ylabel('')

        fig.text(x=ylabel_x, y=0.5, s=ylabel, va='center', rotation=90)
    else:
        for ax in axes[1:]:
            ax.set_ylabel('')


def idxstats(file_path):
    """Returns idxstats output as a DataFrame."""

    buf = StringIO()
    buf.write(pysam.idxstats(file_path))
    buf.seek(0)

    names = ['seq_name', 'seq_length', 'n_mapped', 'n_unmapped']
    return pd.read_csv(buf, sep='\t', names=names)


def flagstat(file_path):
    output = pysam.flagstat(str(file_path))

    lines = output.split('\n')
    values = [int(line.split()[0]) for line in lines if line]

    index = [
        'total', 'secondary', 'supplementary', 'duplicates', 'mapped',
        'paired', 'read1', 'read2', 'properly_paired', 'both_mapped',
        'singletons', 'mate_diff_chr', 'mate_diff_chr_mapq5'
    ]

    return pd.DataFrame({'count': values}, index=index)
