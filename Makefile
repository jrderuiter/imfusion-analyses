.PHONY: clean env env_tophat lint requirements

################################################################################
# GLOBALS                                                                      #
################################################################################

SNAKEMAKE_ARGS =

ENSEMBL_FASTA = data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
ENSEMBL_GTF =  data/external/ensembl/Mus_musculus.GRCm38.76.gtf
DEXSEQ_GTF = data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf


################################################################################
# COMMANDS                                                                     #
################################################################################

env:
	conda env create -f envs/imfusion_star.yml

env_tophat:
	conda env create -f envs/imfusion_tophat.yml

env_htseq:
	conda env create -f envs/htseq.yml

clean:
	find . -name '*.pyc' -exec rm {} \;

clean_data:
	rm -rf data/external/ensembl
	rm -rf data/interim/*

lint:
	flake8 --exclude=lib/,bin/,docs/conf.py .


################################################################################
# PROJECT RULES                                                                #
################################################################################

download_sanger:
	snakemake -s scripts/download_sanger.snake -p $(SNAKEMAKE_ARGS) \
		--config srdf=data/external/E-ERAD-264.sdrf.txt \

download_sb:
	snakemake -s scripts/download_sb.snake -p $(SNAKEMAKE_ARGS) \
	    --config samples=data/raw/samples.sb.txt

build_reference_star: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	im-fusion build star --reference_seq $(ENSEMBL_FASTA) \
		--reference_gtf $(ENSEMBL_GTF) --transposon_seq $(TRANSPOSON_FASTA) \
        --output_base data/interim/references/star/Mus_musculus.GRCm38.t2onc

build_reference_tophat2: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	im-fusion build tophat2 --reference_seq $(ENSEMBL_FASTA) \
		--reference_gtf $(ENSEMBL_GTF) --transposon_seq $(TRANSPOSON_FASTA) \
        --output_base data/interim/references/tophat2/Mus_musculus.GRCm38.t2onc

imfusion_sb: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	snakemake -s scripts/imfusion.snake -p $(SNAKEMAKE_ARGS) \
		--configfile configs/sb.yml
		--config fastq_dir='data/interim/sb/fastq' \
			     output_dir='data/processed/sb' \
				 paired=False

imfusion_sanger: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	snakemake -s scripts/imfusion.snake -p $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     output_dir='data/processed/sanger'

imfusion_sanger_tophat2: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	snakemake -s scripts/imfusion.snake -p $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     output_dir='data/processed/sanger-tophat2'

imfusion_sanger_subsample: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	snakemake -s scripts/imfusion_subsample.snake -p $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
				 output_dir='data/processed/sanger-subsample' \
				 n_reads='15,25,50,70'

fusion_finder_sanger:
	snakemake -s scripts/fusion_finder.snake -p $(SNAKEMAKE_ARGS) \
		--configfile configs/fusion_finder.yml


################################################################################
# External files                                                               #
################################################################################

data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa:
	wget --quiet -P data/external/ensembl ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	gunzip data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	# samtools faidx data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa

data/external/ensembl/Mus_musculus.GRCm38.76.gtf:
	wget --quiet -P data/external/ensembl ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/Mus_musculus.GRCm38.76.gtf.gz
	gunzip data/external/ensembl/Mus_musculus.GRCm38.76.gtf.gz

data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf: $(ENSEMBL_GTF)
	python scripts/dexseq_prepare_annotation.py --aggregate=no $(ENSEMBL_GTF) data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf
