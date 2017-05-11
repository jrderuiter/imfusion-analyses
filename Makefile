.PHONY: clean clean_data env env_tophat lint requirements

################################################################################
# GLOBALS                                                                      #
################################################################################

SNAKEMAKE = snakemake
SNAKEMAKE_ARGS =

ENSEMBL_FASTA = data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa
ENSEMBL_GTF =  data/external/ensembl/Mus_musculus.GRCm38.76.gtf
DEXSEQ_GTF = data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf

TRANSPOSON_FASTA = src/im-fusion/data/t2onc/t2onc2.sequence.fa
TRANSPOSON_FEATURES = src/im-fusion/data/t2onc/t2onc2.features.txt

FUSION_FINDER_INDEX = data/interim/references/fusion-finder/Mus_musculus.GRCm38.76.t2onc.tophat2

################################################################################
# COMMANDS                                                                     #
################################################################################

env:
	conda env create -f envs/imfusion.yml

env_tophat:
	conda env create -f envs/imfusion-tophat.yml

clean:
	find . -name '*.pyc' -exec rm {} \;

clean_data:
	rm -rf data/external/ensembl
	rm -rf data/interim/*

lint:
	flake8 --exclude=lib/,bin/,docs/conf.py .

doc:
	rm -rf ./docs/_build
	sphinx-autobuild ./docs ./docs/_build

gh-pages:
	git checkout gh-pages
	find ./* -not -path '*/\.*' -prune -exec rm -r "{}" \;
	git checkout develop docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	git reset HEAD
	(cd docs && make html)
	mv -fv docs/_build/html/* ./
	rm -rf docs Makefile src AUTHORS.rst CONTRIBUTING.rst HISTORY.rst README.rst
	touch .nojekyll
	git add -A
	git commit -m "Generated gh-pages for `git log develop -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout develop


################################################################################
# PROJECT RULES                                                                #
################################################################################

build_reference_star: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	imfusion-build star \
		--reference_seq $(ENSEMBL_FASTA) \
		--reference_gtf $(ENSEMBL_GTF) \
		--transposon_seq $(TRANSPOSON_FASTA) \
		--transposon_features $(TRANSPOSON_FEATURES) \
		--output_dir data/interim/references/imfusion/GRCm38.76.t2onc.star \
		--blacklist_genes ENSMUSG00000039095 ENSMUSG00000038402 \
		--overhang 100

build_reference_tophat: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	imfusion-build tophat \
		--reference_seq $(ENSEMBL_FASTA) \
		--reference_gtf $(ENSEMBL_GTF) \
		--transposon_seq $(TRANSPOSON_FASTA) \
		--transposon_features $(TRANSPOSON_FEATURES) \
		--output_dir data/interim/references/imfusion/GRCm38.76.t2onc.tophat \
		--blacklist_genes ENSMUSG00000039095 ENSMUSG00000038402

build_reference_fusion_finder: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	mkdir -p $(dir $(FUSION_FINDER_INDEX))
	bedtools maskfasta \
		-fi $(ENSEMBL_FASTA) \
		-bed data/raw/fusion-finder/mask_regions.bed \
		-fo $(FUSION_FINDER_INDEX).fa
	cat $(TRANSPOSON_FASTA) >> $(FUSION_FINDER_INDEX).fa
	bowtie2-build $(FUSION_FINDER_INDEX).fa $(FUSION_FINDER_INDEX)
	tophat2 -G $(ENSEMBL_GTF) \
		--transcriptome-index=$(FUSION_FINDER_INDEX).transcriptome \
		$(FUSION_FINDER_INDEX)
	rm -rf ./tophat_out

shear_splink_download: data/external/shear_splink

sb_download:
	$(SNAKEMAKE) -s pipelines/prepare-sb.snake -p $(SNAKEMAKE_ARGS) \
	    --config samples=data/raw/sb/samples.txt \
				 fastq_dir=data/interim/sb/fastq \
				 sample_out_path=data/processed/sb/samples.txt

sb_imfusion: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	$(SNAKEMAKE) -s pipelines/imfusion.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sb-imfusion.yml \
		--config fastq_dir='data/interim/sb/fastq' \
				 interim_dir='data/interim/sb/star' \
			     processed_dir='data/processed/sb/star' \
				 log_dir='logs/sb/star' \
				 paired=False

sb_star_fusion:
	$(SNAKEMAKE) -s pipelines/star-fusion.snake $(SNAKEMAKE_ARGS) \
	    --configfile configs/star-fusion.yml \
		--config fastq_dir='data/interim/sb/fastq' \
			     interim_dir='data/interim/sb/star-fusion' \
			     processed_dir='data/processed/sb/star-fusion' \
				 log_dir='logs/sb/star-fusion' \
				 samples='data/processed/sb/samples.txt' \
				 paired=False

sb_gene_expression: # sb_imfusion
	featureCounts -a data/external/ensembl/Mus_musculus.GRCm38.76.gtf \
		-o data/processed/sb/star/expression.fc_gene.txt -T 6 \
		`find data/interim/sb/star -name 'alignment.bam'`

sanger_download:
	$(SNAKEMAKE) -s pipelines/prepare-sanger.snake -p $(SNAKEMAKE_ARGS) \
		--config srdf=data/raw/sanger/E-ERAD-264.sdrf.txt \
		 	     sample_path=data/processed/sanger/samples.txt \
				 download_dir=tmp/download/_sanger
	rm -rf tmp/download/_sanger


sanger_imfusion: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	$(SNAKEMAKE) -s pipelines/imfusion.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger-imfusion.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     interim_dir='data/interim/sanger/star' \
			     processed_dir='data/processed/sanger/star' \
				 log_dir='logs/sanger/star' \
				 paired=True \
				 aligner=star

sanger_imfusion_single: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	$(SNAKEMAKE) -s pipelines/imfusion.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger-imfusion.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     interim_dir='data/interim/sanger/star-single' \
			     processed_dir='data/processed/sanger/star-single' \
				 log_dir='logs/sanger/star-single' \
				 paired=False \
				 aligner=star

sanger_imfusion_subsample: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	snakemake -s pipelines/imfusion-subsample.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger-imfusion.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     interim_dir='data/interim/sanger/star-subsample' \
			     processed_dir='data/processed/sanger/star-subsample' \
				 log_dir='logs/sanger/star-subsample' \
				 paired=True \
				 aligner=star \
				 n_reads='15,30,50,70'

sanger_imfusion_tophat: $(ENSEMBL_FASTA) $(ENSEMBL_GTF)
	$(SNAKEMAKE) -s pipelines/imfusion.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger-imfusion-tophat.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     interim_dir='data/interim/sanger/tophat' \
			     processed_dir='data/processed/sanger/tophat' \
				 log_dir='logs/sanger/tophat' \
				 paired=True \
				 aligner=tophat

sanger_fusion_finder:
	$(SNAKEMAKE) -s pipelines/fusion-finder.snake $(SNAKEMAKE_ARGS) \
		--configfile configs/sanger-fusion-finder.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
				 interim_dir='data/interim/sanger/fusion-finder' \
				 processed_dir='data/processed/sanger/fusion-finder' \
				 log_dir='logs/sanger/fusion-finder'

sanger_star_fusion:
	$(SNAKEMAKE) -s pipelines/star-fusion.snake $(SNAKEMAKE_ARGS) \
	    --configfile configs/star-fusion.yml \
		--config fastq_dir='data/interim/sanger/fastq' \
			     interim_dir='data/interim/sanger/star-fusion' \
			     processed_dir='data/processed/sanger/star-fusion' \
				 log_dir='logs/sanger/star-fusion' \
				 samples='data/processed/sanger/samples.txt'

################################################################################
# External files                                                               #
################################################################################

data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa:
	wget --quiet -P data/external/ensembl ftp://ftp.ensembl.org/pub/release-76/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	gunzip data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
	samtools faidx data/external/ensembl/Mus_musculus.GRCm38.dna.primary_assembly.fa

data/external/ensembl/Mus_musculus.GRCm38.76.gtf:
	wget --quiet -P data/external/ensembl ftp://ftp.ensembl.org/pub/release-76/gtf/mus_musculus/Mus_musculus.GRCm38.76.gtf.gz
	gunzip data/external/ensembl/Mus_musculus.GRCm38.76.gtf.gz

data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf: $(ENSEMBL_GTF)
	python scripts/dexseq_prepare_annotation.py --aggregate=no $(ENSEMBL_GTF) data/external/ensembl/Mus_musculus.GRCm38.76.DEXSeq.gtf

data/external/shear_splink:
	wget -O tmp/processed_freeze.tar.gz https://ndownloader.figshare.com/files/8297795?private_link=dd2515b13a5d022eba4d
	mkdir -p tmp/sb_screen
	tar -xf tmp/processed_freeze.tar.gz --directory tmp/sb_screen
	mv tmp/sb_screen/sb/shear_splink/full/all data/external/shear_splink
	rm -rf tmp/processed_freeze.tar.gz tmp/sb_screen
