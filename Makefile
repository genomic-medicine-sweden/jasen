# ==============================================================================
# Brief explanation of Makefile syntax
# ==============================================================================
# Since Makefiles are a bit special, and there are not a lot of easy to read
# tutorials out there, here follows a very brief introduction to the syntax in
# Makefiles.
#
# Basics
# ------------------------------------------------------------------------------
#
# Makefiles are in many ways similar to bash-scripts, but unlike bash-scripts
# they are not written in a linear sequential fashion, but rather divided into
# so called rules, that are typically tightly connected to specific output file
# paths or file path pattern.
#
# Firstly, Makefiles should be named `Makefile` and be put inside a folder
# where one wants to run the make command.
#
# The rule syntax
# ------------------------------------------------------------------------------
#
# The syntax of a rule is at its core very simple:
#
# output file(s)/target(s) : any dependent input files
#     commands to produce the output files from input files
#     possibly more commands
#
# So, the significant part of the rule syntax is the ':'-character, as well as
# the indentation of the rows following the head of the rule, to indicate the
# "recipe" or the commands to produce the outputs of the rule.
#
# A rule can also have just a name, and no concrete output file. That is, it
# would have the form:
#
# name_of_the_rule : any dependent input files
#     one commands
#     more commands
#
# Now, there is one big caveat here, related to our scripts: Make will rebuild
# a target as soon as any of its inputs are updated, or have a newer timestamp
# than the target. This is typically not desired for us, since we might have
# files unpacked with zip with all kinds of different dates, and for a one-run
# installation, we are mostly interested in wheter an output file already
# exists or not, not so much about timestamps.
#
# To change so that Make only cares about whether files exist, and not timestamps,
# one can add a | character before those input files, like so:
#
# name_of_the_rule : input files where timestamp matters | input files where only existence matters
#     one commands
#     more commands
#
# Of course, one can have everything on the right side of the | character, so:
#
# name_of_the_rule : | input files where only existence matters
#     one commands
#     more commands
#
# Running Makefiles
# ------------------------------------------------------------------------------
#
# To run a named rule, you would then do:
#
# $ make name_of_the_rule
#
# For normal rules, you would hit:
#
# $ make <outputfile>
#
# Tip: Type "make" and hit tab twice in the shell, to show available targets to
# run in the current makefile.
#
# Special variables and more
# ------------------------------------------------------------------------------
#
# Inside the commands of the rules, one can write pretty much bash, especially
# if setting `SHELL := /bin/bash` in the beginning of the Makefile which we have
# done below.
#
# There are a some differences though:
#
# 1. Variable syntax using only a single $-character refer to MAKE-variables.
#    Thus, to use variables set in bash, you have to use double-$.
#    So:
#    echo $(this-is-a-make-variable) and $$(this-is-a-bash-variable)
#
#    This also goes for command substitution, so rather than:
#
#    	echo "Lines in file:" $(wc -l somefile.txt)
#
#    ... you would write:
#
#    	echo "Lines in file:" $$(wc -l somefile.txt)
#
# 2. Makefiles also use some special variables with special meaning, to access
#    things like the output and input files:
#
#    $@   The output target (In case of multiple targets, the same command will
#         be run once for every output file, feeding this variable with only a
#         single one at a time)
#
#    $<   The first input file or dependency
#
#    $^   ALL input files / dependencies.
#
# 3. Another speciality is that adding a '@' character before any command, will
#    stop this command from being printed out when it is executed (only its output
#    will be printed).
#
# That's all for now. For the official docs, see:
# https://www.gnu.org/software/make/manual/make.html
#
# Some more caveats
# ------------------------------------------------------------------------------
#  Here are some things that might not be obvious at first:
#
#  To create a rule with a single output but many dependencies, you can add
#  these dependencies on multiple lines by using the continuation character \
#  like so:
#
#  name_of_rule: dependency1 \
#  		dependency2 \
#  		dependency3 \
#  		dependency4 \

# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash

# Define path variables
SCRIPT_DIR := $(shell pwd)
ASSETS_DIR := $(shell realpath $(SCRIPT_DIR)/assets/)
CONTAINER_DIR := $(realpath $(SCRIPT_DIR)/container/)
PRODIGAL_TRAINING_DIR := $(ASSETS_DIR)/prodigal_training_files
# The root folder where the pipeline is currently located. To be mounted into
# the Singularity containers below.
MNT_ROOT := /$(shell readlink -f . | cut -d"/" -f2)
INSTALL_LOG := "$(SCRIPT_DIR)/.install.log"

define log_message
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
	@echo "$$(date "+%Y-%m-%d %H:%M:%S"): $1" | tee -a $(INSTALL_LOG);
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
endef

# ==============================================================================
# Main rules
# ==============================================================================

print_paths:
	@echo "SCRIPT_DIR:" $(SCRIPT_DIR)
	@echo "ASSETS_DIR:" $(ASSETS_DIR)
	@echo "CONTAINER_DIR:" $(CONTAINER_DIR)
	@echo "MNT_ROOT:" $(MNT_ROOT)

install: download_or_build_containers \
	update_databases \
	update_organisms

update_databases: update_amrfinderplus \
	update_mlst_db \
	update_blast_db \
	update_finder_dbs \
	update_shigapass_db

update_organisms: staphylococcus_aureus_all \
	ecoli_all \
	kpneumoniae_all \
	mtuberculosis_all \
	spyogenes_all \
	streptococcus_all

check:	check_chewbbaca \
	check_bwa \
	check_blastdb

# ==============================================================================
# Check and update git submodules
# ==============================================================================
check-and-reinit-git-submodules:
	@if git submodule status | egrep -q '^[-]|^[+]' ; then \
		echo "Updating / re-initializing git submodules ..." \
		&& git submodule update --init; \
	fi

# ==============================================================================
# Build containers
# ==============================================================================

download_or_build_containers:
	$(call log_message,"Checking if any containers need to be built ...")
	@set -euo \
	&& cd $(CONTAINER_DIR) \
	&& make all; \
	cd -

# ==============================================================================
# Update databases
# ==============================================================================

# -----------------------------
# Update ShigaPass database
# -----------------------------
SHIGAPASS_DIR := $(ASSETS_DIR)/ShigaPass
update_shigapass_db: $(SHIGAPASS_DIR)/SCRIPT/ShigaPass_DataBases/IPAH/ipaH_150-mers.fasta.ndb

$(SHIGAPASS_DIR)/SCRIPT/ShigaPass_DataBases/IPAH/ipaH_150-mers.fasta.ndb:
	$(call log_message,"Starting update of ShigaPass database")
	cd $(SHIGAPASS_DIR) \
	&& chmod +x SCRIPT/ShigaPass.sh \
	&& gunzip $(SHIGAPASS_DIR)/Example/Input/*.fasta.gz \
	&& singularity exec \
		--bind $(MNT_ROOT) \
		$(CONTAINER_DIR)/shigapass.sif \
		bash ShigaPass.sh -u \
		-p $(SHIGAPASS_DIR)/SCRIPT/ShigaPass_DataBases \
		-l $(SHIGAPASS_DIR)/Example/Input/ShigaPass_test.txt \
		-o $(SHIGAPASS_DIR)/Example/ShigaPass_Results |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update AMRFinderPlus database
# -----------------------------
# For more info, see:
# https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/
AMRFINDERDB_DIR := $(ASSETS_DIR)/amrfinder_db

update_amrfinderplus: $(AMRFINDERDB_DIR)/latest

$(AMRFINDERDB_DIR)/latest:
	$(call log_message,"Starting update of AMRFinderPlus database ...")
	singularity exec \
		--bind $(MNT_ROOT) \
		$(CONTAINER_DIR)/ncbi-amrfinderplus.sif \
		amrfinder_update \
		--database $(AMRFINDERDB_DIR) |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update MLST database
# -----------------------------
MLSTDB_DIR := $(ASSETS_DIR)/mlst_db

update_mlst_db: $(MLSTDB_DIR)/pubmlst/dbases.xml

$(ASSETS_DIR)/mlst_db/pubmlst/dbases.xml:
	$(call log_message,"Starting update of MLST database ...")
	cd $(MLSTDB_DIR) \
	&& bash mlst-download_pub_mlst.sh |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update Blast database
# -----------------------------
update_blast_db: $(MLSTDB_DIR)/blast

$(MLSTDB_DIR)/blast:
	$(call log_message,"Starting update of Blast database")
	cd $(MLSTDB_DIR) \
	&& singularity exec \
		--bind $(MNT_ROOT) \
		$(CONTAINER_DIR)/blast.sif \
		bash $(MLSTDB_DIR)/mlst-make_blast_db.sh |& tee -a $(INSTALL_LOG)

# -----------------------------
# Update VirulenceFinder database
# -----------------------------
VIRULENCEFINDERDB_DIR := $(ASSETS_DIR)/virulencefinder_db
KMA_DIR := $(ASSETS_DIR)/kma

update_finder_dbs: $(VIRULENCEFINDERDB_DIR)/stx.name

$(ASSETS_DIR)/virulencefinder_db/stx.name: | check-and-reinit-git-submodules
	$(call log_message,"Starting update of VirulenceFinder database")
	cd $(ASSETS_DIR)/kma \
	&& make \
	&& cd $(VIRULENCEFINDERDB_DIR) \
	&& export PATH=$(ASSETS_DIR)/kma:$$PATH \
	&& singularity exec	--bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG) \
	&& cd $(ASSETS_DIR)/resfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG) \
	&& cd $(ASSETS_DIR)/pointfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG) \
	&& cd $(ASSETS_DIR)/serotypefinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index |& tee -a $(INSTALL_LOG)

# ==============================================================================
# Download, index and prep reference genomes for organisms
# ==============================================================================

# -----------------------------
# S. aureus
# -----------------------------
staphylococcus_aureus_all: staphylococcus_aureus_download_reference \
	staphylococcus_aureus_faidx_reference \
	staphylococcus_aureus_bwaidx_reference \
	staphylococcus_aureus_download_prodigal_training_file \
	staphylococcus_aureus_download_cgmlst_schema \
	staphylococcus_aureus_unpack_cgmlst_schema \
	staphylococcus_aureus_prep_cgmlst_schema

SAUR_GENOMES_DIR := $(ASSETS_DIR)/genomes/staphylococcus_aureus
SAUR_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/staphylococcus_aureus
SAUR_REFSEQ_ACC := GCF_000012045.1


staphylococcus_aureus_download_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading S. aureus reference genome ...")
	mkdir -p $(SAUR_GENOMES_DIR) \
	&& cd $(SCRIPT_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 bin/download_ncbi.py \
		-i $(SAUR_REFSEQ_ACC) \
		-o $(SAUR_GENOMES_DIR) |& tee -a $(INSTALL_LOG) \


staphylococcus_aureus_faidx_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.fai

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.fai: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. aureus reference genome using samtools...")
	cd $(SAUR_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/samtools.sif \
		samtools faidx $$(basename $<) |& tee -a $(INSTALL_LOG)


staphylococcus_aureus_bwaidx_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. aureus reference genome using bwa...")
	cd $(SAUR_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


staphylococcus_aureus_download_prodigal_training_file: $(PRODIGAL_TRAINING_DIR)/Staphylococcus_aureus.trn

$(PRODIGAL_TRAINING_DIR)/Staphylococcus_aureus.trn:
	$(call log_message,"Downloading S. aureus prodigal training file ...")
	mkdir -p $(PRODIGAL_TRAINING_DIR) \
	&& cd $(PRODIGAL_TRAINING_DIR) \
	&& wget https://raw.githubusercontent.com/B-UMMI/chewBBACA/master/CHEWBBACA/prodigal_training_files/Staphylococcus_aureus.trn \
		-O $@ \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


staphylococcus_aureus_download_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip

$(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip:
	$(call log_message,"Downloading S. aureus cgMLST schema ...")
	mkdir -p $(SAUR_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(SAUR_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/141106/alleles/ \
		-O $@ \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


staphylococcus_aureus_unpack_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/unpacking.done

$(SAUR_CGMLST_DIR)/alleles/unpacking.done: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip
	$(call log_message,"Unpacking S. aureus cgMLST schema ...")
	cd $$(dirname $<) \
		&& unzip -DDq $$(basename $<) |& tee -a $(INSTALL_LOG) \
		&& echo $$(date "+%Y%m%d %H:%M:%S")": Done unpacking zip file: " $< > $@

staphylococcus_aureus_prep_cgmlst_schema: | $(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn

$(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn: $(SAUR_CGMLST_DIR)/alleles_rereffed

$(SAUR_CGMLST_DIR)/alleles_rereffed: | $(SAUR_CGMLST_DIR)/alleles/unpacking.done check-and-reinit-git-submodules
	$(call log_message,"Prepping S. aureus cgMLST schema ...")
	cd $(SAUR_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-g $(SAUR_CGMLST_DIR)/alleles \
		-o $(SAUR_CGMLST_DIR)/alleles_rereffed \
		--cpu 2 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Staphylococcus_aureus.trn |& tee -a $(INSTALL_LOG)

# -----------------------------
# E. coli
# -----------------------------

ecoli_all: ecoli_download_reference \
	ecoli_faidx_reference \
	ecoli_bwaidx_reference \
	ecoli_download_prodigal_training_file \
	ecoli_download_wgmlst_schema \
	ecoli_prep_ecoli_cgmlst_schema

ECOLI_GENOMES_DIR := $(ASSETS_DIR)/genomes/escherichia_coli
ECOLI_WGMLST_DIR := $(ASSETS_DIR)/wgmlst/escherichia_coli
ECOLI_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/escherichia_coli
ECOLI_REFSEQ_ACC := GCF_000005845.2


ecoli_download_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading E. coli genome ...")
	cd $(SCRIPT_DIR) \
	&& mkdir -p $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 bin/download_ncbi.py \
		-i $(ECOLI_REFSEQ_ACC) \
		-o $(ECOLI_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


ecoli_faidx_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.fai

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.fai: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta
	$(call log_message,"Indexing E. coli genome using samtools...")
	cd $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/samtools.sif \
		samtools faidx $$(basename $<) |& tee -a $(INSTALL_LOG)


ecoli_bwaidx_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta
	$(call log_message,"Indexing E. coli genome using bwa...")
	cd $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


ecoli_download_prodigal_training_file: $(PRODIGAL_TRAINING_DIR)/Escherichia_coli.trn

$(PRODIGAL_TRAINING_DIR)/Escherichia_coli.trn:
	$(call log_message,"Downloading E. coli prodigal training file ...")
	mkdir -p $(PRODIGAL_TRAINING_DIR) \
	&& cd $(PRODIGAL_TRAINING_DIR) \
	&& wget https://raw.githubusercontent.com/B-UMMI/chewBBACA/master/CHEWBBACA/prodigal_training_files/Escherichia_coli.trn \
		-O $@ \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


# Download Ecoli wgmlst INNUENDO schema
ecoli_download_wgmlst_schema: $(ECOLI_WGMLST_DIR)/alleles/ecoli_INNUENDO_wgMLST/Escherichia_coli.trn

$(ECOLI_WGMLST_DIR)/alleles/ecoli_INNUENDO_wgMLST/Escherichia_coli.trn:
	$(call log_message,"Downloading E. coli wgMLST schema ...")
	mkdir -p $(ECOLI_WGMLST_DIR)/alleles &> /dev/null \
	&& cd $(ECOLI_WGMLST_DIR)/alleles \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie DownloadSchema \
		-sp 5 \
		-sc 1 \
		-o $(ECOLI_WGMLST_DIR)/alleles \
		--latest |& tee -a $(INSTALL_LOG)


# Download Ecoli cgmlst cgmlst.org schema
ecoli_download_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip

$(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip:
	$(call log_message,"Downloading E. coli cgMLST schema ...")
	mkdir -p $(ECOLI_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(ECOLI_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/5064703/alleles/ \
		-O $$(basename $@) \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


# Unpack Ecoli cgmlst schema
ecoli_unpack_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/unpacking.done

$(ECOLI_CGMLST_DIR)/alleles/unpacking.done: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip
	$(call log_message,"Unpacking E. coli cgMLST schema ...")
	cd $(ECOLI_CGMLST_DIR)/alleles \
	&& unzip -DDq $$(basename $<) |& tee -a $(INSTALL_LOG) \
	&& echo $$(date +"%Y%m%d %H:%M:%S")": Done unpacking zip file: " $< > $@


# Prepping Ecoli cgmlst cgmlst.org schema
ecoli_prep_ecoli_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn

$(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn: $(ECOLI_CGMLST_DIR)/alleles_rereffed

$(ECOLI_CGMLST_DIR)/alleles_rereffed: | $(ECOLI_CGMLST_DIR)/alleles/unpacking.done check-and-reinit-git-submodules
	$(call log_message,"Prepping E. coli cgMLST schema ... WARNING: This takes a looong time. Put on some coffee")
	cd $(ECOLI_CGMLST_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-g $(ECOLI_CGMLST_DIR)/alleles \
		-o $(ECOLI_CGMLST_DIR)/alleles_rereffed \
		--cpu 2 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Escherichia_coli.trn |& tee -a $(INSTALL_LOG)


# -----------------------------
# K. pneumoniae
# -----------------------------

kpneumoniae_all: kpneumoniae_download_reference \
	kpneumoniae_faidx_reference \
	kpneumoniae_bwaidx_reference \
	kpneumoniae_download_prodigal_training_file \
	kpneumoniae_download_cgmlst_schema \
	kpneumoniae_unpack_cgmlst_schema \
	kpneumoniae_prep_cgmlst_schema


KPNEU_GENOMES_DIR := $(ASSETS_DIR)/genomes/klebsiella_pneumoniae
KPNEU_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/klebsiella_pneumoniae
KPNEU_REFSEQ_ACC := GCF_000240185.1


kpneumoniae_download_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading K. pneumoniae genome ...")
	cd $(SCRIPT_DIR) \
	&& mkdir -p $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 bin/download_ncbi.py \
		-i $(KPNEU_REFSEQ_ACC) \
		-o $(KPNEU_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


kpneumoniae_faidx_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.fai

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.fai: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta
	$(call log_message,"Indexing K. pneumoniae genome using samtools...")
	cd $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/samtools.sif \
		samtools faidx $$(basename $<) |& tee -a $(INSTALL_LOG)


kpneumoniae_bwaidx_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta
	$(call log_message,"Indexing K. pneumoniae genome using bwa...")
	cd $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


kpneumoniae_download_prodigal_training_file: $(PRODIGAL_TRAINING_DIR)/Klebsiella_pneumoniae.trn

$(PRODIGAL_TRAINING_DIR)/Klebsiella_pneumoniae.trn:
	$(call log_message,"Downloading K. pneumonia prodigal training file ...")
	mkdir -p $(PRODIGAL_TRAINING_DIR) \
	&& cd $(PRODIGAL_TRAINING_DIR) \
	&& wget https://raw.githubusercontent.com/B-UMMI/chewBBACA/master/CHEWBBACA/prodigal_training_files/Klebsiella_pneumoniae.trn \
		-O $@ \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


# Download Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_download_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip

$(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip:
	$(call log_message,"Downloading K. pneumoniae cgMLST schema ...")
	mkdir -p $(KPNEU_CGMLST_DIR)/alleles \
	&& cd $(KPNEU_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/2187931/alleles/ \
		-O $$(basename $@) \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


kpneumoniae_unpack_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/unpacking.done

$(KPNEU_CGMLST_DIR)/alleles/unpacking.done: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip
	$(call log_message,"Unpacking K. pneumoniae cgMLST schema ...")
	cd $(KPNEU_CGMLST_DIR)/alleles \
	&& unzip -DDq $$(basename $<) |& tee -a $(INSTALL_LOG) \
	&& echo $$(date "+%Y%m%d %H:%M:%S")": Done unpacking zip file: " $< > $@


# Prep Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_prep_cgmlst_schema: | $(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn

$(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn: $(KPNEU_CGMLST_DIR)/alleles_rereffed

$(KPNEU_CGMLST_DIR)/alleles_rereffed: | $(KPNEU_CGMLST_DIR)/alleles/unpacking.done check-and-reinit-git-submodules
	$(call log_message,"Prepping K. pneumoniae cgMLST schema ... Warning: This takes a looong time. Put on some coffee!")
	cd $(KPNEU_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-g $(KPNEU_CGMLST_DIR)/alleles \
		-o $(KPNEU_CGMLST_DIR)/alleles_rereffed \
		--cpu 2 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Klebsiella_pneumoniae.trn |& tee -a $(INSTALL_LOG)

# -----------------------------
# Streptococcus pyogenes
# -----------------------------

spyogenes_all: spyogenes_download_reference \
	spyogenes_faidx_reference \
	spyogenes_bwaidx_reference \
	spyogenes_download_prodigal_training_file \
	spyogenes_download_cgmlst_schema \
	spyogenes_unpack_cgmlst_schema \
	spyogenes_prep_cgmlst_schema


SPYO_GENOMES_DIR := $(ASSETS_DIR)/genomes/streptococcus_pyogenes
SPYO_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/streptococcus_pyogenes
SPYO_REFSEQ_ACC := GCF_000006785.2


spyogenes_download_reference: $(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta

$(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading S. pyogenes genome ...")
	cd $(SCRIPT_DIR) \
	&& mkdir -p $(SPYO_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 bin/download_ncbi.py \
		-i $(SPYO_REFSEQ_ACC) \
		-o $(SPYO_GENOMES_DIR) |& tee -a $(INSTALL_LOG)


spyogenes_faidx_reference: $(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta.fai

$(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta.fai: $(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. pyogenes genome using samtools...")
	cd $(SPYO_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/samtools.sif \
		samtools faidx $$(basename $<) |& tee -a $(INSTALL_LOG)


spyogenes_bwaidx_reference: $(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta.bwt

$(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta.bwt: $(SPYO_GENOMES_DIR)/$(SPYO_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. pyogenes genome using bwa...")
	cd $(SPYO_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)


spyogenes_download_prodigal_training_file: $(PRODIGAL_TRAINING_DIR)/Streptococcus_pyogenes.trn

$(PRODIGAL_TRAINING_DIR)/Streptococcus_pyogenes.trn:
	$(call log_message,"Downloading K. pneumonia prodigal training file ...")
	mkdir -p $(PRODIGAL_TRAINING_DIR) \
	&& cd $(PRODIGAL_TRAINING_DIR) \
	&& wget https://raw.githubusercontent.com/B-UMMI/chewBBACA/master/CHEWBBACA/prodigal_training_files/Streptococcus_pyogenes.trn \
		-O $@ \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


# Download Streptococcus pyogenes cgmlst cgmlst.org schema
spyogenes_download_cgmlst_schema: $(SPYO_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip

$(SPYO_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip:
	$(call log_message,"Downloading S. pyogenes cgMLST schema ...")
	mkdir -p $(SPYO_CGMLST_DIR)/alleles \
	&& cd $(SPYO_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/2187931/alleles/ \
		-O $$(basename $@) \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)


spyogenes_unpack_cgmlst_schema: $(SPYO_CGMLST_DIR)/alleles/unpacking.done

$(SPYO_CGMLST_DIR)/alleles/unpacking.done: $(SPYO_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip
	$(call log_message,"Unpacking S. pyogenes cgMLST schema ...")
	cd $(SPYO_CGMLST_DIR)/alleles \
	&& unzip -DDq $$(basename $<) |& tee -a $(INSTALL_LOG) \
	&& echo $$(date "+%Y%m%d %H:%M:%S")": Done unpacking zip file: " $< > $@


# Prep Streptococcus cgmlst cgmlst.org schema
spyogenes_prep_cgmlst_schema: | $(SPYO_CGMLST_DIR)/alleles_rereffed/Streptococcus_pyogenes.trn

$(SPYO_CGMLST_DIR)/alleles_rereffed/Streptococcus_pyogenes.trn: $(SPYO_CGMLST_DIR)/alleles_rereffed

$(SPYO_CGMLST_DIR)/alleles_rereffed: | $(SPYO_CGMLST_DIR)/alleles/unpacking.done check-and-reinit-git-submodules
	$(call log_message,"Prepping S. pyogenes cgMLST schema ... Warning: This takes a looong time. Put on some coffee!")
	cd $(SPYO_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-g $(SPYO_CGMLST_DIR)/alleles \
		-o $(SPYO_CGMLST_DIR)/alleles_rereffed \
		--cpu 2 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Streptococcus_pyogenes.trn |& tee -a $(INSTALL_LOG)

# -----------------------------
# Streptococcus
# -----------------------------

streptococcus_all: streptococcus_download_cgmlst_schema \
	streptococcus_unpack_cgmlst_schema \
	streptococcus_prep_cgmlst_schema


STREP_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/streptococcus

# Download Streptococcus cgmlst cgmlst.org schema
streptococcus_download_cgmlst_schema: $(STREP_CGMLST_DIR)/alleles/index.html

$(STREP_CGMLST_DIR)/alleles/index.html:
	$(call log_message,"Downloading Streptococcus cgMLST schema ...")
	mkdir -p $(STREP_CGMLST_DIR)/alleles \
	&& cd $(STREP_CGMLST_DIR)/alleles \
	&& wget -r https://enterobase.warwick.ac.uk/schemes/Streptococcus.cgMLSTv1/ \
		-np \
		-nH \
		--cut-dirs 2 \
		--no-verbose \
		--no-check-certificate |& tee -a $(INSTALL_LOG)

streptococcus_unpack_cgmlst_schema: $(STREP_CGMLST_DIR)/alleles/unpacking.done

$(STREP_CGMLST_DIR)/alleles/unpacking.done: $(STREP_CGMLST_DIR)/alleles/profiles.list
	$(call log_message,"Unpacking S. pyogenes cgMLST schema ...")
	cd $(STREP_CGMLST_DIR)/alleles \
	&& gunzip *.gz |& tee -a $(INSTALL_LOG) \
	&& echo $$(date "+%Y%m%d %H:%M:%S")": Done unpacking gz files: " $< > $@

# Prep Streptococcus cgmlst cgmlst.org schema
streptococcus_prep_cgmlst_schema: $(STREP_CGMLST_DIR)/alleles_rereffed_summary_stats.tsv

$(STREP_CGMLST_DIR)/alleles_rereffed_summary_stats.tsv: $(STREP_CGMLST_DIR)/alleles_rereffed

$(STREP_CGMLST_DIR)/alleles_rereffed: | $(STREP_CGMLST_DIR)/alleles/index.html
	$(call log_message,"Prepping Streptococcus cgMLST schema ... Warning: This takes a looong time. Put on some coffee!")
	cd $(STREP_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-g $(STREP_CGMLST_DIR)/alleles \
		-o $(STREP_CGMLST_DIR)/alleles_rereffed \
		--cpu 2 |& tee -a $(INSTALL_LOG)

# -----------------------------
# M. tuberculosis
# -----------------------------

MTUBE_GENOMES_DIR := $(ASSETS_DIR)/genomes/mycobacterium_tuberculosis
MTUBE_TBDB_DIR := $(ASSETS_DIR)/tbdb
MTUBE_TBPROFILER_DBS_DIR := $(ASSETS_DIR)/tbprofiler_dbs
MTUBE_REFSEQ_ACC := GCF_000195955.2

mtuberculosis_all: mtuberculosis_download_reference mtuberculosis_faidx_reference mtuberculosis_bwaidx_reference mtuberculosis_converged_who_fohm_tbdb mtuberculosis_bgzip_bed mtuberculosis_index_bed

mtuberculosis_download_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading M. tuberculosis genome ...")
	mkdir -p $(MTUBE_GENOMES_DIR) \
	&& cd $(SCRIPT_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bonsai-prp.sif \
		python3 bin/download_ncbi.py \
		-i $(MTUBE_REFSEQ_ACC) \
		-o $(MTUBE_GENOMES_DIR) |& tee -a $(INSTALL_LOG)

mtuberculosis_faidx_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.fai

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.fai: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta
	$(call log_message,"Indexing M. tuberculosis genome using samtools...")
	cd $(MTUBE_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/samtools.sif \
		samtools faidx $$(basename $<) |& tee -a $(INSTALL_LOG)

mtuberculosis_bwaidx_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta
	$(call log_message,"Indexing M. tuberculosis genome using bwa...")
	cd $(MTUBE_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<) |& tee -a $(INSTALL_LOG)

mtuberculosis_converged_who_fohm_tbdb: $(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.variables.json

$(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.variables.json: $(MTUBE_TBPROFILER_DBS_DIR)/converged_who_fohm_tbdb.csv
	$(call log_message,"Creating WHO FoHM TBDB ...")
	cd $(MTUBE_TBDB_DIR) \
	&& cp $(MTUBE_TBPROFILER_DBS_DIR)/converged_who_fohm_tbdb.csv $(MTUBE_TBDB_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/tb-profiler.sif \
		tb-profiler create_db --prefix converged_who_fohm_tbdb --dir $(MTUBE_TBDB_DIR) \
		--match_ref $(MTUBE_GENOMES_DIR)/GCF_000195955.2.fasta --csv converged_who_fohm_tbdb.csv \
	&& tb-profiler load_library converged_who_fohm_tbdb --dir $(MTUBE_TBDB_DIR) |& tee -a $(INSTALL_LOG)

mtuberculosis_bgzip_bed: $(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed.gz

$(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed.gz: $(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed
	$(call log_message,"Bgzipping converged WHO + FoHM + TBDB bed file ...")
	cd $(MTUBE_TBDB_DIR) \
	&& bgzip $$(basename $<) |& tee -a $(INSTALL_LOG)

mtuberculosis_index_bed: $(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed.gz.tbi

$(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed.gz.tbi: $(MTUBE_TBDB_DIR)/converged_who_fohm_tbdb.bed.gz
	$(call log_message,"Indexing converged WHO + FoHM + TBDB bgzipped bed file ...")
	cd $(MTUBE_TBDB_DIR) \
	&& tabix -p bed $$(basename $<) |& tee -a $(INSTALL_LOG)

# ==============================================================================
# Perform checks
# ==============================================================================
# TODO: Integrate with main installation rules?

# -----------------------------
# Check chewBBACA
# -----------------------------
check_chewbbaca:
	@cd $(SCRIPT_DIR) \
	&& saureus=$(SAUR_CGMLST_DIR)/alleles_rereffed \
	&& ecoli=$(ECOLI_CGMLST_DIR)/alleles_rereffed \
	&& kpneumoniae=$(KPNEU_CGMLST_DIR)/alleles_rereffed \
	&& spyogenes=$(SPYO_CGMLST_DIR)/alleles_rereffed \
	&& streptococcus=$(STREP_CGMLST_DIR)/alleles_rereffed \
	&& if [[ -d "$$saureus" && -d "$$ecoli" && -d "$$kpneumoniae" && -d "$$spyogenes" && -d "$$streptococcus" ]]; then \
		echo "[✓] PASSED check for chewBBACA: Directories exist:" |& tee -a $(INSTALL_LOG) \
		&& echo "- $$saureus" |& tee -a $(INSTALL_LOG) \
		&& echo "- $$ecoli" |& tee -a $(INSTALL_LOG) \
		&& echo "- $$kpneumoniae" |& tee -a $(INSTALL_LOG) \
		&& echo "- $$spyogenes" |& tee -a $(INSTALL_LOG) \
		&& echo "- $$streptococcus" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for chewBBACA: Some directories do not exist:"; \
		if [[ ! -d $$saureus ]]; then \
			echo "    Missing directory: $$saureus" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$ecoli ]]; then \
			echo "    Missing directory: $$ecoli" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$kpneumoniae ]]; then \
			echo "    Missing directory: $$kpneumoniae" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$spyogenes ]]; then \
			echo "    Missing directory: $$spyogenes" |& tee -a $(INSTALL_LOG);  \
		fi; \
		if [[ ! -d $$streptococcus ]]; then \
			echo "    Missing directory: $$streptococcus" |& tee -a $(INSTALL_LOG);  \
		fi; \
		echo "    Please report this in an issue on the JASEN repo: https://github.com/genomic-medicine-sweden/JASEN/issues" |& tee -a $(INSTALL_LOG); \
	fi

# -----------------------------
# Check BWA
# -----------------------------
check_bwa:
	@cd $(SCRIPT_DIR) \
	&& ref=$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta; \
		refamb=$${ref}.amb; \
		refann=$${ref}.ann; \
		refbwt=$${ref}.bwt; \
		refpac=$${ref}.pac; \
		refsa=$${ref}.sa \
	&& if [[ -f $$ref \
		&& -f $${refamb} \
		&& -f $${refann} \
		&& -f $${refbwt} \
		&& -f $${refpac} \
		&& -f $${refsa} \
	]]; then \
		echo "[✓] PASSED check for bwa: Indexes exist in $(SAUR_GENOMES_DIR)" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for bwa: Indexes do not exist in $(SAUR_GENOMES_DIR)" |& tee -a $(INSTALL_LOG); \
	fi

# -----------------------------
# Check BlastDB
# -----------------------------
MLST_BLAST_DIR := $(ASSETS_DIR)/mlst_db/blast
check_blastdb:
	@cd $(SCRIPT_DIR) \
	&& mlst=$(MLST_BLAST_DIR)/mlst.fa; \
	 mlstndb=$${mlst}.ndb; \
	 mlstnhd=$${mlst}.nhd; \
	 mlstnhi=$${mlst}.nhi; \
	 mlstnhr=$${mlst}.nhr; \
	 mlstnin=$${mlst}.nin; \
	 mlstnog=$${mlst}.nog; \
	 mlstnos=$${mlst}.nos; \
	 mlstnot=$${mlst}.not; \
	 mlstnsq=$${mlst}.nsq; \
	 mlstntf=$${mlst}.ntf; \
	 mlstnto=$${mlst}.nto \
	&& if [[ -f $${mlst} \
		&& -f $${mlstndb} \
		&& -f $${mlstnhd} \
		&& -f $${mlstnhi} \
		&& -f $${mlstnhr} \
		&& -f $${mlstnin} \
		&& -f $${mlstnog} \
		&& -f $${mlstnos} \
		&& -f $${mlstnot} \
		&& -f $${mlstnsq} \
		&& -f $${mlstntf} \
		&& -f $${mlstnto} \
	]]; then \
		echo "[✓] PASSED check for blast: Indexes exist in $(MLST_BLAST_DIR)" |& tee -a $(INSTALL_LOG); \
	else \
		echo "[!] FAILED check for blast: Indexes do not exist in $(MLST_BLAST_DIR)" |& tee -a $(INSTALL_LOG); \
	fi
