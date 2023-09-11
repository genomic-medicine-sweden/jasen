# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash

# Define path variables
SCRIPT_DIR := $(shell pwd)
ASSETS_DIR := $(shell realpath $(SCRIPT_DIR)/../assets/)
CONTAINER_DIR := $(realpath $(SCRIPT_DIR)/../container/)
PRODIGAL_TRAINING_DIR := $(ASSETS_DIR)/prodigal_training_files
# The root folder where the pipeline is currently located. To be mounted into
# the Singularity containers below.
MNT_ROOT := /$(shell readlink -f . | cut -d"/" -f2)
INSTALL_LOG := "$(SCRIPT_DIR)/.install.log"

# Declare that these commands don't produce outputs
.PHONY: print_vars install update_amrfinderplus update_mlst_db update_blast_db update_virulencefinder_db

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

# NOTE: We don't run *_all for kpneumoniae, since prodigal training file is missing
install: build_containers \
	update_amrfinderplus \
	update_mlst_db \
	update_blast_db \
	update_virulencefinder_db \
	saureus_all \
	ecoli_all \
	kpneumoniae_unpack_cgmlst_schema \
	mtuberculosis_all

check:	check_chewbbaca \
	check_bwa \
	check_blastdb


# ==============================================================================
# Build containers
# ==============================================================================

build_containers:
	cd ../container \
	&& make all

# ==============================================================================
# Update databases
# ==============================================================================

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
		$(CONTAINER_DIR)/amrfinderplus.sif \
		amrfinder_update \
		--database $(AMRFINDERDB_DIR) &>> $(INSTALL_LOG)

# -----------------------------
# Update MLST database
# -----------------------------
MLSTDB_DIR := $(ASSETS_DIR)/mlst_db

update_mlst_db: $(MLSTDB_DIR)/pubmlst/dbases.xml

$(ASSETS_DIR)/mlst_db/pubmlst/dbases.xml:
	$(call log_message,"Starting update of MLST database ...")
	cd $(MLSTDB_DIR) \
	&& bash mlst-download_pub_mlst.sh &>> $(INSTALL_LOG)

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
		bash $(MLSTDB_DIR)/mlst-make_blast_db.sh &>> $(INSTALL_LOG)

# -----------------------------
# Update VirulenceFinder database
# -----------------------------
VIRULENCEFINDERDB_DIR := $(ASSETS_DIR)/virulencefinder_db
KMA_DIR := $(ASSETS_DIR)/kma

update_virulencefinder_db: $(VIRULENCEFINDERDB_DIR)/stx.name

$(ASSETS_DIR)/virulencefinder_db/stx.name:
	$(call log_message,"Starting update of VirulenceFinder database")
	cd $(ASSETS_DIR)/kma \
	&& make \
	&& cd $(VIRULENCEFINDERDB_DIR) \
	&& export PATH=$(ASSETS_DIR)/kma:$$PATH \
	&& singularity exec	--bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index \
	&& cd $(ASSETS_DIR)/resfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index \
	&& cd $(ASSETS_DIR)/pointfinder_db \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 INSTALL.py \
		$(KMA_DIR)/kma_index

# ==============================================================================
# Download, index and prep reference genomes for organisms
# ==============================================================================

# -----------------------------
# S. Aureus
# -----------------------------
.PHONY: saureus_all \
	saureus_download_reference \
	saureus_index_reference \
	saureus_download_cgmlst_schema \
	saureus_unpack_cgmlst_schema \
	saureus_prep_cgmlst_schema

saureus_all: saureus_download_reference \
	saureus_index_reference \
	saureus_download_cgmlst_schema \
	saureus_unpack_cgmlst_schema \
	saureus_prep_cgmlst_schema

SAUR_GENOMES_DIR := $(ASSETS_DIR)/genomes/staphylococcus_aureus
SAUR_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/staphylococcus_aureus
SAUR_REFSEQ_ACC := NC_002951.2


saureus_download_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading S. Aureus reference genome ...")
	mkdir -p $(SAUR_GENOMES_DIR) \
	&& cd $(ASSETS_DIR)/.. \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(SAUR_REFSEQ_ACC) \
		-o $(SAUR_GENOMES_DIR)


saureus_index_reference: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt

$(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta.bwt: $(SAUR_GENOMES_DIR)/$(SAUR_REFSEQ_ACC).fasta
	$(call log_message,"Indexing S. Aureus reference genome ...")
	cd $(SAUR_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<)


saureus_download_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip

$(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip:
	$(call log_message,"Downloading S. Aureus cgMLST schema ...")
	mkdir -p $(SAUR_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(SAUR_CGMLST_DIR)/alleles \
	&& wget https://www.cgmlst.org/ncs/schema/141106/alleles/ \
		-O $@ \
		--no-check-certificate &> /dev/null \


saureus_unpack_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta

$(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta: $(SAUR_CGMLST_DIR)/alleles/cgmlst_141106.zip
	$(call log_message,"Unpacking S. Aureus cgMLST schema ...")
	cd $$(dirname $<) \
		&& unzip -Duq $$(basename $<)


saureus_prep_cgmlst_schema: $(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn

$(SAUR_CGMLST_DIR)/alleles_rereffed/Staphylococcus_aureus.trn: $(SAUR_CGMLST_DIR)/alleles/SACOL2694.fasta
	$(call log_message,"Prepping S. Aureus cgMLST schema ...")
	cd $(SAUR_CGMLST_DIR) \
	&& echo "WARNING! Prepping cgMLST schema. This takes a looong time. Put on some coffee" \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(SAUR_CGMLST_DIR)/alleles \
		-o $(SAUR_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Staphylococcus_aureus.trn

# -----------------------------
# E. coli
# -----------------------------

.PHONY: ecoli_all \
	ecoli_download_reference \
	ecoli_index_reference \
	ecoli_download_wgmlst_schema \
	ecoli_download_cgmlst_schema \
	ecoli_prep_ecoli_cgmlst_schema

ecoli_all: ecoli_index_reference \
	ecoli_download_wgmlst_schema \
	ecoli_prep_ecoli_cgmlst_schema

ECOLI_GENOMES_DIR := $(ASSETS_DIR)/genomes/escherichia_coli
ECOLI_WGMLST_DIR := $(ASSETS_DIR)/wgmlst/escherichia_coli
ECOLI_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/escherichia_coli
ECOLI_REFSEQ_ACC := NC_000913.3


ecoli_download_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading E. coli genome ...")
	cd $(ASSETS_DIR)/.. \
	&& mkdir -p $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(ECOLI_REFSEQ_ACC) \
		-o $(ECOLI_GENOMES_DIR)


ecoli_index_reference: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt

$(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta.bwt: $(ECOLI_GENOMES_DIR)/$(ECOLI_REFSEQ_ACC).fasta
	$(call log_message,"Indexing E. coli genome ...")
	cd $(ECOLI_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<)


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
		--latest


# Download Ecoli cgmlst cgmlst.org schema
ecoli_download_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip

$(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip:
	$(call log_message,"Downloading E. coli cgMLST schema ...")
	mkdir -p $(ECOLI_CGMLST_DIR)/alleles &> /dev/null \
	&& cd $(ECOLI_CGMLST_DIR)/alleles \
	&& wget -O $$(basename $@) https://www.cgmlst.org/ncs/schema/5064703/alleles/ --no-check-certificate &> /dev/null


# Unpack Ecoli cgmlst schema
ecoli_unpack_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles/b4383.fasta

$(ECOLI_CGMLST_DIR)/alleles/b4383.fasta: $(ECOLI_CGMLST_DIR)/alleles/ecoli_cgmlst_alleles_5064703.zip
	$(call log_message,"Unpacking E. coli cgMLST schema ...")
	cd $(ECOLI_CGMLST_DIR)/alleles \
	&& unzip -Duq $$(basename $<)


# Prepping Ecoli cgmlst cgmlst.org schema
ecoli_prep_ecoli_cgmlst_schema: $(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn

$(ECOLI_CGMLST_DIR)/alleles_rereffed/Escherichia_coli.trn: $(ECOLI_CGMLST_DIR)/alleles/b4383.fasta
	$(call log_message,"Prepping E. coli cgMLST schema ... WARNING: This takes a looong time. Put on some coffee")
	cd $(ECOLI_CGMLST_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(ECOLI_CGMLST_DIR)/alleles \
		-o $(ECOLI_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Escherichia_coli.trn


# -----------------------------
# K. pneumoniae
# -----------------------------

.PHONY: kpneumoniae_all \
	kpneumoniae_download_reference \
	kpneumoniae_index_reference \
	kpneumoniae_download_cgmlst_schema \
	kpneumoniae_prep_cgmlst_schema

kpneumoniae_all: kpneumoniae_download_reference \
	kpneumoniae_index_reference \
	kpneumoniae_download_cgmlst_schema \
	kpneumoniae_prep_cgmlst_schema


KPNEU_GENOMES_DIR := $(ASSETS_DIR)/genomes/klebsiella_pneumoniae
KPNEU_CGMLST_DIR := $(ASSETS_DIR)/cgmlst/klebsiella_pneumoniae
KPNEU_REFSEQ_ACC := NC_016845.1


kpneumoniae_download_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading K pneumoniae genome ...")
	cd $(ASSETS_DIR)/.. \
	&& mkdir -p $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(KPNEU_REFSEQ_ACC) \
		-o $(KPNEU_GENOMES_DIR)


kpneumoniae_index_reference: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt

$(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta.bwt: $(KPNEU_GENOMES_DIR)/$(KPNEU_REFSEQ_ACC).fasta
	$(call log_message,"Indexing K pneumoniae genome ...")
	cd $(KPNEU_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<)


# Download Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_download_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip

$(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip:
	$(call log_message,"Downloading K. pneumoniae cgMLST schema ...")
	mkdir -p $(KPNEU_CGMLST_DIR)/alleles \
	&& cd $(KPNEU_CGMLST_DIR)/alleles \
	&& wget -O $$(basename $@) https://www.cgmlst.org/ncs/schema/2187931/alleles/ --no-check-certificate &> /dev/null


kpneumoniae_unpack_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles/KP1_RS24625.fasta

$(KPNEU_CGMLST_DIR)/alleles/KP1_RS24625.fasta: $(KPNEU_CGMLST_DIR)/alleles/cgmlst_schema_2187931.zip
	$(call log_message,"Unpacking K. pneumoniae cgMLST schema ...")
	cd $(KPNEU_CGMLST_DIR)/alleles \
	&& unzip -Duq $$(basename $<)


# Prep Kpneumoniae cgmlst cgmlst.org schema
kpneumoniae_prep_cgmlst_schema: $(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn

$(KPNEU_CGMLST_DIR)/alleles_rereffed/Klebsiella_pneumoniae.trn:
	$(call log_message,"Prepping K. pneumoniae cgMLST schema ... Warning: This takes a looong time. Put on some coffee!")
	mkdir -p $(KPNEU_CGMLST_DIR)/alleles_rereffed \
	&& cd $(KPNEU_CGMLST_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/chewbbaca.sif \
		chewie PrepExternalSchema \
		-i $(KPNEU_CGMLST_DIR)/alleles \
		-o $(KPNEU_CGMLST_DIR)/alleles_rereffed \
		--cpu 1 \
		--ptf $(PRODIGAL_TRAINING_DIR)/Klebsiella_pneumoniae.trn #TODO: create Kpneumoniae training file

# -----------------------------
# M. tuberculosis
# -----------------------------
MTUBE_GENOMES_DIR := $(ASSETS_DIR)/genomes/mycobacterium_tuberculosis
MTUBE_REFSEQ_ACC := NC_000962.3

.PHONY: mtuberculosis_all \
	mtuberculosis_download_reference\
	mtuberculosis_index_reference

mtuberculosis_all: mtuberculosis_download_reference mtuberculosis_index_reference

mtuberculosis_download_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta:
	$(call log_message,"Downloading M. tuberculosis genome ...")
	mkdir -p $(MTUBE_GENOMES_DIR) \
	&& cd $(ASSETS_DIR)/.. \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/pythonScripts.sif \
		python3 bin/download_ncbi.py \
		-i $(MTUBE_REFSEQ_ACC) \
		-o $(MTUBE_GENOMES_DIR)

mtuberculosis_index_reference: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt

$(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta.bwt: $(MTUBE_GENOMES_DIR)/$(MTUBE_REFSEQ_ACC).fasta
	$(call log_message,"Indexing M. tuberculosis genome ...")
	cd $(MTUBE_GENOMES_DIR) \
	&& singularity exec --bind $(MNT_ROOT) $(CONTAINER_DIR)/bwakit.sif \
		bwa index $$(basename $<)

# ==============================================================================
# Perform checks
# ==============================================================================
# TODO: Integrate with main installation rules?

# -----------------------------
# Check chewBBACA
# -----------------------------
check_chewbbaca:
	@cd $(ASSETS_DIR)/.. \
	&& saureus=$(SAUR_CGMLST_DIR)/alleles_rereffed \
	&& ecoli=$(ECOLI_CGMLST_DIR)/alleles_rereffed \
	&& kpneumoniae=$(KPNEU_CGMLST_DIR)/alleles_rereffed \
	&& if [[ -d "$$saureus" && -d "$$ecoli" && -d "$$kpneumoniae" ]]; then \
		echo "[✓] PASSED check for chewBBACA: Directories exist:" \
		&& echo "- "$$saureus \
		&& echo "- "$$ecoli \
		&& echo "- $$kpneumoniae"; \
	else \
		echo "[!] FAILED check for chewBBACA: Some directories do not exist:"; \
		if [[ ! -d $$saureus ]]; then \
			echo "    Missing directory: $$saureus";  \
		fi; \
		if [[ ! -d $$ecoli ]]; then \
			echo "    Missing directory: $$ecoli";  \
		fi; \
		if [[ ! -d $$kpneumoniae ]]; then \
			echo "    Missing directory: $$kpneumoniae";  \
		fi; \
		echo "    Please report this in an issue on the JASEN repo: https://github.com/genomic-medicine-sweden/JASEN/issues"; \
	fi

# -----------------------------
# Check BWA
# -----------------------------
check_bwa:
	@cd $(ASSETS_DIR)/.. \
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
		echo "[✓] PASSED check for bwa: Indexes exist in $(SAUR_GENOMES_DIR)"; \
	else \
		echo "[!] FAILED check for bwa: Indexes do not exist in $(SAUR_GENOMES_DIR)"; \
	fi

# -----------------------------
# Check BlastDB
# -----------------------------
MLST_BLAST_DIR := $(ASSETS_DIR)/mlst_db/blast
check_blastdb:
	@cd $(ASSETS_DIR)/.. \
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
		echo "[✓] PASSED check for blast: Indexes exist in $(MLST_BLAST_DIR)"; \
	else \
		echo "[!] FAILED check for blast: Indexes do not exist in $(MLST_BLAST_DIR)"; \
	fi
