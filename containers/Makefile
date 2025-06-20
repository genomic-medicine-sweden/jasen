# N.B: See the Makefile in the main JASEN directory for a brief explanation of
# the Makefile syntax(!)
# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash
INSTALL_LOG := .install.log

all: \
	build_local_containers \
	download_remote_containers

# These are containers that will be built locally using the definition files in
# the current folder
local_containers := bonsai-prp.sif \
					gambitcore.sif

# These containers will be downloaded remotely
remote_containers := ncbi-amrfinderplus.sif \
					blast.sif \
					bracken.sif \
					bwakit.sif \
					chewbbaca.sif \
					emmtyper.sif \
					fastqc.sif \
					flye.sif \
					freebayes.sif \
					hostile.sif \
					htslib.sif \
					kraken2.sif \
					medaka.sif \
					minimap2.sif \
					mlst.sif \
					mykrobe.sif \
					nanoplot.sif \
					perl-json.sif \
					prodigal.sif \
					quast.sif \
					resfinder.sif \
					samtools.sif \
					sccmec.sif \
					seqtk.sif \
					serotypefinder.sif \
					shigapass.sif \
					ska2.sif \
					skesa.sif \
					snippy.sif \
					sourmash.sif \
					spades.sif \
					spatyper.sif \
					tb-profiler.sif \
					virulencefinder.sif

# URLs to Docker containers
DOCKER_bonsai-prp := docker://clinicalgenomicslund/bonsai-prp:1.3.1
DOCKER_gambitcore := docker://clinicalgenomicslund/gambitcore:0.0.2

# URLs to remote containers
URL_ncbi-amrfinderplus := https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.11.11--h6e70893_0
URL_blast := https://depot.galaxyproject.org/singularity/blast:2.14.0--hf3cf87c_0
URL_bracken := https://depot.galaxyproject.org/singularity/bracken:2.8--py39hc16433a_0
URL_bwakit := https://depot.galaxyproject.org/singularity/bwakit:0.7.17.dev1--hdfd78af_1
URL_chewbbaca := https://depot.galaxyproject.org/singularity/chewbbaca:3.3.2--pyhdfd78af_0
URL_emmtyper := https://depot.galaxyproject.org/singularity/emmtyper:0.2.0--py_0
URL_fastqc := https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0
URL_flye := https://depot.galaxyproject.org/singularity/flye:2.9.3--py39hd65a603_0
URL_freebayes := https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2
URL_hostile := https://depot.galaxyproject.org/singularity/hostile:2.0.0--pyhdfd78af_0
URL_htslib := https://depot.galaxyproject.org/singularity/htslib:1.21--h566b1c6_1
URL_kraken2 := https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_3
URL_medaka := https://depot.galaxyproject.org/singularity/medaka:2.0.1--py39hf77f13f_0
URL_minimap2 := https://depot.galaxyproject.org/singularity/minimap2:2.28--h577a1d6_4
URL_mlst := https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_1
URL_mykrobe := https://depot.galaxyproject.org/singularity/mykrobe:0.12.2--py39h70e0db4_0
URL_nanoplot := https://depot.galaxyproject.org/singularity/nanoplot:1.43.0--pyhdfd78af_1
URL_perl-json := https://depot.galaxyproject.org/singularity/perl-json:4.10--pl5321hdfd78af_0
URL_prodigal := https://depot.galaxyproject.org/singularity/prodigal:2.6.3--h7b50bb2_10
URL_quast := https://depot.galaxyproject.org/singularity/quast:5.2.0--py310pl5321hc8f18ef_2
URL_resfinder := https://depot.galaxyproject.org/singularity/resfinder:4.7.2--pyhdfd78af_0
URL_samtools := https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0
URL_sccmec := https://depot.galaxyproject.org/singularity/sccmec:1.2.0--hdfd78af_0
URL_seqtk := https://depot.galaxyproject.org/singularity/seqtk:1.4--he4a0461_2
URL_serotypefinder := https://depot.galaxyproject.org/singularity/serotypefinder:2.0.2--py312hdfd78af_1
URL_shigapass := https://depot.galaxyproject.org/singularity/shigapass:1.5.0--hdfd78af_0
URL_ska2 := https://depot.galaxyproject.org/singularity/ska2:0.3.10--h4349ce8_0
URL_skesa := https://depot.galaxyproject.org/singularity/skesa:2.5.1--hdcf5f25_0
URL_snippy := https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2
URL_sourmash := https://depot.galaxyproject.org/singularity/sourmash:4.8.2--hdfd78af_0
URL_spades := https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1
URL_spatyper := https://depot.galaxyproject.org/singularity/spatyper:0.3.3--pyhdfd78af_3
URL_tb-profiler := https://depot.galaxyproject.org/singularity/tb-profiler:6.3.0--pyhdfd78af_0
URL_virulencefinder := https://depot.galaxyproject.org/singularity/virulencefinder:2.0.4--hdfd78af_1

define log_message
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
	@echo "$$(date "+%Y-%m-%d %H:%M:%S"): $1" | tee -a $(INSTALL_LOG);
	@echo "--------------------------------------------------------------------------------" | tee -a $(INSTALL_LOG);
endef

# ==============================================================================
# Main rules
# ==============================================================================

print_local_containers:
	@for def in $(local_containers); do echo $$def; done

print_remote_containers:
	@for cont in $(remote_containers); do echo $$cont; done


# ----------------------------------------------------------------
# Download containers stored remotely
# ----------------------------------------------------------------
download_remote_containers: $(remote_containers)

$(remote_containers):
	$(call log_message,"Downloading remote container $@ ...")
	$(eval tool := $(basename $@))
	$(eval container_url := ${URL_${tool}})
	$(eval sanitized_url := $(shell echo ${container_url} | sed -e 's/%3A/:/g'))
	output_file=$$(echo ${sanitized_url} | sed -e 's#https://##' -e 's#/#-#g' -e 's#:#-#g').img \
	&& wget -O $${output_file} ${sanitized_url} --no-check-certificate \
	&& ln -sf $${output_file} $@ 


# ----------------------------------------------------------------
# Build some containers locally (which are not available remotely)
# ----------------------------------------------------------------
build_local_containers: $(local_containers)

$(local_containers):
	$(call log_message,"Building local container $@ ...")
	$(eval tool := $(basename $@))
	$(eval container_url := ${DOCKER_${tool}})
	@chmod o+w $(INSTALL_LOG) \
	&& output_file=$$(echo ${container_url} | sed -e 's#docker://##' -e 's#/#-#g' -e 's#:#-#g').img \
	&& apptainer build --force $${output_file} ${container_url} \
	&& ln -sf $${output_file} $@ 
