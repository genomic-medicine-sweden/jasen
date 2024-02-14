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
local_containers := bonsai-prp.sif

# These containers will be downloaded remotely
remote_containers := ncbi-amrfinderplus.sif \
					blast.sif \
					bracken.sif \
					bwakit.sif \
					chewbbaca.sif \
					freebayes.sif \
					kraken2.sif \
					mlst.sif \
					mykrobe.sif \
					perl-json.sif \
					quast.sif \
					resfinder.sif \
					samtools.sif \
					skesa.sif \
					snippy.sif \
					sourmash.sif \
					spades.sif \
					tb-profiler.sif \
					virulencefinder.sif

# URLs to Docker containers
DOCKER_bonsai-prp := docker://clinicalgenomicslund/bonsai-prp:0.4.0

# URLs to remote containers
URL_ncbi-amrfinderplus := https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.11.11--h6e70893_0
URL_blast := https://depot.galaxyproject.org/singularity/blast:2.14.0--hf3cf87c_0
URL_bracken := https://depot.galaxyproject.org/singularity/bracken:2.8--py39hc16433a_0
URL_bwakit := https://depot.galaxyproject.org/singularity/bwakit:0.7.17.dev1--hdfd78af_1
URL_chewbbaca := https://depot.galaxyproject.org/singularity/chewbbaca:3.3.2--pyhdfd78af_0
URL_freebayes := https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2
URL_kraken2 := https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_3
URL_mlst := https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_1
URL_mykrobe := https://depot.galaxyproject.org/singularity/mykrobe:0.12.2--py39h70e0db4_0
URL_perl-json := https://depot.galaxyproject.org/singularity/perl-json:4.10--pl5321hdfd78af_0
URL_quast := https://depot.galaxyproject.org/singularity/quast:5.2.0--py310pl5321hc8f18ef_2
URL_resfinder := https://depot.galaxyproject.org/singularity/resfinder:4.4.2--pyhdfd78af_1
URL_samtools := https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0
URL_skesa := https://depot.galaxyproject.org/singularity/skesa:2.5.1--hdcf5f25_0
URL_snippy := https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2
URL_sourmash := https://depot.galaxyproject.org/singularity/sourmash:4.8.2--hdfd78af_0
URL_spades := https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1
URL_tb-profiler := https://depot.galaxyproject.org/singularity/tb-profiler:5.0.1--pyhdfd78af_1
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
	output_file=$$(echo ${container_url} | sed -e 's#https://##' -e 's#/#-#g' -e 's#:#-#g').img \
	&& wget -O $${output_file} ${container_url} --no-check-certificate \
	&& ln -sf $${output_file} $@ 


# ----------------------------------------------------------------
# Build some containers locally (which are not available remotely)
# ----------------------------------------------------------------
build_local_containers: $(local_containers)

$(local_containers):
ifneq ($(shell id -u), 0)
	$(error ERROR: You must run the build_local_containers rule as sudo! Please check the README.md file or cd into the container directory and run: `sudo make build_local_containers` followed by `make download_remote_containers`. Current UID: ${curr_uid});
else
	$(call log_message,"Building local container $@ ...")
	$(eval tool := $(basename $@))
	$(eval container_url := ${DOCKER_${tool}})
	@chmod o+w $(INSTALL_LOG) \
	&& output_file=$$(echo ${container_url} | sed -e 's#https://##' -e 's#/#-#g' -e 's#:#-#g').img \
	&& sudo singularity build --force $${output_file} ${container_url} \
	&& ln -sf $${output_file} $@ 
endif
