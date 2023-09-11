# ==============================================================================
# Various definitions
# ==============================================================================
# Make sure bash is used as shell, for consistency and to make some more
# advanced scripting possible than with /bin/sh
SHELL := /bin/bash

all: \
	download_remote_containers \
	build_local_containers

# These are containers that will be built locally using the definition files in
# the current folder
local_containers := chewbbaca.sif \
					postAlignQc.sif \
					resfinder.sif \
					virulencefinder.sif \
					pythonScripts.sif

# These containers will be downloaded remotely
remote_containers := abritamr.sif \
					amrfinderplus.sif \
					blast.sif \
					bracken.sif \
					bwakit.sif \
					freebayes.sif \
					kraken2.sif \
					mlst.sif \
					mykrobe.sif \
					perl_json.sif \
					quast.sif \
					samtools.sif \
					skesa.sif \
					snippy.sif \
					sourmash.sif \
					spades.sif \
					tbprofiler.sif

# URLs to remote containers
URL_abritamr := https://depot.galaxyproject.org/singularity/abritamr:1.0.14--pyhdfd78af_0
URL_amrfinderplus := https://depot.galaxyproject.org/singularity/ncbi-amrfinderplus:3.11.11--h6e70893_0
URL_blast := https://depot.galaxyproject.org/singularity/blast:2.14.0--hf3cf87c_0
URL_bracken := https://depot.galaxyproject.org/singularity/bracken:2.8--py39hc16433a_0
URL_bwakit := https://depot.galaxyproject.org/singularity/bwakit:0.7.17.dev1--hdfd78af_1
URL_freebayes := https://depot.galaxyproject.org/singularity/freebayes:1.3.6--hbfe0e7f_2
URL_kraken2 := https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_3
URL_mlst := https://depot.galaxyproject.org/singularity/mlst:2.23.0--hdfd78af_1
URL_mykrobe := https://depot.galaxyproject.org/singularity/mykrobe:0.12.2--py39h70e0db4_0
URL_perl_json := https://depot.galaxyproject.org/singularity/perl-json:4.10--pl5321hdfd78af_0
URL_quast := https://depot.galaxyproject.org/singularity/quast:5.2.0--py310pl5321hc8f18ef_2
URL_samtools := https://depot.galaxyproject.org/singularity/samtools:1.17--h00cdaf9_0
URL_skesa := https://depot.galaxyproject.org/singularity/skesa:2.4.0--he1c1bb9_0
URL_snippy := https://depot.galaxyproject.org/singularity/snippy:4.6.0--hdfd78af_2
URL_sourmash := https://depot.galaxyproject.org/singularity/sourmash:4.8.2--hdfd78af_0
URL_spades := https://depot.galaxyproject.org/singularity/spades:3.15.5--h95f258a_1
URL_tbprofiler := https://depot.galaxyproject.org/singularity/tb-profiler:4.4.2--pyh7cba7a3_0

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
	version=$$(echo ${container_url} | grep -oP "([0-9]+\.)+[0-9]+") \
	&& output_file="${tool}_$${version}.sif" \
	&& echo "Downloading $$tool as $$output_file ..." \
	&& wget -O $$output_file ${container_url} --no-check-certificate \
	&& ln -sf $$output_file $@


# ----------------------------------------------------------------
# Build some containers locally (which are not available remotely)
# ----------------------------------------------------------------
build_local_containers: $(local_containers)

$(local_containers):
	$(call log_message,"Building local container $@ ...")
	$(eval tool := $(basename $@))
	@tool=${tool} \
	&& if [ $$(grep -q "$${tool^^}_VERSION=" $${tool} ) ]; then \
		version=$$(grep $${tool^^}"_VERSION=" $${tool} | sed "s/.*=//"); \
	else \
		version=$$(grep "VERSION " "${tool}" | sed "s/.* //"); \
	fi \
    && output_file="$${tool}_$${version}.sif" \
	&& singularity build --fakeroot --force "$${output_file}" "$${tool}" \
	&& ln -sf $${output_file} $@
