import pandas as pd
import os

# dataset and samples
metadata = pd.read_table(config["metadata"], dtype={"sample": str}, sep = " ").set_index(["sample"], drop=False)
sample_set = metadata["sample"].tolist()

# target rule
rule all:
	input:
		count_summary = config['result_dir'] + "count.txt",
		alignment_summary = expand(config['result_dir'] + "{sample}_alignment_summary.out", sample = sample_set),
		contamination_summary = config['result_dir'] + "all_contamination_summary.out" if config["remove_contamination"] else [],
		dedup_summary = config['result_dir'] + "all_dedup_summary.out" if config["remove_duplicates"] else []
		
### 
# Duplicate Removal with cd-hit-dup (Optional) ####
###

if "remove_duplicates_method" not in config:
  config["remove_duplicates_method"] = "cd-hit-dup"

rule cd_hit_dup:
	input:
		read1 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read1"],
		read2 = lambda wildcards: config["fastq_dir"] + metadata.loc[wildcards.sample, "read2"]
	params:
		qsub = config['qsub_duplicateRemoval'],
		method = config["remove_duplicates_method"],
		result_dir = config['result_dir'] 
	output:
		read1_nodup= config['result_dir'] + "fastq_files/{sample}_1_nodup.fastq",
		read2_nodup= config['result_dir'] + "fastq_files/{sample}_2_nodup.fastq",
		dedup_summary = temp(config['result_dir'] + "{sample}_dedup_summary.out")
	shell:
		'''
		#!/bin/bash
	
		# Dedup search
		mkdir -p {params.result_dir}/fastq_files
		zcat {input.read1} > {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq
		zcat {input.read2} > {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq
		
		if [ "{params.method}" == "cd-hit-dup" ]
		then
			cd-hit-dup -i {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq -i2 {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq -o {output.read1_nodup} -o2 {output.read2_nodup} 
		else
			fastp -i {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq -I {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq -o {output.read1_nodup} -O {output.read2_nodup} --dedup
		fi

		# Dedup summary
		mkdir -p {params.result_dir}
		total_reads=$(wc -l < {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq)
		dedup_reads=$(wc -l < {output.read1_nodup})
		echo $((total_reads/4)) $((dedup_reads/4)) >> {output.dedup_summary}

		# remove temporary fastq files
		rm {params.result_dir}/fastq_files/{wildcards.sample}_1.fastq {params.result_dir}/fastq_files/{wildcards.sample}_2.fastq 
		rm -f {params.result_dir}/fastq_files/{wildcards.sample}_*.clstr
		'''
  
rule cd_hit_dup_summary:
	input:
		dedup_summary = expand(config['result_dir'] + "{sample}_dedup_summary.out", sample = sample_set) if config["remove_duplicates"] else [],
	params:
		qsub="",
		result_dir = config['result_dir']
	output:
		config['result_dir'] + "all_dedup_summary.out"
	script: 
		"scripts/dedup_summary.py"
	
### 
# Remove Contamination (Optional) ####
###

# select fastq files
def get_fastq_files_afterduplicateremoval(sample, readn):
        if config["remove_duplicates"]:
                return config['result_dir'] + "fastq_files/" + sample + "_" + readn + "_nodup.fastq"
        else:
                return config["fastq_dir"] + metadata.loc[sample, "read" + readn]

rule removeContamination:
	input:
		read1 = lambda wildcards: get_fastq_files_afterduplicateremoval(wildcards.sample, "1"),
		read2 = lambda wildcards: get_fastq_files_afterduplicateremoval(wildcards.sample, "2")
	params:
		qsub = config['qsub_removeContamination'],
		result_dir = config['result_dir'],
		contamination_reference = config['contamination_reference'],
		bowtie2_param = config["bowtie2_param"]
	output:
	  read1_norrna= config['result_dir'] + "fastq_files/{sample}_1_norrna.fastq",
	  read2_norrna= config['result_dir'] + "fastq_files/{sample}_2_norrna.fastq",
	  contReads_summary = temp(config['result_dir'] + "{sample}_contamination_summary.out")
	shell:
		'''
		mkdir -p {params.result_dir}/fastq_files
	  bowtie2 {params.bowtie2_param} --un-conc {params.result_dir}/fastq_files/{wildcards.sample}_%_norrna.fastq --al-conc {params.result_dir}/fastq_files/{wildcards.sample}_%_rrna.fastq -x {params.contamination_reference} -1 {input.read1} -2 {input.read2} 2>> {params.result_dir}/{wildcards.sample}_hisat2_summary.out
	
		# Contamination summary
		total_reads=$(grep "reads; of these:" {params.result_dir}/{wildcards.sample}_hisat2_summary.out | awk '{{print $1}}')
		contam_reads=$(grep "aligned concordantly exactly 1 time" {params.result_dir}/{wildcards.sample}_hisat2_summary.out | awk '{{print $1}}')
		echo $((total_reads)) $((contam_reads)) >> {output.contReads_summary}
		rm {params.result_dir}/{wildcards.sample}_hisat2_summary.out
		'''
	
rule contam_summary:
	input:
		contamination_summary = expand(config['result_dir'] + "{sample}_contamination_summary.out", sample = sample_set)
	params:
		qsub="",
		result_dir = config['result_dir']
	output:
	  config['result_dir'] + "all_contamination_summary.out"
	script: 
		"scripts/contamination_summary.py"
		
### 
# Align Reads ####
###

# select fastq files
def get_fastq_files_foralignment(sample, readn):
	if config["remove_contamination"]:
		return config['result_dir'] + "fastq_files/" + sample + "_" + readn + "_norrna.fastq"
	else:
		return get_fastq_files_afterduplicateremoval(sample, readn)

rule run_Hisat2:
	input:
		read1 = lambda wildcards: get_fastq_files_foralignment(wildcards.sample, "1"),
		read2 = lambda wildcards: get_fastq_files_foralignment(wildcards.sample, "2")
	params:
		qsub = config['qsub_hisat2'],
		result_dir = config['result_dir'],
		reference = config['reference'], 
		hisat2_param = config['hisat2_param'],
		hisat2_splice_sites = config['hisat2_splice_sites']
	output:
		sortedreads = config['result_dir'] + "bam_files/{sample}_sorted.bam", 
		alignment_summary = config['result_dir'] + "{sample}_alignment_summary.out"
	shell:
		'''
		mkdir -p {params.result_dir}/bam_files
		if [ "{params.hisat2_splice_sites}" == "" ]
		then
			hisat2-align-s {params.hisat2_param} -x {params.reference} -1 {input.read1} -2 {input.read2} -S {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam 2>> {output.alignment_summary}
		else
			hisat2-align-s {params.hisat2_param} --known-splicesite-infile {params.hisat2_splice_sites} -x {params.reference} -1 {input.read1} -2 {input.read2} -S {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam 2>> {output.alignment_summary}
		fi
		samtools sort {params.result_dir}/bam_files/{wildcards.sample}_aligned.sam -o {output.sortedreads}
		samtools index {output.sortedreads}
		'''
		
rule align_summary:
	input:
		alignedReads_summary = expand(config['result_dir'] + "{sample}_alignment_summary.out", sample = sample_set)
	params:
		qsub="",
		result_dir = config['result_dir']
	output:
	  temp(config['result_dir'] + "all_alignment_summary.out")
	script: 
		"scripts/align_summary.py"

### 
# FeatureCounts ####
###

if "FeatureCount_param" not in config:
  config["FeatureCount_param"] = "-p -g gene_name"

rule FeatureCount:
	input:
		bam_files = expand(config['result_dir'] + "bam_files/{sample}_sorted.bam", sample = sample_set)
	params:
		qsub = config['qsub_featureCounts'],
		FeatureCount_param = config["FeatureCount_param"],
		result_dir = config['result_dir'],
		annotation = config['annotation'],
		script_featurecounts = config['featurecounts']
	output:
		count_summary = config['result_dir'] + "count.txt"
	shell:
		'''
		featureCounts {params.FeatureCount_param} -a {params.annotation} -o {output.count_summary} {input.bam_files}
		'''
