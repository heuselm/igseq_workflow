# Targeted IG Variable domain sequencing analysis pipeline
# heuselm@IMP 03.09.2020
##########################################################

# Step 1: Define variables
fastqpairedreads, = glob_wildcards("00_input/{reads}.fastq.gz")
# print(fastqpairedreads)
import pandas as pd
barcodes = pd.read_table('00_input/barcodes.txt', names=['sample_id','master_barcode', 'slave_barcode', 'species'])
SAMPLES = barcodes['sample_id']
SPECIES = barcodes['species']
read_num = ['R1', 'R2']


rule all:
	input: expand("05_finalclones/{sample}_clones_all.txt", sample=SAMPLES)


rule step1_checkout:
	input: 	expand("00_input/{reads}.fastq.gz", reads = fastqpairedreads)
	threads: 4
	singularity: "../../bin/containers/igseq_container.simg"
	output:
		#"01_checkout/complete_checkout.txt"
		expand("01_checkout/{sample}_{read}.fastq.gz", sample=SAMPLES, read=read_num)
	shell: 	"""
			# mkdir 01_checkout
			migec Checkout -cute 00_input/barcodes.txt {input} 01_checkout/
			# touch 01_checkout/complete_checkout.txt
			"""


rule step2_histogram:
	input: expand("01_checkout/{sample}_{read}.fastq.gz", sample=SAMPLES, read=read_num)
	output: "02_histogram/overseq.txt"
	singularity: "../../bin/containers/igseq_container.simg"
	shell: 	"""
			# mkdir 02_histogram
			migec Histogram 01_checkout/ 02_histogram/
			"""


rule step3_assembleReadsByUmi:
	input: "02_histogram/overseq.txt"
	output: expand("03_assemble/{sample}_{read}.t4.cf.fastq", sample=SAMPLES, read=read_num)
	threads: 6
	resources:
		mem_mb=24576
	singularity: "../../bin/containers/igseq_container.simg"
	shell:	"""
			ls -lisa 01_checkout/
			ls -lisa 02_histogram/
			java -Xmx24G -jar /usr/bin/migec/migec-1.2.9.jar AssembleBatch \
			--force-collision-filter --force-overseq 4 01_checkout/ 02_histogram/ 03_assemble/
			touch 03_assemble/assembleReadsByUmi.done
			ls -lisa 03_assemble/*
			"""


rule step4_alignToGermline:
	input: 
		read_1="03_assemble/{sample}_R1.t4.cf.fastq",
		read_2="03_assemble/{sample}_R2.t4.cf.fastq",
	output: 
		aligned="04_aligned/{sample}_R12.t4.cf.fastq.aligned.vdjca"
	params:
		spc = lambda wildcards: list(SPECIES[SAMPLES==wildcards.sample])
	threads: 6
	resources:
		mem_mb=24576
	singularity: "../../bin/containers/igseq_container.simg"
	shell: 	"""
			echo 'analyzing sample: {input}'
			echo 'current species: {params.spc}'
			java -Xmx24G -jar /usr/bin/mixcr/mixcr.jar align -p kaligner2 \
			-s {params.spc} -OreadsLayout=Collinear --trimming-window-size 4 \
			--trimming-quality-threshold 15 -r 04_aligned/alignmentReport_{wildcards.sample}.txt \
			{input.read_1} {input.read_2} {output.aligned}
			"""


rule step5_assembleFinalClones:
	input:
		expand("04_aligned/{{sample}}_R12.t4.cf.fastq.aligned.vdjca", sample=SAMPLES)
	output:
		expand("05_finalclones/{{sample}}.clns", sample=SAMPLES)
	threads: 6
	singularity: "../../bin/containers/igseq_container.simg"
	shell:	"""
			echo 'assembling final clones for: {input}'
			java -Xmx24G -jar /usr/bin/mixcr/mixcr.jar assemble \
			-r 05_finalclones/cloneassemblyReport{wildcards.sample}.txt \
			-OassemblingFeatures=VDJRegion \
			-OseparateByC=true -f \
			{input} {output}
			"""


rule step6_exports:
	input:
		"05_finalclones/{sample}.clns"
	output:
		"05_finalclones/{sample}_clones_all.txt"
	threads: 6
	singularity: "../../bin/containers/igseq_container.simg"
	shell:	"java -jar /usr/bin/mixcr/mixcr.jar exportClones -ot {input} {output}"

