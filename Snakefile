configfile: "config.yaml"

samples=glob_wildcards(config["read_dir"]+"/{sample}_R1_001.fastq.gz")[0]
#print(f"Using samples {samples}.")

rule all:
	input:
	        "fastqc/multiqc_report.html",
	        "fastqc_trimmed/multiqc_report.html",
		"analyses/gtdbtk/",
		"analyses/coverm_relative_abundance_minimap2.txt",
		"analyses/dram/"

rule fastqc:
	input:
	        R1=config["read_dir"]+"/{sample}_R1_001.fastq.gz",
	        R2=config["read_dir"]+"/{sample}_R2_001.fastq.gz"
	output:
	        "fastqc/{sample}/"
	threads:
	        6
	log:
	        "logs/fastqc/{sample}.fastqc.log"
	shell:
	        "(fastqc -o {output} -t {threads} {input.R1} {input.R2}) 2> {log}"

rule multiqc:
	input:
	        expand("fastqc/{sample}/",sample=samples)
	output:
	        "fastqc/multiqc_report.html"
	params:
	        exec=config["path_to_multiqc"]
	shell:
	        '''cd fastqc; {params.exec} */*'''

rule trimmomatic:
	input:
	        R1=config["read_dir"]+"/{sample}_R1_001.fastq.gz",
	        R2=config["read_dir"]+"/{sample}_R2_001.fastq.gz"
	output:
	        paired_R1="trimmomatic/{sample}_paired_1.fastq.gz",
	        paired_R2="trimmomatic/{sample}_paired_2.fastq.gz",
	        unpaired_R1="trimmomatic/{sample}_unpaired_1.fastq.gz",
	        unpaired_R2="trimmomatic/{sample}_unpaired_2.fastq.gz",
	        summary="trimmomatic/{sample}_stats_summary.txt"
	params:
	        adapters=config["adapters_file"]
	threads:
	        36
	log:
	        "logs/trimmomatic/{sample}.trimmomatic.log"
	shell:
	        "(trimmomatic PE -threads {threads} -trimlog {log} -summary {output.summary} {input.R1} {input.R2} {output.paired_R1} {output.unpaired_R1} {output.paired_R2} {output.unpaired_R2} "
	        " ILLUMINACLIP:{params.adapters}:2:30:10 LEADING:3 TRAILING:20 SLIDINGWINDOW:06:20  MINLEN:100 -phred33) 2> {log}"

rule fastqc_trimmed:
	input:
	        paired_R1="trimmomatic/{sample}_paired_1.fastq.gz",
	        paired_R2="trimmomatic/{sample}_paired_2.fastq.gz",
	        unpaired_R1="trimmomatic/{sample}_unpaired_1.fastq.gz",
	        unpaired_R2="trimmomatic/{sample}_unpaired_2.fastq.gz",
	output:
	        "fastqc_trimmed/{sample}/"
	threads:
	        6
	log:
	        "logs/fastqc_trimmed/{sample}.fastqc.log"
	shell:
	        "(fastqc -o {output} -t {threads} {input.paired_R1} {input.paired_R2} {input.unpaired_R1} {input.unpaired_R2}) 2> {log}"

rule multiqc_trimmed:
	input:
	        expand("fastqc_trimmed/{sample}/",sample=samples)
	output:
	        "fastqc_trimmed/multiqc_report.html"
	params:
	        exec=config["path_to_multiqc"]
	shell:
	        '''cd fastqc_trimmed; {params.exec} */*'''

rule megahit:
	input:
	        R1="trimmomatic/{sample}_paired_1.fastq.gz",
	        R2="trimmomatic/{sample}_paired_2.fastq.gz"
	output:
	        "assemblies/megahit/{sample}/final.contigs.fa"
	threads:
	        36
	resources:
	        mem_gb=config["megahit_memory"]
	log:
	        "logs/assemblies/megahit/{sample}_megahit.log"
	params:
	        min_contig_length=config["min_contig_length"],
	        folder="assemblies/megahit/{sample}/"
	shell:
	        "(rm -rf {params.folder}; megahit -1 {input.R1} -2 {input.R2} -t {threads} -o {params.folder} ) 2> {log}"

rule spades:
	input:
	        R1="trimmomatic/{sample}_paired_1.fastq.gz",
	        R2="trimmomatic/{sample}_paired_2.fastq.gz"
	output:
	        "assemblies/spades/{sample}/contigs.fasta"
	threads:
	        36
	resources:
	        mem_gb=config["spades_memory"]
	log:
	        "logs/assemblies/spades/{sample}_spades.log"
	params:
	        folder="assemblies/spades/{sample}/"
	shell:
	        "(spades.py --meta -1 {input.R1} -2 {input.R2} -t {threads} -m {resources.mem_gb} -o {params.folder}) 2> {log}"

rule gunzip_reads:
	input:
	        R1="trimmomatic/{sample}_paired_1.fastq.gz",
	        R2="trimmomatic/{sample}_paired_2.fastq.gz"
	output:
	        R1=temp("trimmomatic/{sample}_paired_1.fastq"),
	        R2=temp("trimmomatic/{sample}_paired_2.fastq")
	shell:
	        "gunzip -c {input.R1} > {output.R1}; gunzip -c {input.R2} > {output.R2}"

rule megahit_metawrap:
	input:
	        assembly="assemblies/megahit/{sample}/final.contigs.fa",
	        R1="trimmomatic/{sample}_paired_1.fastq",
	        R2="trimmomatic/{sample}_paired_2.fastq"
	output:
	        a="INITIAL_BINNING_1k/megahit/{sample}/metabat2_bins/",
	        b="INITIAL_BINNING_1k/megahit/{sample}/maxbin2_bins/",
	        c="INITIAL_BINNING_1k/megahit/{sample}/concoct_bins/"
	params:
	        output_folder="INITIAL_BINNING_1k/megahit/{sample}/"
	threads:
	        36
	resources:
	        mem_gb=config["metawrap_memory"]
	log:
	        "logs/metawrap/megahit_{sample}.log"
	shell:
	        "rm -rf {params.output_folder}; (metawrap binning -o {params.output_folder} -t {threads} -m {resources.mem_gb} -a {input.assembly} "
	        "  --universal --metabat2 --maxbin2 --concoct {input.R1} {input.R2} ) 2> {log}"

rule spades_metawrap:
	input:
	        assembly="assemblies/spades/{sample}/contigs.fasta",
	        R1="trimmomatic/{sample}_paired_1.fastq",
	        R2="trimmomatic/{sample}_paired_2.fastq"
	output:
	        a="INITIAL_BINNING_1k/spades/{sample}/metabat2_bins/",
	        b="INITIAL_BINNING_1k/spades/{sample}/maxbin2_bins/",
	        c="INITIAL_BINNING_1k/spades/{sample}/concoct_bins/"
	params:
	        output_folder="INITIAL_BINNING_1k/spades/{sample}/"
	threads:
	        36
	resources:
	        mem_gb=config["metawrap_memory"]
	log:
	        "logs/metawrap/spades_{sample}.log"
	shell:
	        "rm -rf {params.output_folder}; (metawrap binning -o {params.output_folder} -t {threads} -m {resources.mem_gb} -a {input.assembly} "
	        "  --universal --metabat2 --maxbin2 --concoct {input.R1} {input.R2} ) 2> {log}"

rule megahit_metawrap_refinement:
	input:
	        a="INITIAL_BINNING_1k/megahit/{sample}/metabat2_bins/",
	        b="INITIAL_BINNING_1k/megahit/{sample}/maxbin2_bins/",
	        c="INITIAL_BINNING_1k/megahit/{sample}/concoct_bins/"
	output:
	        "BIN_REFINEMENT_1k_50_5/megahit/{sample}/metawrap_50_5_bins/"
	params:
	        folder="INITIAL_BINNING_1k/megahit/{sample}/",
	        output_folder="BIN_REFINEMENT_1k_50_5/megahit/{sample}/"
	threads:
	        36
	resources:
	        mem_gb=config["metawrap_memory"]
	log:
	        "logs/metawrap/{sample}_megahit_refine.log"
	shell:
	        "(metawrap bin_refinement -o {params.output_folder} -A {input.a} -B {input.b} -C {input.c} -t {threads} "
	        " -c 50 -x 5 ) 2> {log}"

rule spades_metawrap_refinement:
	input:
	        a="INITIAL_BINNING_1k/spades/{sample}/metabat2_bins/",
	        b="INITIAL_BINNING_1k/spades/{sample}/maxbin2_bins/",
	        c="INITIAL_BINNING_1k/spades/{sample}/concoct_bins/"
	output:
	        "BIN_REFINEMENT_1k_50_5/spades/{sample}/metawrap_50_5_bins/"
	params:
	        folder="INITIAL_BINNING_1k/spades/{sample}/",
	        output_folder="BIN_REFINEMENT_1k_50_5/spades/{sample}/"
	threads:
	        36
	resources:
	        mem_gb=config["metawrap_memory"]
	log:
	        "logs/metawrap/{sample}_spades_refine.log"
	shell:
	        "(metawrap bin_refinement -o {params.output_folder} -A {input.a} -B {input.b} -C {input.c} -t {threads} "
	        " -c 50 -x 5 ) 2> {log}"

rule symlink_refined_bins:
	input:
	        spades=expand("BIN_REFINEMENT_1k_50_5/spades/{sample}/metawrap_50_5_bins/",sample=samples),
	        megahit=expand("BIN_REFINEMENT_1k_50_5/megahit/{sample}/metawrap_50_5_bins/",sample=samples)
	output:
	        folder= "Bins_metawrap_50_5_all/"
	shell:
	        "cd {output.folder};"
	        "for folder in {input.spades}; do for file in ../${{folder}}/*.fa; do if [ -f $file ]; then  temp=$(basename $file); sample_name=$(echo $file | cut -d '/' -f4); ln -s $file ${{sample_name}}.spades.${{temp}}; fi; done; done;"
	        "for folder in {input.megahit}; do for file in ../${{folder}}/*.fa; do if [ -f $file ]; then temp=$(basename $file); sample_name=$(echo $file | cut -d '/' -f4); ln -s $file ${{sample_name}}.megahit.${{temp}}; fi; done; done;"

rule checkm2:
	input:
	        "Bins_metawrap_50_5_all/"
	output:
	        "checkm2/quality_report.tsv"
	log:
	        "logs/checkm2.log"
	params:
	        env=config["checkm2_env"],
		folder="checkm2"
	shell:
	        '''eval "$(conda shell.bash hook)"; conda activate {params.env};'''
	        '''(checkm2 predict -i {input} -o {params.folder} -x .fa --force --threads 32) 2> {log}'''

rule condense_checkm2_stats:
	input:
	        "checkm2/quality_report.tsv"
	output:
	        "checkm2_stats.csv"
	shell:
	        '''echo genome,completeness,contamination,N50,length > {output}; '''
	        '''cat {input} | cut -d$'\t' -f1,2,3,6,8 | grep bin | sed 's/_checkm2//g' | sed 's/\t/.fa,/' | sed "s/\t/,/g" >> {output}'''

rule drep:
	input:
	        folder="Bins_metawrap_50_5_all/",
	        stats="checkm2_stats.csv"
	output:
	        "analyses/dRep/dereplicated_genomes/"
	log:
	        "logs/drep.log"
	params:
	        output_folder="analyses/dRep/"
	shell:
	        "(dRep dereplicate -p 36 -comp 50 -con 5 -strW 0 {params.output_folder} -g {input.folder}* --genomeInfo {input.stats} ) 2> {log}"

rule gtdbtk:
	input:
		"analyses/dRep/dereplicated_genomes/"
	output:
		"analyses/gtdbtk/"
	log:
		"logs/gtdbtk.log"
	run:
		shell("(gtdbtk classify_wf --genome_dir {input} --extension fa  --force --out_dir {output}) 2> {log}")


rule coverm:
	input:
		assemblies="analyses/dRep/dereplicated_genomes/",
		R1=expand("trimmomatic/{sample}_paired_1.fastq.gz",sample=samples),
		R2=expand("trimmomatic/{sample}_paired_2.fastq.gz",sample=samples)
	output:
		"analyses/coverm_relative_abundance_minimap2.txt"
		
	params:
		read_folder="trimmomatic/"
	log:
		"logs/coverm.log"
	threads:
		36
	shell:
		"(coverm genome --coupled trimmomatic/*_paired*.fastq.gz --genome-fasta-files {input.assemblies}*.fa -x fa -p  minimap2-sr "
		"--proper-pairs-only --min-read-percent-identity 95 --min-read-aligned-length 75 -m relative_abundance "
		" -t {threads} > {output} ) 2> {log}"

rule dram:
	input:
		"analyses/dRep/dereplicated_genomes/"
	output:
		"analyses/dram/"
	threads:
		36
	params:
		dram_env=config["dram_env"]
	log:
		"logs/dram.log"
	shell:
		'''eval "$(conda shell.bash hook)"; conda activate {params.dram_env}; rm -rf {output}; (DRAM.py annotate -i '{input}*.fa' --min_contig_size 1000 --bit_score_threshold '''
		''' 50 --rbh_bit_score_threshold 250 --thread {threads} -o {output} ) 2> {log} '''


