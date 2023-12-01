#libraries, packages
import yaml
import os
import time
from collections import OrderedDict #

timestamp = time.strftime("%d-%b-%Y-%H:%M", time.localtime())
# Define result folder with timestamp
result_folder = f"result{timestamp}"

# Load the samples from the YAML file
config_file = yaml.load(open("config.yaml"), Loader=yaml.FullLoader) #for config file

#config file parameters
working_dir = config_file['working_directory']
samples = config_file['samples']
trimmomatic_params = config_file.get('trimmomatic_params', {})
#penguin_flag = config_file['penguin_plot_enabled']

flag_penguin = True

rule all:
    input:
        #fastqc
        expand("{wd}/output/{rs}/{sample}/fastqc/{sample}_R1_fastqc.html", sample=samples, wd=working_dir, rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/fastqc/{sample}_R1_fastqc.zip", sample=samples, wd=working_dir,rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/fastqc/{sample}_R2_fastqc.html", sample=samples, wd=working_dir,rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/fastqc/{sample}_R2_fastqc.zip", sample=samples, wd=working_dir,rs=result_folder),
        #trimmomatic
        expand("{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_paired_forward.fastq", sample=samples, wd=working_dir, rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_paired_reverse.fastq", sample=samples, wd=working_dir, rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_unpaired_forward.fastq", sample=samples, wd=working_dir, rs=result_folder),
        expand("{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_unpaired_reverse.fastq", sample=samples, wd=working_dir, rs=result_folder),
        #generate_config_file
        expand("{wd}/output/{rs}/config_{ts}.yaml", wd=working_dir, ts=timestamp, rs=result_folder),
        # R - penguins
        expand("{wd}/output/{rs}/penguin_plot.pdf", wd=working_dir, rs=result_folder)


#This rule generates a copy of the config file in the output folder - for tracking of the parameters 
rule generate_config_file:
    output: "{wd}/output/{rs}/config_{ts}.yaml"
    run:
        with open(output[0], 'w') as f:
            ts = time.strftime("%d%b%Y", time.gmtime())
            yaml.dump(config_file, f, default_flow_style=False, Dumper=yaml.Dumper, sort_keys=False)
        print(f"Generated config file: {output[0]}")

#This rule takes the R1 and R2 fastq files runs quality control
#Outputs html and zip for each sample
rule fastqc:
    input:
        r1="{wd}/data/{sample}_R1.fastq",
        r2="{wd}/data/{sample}_R2.fastq"
    output:
        html_r1="{wd}/output/{rs}/{sample}/fastqc/{sample}_R1_fastqc.html",
        zip_r1="{wd}/output/{rs}/{sample}/fastqc/{sample}_R1_fastqc.zip",
        html_r2="{wd}/output/{rs}/{sample}/fastqc/{sample}_R2_fastqc.html",
        zip_r2="{wd}/output/{rs}/{sample}/fastqc/{sample}_R2_fastqc.zip"
    shell:
        """
        mkdir -p output/{wildcards.rs}/{wildcards.sample}/fastqc &&
        fastqc {input.r1} --outdir output/{wildcards.rs}/{wildcards.sample}/fastqc &&
        fastqc {input.r2} --outdir output/{wildcards.rs}/{wildcards.sample}/fastqc &&
        echo 'FastQC rule ran successfully for {wildcards.sample}'
        """

#This rule runs trimmomatic on the samples
rule trimmomatic:
    input:
        forward1="{wd}/data/{sample}_R1.fastq",
        reverse1="{wd}/data/{sample}_R2.fastq"
    output: 
        paired_forward="{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_paired_forward.fastq",
        unpaired_forward="{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_unpaired_forward.fastq", 
        paired_reverse="{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_paired_reverse.fastq",
        unpaired_reverse="{wd}/output/{rs}/{sample}/trimmomatic/{sample}_trimmed_data_unpaired_reverse.fastq"
    params:
        ADAPTERS= trimmomatic_params.get("ADAPTERS", "TruSeq2-SE.fa"),
        LEADING = trimmomatic_params.get("LEADING",3),
        TRAILING = trimmomatic_params.get("TRAILING", 3),
        SLIDINGWINDOW = trimmomatic_params.get("SLIDINGWINDOW", "4:15"),
        MINLEN = trimmomatic_params.get("MINLEN", 36)
    shell:
        """
        mkdir -p output/{wildcards.rs}/{wildcards.sample}/trimmomatic &&
        trimmomatic PE -threads 5 \
            "{input.forward1}" "{input.reverse1}" \
            "{output.paired_forward}" "{output.unpaired_forward}" \
            "{output.paired_reverse}" "{output.unpaired_reverse}" \
            ILLUMINACLIP:"{params.ADAPTERS}":2:30:10 \
            LEADING:{params.LEADING} TRAILING:{params.TRAILING} SLIDINGWINDOW:{params.SLIDINGWINDOW} MINLEN:{params.MINLEN}
        echo 'Trimmomatic rule ran successfully for {wildcards.sample}'
        """

#This rule plots some random stuff from the penguin dataset
rule penguin_plot: 
    input:
        peng_dataset ="{wd}/data/penguins.csv"
    output: "{wd}/output/{rs}/penguin_plot.pdf"
    shell:  "Rscript ./scripts/plot_things.r {input} {output}"
