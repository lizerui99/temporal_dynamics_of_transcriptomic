"""
1.合并不同lane
2.质控
3.比对
4.过滤q20
5.组装转录本、转录本定量
6.获取表，获取exp矩阵
7.修改id
"""
import os
path = './sample'
dirs = os.listdir(path)

SAMPLES = []
for file in dirs:
    if file.strip().split(".")[-1] == "gz":
        new_file = '_'.join(file.split('_')[0:2])
        SAMPLES.append(new_file)
    else:
        continue
SAMPLES = list(set(SAMPLES))

rule all:
    input:
        expand("06_transcripts_quant/{sample}/{sample}.gtf",sample=SAMPLES),
        expand("06_transcripts_quant/{sample}/{sample}.tab",sample=SAMPLES)

rule fastp:
    input:
        fwd= "sample/{sample}_R1.fq.gz",
        rev= "sample/{sample}_R2.fq.gz"
    output:
        fwd= "01_fastp/{sample}_R1.clean.fq.gz",
        rev= "01_fastp/{sample}_R2.clean.fq.gz"
    threads: 4
    priority: 0
    params:
        json = "01_fastp/{sample}.clean.json",
        html = "01_fastp/{sample}.clean.html"
    shell:
        "fastp -w {threads} -i {input.fwd} -I {input.rev} -o {output.fwd} -O {output.rev} -h {params.html} -j {params.json}"

rule mapping:
    input:
        fwd = "01_fastp/{sample}_R1.clean.fq.gz",
        rev = "01_fastp/{sample}_R2.clean.fq.gz"
    output:
        "02_mapping/{sample}.bam"
    priority: 1
    params:
        index_ref = "/mnt/maxuxu/data/genome/maize_v4_new/Zm-B73-REFERENCE-GRAMENE-4.0.fa.hisat2"
    log:
        "02_mapping/log/{sample}.log"
    shell:
        "hisat2 -p 4 --dta -x {params.index_ref} -1 {input.fwd} -2 {input.rev} --summary-file {log} |samtools view -bS > {output}"

rule rm_low_quality_mapping:
    input:
        "02_mapping/{sample}.bam"
    priority: 2
    output:
        "03_q20_mapping/{sample}.q20.sorted.bam"
    shell:
        "samtools view -q 20 -b {input} | samtools sort -o {output}"

rule assembly_based_ref:
    input:
        "03_q20_mapping/{sample}.q20.sorted.bam"
    output:
        gtf = "04_assembly_ref/{sample}.gtf"
    params:
        gtf = "/mnt/maxuxu/data/genome/maize_v4_new/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3",
        lable = "{sample}"
    priority: 3
    threads:1
    shell:
        "stringtie -p {threads} -G {params.gtf} -o {output.gtf} -l {params.lable} {input}"

rule combination_trans_assmblies:
    input:
        expand("04_assembly_ref/{sample}.gtf",sample=SAMPLES)
    output:
        "05_stringtie_merged/stringtie_merged.gtf"
    params:
        gtf = "/mnt/maxuxu/data/genome/maize_v4_new/Zm-B73-REFERENCE-GRAMENE-4.0_Zm00001d.2.gff3"
    shell:
        "stringtie --merge -G {params.gtf} -o {output} {input}"

rule transcripts_quantification:
    input:
        merge_gtf = "05_stringtie_merged/stringtie_merged.gtf",
        bam = "03_q20_mapping/{sample}.q20.sorted.bam",
    output:
        tab = "06_transcripts_quant/{sample}/{sample}.tab",
        gtf = "06_transcripts_quant/{sample}/{sample}.gtf",
    params:
        lable = "{sample}"
    shell:
        "stringtie -e -B -G {input.merge_gtf} -A {output.tab} -o {output.gtf} -l {params.lable} {input.bam}"

