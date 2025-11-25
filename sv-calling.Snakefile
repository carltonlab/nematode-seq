# --- Snakemake workflow: SV callers (Manta, GRIDSS, smoove/LUMPY) + SURVIVOR merge ---
# Inputs: duplicate-marked, worm-only BAMs from the prep workflow:
#         results/{sample}.worm.manta.bam (+ .bai)
# Config: set WORM_FA to your worm-only reference FASTA
# Usage:  snakemake -j 16

WORM_FA = "/home/pcarlton/mapping/worm.fa"   # e.g., C_elegans.CGC1.fa
THREADS = 16
CONFIGMANTA = "/home/pcarlton/bin/manta-1.6.0.centos6_x86_64/bin/configManta.py"

from snakemake.io import glob_wildcards
SAMPLES = sorted(set(glob_wildcards("results/{sample}.worm.manta.bam").sample))

rule all:
    input:
        expand("sv/{sample}/manta/results/variants/diploidSV.vcf.gz.tbi", sample=SAMPLES),
        expand("sv/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz.csi", sample=SAMPLES),
        expand("sv/{sample}/gridss/{sample}.gridss.vcf.gz.tbi", sample=SAMPLES),
        expand("sv/{sample}/merge/{sample}.SURVIVOR.vcf", sample=SAMPLES)

# Ensure FASTA index exists (required by all callers) and BWA index for GRIDSS
rule faidx_worm:
    input:  WORM_FA
    output: WORM_FA + ".fai"
    shell:  "samtools faidx {input}"

rule bwa_index_worm:
    input:  WORM_FA
    output:
        expand(WORM_FA + ".{ext}", ext=["amb","ann","pac","0123","bwt.2bit.64"])  # bwa-mem2 index
    shell:
        "bwa-mem2 index {input}"

# -------------------- MANTA --------------------
rule manta_config:
    input:
        bam="results/{sample}.worm.manta.bam",
        bai="results/{sample}.worm.manta.bam.bai",
        fai=WORM_FA + ".fai"
    output:
        cfg="sv/{sample}/manta/_configured.ok"
    threads: THREADS
    params:
        rundir=lambda wc: f"sv/{wc.sample}/manta"
    shell:
        (
        "mkdir -p {params.rundir}; "
        "{CONFIGMANTA} --bam {input.bam} --referenceFasta {WORM_FA} --runDir {params.rundir}; "
        "touch {output.cfg}"
        )

rule manta_run:
    input:
        cfg="sv/{sample}/manta/_configured.ok"
    output:
        vcf="sv/{sample}/manta/results/variants/diploidSV.vcf.gz",
        tbi="sv/{sample}/manta/results/variants/diploidSV.vcf.gz.tbi"
    threads: THREADS
    params:
        rundir=lambda wc: f"sv/{wc.sample}/manta"
    shell:
        "{params.rundir}/runWorkflow.py -m local -j {threads}"
        #"python {params.rundir}/runWorkflow.py -m local -j {threads}"

# -------------------- SMOOVE / LUMPY --------------------
rule smoove_call:
    conda: "envs/smoove.yml"
    input:
        bam="results/{sample}.worm.manta.bam",
        bai="results/{sample}.worm.manta.bam.bai",
        fai=WORM_FA + ".fai"
    output:
        vcf="sv/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz",
        tbi="sv/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz.csi"
    threads: THREADS
    shell:
        (
        "mkdir -p sv/{wildcards.sample}/smoove; "
        "smoove call --outdir sv/{wildcards.sample}/smoove --name {wildcards.sample} "
        "--fasta {WORM_FA} --genotype {input.bam} --processes {threads} "
        "&& bcftools index -f -c {output.vcf}"
        )
# -------------------- GRIDSS --------------------
rule gridss_call:
    input:
        bam="results/{sample}.worm.manta.bam",
        bai="results/{sample}.worm.manta.bam.bai",
        fai=WORM_FA + ".fai",
        bwaidx=WORM_FA + ".bwt.2bit.64"
    output:
        vcf="sv/{sample}/gridss/{sample}.gridss.vcf.gz",
        tbi="sv/{sample}/gridss/{sample}.gridss.vcf.gz.tbi",
        asm="sv/{sample}/gridss/{sample}.gridss.assembly.bam"
    threads: THREADS
    params:
        workdir=lambda wc: f"sv/{wc.sample}/gridss/work"
    shell:
        (
        "mkdir -p sv/{wildcards.sample}/gridss; "
        "gridss --reference {WORM_FA} --output {output.vcf} --assembly {output.asm} "
        "--workingdir {params.workdir} --threads {threads} --jvmheap 31g -w {params.workdir}/tmp "
        "{input.bam}; "
        "tabix -f -p vcf {output.vcf}"
        )

# -------------------- SURVIVOR MERGE (use UNZIPPED VCFs) --------------------

# Unzip each caller's VCF (SURVIVOR is more reliable on plain .vcf)
rule manta_unzip:
    input:
        "sv/{sample}/manta/results/variants/diploidSV.vcf.gz"
    output:
        "sv/{sample}/manta/results/variants/diploidSV.vcf"
    shell:
        "zcat {input} > {output}"

rule gridss_unzip:
    input:
        "sv/{sample}/gridss/{sample}.gridss.vcf.gz"
    output:
        "sv/{sample}/gridss/{sample}.gridss.vcf"
    shell:
        "zcat {input} > {output}"

rule smoove_unzip:
    input:
        "sv/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz"
    output:
        "sv/{sample}/smoove/{sample}-smoove.genotyped.vcf"
    shell:
        "zcat {input} > {output}"

# Prepare per-sample list of UNZIPPED VCFs to merge
rule survivor_list:
    input:
        manta="sv/{sample}/manta/results/variants/diploidSV.vcf",
        gridss="sv/{sample}/gridss/{sample}.gridss.vcf",
        smoove="sv/{sample}/smoove/{sample}-smoove.genotyped.vcf"
    output:
        lst="sv/{sample}/merge/vcf_list.txt"
    shell:
        (
        "mkdir -p sv/{wildcards.sample}/merge; "
        "printf '%s\n%s\n%s\n' {input.manta} {input.gridss} {input.smoove} > {output.lst}"
        )

rule survivor_merge:
    input:
        lst="sv/{sample}/merge/vcf_list.txt"
    output:
        vcf="sv/{sample}/merge/{sample}.SURVIVOR.vcf"
    threads: 4
    params:
        maxdist=1000,
        minsupp=2,
        taketype=1,
        takestrand=1,
        scale=0,
        minsize=50
    shell:
        (
        "SURVIVOR merge {input.lst} {params.maxdist} {params.minsupp} "
        "{params.taketype} {params.takestrand} {params.scale} {params.minsize} {output.vcf}"
        )

#
# -------------------- SURVIVOR MERGE --------------------
## Prepare per-sample list of VCFs to merge, then run SURVIVOR
#rule survivor_list:
#    input:
#        manta="sv/{sample}/manta/results/variants/diploidSV.vcf.gz",
#        gridss="sv/{sample}/gridss/{sample}.gridss.vcf.gz",
#        smoove="sv/{sample}/smoove/{sample}-smoove.genotyped.vcf.gz"
#    output:
#        lst="sv/{sample}/merge/vcf_list.txt"
#    shell:
#        (
#        "mkdir -p sv/{wildcards.sample}/merge; "
#        "printf '%s\n%s\n%s\n' {input.manta} {input.gridss} {input.smoove} > {output.lst}"
#        )
#
#rule survivor_merge:
#    input:
#        lst="sv/{sample}/merge/vcf_list.txt"
#    output:
#        vcf="sv/{sample}/merge/{sample}.SURVIVOR.vcf"
#    threads: 4
#    params:
#        maxdist=1000,
#        minsupp=2,
#        taketype=1,
#        takestrand=1,
#        scale=0,
#        minsize=50
#    shell:
#        (
#        # Merge with 1kb proximity, require >=2 callers, same type/strand, min size 50bp
#        "SURVIVOR merge {input.lst} {params.maxdist} {params.minsupp} {params.taketype} {params.takestrand} {params.scale} {params.minsize} {output.vcf}"
#        )
#
