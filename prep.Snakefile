# --- Snakemake workflow: worm WGS early prep (trim → map to combined → mkdup → worm-only BAM) ---
# Produces: results/{sample}.worm.manta.bam (+ .bai)
# Rationale: map to COMBINED reference (worm + contaminants), mark dups, then
# filter to worm contigs and strip SA tags so downstream tools (e.g., Manta) don’t
# choke on supplementary alignments pointing to non-worm contigs.
#
# Edit the two paths below to your FASTA files. Place your FASTQs in this directory
# named as {sample}_1.fq.gz and {sample}_2.fq.gz.
#
# Usage:
#   snakemake -j 16 --use-conda     # if you add `conda:` blocks (here we assume tools on PATH)
# or
#   snakemake -j 16

# -------------------- CONFIG --------------------
COMBINED_FA = "/home/pcarlton/mapping/combined.fa"   # worm + E. coli + PhiX (etc.)
WORM_FA     = "/home/pcarlton/mapping/worm.fa"       # worm-only FASTA (e.g., CGC1)
THREADS     = 16

# Tools expected on PATH: fastp, bwa-mem2 (or bwa), samtools, sambamba, bedtools
# If you prefer minimap2, see the mapping rule comment.

# -------------------- SAMPLE DISCOVERY --------------------
# Samples are inferred from files like: {sample}_1.fq.gz and {sample}_2.fq.gz
from snakemake.io import glob_wildcards
SAMPLES = sorted(set(glob_wildcards("raw/{sample}_1.fq.gz").sample))

# -------------------- TARGETS --------------------
rule all:
    input:
        expand("results/{sample}.worm.manta.bam.bai", sample=SAMPLES),
        expand("qc/{sample}.combined.flagstat.txt", sample=SAMPLES),
        expand("qc/{sample}.combined.idxstats.txt", sample=SAMPLES)

# -------------------- REFERENCE INDEXING --------------------
rule faidx_worm:
    input:  WORM_FA
    output: WORM_FA + ".fai"
    shell:  "samtools faidx {input}"

rule bwa_index_combined:
    input:  COMBINED_FA
    output:
        expand(COMBINED_FA + ".{ext}", ext=["amb","ann","bwt.2bit.64","pac","0123"])  # bwa-mem2 index files
    shell:
        # For classic bwa, the index set is different; this is fine for bwa-mem2
        "bwa-mem2 index {input}"

# Build BED and contig list for worm to enable region filtering and header fixes
rule worm_regions:
    input:  WORM_FA + ".fai"
    output:
        bed="ref/worm.contigs.bed",
        lst="ref/worm.contigs.list"
    shell:
        (
        "mkdir -p ref; "
        # BED with [0, length) for each contig
        "awk -v OFS='\\t' '{{print $1,0,$2}}' {input} > {output.bed}; "
        # Plain list of contig names
        "cut -f1 {input} > {output.lst}"
        )

# -------------------- QC & TRIMMING --------------------
rule fastp_trim:
    input:
        r1="raw/{sample}_1.fq.gz",
        r2="raw/{sample}_2.fq.gz"
    output:
        r1="work/trim/{sample}_1.trim.fq.gz",
        r2="work/trim/{sample}_2.trim.fq.gz",
        html="qc/fastp/{sample}.fastp.html",
        json="qc/fastp/{sample}.fastp.json"
    threads: THREADS
    shell:
        (
        "mkdir -p work/trim qc/fastp; "
        "fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} "
        "--detect_adapter_for_pe --cut_right --cut_right_window_size 4 --cut_right_mean_quality 20 "
        "--length_required 30 -h {output.html} -j {output.json} --thread {threads}"
        )

# -------------------- MAPPING TO COMBINED --------------------
rule map_combined_bwa:
    input:
        idx=COMBINED_FA + ".bwt.2bit.64",
        r1="work/trim/{sample}_1.trim.fq.gz",
        r2="work/trim/{sample}_2.trim.fq.gz"
    output:
        bam="work/map/{sample}.combined.bam",
        bai="work/map/{sample}.combined.bam.bai"
    threads: THREADS
    params:
        rg=lambda wildcards: f"@RG\tID:{wildcards.sample}\tSM:{wildcards.sample}\tLB:lib1\tPL:ILLUMINA\tPU:lane"
    shell:
        (
        "mkdir -p work/map; "
        # If you prefer minimap2: minimap2 -ax sr -t {threads} {COMBINED_FA} ...
        "bwa-mem2 mem -t {threads} -R '{params.rg}' {COMBINED_FA} {input.r1} {input.r2} "
        "| samtools sort -@ {threads} -o {output.bam}; "
        "samtools index {output.bam}"
        )

# -------------------- MARK DUPLICATES (combined) --------------------
rule markdup_combined:
    input:  bam="work/map/{sample}.combined.bam"
    output: bam="work/map/{sample}.combined.mkdup.bam", bai="work/map/{sample}.combined.mkdup.bam.bai"
    threads: THREADS
    shell:
        (
        "sambamba markdup -t {threads} {input.bam} {output.bam}; "
        "samtools index {output.bam}"
        )

# -------------------- QC STATS ON COMBINED --------------------
rule flagstat_combined:
    input:  "work/map/{sample}.combined.mkdup.bam"
    output: "qc/{sample}.combined.flagstat.txt"
    threads: 4
    shell:  "mkdir -p qc; samtools flagstat {input} > {output}"

rule idxstats_combined:
    input:  "work/map/{sample}.combined.mkdup.bam"
    output: "qc/{sample}.combined.idxstats.txt"
    threads: 4
    shell:  "samtools idxstats {input} > {output}"

# -------------------- FILTER TO WORM + STRIP CROSS-SPECIES SA TAGS --------------------
# We first restrict alignments to worm contigs via -L, then drop any SA:Z tag entirely to avoid
# supplementary references to non-worm contigs confusing SV callers like Manta. If you want to
# preserve SA tags that point exclusively to worm contigs, replace the sed with a small Python filter.

rule worm_only_bam:
    input:
        bam="work/map/{sample}.combined.mkdup.bam",
        bed="ref/worm.contigs.bed",
        fai=WORM_FA + ".fai",
        contigs_list="ref/worm.contigs.list",
        sa_filter="scripts/filter_SA.py"
    output:
        bam="results/{sample}.worm.manta.bam",
        bai="results/{sample}.worm.manta.bam.bai"
    threads: THREADS
    params:
        tmp="work/map/{sample}.worm.tmp",
        hdr="work/map/{sample}.worm.header.sam"
    shell:
        (
        "mkdir -p work/map results; "

        # 1) keep only worm alignments (header preserved)
        "echo in step 1;"
        "samtools view -@ {threads} -h -b -L {input.bed} {input.bam} -o {params.tmp}.Lfilter.bam; "
        "echo finished step 1;"

        # 2) selectively filter SA: keep SA entries pointing to worm contigs; drop others
        "echo in step 2;"
        "python {input.sa_filter} {input.contigs_list} {params.tmp}.Lfilter.bam {params.tmp}.SAfiltered.bam; "
        "echo finished step 2;"

        # 3) fix mate fields after filtering (name-sort → fixmate → coord-sort)
        "echo in step 3.1;"
        "samtools sort -n -@ {threads} -o {params.tmp}.namesort.bam {params.tmp}.SAfiltered.bam; "

        "echo in step 3.2;"
        "samtools fixmate -r -@ {threads} {params.tmp}.namesort.bam {params.tmp}.fixmate.bam; "

        "echo in step 3.3;"
        "samtools sort -@ {threads} -o {params.tmp}.coord.bam {params.tmp}.fixmate.bam; "

        #"echo in step 3.4;"
        #"samtools index {params.tmp}.coord.bam; "

        "echo in step 3.5;"
        # 3b) enforce valid mate fields post-fixmate (handles PNEXT=0 with RNEXT='*')
        "python scripts/sanitize_mates.py {input.contigs_list} {params.tmp}.coord.bam {params.tmp}.coord.sane.bam; "
        "mv {params.tmp}.coord.sane.bam {params.tmp}.coord.bam; "

        "echo finished step 3;"

        # 4) rebuild header: drop old @SQ and regenerate from worm.fa.fai (preserve @HD/@RG/@PG)
        "echo in step 4;"

	# rebuild header text
	"(samtools view -H {params.tmp}.coord.bam | grep -v '^@SQ'; "
	" awk '{{print \"@SQ\\tSN:\"$1\"\\tLN:\"$2}}' {input.fai}) > {params.hdr}; "

	# rebuild BAM from SAM so TIDs are remapped safely
	"(cat {params.hdr}; samtools view {params.tmp}.coord.bam) "
	"| samtools view -@ {threads} -b -o {output.bam} -; "

	# index
	"samtools index -@ {threads} {output.bam}; "

        #"(samtools view -H {params.tmp}.coord.bam | grep -v '^@SQ'; "
        #" awk '{{print \"@SQ\\tSN:\"$1\"\\tLN:\"$2}}' {input.fai}) > {params.hdr}; "
        #"samtools reheader {params.hdr} {params.tmp}.coord.bam > {output.bam}; "
        #"cp {output.bam} bak.bam;"
        #"samtools index {output.bam}; "

        "echo finished step 4;"

        # cleanup
        #"rm -f {params.tmp}.* {params.hdr}"
        )


#rule worm_only_bam:
#    input:
#        bam="work/map/{sample}.combined.mkdup.bam",
#        bed="ref/worm.contigs.bed"
#    output:
#        bam="results/{sample}.worm.manta.bam",
#        bai="results/{sample}.worm.manta.bam.bai"
#    threads: THREADS
#    shell:
#        (
#        "mkdir -p results; "
#        # keep only alignments overlapping worm contigs, keep header (-h), then strip SA:Z tags
#        "samtools view -@ {threads} -h -b -L {input.bed} {input.bam} "
#        "| samtools view -h - "
#        "| sed 's/\tSA:Z:[^\t]*//g' "
#        "| samtools view -@ {threads} -b -o {output.bam}; "
#        "samtools index {output.bam}"
#        )

# -------------------- OPTIONAL: MOSDEPTH COVERAGE BINS (uncomment to use) --------------------
# rule mosdepth_worm:
#     input:  bam="results/{sample}.worm.manta.bam"
#     output: "qc/mosdepth/{sample}.regions.bed.gz"
#     threads: THREADS
#     shell:
#         (
#         "mkdir -p qc/mosdepth; "
#         "mosdepth -t {threads} --by 10000 qc/mosdepth/{wildcards.sample} {input.bam}"
#         )
