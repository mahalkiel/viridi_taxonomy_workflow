"""
Reads Based Metagenomics Workflow Snakefile
Before running this workflow:

Useful commands:
    I advise you consult the Snakemake docs.
    It's helpful to set up a default profile and run that
    rather than putting in many commands each run.
$conda activate snakemake
unlock in case there was an issue last run: 
$snakemake --unlock
dry run to check what files are missing/will be created if all goes well:
$snakemake -n
run: (modify the number of cores for your environment)
$nohup snakemake --use-conda --cores 200 &
@ Author: Michael Snyder
@ The Mackelprang Lab, California State University, Northridge
@ Created: 10th October 2022
"""

# python dependencies (in conda snakemake environment):
import os
import pickle

# global variables:
configfile: "config/config.yml"
DIRECTORY = config['DIRECTORY']
KEGG_DB = config['KEGG_DB']
SRCDIR = config['SRCDIR']
# make a symlink$ ln -s source_file myfile
DAA_FILES_DIR = config['DAA_FILES_DIR']
PROJECT = config['PROJECT']
REFSEQ_DB = config['REFSEQ_DB']
MEGAN_DB = config['MEGAN_DB']
DBCAN_OUTS = config['DBCAN_OUTS'].split()
DBCAN_DB = config['DBCAN_DB']
BT2_IDX = config['BT2_IDX']
READ = ["R1", "R2"]



def samplenamemaker(DIRECTORY):
    """
    Return the basenames for R1 files in DIRECTORY path as a list
    """
    my_file = [f for f in os.listdir(DAA_FILES_DIR) if os.path.isfile(os.path.join(DAA_FILES_DIR, f))]
    basename = []
    for fn in my_file:
        if fn.endswith('.fasta.daa'):
            index = fn.find('.fasta.daa')
            base = fn[:index]
            basename.append(base)
    return(basename)

def srcfile_finder(SRCDIR):
    """
    Return the basenames for R1 files in DIRECTORY path as a listresults/sam/{srcfile}.sam
    """
    my_file = [f for f in os.listdir(SRCDIR) if os.path.isfile(os.path.join(SRCDIR, f))]
    print(my_file)
    basename = []
    for fn in my_file:
        if fn.endswith('_R1.fastq'):
            index = fn.find('_R1.fastq')
            base = fn[:index]
            basename.append(base)
    return(basename)

# the basename of each sample is the {sample} wildcard
#sampleroots = samplenamemaker(DIRECTORY)
SRCFILES = srcfile_finder(SRCDIR)
print("pumpkins")
print(SRCFILES)

# path to outputs
# uncomment only the target files
# snakemake will run all needed steps to output them
rule all:
    input:
        ## DBCAN:
        #expand("results/dbcan_out/{sample}_{dbcanouts}", sample = sampleroots, dbcanouts = DBCAN_OUTS),
        ## DIAMOND VS REFSEQ:
        #expand("results/taxonomy_daa/{sample}.daa", sample = sampleroots),
        #expand("results/sam/{sample}.sam", sample = sampleroots),
        #expand("results/meganized_done/{sample}.done", sample = sampleroots),
        #"results/feature_tables/all_otu.csv",
        #expand("results/meganized_naive_done/{sample}.done", sample = sampleroots),
        #expand("results/tax_files/{sample}_lca80.csv", sample = sampleroots),
        #expand("results/tax_files/{sample}_lca50.csv", sample = sampleroots),
        #expand("results/tax_files/{sample}_ms100_tp1.csv", sample = sampleroots),
        #expand("results/tax_files/{sample}_default.csv", sample = sampleroots),
        #expand("results/tax_files/{sample}_me000001_sup1_tp01.csv", sample = sampleroots),
        expand("results/sam/{srcfile}.sam", srcfile = SRCFILES),
        #expand("results/taxonomy_m6/{sample}.m6", sample = sampleroots),
        #expand("results/taxonomy_out/{sample}.tax", sample = sampleroots)
        ## DIAMOND VS CAZYME DB ONLY:
        #expand("results/cazy_annotated/{sample}.cazy.m8", sample = sampleroots),q
        #expand("results/cazy_annotated/{sample}.outu.cazy.m8", sample = sampleroots)
        #expand("results/cazy_annotated/{sample}.filt.cazy.m8", sample = sampleroots),
        ## DIAMOND VS KEGG DATABASE:
        #expand("results/kegg_out/{sample}.m8", sample = sampleroots),
        #expand("results/kegg_out/{sample}.outu.m8", sample = sampleroots)
        ## CREATE FEATURETABLES:
        #"results/featuretables/tax_counts.csv"
        #"results/featuretables/kegg_outm_counts.txt",
        #"results/featuretables/kegg_outu_counts.txt",
        #"results/featuretables/cazy_u_counts.csv",
        #"results/featuretables/cazy_m_counts.csv"
        ## MEAGANIZE
        #expand("results/rma_files/{sample}.rma6", sample = sampleroots)

        
# there must be a better way to acheive this than to mv each file - with prefixes? with setting the output file?
# --out_dir dbcanouts 


rule diamond_cazy:
    output:
        out="results/cazy_annotated/{sample}.cazy.m8"
    input:
        fname="results/fasta_files/{sample}.cleaned.outm.fa"
    params:
        command="blastx",
        db=DBCAN_DB,
        max_target_seqs="1", #only report the best hit
        sensitivity="--fast",
        blocksize="20",
        index_chunks="1",
        minSupportPercent="0"
    log:
        outm="results/diamond_logs/cazy/{sample}_outm2.log"
    shell:
        """
        diamond {params.command} \
        -q {input.fname} \
        -d {params.db} -o {output} \
        --max-target-seqs {params.max_target_seqs} \
        -b {params.blocksize} \
        -c {params.index_chunks} \
        {params.sensitivity} >> {log.outm} 2>&1
        """


rule diamond_taxonomy:
    output:
        out="results/taxonomy_daa/{sample}.daa"
    input:
        fname=os.path.join(DIRECTORY, "{sample}.fasta")
    params:
        command="blastx",
        db=REFSEQ_DB,
        # sensitivity="--fast",
        blocksize="20",
        index_chunks="1",
        out_format="100"
    log:
        outu="results/diamond_logs/taxonomy/{sample}_daa.log",
    shell:
        """
        diamond {params.command} \
        -q {input.fname} \
        -d {params.db} -o {output} \
        -b {params.blocksize} \
        -c {params.index_chunks} \
        -f {params.out_format} \
        {params.sensitivity} >> {log.outu} 2>&1
        """


rule diamond_taxonomy_taxformat:
    output:
        out="results/taxonomy_out/{sample}.tax"
    input:
        fname=os.path.join(DIRECTORY, "{sample}.fasta")
    params:
        command="blastx",
        db=REFSEQ_DB,
        max_target_seqs="5", #only LCA 5 best hits
        sensitivity="--fast",
        blocksize="20",
        index_chunks="1",
        out_format="102"
    log:
        outu="results/diamond_logs/taxonomy/{sample}_tax.log",
    shell:
        """
        diamond {params.command} \
        -q {input.fname} \
        -d {params.db} -o {output} \
        --max-target-seqs {params.max_target_seqs} \
        -b {params.blocksize} \
        -c {params.index_chunks} \
        -f {params.out_format} \
        {params.sensitivity} >> {log.outu} 2>&1
        """


rule diamond_v_kegg:
    output:
        out="results/kegg_out/{sample}.m8"
    input:
        fname="results/fasta_files/{sample}.cleaned.outm.fa"
    params:
        command="blastx",
        db=KEGG_DB,
        max_target_seqs="1", #only report the best hit
        sensitivity="--fast",
        blocksize="20",
        index_chunks="1"
    threads: workflow.cores
    log:
        outm="results/diamond_logs/kegg/{sample}_outm.log"
    shell:
        """
        diamond {params.command} \
        -q {input} \
        -d {params.db} \
        -o {output} \
        -k {params.max_target_seqs} \
        -b {params.blocksize} \
        -c {params.index_chunks} \
        {params.sensitivity} >> {log.outm} 2>&1
        """


rule diamond_v_refseq:
    output:
        out="results/kegg_out/{sample}.daa",
    input:
        fname="results/fasta_files/{sample}.cleaned.outu.fa"
    params:
        command="blastx",
        db=KEGG_DB,
        max_target_seqs="1", #only report the best hit
        sensitivity="--fast",
        blocksize="20",
        index_chunks="1"
    threads: workflow.cores
    log:
        outu="results/diamond_logs/kegg/{sample}_outu.log"    
    shell:
        """
        diamond {params.command} \
        -q {input.fname} \
        -d {params.db} \
        -o {output} \
        -k {params.max_target_seqs} \
        -b {params.blocksize} \
        -c {params.index_chunks} \
        {params.sensitivity} >> {log.outu} 2>&1
        """


rule meganize_naive:
    output:
        out=touch("results/meganized_naive_done/{sample}.done")
    input:
        infile=os.path.join(DAA_FILES_DIR, "{sample}.fasta.daa")
    conda:
        "workflow/envs/megan_env.yml"
    params:
        lcaAlgorithm="naive",
        mapDB=MEGAN_DB,
        threads="16",
    shell:
        """
        daa-meganizer \
        -i {input.infile} \
        -alg {params.lcaAlgorithm}
        """


# changed minScore from 10 -> 40
# changed minScore from 40 -> 50
# changed topPercent 1 -> 0.1
# changed maxExpected 0.01 -> 0

rule meganize:
    output:
        out=touch("results/meganized_done/{sample}.done")
    input:
        infile=os.path.join(DAA_FILES_DIR, "{sample}.fasta.daa")
    conda:
        "workflow/envs/megan_env.yml"
    params:
        lcaAlgorithm="weighted",
        mapDB=MEGAN_DB,
        threads="16",
        maxExpected="0.000001",
        minSupportPercent="0",
        minSupport="1",
        # minScore="100",
        topPercent="0.01",
        # lcaCoveragePercent="100"
    shell:
        """
        daa-meganizer \
        -i {input.infile} \
        -alg {params.lcaAlgorithm} \
        -mdb {params.mapDB} \
        -me {params.maxExpected} \
        -sup {params.minSupport} \
        -top {params.topPercent} \
        -t {params.threads} 
        """
        # -supp {params.minSupportPercent} \
        # -ms {params.minScore} \
        # -lcp {params.lcaCoveragePercent}
        


rule daa2info:
    output:
        "results/tax_files/{sample}_me000001_sup1_tp01.csv"
    input:
        "results/meganized_done/{sample}.done"
    threads: 16
    params:
        class_to_count="Taxonomy",
        infile=os.path.join(DAA_FILES_DIR, "{sample}.fasta.daa")
    conda:
        "workflow/envs/megan_env.yml"
    shell:
        """
        daa2info -i {params.infile} -o {output} -c2c {params.class_to_count}
        """ 


rule bowtie_alignment:
#bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r> | --interleaved <i> | --sra-acc <acc> | b <bam>} -S [<sam>]
    input:
        in1=os.path.join(SRCDIR, "{srcfile}_R1.fastq"),
        in2=os.path.join(SRCDIR, "{srcfile}_R2.fastq"),
    output:
        "results/sam/{srcfile}.sam"
    conda:
        "workflow/envs/bowtie2-env.yml"
    params:
        bt2_idx = BT2_IDX
    shell:
        """
        bowtie2 -x {params.bt2_idx:q} -1 {input.in1:q} -2 {input.in2:q} -S {output:q}
        """

