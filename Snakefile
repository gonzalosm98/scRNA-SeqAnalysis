configfile: "config.yaml"

import sys
from snakemake.utils import makedirs
from pathlib import Path
pipeline = "single-cell-data-processing"

include: "rules/create_file_log.smk"

###### SET OUTDIR AND WORKDIR ######

if "OUTDIR" in config:
    OUTDIR = config["OUTDIR"]
    workdir: config["OUTDIR"]
else:
    OUTDIR = os.getcwd()


###### CREATE LOGS SLURM DIR ######

makedirs("logs_slurm")


###### LOAD DATA FROM CONFIG FILE ######
# Basic options
DATA_DIR = config["DATA"]
CELLRANGER_PATH = config["CELLRANGER_PATH"]
PREFIX = config["PREFIX"]


# Extra settings for cellranger count
CR_COUNT_EXTRA = config["CR_COUNT_extra"]

# QC parameters
MITO_PERCENTAGE = config["MITO_PERCENTAGE"] # keep cells with less than X% mitochondrial read fraction
NUMBER_GENES_PER_CELL = config["NUMBER_GENES_PER_CELL"] # keep cells with more than X genes
NUMBER_UMI_PER_CELL = config["NUMBER_UMI_PER_CELL"] # keep cells with more than X UMIs
ENSEMBLE_BIOMART_SPECIES = config["ENSEMBLE_BIOMART_SPECIES"] # ensembl biomart species used to get the mitochondrial genes for that species

# Doublet removal - threshold doublet score
SCRUB_THRESHOLD = config['SCRUB_THRESHOLD']


###### GET SAMPLE NAMES FROM INPUT FILES ######

    SAMPLES, REST, = glob_wildcards(os.path.join(DATA_DIR, "{samp}_S1_L00{1,2}_R{1,2}_001.fastq.gz"))


# Create dictionary with sample names 
# If strings in SAMPLES are in a <sample>/<sample>.fastq format, the name of the directory ("sample") will be used as the sample name
samples_dict = {}
for el in SAMPLES:
    if "/" in el: 
    # if the string corresponds to directory/file.fastq.gz, split the string at the separator "/"
        op = el.split("/")
        # if the sample name op[0], the first element after the split, is already in the dictionary, add the 
        if op[0] in samples_dict:
            samples_dict[op[0]].append(op[1])
        else:
            samples_dict[op[0]] = [op[1]]
    else:
        samples_dict[el] = [el]


###### Create lists with cellranger count output ######

cellranger_count_outfiles = ["web_summary.html",
                            "metrics_summary.csv", 
                            "possorted_genome_bam.bam",
                            "possorted_genome_bam.bam.bai",
                            "filtered_feature_bc_matrix.h5",
                            "raw_feature_bc_matrix.h5",
                            "molecule_info.h5",
                            "cloupe.cloupe"]
cellranger_count_outdirs = ["filtered_feature_bc_matrix",
                            "raw_feature_bc_matrix",
                            "analysis"]


###### FUNCTIONS ######

def set_scrub_treshold(wildcards):
    """
    If SCRUB_THRESHOLD is set, return the treshold per sample
    """
    if not SCRUB_THRESHOLD:
        return("")
    else:
        return(SCRUB_THRESHOLD[wildcards.samples])


###### TARGETS ######

# Rules that don't need to run in the cluster
localrules:  create_file_log, remove_ambient_RNA, combine_cellranger_counter_metrics, get_mito_genes, edit_gtf, filter_GTF, QC, remove_doublets

# Target rules (desired output)
rule all:
    input:
        files_log,
        'cellranger_count_metrics_allsamples.tsv',
        expand('4_Doublets/{samples}_QC_doublets.h5ad', samples = samples_dict.keys()),


###### MAIN PIPELINE ######


"""
Cellranger count:
Can't specify output of the rule with the actual output files since cellranger count gives the error:
RuntimeError: <directory> is not a pipestance directory
if the cellranger count output directory was not created by cellranger count. 
By adding the output files to the output of the rule, snakemake would create the necessary
directories, and therefore trigger the error
One option would be to keep the complete snakemake output and do "rm -r {{samples}}" before the
cellranger_count rule
"""

rule cellranger_count:
    input:
        DATA_DIR
    output:
        expand("{{samples}}/outs/{counts_out}",counts_out = cellranger_count_outfiles),
        directory(expand("{{samples}}/outs/{counts_outdirs}",counts_outdirs = cellranger_count_outdirs)),
        # touch("steps_done/cellranger_count_{samples}.done")
    message:
        'Rule {rule} processing'
    params:
        extra = CR_COUNT_EXTRA,
        transcriptome = os.path.join(workflow.basedir, f"{PREFIX}_genome"),
        fastqs = FASTQS_DIR,
        samples = lambda wildcards: ",".join(samples_dict[wildcards.samples]),
        cr_path = CELLRANGER_PATH
    shell:
        """
rm -r {wildcards.samples} 
{params.cr_path}/cellranger count \
{params.extra} \
--id={wildcards.samples} \
--transcriptome={params.transcriptome} \
--fastqs={params.fastqs} \
--sample={params.samples}
        """


rule remove_ambient_RNA:
    """
    For each sample, creates a notebook, with the analysis steps and results. 
    Results are saved to SoupX/<sample>
    These files will be created: matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
    """
    input:
        rules.cellranger_count.output
    output:
        '2_ambient_RNA_correction/Ambient_RNA_correction_{samples}.html',
        # expand("SoupX/{{samples}}/{soupxfile}", soupxfile = ["matrix.mtx.gz", "features.tsv.gz", "barcodes.tsv.gz"]), 
    message:
        'Rule {rule} processing'
    params:
        wd = os.getcwd(),
        input = "{samples}/outs/",
        sample = "{samples}"
    script:
        'scripts/remove_ambient_RNA.Rmd'

rule combine_cellranger_counter_metrics:
    """
    Combines cellranger count metrics for all samples in a single table
    """
    input:
        expand("{samples}/outs/{counts_out}",samples = samples_dict.keys(),counts_out = cellranger_count_outfiles),
        # expand("steps_done/cellranger_count_{samples}.done", samples = samples_dict.keys())
    output:
        'cellranger_count_metrics_allsamples.tsv'
    message:
        'Rule {rule} processing'
    script:
        'scripts/combine_cellrange_counter_metrics.R'

rule get_mito_genes:
    """
    Extract mitochondrial gene names for your species from Ensembl
    Will create two files: one with ensembl gene symbols and the other with gene IDs
    """
    output:
        ensembl = f"{ENSEMBLE_BIOMART_SPECIES}_mito_genes_ensembl.csv",
        symbol = f'{ENSEMBLE_BIOMART_SPECIES}_mito_genes_symbol.csv'
    message:
        'Rule {rule} processing'
    params:
        species =  ENSEMBLE_BIOMART_SPECIES
    run:
        import pandas as pd
        import scanpy as sc
        mito_ensembl_ids = sc.queries.mitochondrial_genes(params.species, attrname="ensembl_gene_id")
        mito_gene_ids = sc.queries.mitochondrial_genes(params.species, attrname="external_gene_name")
        mito_ensembl_ids.to_csv(output[0])  
        mito_gene_ids.to_csv(output[1])


rule QC:
    """
    For each sample, creates a notebook with QC steps and results.
    It will also create a data object file for each sample (_QC.h5ad) and
    a directory with the plots created
    """
    input:
        ambient_RNA = rules.remove_ambient_RNA.output,
        mito_genes_ensembl = rules.get_mito_genes.output.ensembl,
        mito_genes_symbol = rules.get_mito_genes.output.symbol
    output:
        "3_QC/{samples}_QC.h5ad"
    log:
        notebook = "3_QC/processed_notebook_{samples}.ipynb"
    params:
        mito_percentage = MITO_PERCENTAGE,
        number_genes_per_cell = NUMBER_GENES_PER_CELL,
        number_UMI_per_cell = NUMBER_UMI_PER_CELL,
        ensemble_biomart_species = ENSEMBLE_BIOMART_SPECIES,
        sample = "{samples}",
    message:
        'Rule {rule} processing'
    notebook:
        'scripts/QC_Scanpy.py.ipynb'


rule remove_doublets:
    """
    For each sample creates a notebook with steps to remove doublets and results.
    It will create a data object file (_doublets.h5ad) and a directory with the plots created
    """
    input:
        rules.QC.output
    output:
        '4_Doublets/{samples}_QC_doublets.h5ad'
    log:
        notebook = "4_Doublets/processed_notebook_{samples}.ipynb"
    message:
        'Rule {rule} processing'
    params:
        sample = "{samples}",
        scrub_threshold = lambda wildcards: set_scrub_treshold(wildcards) #treshold per sample defined in the configfile
    notebook:
        'scripts/Doublet_removal.py.ipynb'

onsuccess:
    print("Workflow finished")
