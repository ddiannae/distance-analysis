## Snakefile for ASCAT2 files from GDC
##
## Tissue type just like in GDC, lowercase is fine
### DONE_IMAGES = ["esophagus", "skin", "breast", "kidney","lung", "liver", "prostate", "pancreas", "bladder", "brain", "colorectal", "uterus", "thyroid"]  
#
#TISSUES = ["breast", "prostate", "pancreas", "bladder", "skin", "brain", "liver", "esophagus",  "lung", "kidney", "colorectal", "uterus", "thyroid"]
#"bladder", "brain", "breast", "colorectal", "esophagus", "kidney", "liver", "lung", "pancreas", "prostate", "skin", "thyroid", "uterus"
import glob
TISSUES = ["bladder"]
DATADIR ="/datos/ot/diana/regulacion-trans"
INTERCUT=18168
# Adjust the umber of cores according to the machine and number of tissues
MCCORES = 75
files = [] 
for t in TISSUES:
  ## Example:data/breast/manifests/breast-cancer-rna_counts.txt"
  files.append(DATADIR+ "/" + t + "/network-tables/normal-interactions.tsv")
  files.append(DATADIR+ "/" + t + "/network-tables/normal-vertices.tsv")
  files.append(DATADIR+ "/" + t + "/network-tables/cancer-interactions.tsv")
  files.append(DATADIR+ "/" + t + "/network-tables/cancer-vertices.tsv")
  #files.append(DATADIR+"/" + t+ "/manifests/"+t+"-cancer-rna_counts.txt")
  #files.append(DATADIR+"/" + t+ "/manifests/"+t+"-normal-rna_counts.txt")
  #files.append(DATADIR+"/" + t + "/plots/normalization_plots.pdf")
  #files.apped(DATADIR+ "/" + t + "/rdata/normalization_results.tsv")

#network=glob_wildcards(DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_network_{pval}.sif"),
def getNetworkSif(wildcards):
  return [file for file in glob.glob(DATADIR+"/" + wildcards["tissue"] + "/*_*_*_networks/" + wildcards["type"] + "_network_*.sif")]

rule all:
  input:
    files

rule get_min_networks:
  input:
    DATADIR+"/{tissue}/network-tables/{type}-interactions.tsv"
  output:
    DATADIR+"/{tissue}/network-tables/{type}-interactions-"+INTERCUT+".tsv"
  shell:
    """
    head -n INTERCUT {wildcards.input} > {wildcards.output}
    """

rule get_network_tables:
  input:
    network=getNetworkSif,
    annot=DATADIR+"/{tissue}/rdata/annot.RData"
  output:
    interactions=DATADIR+"/{tissue}/network-tables/{type}-interactions.tsv",
    vertices=DATADIR+"/{tissue}/network-tables/{type}-vertices.tsv"
  shell:
    """
    mkdir -p DATADIR/{wildcards.tissue}/network-tables
    Rscript getNetworkTables.R {input.annot} {input.network} {output.interactions} {output.vertices}
    """
    
rule run_aracne_all:
  input: 
    matrix=DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_{type}.tsv",
    genelist=DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
  output:
    directory(DATADIR+"/{tissue}/{step1}_{step2}_{step3}_networks/{type}_adj_{pval}"),
  shell:
    """
    export ARACNEHOME={ARACNEHOME}
    mkdir -p {output}
    python src/python/aracne-par.py {input.matrix} {input.genelist} {MCCORES} {wildcards.pval} {output} > {DATADIR}/{wildcards.tissue}/log/aracne_{wildcards.type}_{wildcards.pval}.log
    """

rule user_normalization:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData"
  output:
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn.RData",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_cancer.tsv",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_arsyn_normal.tsv",
    DATADIR+"/{tissue}/rdata/{step1}_{step2}_{step3}_norm_cpm10_genelist.txt"
  shell:
    "Rscript src/userNormalization.R {wildcards.tissue} {DATADIR} {wildcards.step1} {wildcards.step2} {wildcards.step3} > {DATADIR}/{wildcards.tissue}/log/{wildcards.step1}_{wildcards.step2}_{wildcards.step3}normalization_plots.log" 

rule normalization_plots:
  input:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  output:
    DATADIR+"/{tissue}/plots/normalization_plots.pdf"
  shell:
    "Rscript src/normalizationPlots.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/log/normalization_plots.log" 

rule normalization_test:
  input:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  output:
    DATADIR+"/{tissue}/rdata/normalization_results.tsv"
  shell:
    "Rscript src/normalizationTest.R {wildcards.tissue} {DATADIR} {MCCORES} > {DATADIR}/{wildcards.tissue}/log/normalization_test.log"

rule filter_low_expression:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData",
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  output:
    DATADIR+"/{tissue}/rdata/mean10_proteinCoding.RData",
  shell:
    "Rscript src/filterLowExpression.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/filter_low.log"

rule qc:
  input:
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  output:
    DATADIR+"/{tissue}/plots/pca_score_raw.png"
  shell:
    "Rscript src/QC.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/qc.log"
    
## We need to run these two together because the output of the download_files
## tasks depends on the manifest and there is no easy way to specify this on 
## snakemake
rule download_files_and_get_ascat_matrix:
  input:
    ## Manifest file
    normal=DATADIR+"/{tissue}/manifests/{tissue}-normal-rna_counts.txt",
    cancer=DATADIR+"/{tissue}/manifests/{tissue}-cancer-rna_counts.txt",
  output: 
    DATADIR+"/{tissue}/{tissue}-matrix.tsv",
    DATADIR+"/{tissue}/{tissue}-samples.tsv",
    DATADIR+"/{tissue}/rdata/raw_full.RData"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-normal-rna
    mkdir -p {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-cancer-rna
    ./bin/gdc-client download -d {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-normal-rna -m {input.normal} --retry-amount 3
    ./bin/gdc-client download -d {DATADIR}/{wildcards.tissue}/raw/{wildcards.tissue}-cancer-rna -m {input.cancer} --retry-amount 3
    Rscript src/getRawMatrices.R {wildcards.tissue} {DATADIR} > {DATADIR}/{wildcards.tissue}/log/get-matrix.log
    """
rule get_manifest:
  output:
    ## Example: data/breast/manifests/breast-cancer-rna_counts.txt"
    DATADIR+"/{tissue}/manifests/{tissue}-{type}-rna_counts.txt"
  shell:
    """
    mkdir -p {DATADIR}/{wildcards.tissue}/log
    mkdir -p {DATADIR}/{wildcards.tissue}/manifests
    python src/queryGDC.py {wildcards.tissue} {wildcards.type} {DATADIR}/{wildcards.tissue}/manifests false 
    """
