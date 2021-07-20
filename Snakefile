## Configfile
configfile: "config.yaml"

import glob
 
files = [] 
for t in config["tissues"]:
  #files.append(config["datadir"]+ "/" + t + "/distance_plots/mi_boxplot_network_1e-8.png")
  files.append(f'{config["datadir"]}/{t}/{config["figdir"]}/all-bins-fixed-size.png')
  #files.append(f'{config["datadir"]}/{t}/network_tables_{config["algorithm"]}/cancer_vertices_{config["intercut"]}.tsv')
  #files.append(config["datadir"]+ "/" + t + "/distance_plots/comm_diameter_histogram_network_"+str(config["intercut"])+".png")
#   files.append(config["datadir"]+ "/" + t + "/network_tables_{config['algorithm']}/normal_interactions_"+str(config["intercut"])+".tsv")
#   files.append(config["datadir"]+ "/" + t + "/network_tables_{config['algorithm']}/cancer_interactions_"+str(config["intercut"])+".tsv")

def getMIMatrix(wildcards):
  if config["algorithm"] == "aracne":
  	return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] + "/*_*_*_networks/" + wildcards["type"] + "_network_1.adj")]
  elif config["algorithm"] == "infotheo":
  	return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] + "/*_*_*_mi/" + wildcards["type"] + "_mi_matrix.adj")]

rule all:
  input:
    files

rule get_bin_plots:
  input:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-fixed-size-all-"+str(config["sizebin"])+".tsv",
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-fixed-size-all-"+str(config["sizebin"])+".tsv",
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-fixed-distance-all-"+str(config["distbin"])+".tsv",
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-fixed-distance-all-"+str(config["distbin"])+".tsv"
  output:
    config["datadir"]+"/{tissue}/"+config["figdir"]+"/all-bins-fixed-distance.png",
    config["datadir"]+"/{tissue}/"+config["figdir"]+"/all-bins-fixed-size.png"
  run:
    shell(f'Rscript getBinDistancePlots.R {config["datadir"]}/{wildcards.tissue}/{config["distdir"]}/ {config["datadir"]}/{wildcards.tissue}/{config["figdir"]}/ {config["distbin"]} {config["sizebin"]}')

rule get_size_bin:
  input:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
  output:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-size-bychr-"+str(config["sizebin"])+".tsv",
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-size-all-"+str(config["sizebin"])+".tsv"
  run:
    shell(f'Rscript getBinStats.R {input} {config["datadir"]}/{wildcards.tissue}/{config["distdir"]}/ {wildcards.type} {config["sizebin"]} "size" {config["mccores"]}')

rule get_dist_bins:
  input:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
  output:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-distance-bychr-"+str(config["distbin"])+".tsv",
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-distance-all-"+str(config["distbin"])+".tsv"
  run:
    shell(f'Rscript getBinStats.R {input} {config["datadir"]}/{wildcards.tissue}/{config["distdir"]}/ {wildcards.type} {config["distbin"]} "distance" {config["mccores"]}')

rule get_intra_inter:
  input:
    network=getMIMatrix,
    annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
  output:
    config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
  run:
    shell(f'Rscript getIntraInteractions.R {input.annot} {input.network} {config["datadir"]}/{wildcards.tissue}/{config["distdir"]} {config["mccores"]}')

rule get_intra_comms_distance_plots:
  input:
    inter_normal=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/normal_interactions_{cutoff}.tsv",
    inter_cancer=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/cancer_interactions_{cutoff}.tsv",
    ver_normal=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/normal_vertices_{cutoff}.tsv",
    ver_cancer=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/cancer_vertices_{cutoff}.tsv",
    comm_normal=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/normal_comm_{cutoff}.tsv",
    comm_cancer=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/cancer_comm_{cutoff}.tsv",
  output:
    expand(config["datadir"]+"/{{tissue}}/distance_plots/comm_{plotfactor}_{plottype}_network_{{cutoff}}.png", plotfactor=["diameter","meandistance"],plottype=["boxplot","histogram"])
  shell:
    """
    Rscript getDistancePlots.R {input.inter_normal} {input.ver_normal} {input.comm_normal} {input.inter_cancer} {input.ver_cancer} {input.comm_cancer} {config["datadir"]}/{wildcards.tissue}/distance_plots {wildcards.cutoff}
    """

rule get_intra_comms_plots:
  input:
    ci_normal=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/normal_comm-info_{cutoff}.tsv",
    ci_cancer=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/cancer_comm-info_{cutoff}.tsv",
  output:
    expand(config["datadir"]+"/{{tissue}}/distance_plots/comm_{plotfactor}_{plottype}_network_{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
  shell:
    """
    Rscript getIntraCommunitiesPlots.R {input.ci_normal} {input.ci_cancer} {config["datadir"]}/{wildcards.tissue}/distance_plots {wildcards.cutoff}
    """
   
rule get_intra_comms:
  input:
    interactions=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/{type}_interactions_{cutoff}.tsv",
    vertices=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/{type}_vertices_{cutoff}.tsv"
  output:
    comms=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/{type}_comm_{cutoff}.tsv",
    comm_info=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/{type}_comm-info_{cutoff}.tsv",
  shell:
    """
    Rscript getIntraCommunities.R {input.interactions} {input.vertices} {output.comms} {output.comm_info}
    """

rule get_distribution_plots:
  input:
    normal=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/normal_interactions_{cutoff}.tsv",
    cancer=config["datadir"]+"/{tissue}/network_tables_{config['algorithm']}/cancer_interactions_{cutoff}.tsv"
  output:
    config["datadir"]+"/{tissue}/distance_plots/mi_boxplot_network_{cutoff}.png",
    config["datadir"]+"/{tissue}/distance_plots/mi_density_network_{cutoff}.png"
  shell:
    """
    mkdir -p {config["datadir"]}/{wildcards.tissue}/distance_plots
    Rscript getMIDistributionPlots.R {input.normal} {input.cancer} {config["datadir"]}/{wildcards.tissue}/distance_plots {wildcards.cutoff}
    """

rule get_network_tables:
  input:
    # It is not ok to use python code here because of snakemake's wildcards pattern matching 
    network=getMIMatrix,
    annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
  output:
    config["datadir"] + '/{tissue}/network_tables_' + config["algorithm"] + '/{type}_interactions_' + str(config["intercut"]) + '.tsv',
    config["datadir"] + '/{tissue}/network_tables_' + config["algorithm"] + '/{type}_vertices_' + str(config["intercut"]) + '.tsv'
  run:
    # It is ok to use python code here because all variables exist
    shell(f'mkdir -p {config["datadir"]}/{wildcards.tissue}/network_tables_{config["algorithm"]}')
    shell(f'Rscript getNetworkTables.R {input.annot} {input.network} {config["intercut"]} {config["datadir"]}/{wildcards.tissue}/network_tables_{config["algorithm"]} {config["mccores"]}')
    
