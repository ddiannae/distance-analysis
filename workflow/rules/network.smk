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
        mi_matrix=getMIMatrix
    output:
        interactions=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}_interactions_" + str(config["intercut"]) + ".tsv",
        vertices=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}_vertices_" + str(config["intercut"]) + ".tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cut=config["intercut"],
    threads: 8
    log
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_network_table.log"
    script:
        "../scripts/networkTables.R"
# It is ok to use python code here because all variables exist
#shell(f'mkdir -p {config["datadir"]}/{wildcards.tissue}/network_tables_{config["algorithm"]}')
    
