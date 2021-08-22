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
        comm_info_normal=config["datadir"]+"/{tissue}/network_tables_" + config["algorithm"] + "/normal_comm-info_{cutoff}.tsv",
        comm_info_cancer=config["datadir"]+"/{tissue}/network_tables_" + config["algorithm"] + "/cancer_comm-info_{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/distance_plots/comm_{plotfactor}_{plottype}_network_{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/intra_communities_{cutoff}_plots.log"
    script:
        "../scripts/intraCommunitiesPlots.R"
   
rule get_intra_comms:
    input:
        interactions=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}_interactions_{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}_vertices_{cutoff}.tsv"
    output:
        comm=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}_comm_{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}_comm-info_{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_intra_communities_{cutoff}.log"
    script:
        "../intraCommunities.R"

rule get_distribution_plots:
    input:
        normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal_interactions_{cutoff}.tsv",
        cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer_interactions_{cutoff}.tsv"
    output:
        boxplot=config["datadir"]+"/{tissue}/distance_plots/mi_boxplot_network_{cutoff}.png",
        desity=config["datadir"]+"/{tissue}/distance_plots/mi_density_network_{cutoff}.png"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/network_{cutoff}_distribution_plots.log"
    script:
        "../scripts/MIDistributionPlots.R"

rule get_network_tables:
    input:
        mi_matrix=getMIMatrix,
        done=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/done.txt"
    output:
        interactions=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}_interactions_{cutoff}.tsv",
        vertices=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}_vertices_{cutoff}.tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cutoff="{cutoff}",
    threads: 8
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_network_table_{cutoff}.log"
    script:
        "../scripts/networkTables.R"
    
