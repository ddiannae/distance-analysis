rule get_intra_comms_distance_plots:
    input:
        inter_normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal-interactions-{cutoff}.tsv",
        inter_cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer-interactions-{cutoff}.tsv",
        ver_normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal-vertices-{cutoff}.tsv",
        ver_cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer-vertices-{cutoff}.tsv",
        comm_normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal-comm-{cutoff}.tsv",
        comm_cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer-comm-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/distance_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["diameter","meandistance"],plottype=["boxplot","histogram"])
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/distance_communities_{cutoff}_plots.log"
    script:
        "../scripts/distancePlots.R" 

rule get_intra_comms_plots:
    input:
        comm_info_normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal-comm-info-{cutoff}.tsv",
        comm_info_cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer-comm-info-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/distance_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/intra_communities_{cutoff}_plots.log"
    script:
        "../scripts/intraCommunitiesPlots.R"
   
rule get_intra_comms:
    input:
        interactions=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}-vertices-{cutoff}.tsv"
    output:
        comm=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}-comm-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/{cond}-comm-info-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_intra_communities_{cutoff}.log"
    script:
        "../scripts/intraCommunities.R"

rule get_distribution_plots:
    input:
        normal=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/normal-interactions-{cutoff}.tsv",
        cancer=config["datadir"]+"/{tissue}/network_" + config["algorithm"] + "/cancer-interactions-{cutoff}.tsv"
    output:
        boxplot=config["datadir"]+"/{tissue}/distance_plots/mi-boxplot-network-{cutoff}.png",
        desity=config["datadir"]+"/{tissue}/distance_plots/mi-density-network-{cutoff}.png"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/network_{cutoff}_distribution_plots.log"
    script:
        "../scripts/MIDistributionPlots.R"

rule get_network_tables:
    input:
        mi_matrix=getMIMatrix,
        done=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/done.txt"
    output:
        interactions=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"] + "/{tissue}/network_" + config["algorithm"] + "/{cond}-vertices-{cutoff}.tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cutoff="{cutoff}",
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_network_table_{cutoff}.log"
    script:
        "../scripts/networkTables.R"
    