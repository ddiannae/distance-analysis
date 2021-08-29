rule get_network_plots:
    input:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/comm-diameter-histogram-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/comm-size-histogram-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/mi-density-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/degree-distribution-{cutoff}.png",
    output:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/network-plots-{cutoff}.txt"
    shell:
        "echo done > {output}"

rule get_intra_comms_distance_plots:
    input:
        inter_normal=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/normal-interactions-{cutoff}.tsv",
        inter_cancer=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/cancer-interactions-{cutoff}.tsv",
        ver_normal=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/normal-vertices-{cutoff}.tsv",
        ver_cancer=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/cancer-vertices-{cutoff}.tsv",
        comm_normal=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/normal-comm-{cutoff}.tsv",
        comm_cancer=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/cancer-comm-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/network_"+config["algorithm"]+"_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["diameter","meandistance"],plottype=["boxplot","histogram"])
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/intra_communities_{cutoff}_distance_plots.log"
    script:
        "../scripts/intraCommunitiesDistancePlots.R"

rule get_intra_comms_plots:
    input:
        comm_info_normal=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/normal-comm-info-{cutoff}.tsv",
        comm_info_cancer=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/cancer-comm-info-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/network_"+config["algorithm"]+"_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/intra_communities_{cutoff}_plots.log"
    script:
        "../scripts/intraCommunitiesStatsPlots.R"
   
rule get_intra_comms:
    input:
        interactions=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        comm=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-comm-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-comm-info-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_intra_communities_{cutoff}.log"
    script:
        "../scripts/intraCommunities.R"

rule get_degree_distribution_plots:
    input:
        network_normal=config["datadir"]+"/{tissue}/rdata/normal_network_"+config["algorithm"]+"_{cutoff}.RData",
        network_cancer=config["datadir"]+"/{tissue}/rdata/cancer_network_"+config["algorithm"]+"_{cutoff}.RData",
    params:
        tissue="{tissue}"
    output:
        dd=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/degree-distribution-{cutoff}.png",
        cdd=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/degree-distribution-cumulative-{cutoff}.png",
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/degree_distribution_{cutoff}_plots.log"
    script:
        "../scripts/degreeDistributionPlots.R"

rule get_network_stats:
    input:
        interactions=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        network=config["datadir"]+"/{tissue}/rdata/{cond}_network_"+config["algorithm"]+"_{cutoff}.RData",
        network_stats=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-network-stats-{cutoff}.tsv",
        node_attributes=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-node-attributes-{cutoff}.tsv",
        edge_attributes=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-edge-attributes-{cutoff}.tsv"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_network_stats_{cutoff}.log"
    script:
        "../scripts/networkStats.R"

rule get_distribution_plots:
    input:
        normal=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/normal-interactions-{cutoff}.tsv",
        cancer=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/cancer-interactions-{cutoff}.tsv"
    output:
        boxplot=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/mi-boxplot-network-{cutoff}.png",
        density=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"_plots/mi-density-network-{cutoff}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/network_{cutoff}_distribution_plots.log"
    script:
        "../scripts/MIDistributionPlots.R"

rule get_network_tables:
    input:
        mi_matrix=getMIMatrix,
        done=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/done.txt"
    output:
        interactions=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cutoff="{cutoff}",
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/{cond}_network_table_{cutoff}.log"
    script:
        "../scripts/networkTables.R"
    
