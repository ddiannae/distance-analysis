rule get_network_plots_output:
    input:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/comm-diameter-histogram-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/comm-size-histogram-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/mi-density-network-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/degree-distribution-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/shared-vertices-{cutoff}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/network-plots-{cutoff}.txt"
    shell:
        "echo done > {output}"


rule get_cancer_normal_interaction:
    input:
        inter_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-interactions-{cutoff}.tsv",
        inter_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-interactions-{cutoff}.tsv",
        ver_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-vertices-{cutoff}.tsv",
        ver_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-vertices-{cutoff}.tsv",
    output:
        inter_normal_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal_only-interactions-{cutoff}.tsv",
        inter_cancer_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer_only-interactions-{cutoff}.tsv",
        inter_shared=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/shared-interactions-{cutoff}.tsv",
        ver_normal_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal_only-vertices-{cutoff}.tsv",
        ver_cancer_only=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer_only-vertices-{cutoff}.tsv",
        ver_shared=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/shared-vertices-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/networks_interaction_{cutoff}.log"
    script:
        "../scripts/networksComparison.R"

rule get_intra_comms_distance_plots:
    input:
        inter_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-interactions-{cutoff}.tsv",
        inter_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-interactions-{cutoff}.tsv",
        ver_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-vertices-{cutoff}.tsv",
        ver_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-vertices-{cutoff}.tsv",
        comm_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-comm-{cutoff}.tsv",
        comm_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-comm-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["diameter","meandistance"],plottype=["boxplot","histogram"])
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/intra_communities_{cutoff}_distance_plots.log"
    script:
        "../scripts/intraCommunitiesDistancePlots.R"

rule get_intra_comms_plots:
    input:
        comm_info_normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-comm-info-{cutoff}.tsv",
        comm_info_cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-comm-info-{cutoff}.tsv"
    output:
        expand(config["datadir"]+"/{{tissue}}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/comm-{plotfactor}-{plottype}-network-{{cutoff}}.png", plotfactor=["order","size", "density"],plottype=["boxplot","histogram"])
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/intra_communities_{cutoff}_plots.log"
    script:
        "../scripts/intraCommunitiesStatsPlots.R"
   
rule get_intra_comms:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        comm=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-comm-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-comm-info-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{cond}_intra_communities_{cutoff}.log"
    script:
        "../scripts/intraCommunities.R"

rule get_degree_distribution_plots:
    input:
        network_normal=config["datadir"]+"/{tissue}/rdata/normal_"+config["netdir"]+"_"+config["algorithm"]+"_{cutoff}.RData",
        network_cancer=config["datadir"]+"/{tissue}/rdata/cancer_"+config["netdir"]+"_"+config["algorithm"]+"_{cutoff}.RData",
    params:
        tissue="{tissue}"
    output:
        dd=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/degree-distribution-{cutoff}.png",
        cdd=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/degree-distribution-cumulative-{cutoff}.png",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/degree_distribution_{cutoff}_plots.log"
    script:
        "../scripts/degreeDistributionPlots.R"

rule get_network_stats:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    output:
        network=config["datadir"]+"/{tissue}/rdata/{cond}_"+config["netdir"]+"_"+config["algorithm"]+"_{cutoff}.RData",
        network_stats=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-network-stats-{cutoff}.tsv",
        node_attributes=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-node-attributes-{cutoff}.tsv",
        edge_attributes=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-edge-attributes-{cutoff}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{cond}_network_stats_{cutoff}.log"
    script:
        "../scripts/networkStats.R"

rule get_distribution_plots:
    input:
        normal=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/normal-interactions-{cutoff}.tsv",
        cancer=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/cancer-interactions-{cutoff}.tsv"
    output:
        boxplot=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/mi-boxplot-network-{cutoff}.png",
        density=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/mi-density-network-{cutoff}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/network_{cutoff}_distribution_plots.log"
    script:
        "../scripts/MIDistributionPlots.R"

rule get_network_tables:
    input:
        mi_matrix=getMIMatrix,
        done=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/done.txt"
    output:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        cutoff="{cutoff}",
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{cond}_network_table_{cutoff}.log"
    script:
        "../scripts/networkTables.R"
    
