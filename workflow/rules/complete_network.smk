rule get_intra_inter:
    input:
        network=getMIMatrix,
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_intra_interactions.log" 
    params:
        distance_dir=get_distance_dir,
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 8
    script:
        "../scripts/getIntraInteractions.R {input.annot} {input.network} {config["datadir"]}/{wildcards.tissue}/{config["distdir"]} {config["mccores"]}')

