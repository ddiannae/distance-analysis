rule get_size_bin:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
    output:
        by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-size-bychr-"+str(config["sizebin"])+".tsv",
        all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-size-all-"+str(config["sizebin"])+".tsv"
    params:
        distance_dir=get_distance_dir,
        sizebin=config["sizebin"],
        type="{wildcards.type}",
        bintype="size"
    threads: 8
    script:
      "../scripts/getBinStats.R"

rule get_dist_bins:
	input:
		config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-all-distance-mi.tsv"
	output:
		by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-distance-bychr-"+str(config["distbin"])+".tsv",
		all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{type}-fixed-distance-all-"+str(config["distbin"])+".tsv"
	params:
        distance_dir=get_distance_dir,
        sizebin=config["distbin"],
        type="{wildcards.type}",
        bintype="distance"
    threads: 8
  	script:
		"../scripts/getBinStats.R"

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

