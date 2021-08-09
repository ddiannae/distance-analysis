rule get_size_bin:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    output:
        by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-size-bychr-"+str(config["sizebin"])+".tsv",
        all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-size-all-"+str(config["sizebin"])+".tsv"
    params:
        distance_dir=get_distance_dir,
        sizebin=config["sizebin"],
        cond="{wildcards.cond}",
        bintype="size"
    threads: 8
    script:
      "../scripts/getBinStats.R"

rule get_dist_bins:
	input:
		config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
	output:
		by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-distance-bychr-"+str(config["distbin"])+".tsv",
		all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-distance-all-"+str(config["distbin"])+".tsv"
	params:
        distance_dir=get_distance_dir,
        sizebin=config["distbin"],
        cond="{wildcards.cond}",
        bintype="distance"
    threads: 8
  	script:
		"../scripts/getBinStats.R"

rule get_intra_inter:
    input:
        network=getMIMatrix,
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_intra_interactions.log" 
    params:
        distance_dir=get_distance_dir,
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 8
    script:
        "../scripts/getIntraInteractions.R"

