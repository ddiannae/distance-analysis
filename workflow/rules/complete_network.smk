rule get_bin_chr_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-fixed-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/all-bins-fixed-{bintype}-{binsize}.png"
    params:
        tissue="{wildcards.tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_plot_all_{binsize}.log" 
    script:
        "../scripts/getBinDistancePlot.R"

rule get_bin_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-fixed-{bintype}-all-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-fixed-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-{bintype}-all-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_fitted_data_all_{binsize}.log" 
    script:
        "../scripts/getFittedData.R"

rule get_bin_chr_plots:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-fixed-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-fixed-{bintype}-bychr-{binsize}.tsv",
        fitted=config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-fixed-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/all-bins-fixed-{bintype}-{binsize}.png",
    params:
        tissue="{wildcards.tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_plot_bychr_{binsize}.log" 
    script:
        "../scripts/getBinDistancByChrPlot.R"

rule get_bin_chr_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-fixed-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-fixed-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-{bintype}-bychr-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_fitted_data_bychr_{binsize}.log" 
    script:
        "../scripts/getFittedDataByChr.R"

rule get_bins:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    output:
        by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-{bintype}-bychr-{binsize}.tsv",
        all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-fixed-{bintype}-all-{binsize}.tsv"
    params:
        distance_dir=get_distance_dir,
        binsize="{wildcards.binsize}",
        cond="{wildcards.cond}",
        bintype="{wildcards.bintype}"
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

