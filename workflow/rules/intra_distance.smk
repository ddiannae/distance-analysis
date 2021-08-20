rule get_bin_distance_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-{bintype}-{binsize}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_plot_{binsize}.log" 
    script:
        "../scripts/getBinDistancePlot.R"

rule get_bin_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-bins-{bintype}-all-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_fitted_bins_{bintype}_all_{binsize}.log" 
    script:
        "../scripts/getFittedData.R"

rule get_bin_chr_plots:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-bins-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-bins-{bintype}-bychr-{binsize}.tsv",
        fitted=config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-{bintype}-bychr-{binsize}.png"
    params:
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_{bintype}_plot_bychr_{binsize}.log" 
    script:
        "../scripts/getBinDistanceByChrPlot.R"

rule get_bin_chr_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-bins-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-bins-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-bychr-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/get_fitted_bins_{bintype}_bychr_{binsize}.log" 
    script:
        "../scripts/getFittedDataByChr.R"

rule get_bins:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    output:
        by_chr=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-bins-{bintype}-bychr-{binsize}.tsv",
        all=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-bins-{bintype}-all-{binsize}.tsv"
    params:
        distance_dir=get_distance_dir,
        binsize="{binsize}",
        cond="{cond}",
        bintype="{bintype}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_get_bins_{bintype}_{binsize}.log" 
    threads: 8
    script:
      "../scripts/getBinStats.R"

rule get_intra_inter:
    input:
        network=getMIMatrix,
        log_file=config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/done.txt"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_get_intra_interactions.log" 
    params:
        distance_dir=get_distance_dir,
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 8
    script:
        "../scripts/getIntraInteractions.R"

