rule get_distance_plots:
    input:
        #config["datadir"]+"/{tissue}/"+config["figdir"]+"/cancer-intra-interactions-count.png",
        #config["datadir"]+"/{tissue}/"+config["figdir"]+"/normal-intra-interactions-count.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-distance-"+str(config["distbin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-distance-bychr-"+str(config["distbin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-size-"+str(config["sizebin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-size-bychr-"+str(config["sizebin"])+".png",
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-plots.txt"
    shell:
        "echo done > {output}"

rule get_intra_windows_cout:
    input:
        intra_mi=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    output:
        tsv=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-interactions-count.tsv",
        plot=config["datadir"]+"/{tissue}/"+config["figdir"]+"/{cond}-intra-interactions-count.png"
    params:
        tissue="{tissue}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 10
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_interactions_count.log"
    script:
        "../scripts/intraWindowsCount.R"

rule get_bin_distance_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-{bintype}-{binsize}.png"
    params:
        tissue = "{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{bintype}_plot_{binsize}.log" 
    script:
        "../scripts/binDistancePlot.R"

rule get_bin_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-bins-{bintype}-all-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-bins-{bintype}-all-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-all-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/fitted_bins_{bintype}_all_{binsize}.log" 
    script:
        "../scripts/fittedData.R"

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
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{bintype}_plot_bychr_{binsize}.log" 
    script:
        "../scripts/binDistanceByChrPlot.R"

rule get_bin_chr_fitted:
    input:
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-bins-{bintype}-bychr-{binsize}.tsv",
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-bins-{bintype}-bychr-{binsize}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/fitted-bins-{bintype}-bychr-{binsize}.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/fitted_bins_{bintype}_bychr_{binsize}.log" 
    script:
        "../scripts/fittedDataByChr.R"

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
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_bins_{bintype}_{binsize}.log" 
    threads: 8
    script:
      "../scripts/binStats.R"

rule get_intra_inter:
    input:
        network=getMIMatrix,
        log_file=config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/done.txt"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_interactions.log" 
    params:
        distance_dir=get_distance_dir,
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 8
    script:
        "../scripts/intraInteractions.R"

