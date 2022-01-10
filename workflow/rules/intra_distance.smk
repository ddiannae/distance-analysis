rule get_distance_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/cancer-intra-interactions-by_window-count.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/cancer-intra-interactions-by_cytoband-count.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/normal-intra-interactions-by_window-count.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/normal-intra-interactions-by_cytoband-count.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-distance-"+str(config["distbin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-distance-bychr-"+str(config["distbin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-size-"+str(config["sizebin"])+".png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/bin-size-bychr-"+str(config["sizebin"])+".png",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-intra-interactions-by_cytoband-tresults.tsv",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-intra-interactions-by_cytoband-tresults.tsv",
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-plots.txt"
    shell:
        "echo done > {output}"

rule get_intra_null_model:
    input:
        mi_matrix=getMIMatrix,
        intra_count=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-interactions-by_cytoband-count.tsv",
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-interactions-by_cytoband-tresults.tsv"
    params:
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        annot_cytoband="input/Biomart_Ensembl80_GRCh38_p2_regions.tsv"
    threads: 36
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_interactions_by_cytoband_tresults.log"
    script:
        "../scripts/intraInteractionsNullModel.R"

rule get_intra_count:
    input:
        mi_matrix=getMIMatrix
    output:
        tsv=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-interactions-{intra_region}-count.tsv",
        plot=config["datadir"]+"/{tissue}/"+config["figdir"]+"/{cond}-intra-interactions-{intra_region}-count.png"
    params:
        tissue="{tissue}",
        region_type="{intra_region}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData",
        annot_cytoband="input/Biomart_Ensembl80_GRCh38_p2_regions.tsv"
    threads: 18 
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_interactions_{intra_region}_count.log"
    script:
        "../scripts/intraInteractionsCount.R"

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
        "../scripts/binFitting.R"

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
        "../scripts/binFittingByChr.R"

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
    threads: 18
    script:
      "../scripts/binStats.R"

rule get_intra_interactions:
    input:
        mi_matrix=getMIMatrix,
        log_file=config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/done.txt"
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-all-distance-mi.tsv"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_interactions.log" 
    params:
        distance_dir=get_distance_dir,
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 18
    script:
        "../scripts/intraInteractions.R"

