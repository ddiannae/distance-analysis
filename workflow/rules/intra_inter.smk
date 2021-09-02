rule get_intra_inter_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-onek-bins.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-log-bins.png"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-plots.txt"
    shell:
        "echo done > {output}"

rule get_intra_plot:
    input:
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-intra-inter-count-{bintype}-bins.tsv",
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-intra-inter-count-{bintype}-bins.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-{bintype}-bins.png"
    params:
        bintype="{bintype}",
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/intra_inter_count_{bintype}_plot.log" 
    script:
        "../scripts/intraInterPlot.R"

rule get_intra_inter_count:
    input:
        mi_matrix=getMIMatrix,
    output:
        onek_bins=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-onek-bins.tsv",
        log_bins=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-log-bins.tsv"
    params:
        cond="{cond}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_get_intra_inter_count.log" 
    script:
        "../scripts/intraInterCount.R"

