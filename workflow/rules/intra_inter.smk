rule get_intra_inter_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-onek-bins.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-log-bins.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-onek-chunks.png",
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-log-chunks.png",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-intra-inter-null_model-log-bins.tsv",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-intra-inter-null_model-log-bins.tsv",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-intra-inter-null_model-onek-bins.tsv",
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-intra-inter-null_model-onek-bins.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-plots.txt"
    shell:
        "echo done > {output}"

rule get_intra_inter_null_model:
    input:
        mi_matrix=getMIMatrix,
        counts=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-{bintype}-bins.tsv",
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-null_model-{bintype}-bins.tsv"
    params:
        cond="{cond}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 18
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_intra_inter_null-model_{bintype}_bins.log" 
    script:
        "../scripts/interInteractionsNullModel.R"

rule get_intra_plot:
    input:
        normal=config["datadir"]+"/{tissue}/"+config["distdir"]+"/cancer-intra-inter-count-{bintype}-{bc}.tsv",
        cancer=config["datadir"]+"/{tissue}/"+config["distdir"]+"/normal-intra-inter-count-{bintype}-{bc}.tsv"
    output:
        config["datadir"]+"/{tissue}/"+config["figdir"]+"/intra-inter-count-{bintype}-{bc}.png"
    params:
        bintype="{bintype}",
        tissue="{tissue}"
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/intra_inter_count_{bintype}_{bc}_plot.log" 
    script:
        "../scripts/intraInterPlot.R"

rule get_intra_inter_count:
    input:
        mi_matrix=getMIMatrix
    output:
        onek_chunks=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-onek-chunks.tsv",
        onek_bins=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-onek-bins.tsv",
        log_bins=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-log-bins.tsv",
        log_chunks=config["datadir"]+"/{tissue}/"+config["distdir"]+"/{cond}-intra-inter-count-log-chunks.tsv"
    params:
        cond="{cond}",
        annot=config["datadir"]+"/{tissue}/rdata/annot.RData"
    threads: 18
    log:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/{cond}_get_intra_inter_count.log" 
    script:
        "../scripts/intraInterCount.R"

