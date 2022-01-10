rule done_assortativity_plots:
    input:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/assortativity/cancer-comm-assort-enrich-0_terms-{cutoff}.png",
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/assortativity/normal-comm-assort-enrich-0_terms-{cutoff}.png",
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/assortativity/assort-{cutoff}.txt",
    shell:
        "echo done > {output}"

rule get_assortativities_plot:
    input:
        chr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-chr-assortativity-{cutoff}.tsv",
        expr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-expr-assortativity-{cutoff}.tsv",
        comm_info=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/communities/{cond}-comm-info-all-louvain-{cutoff}.tsv",
        vertices=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-vertices-{cutoff}.tsv",
        enrich=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/enrichments/go-{cond}-comm-all-{cutoff}.tsv"
    output:
        comm_assortativity_nterms=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-comm-assort-enrich-{terms}_terms_only-{cutoff}.tsv",
        comm_assortativity_tsv=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-comm-assort-enrich-{terms}_terms_all-{cutoff}.tsv",
        comm_assortativity_png=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"_plots/assortativity/{cond}-comm-assort-enrich-{terms}_terms-{cutoff}.png",
    params:
        tissue="{tissue}",
        cond="{cond}",
        terms="{terms}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{cond}_assortativity_plot_{terms}_terms_{cutoff}.log"
    script:
        "../scripts/assortativityEnrichmentPlot.R"

rule get_assortativities:
    input:
        interactions=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/{cond}-interactions-{cutoff}.tsv",
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/communities/{cond}-comm-all-louvain-{cutoff}.tsv",
        expression=getDEGFile
    output:
        chr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-chr-assortativity-{cutoff}.tsv",
        expr_assortativity=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/assortativity/{cond}-expr-assortativity-{cutoff}.tsv",
    params:
        cond="{cond}"
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{cond}_assortativity_{cutoff}.log"
    script:
        "../scripts/communitiesAssortativity.R"
