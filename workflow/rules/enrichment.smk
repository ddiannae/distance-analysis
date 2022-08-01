rule get_all_enrichments:
    input:
        expand(config["datadir"]+"/{{tissue}}/"+config["netdir"]+"_"+config["algorithm"]+"/enrichments/{enrch}-{cond}-comm-all-louvain-{{cutoff}}.tsv",
        enrch=["kegg","go","ncg","onco"], cond=["cancer","normal"])
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/enrichments/all-enrichments-{cutoff}.txt"
    shell:
        "echo done > {output}"

rule get_other_enrichments:
    input:
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/communities/{cond}-comm-{ctype}-louvain-{cutoff}.tsv",
        universe=getGeneUniverse
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/enrichments/{etype}-{cond}-comm-{ctype}-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/{etype}_enrichment_{cond}_{ctype}_{cutoff}.log"
    params:
        enrich_type="{etype}"
    script:
        "../scripts/other_enrichments.R"

rule get_go_enrichments:
    input:
        membership=config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/communities/{cond}-comm-{ctype}-louvain-{cutoff}.tsv",
        universe=getGeneUniverse
    output:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/enrichments/go-{cond}-comm-{ctype}-{cutoff}.tsv",
    log:
        config["datadir"]+"/{tissue}/"+config["netdir"]+"_"+config["algorithm"]+"/log/go_enrichment_{cond}_{ctype}_{cutoff}.log"
    script:
        "../scripts/go_enrichments.R"
