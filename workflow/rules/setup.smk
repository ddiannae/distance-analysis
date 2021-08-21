rule setup_distance_log:
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/done.txt"
    shell:
        """
        touch {output}
        """

rule setup_network_log:
    output:
        config["datadir"]+"/{tissue}/network_"+config["algorithm"]+"/log/done.txt"
    shell:
        """
        touch {output}
        """
