rule setup_log:
    output:
        config["datadir"]+"/{tissue}/"+config["distdir"]+"/log/done.txt"
    shell:
        """
        touch {output}
        """
