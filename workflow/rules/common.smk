import glob

def get_output_files(wildcards):
    files = []
    for t in config["tissues"]:
        files.append(f'{config["datadir"]}/{t}/{config["netdir"]}_{config["algorithm"]}_plots/network-plots-{config["cutoff"]}.txt')
        files.append(f'{config["datadir"]}/{t}/{config["figdir"]}/intra-plots.txt')
        files.append(f'{config["datadir"]}/{t}/{config["figdir"]}/intra-inter-plots.txt')
    return files

def get_distance_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}/{config["distdir"]}' 

def getMIMatrix(wildcards):
    if config["algorithm"] == "aracne":
        return [file for file in glob.glob(config["datadir"]+"/" +
        wildcards["tissue"] + "/correlation/*_*_*_si-arsyn_" + wildcards["cond"] + "_mi.adj")]
    elif config["algorithm"] == "infotheo":
        return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] + "/*_*_*_mi/" + wildcards["type"] + "_mi_matrix.adj")]

