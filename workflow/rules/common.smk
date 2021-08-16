import glob

def get_output_files(wildcards):
    files = []
    for t in config["tissues"]:
        files.append(f'{config["datadir"]}/{t}/{config["figdir"]}/bin-distance-bychr-{config["distbin"]}.png')
        files.append(f'{config["datadir"]}/{t}/{config["figdir"]}/bin-distance-{config["distbin"]}.png')
    return files

def get_distance_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}/{config["distdir"]}' 

def getMIMatrix(wildcards):
    if config["algorithm"] == "aracne":
        return [file for file in glob.glob(config["datadir"]+"/" +
        wildcards["tissue"] + "/networks/" + wildcards["cond"] + "_network.adj")]
    elif config["algorithm"] == "infotheo":
        return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] + "/*_*_*_mi/" + wildcards["type"] + "_mi_matrix.adj")]

