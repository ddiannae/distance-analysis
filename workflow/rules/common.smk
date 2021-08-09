import glob

def get_distance_dir(wildcards):
    return f'{config["datadir"]}/{wildcards.tissue}/{config["distdir"]}' 

def getMIMatrix(wildcards):
  if config["algorithm"] == "aracne":
  	return [file for file in glob.glob(config["datadir"]+"/" +
    wildcards["tissue"] + "/networks/*_" + wildcards["type"] + "_network.adj")]
  elif config["algorithm"] == "infotheo":
  	return [file for file in glob.glob(config["datadir"]+"/" + wildcards["tissue"] + "/*_*_*_mi/" + wildcards["type"] + "_mi_matrix.adj")]

