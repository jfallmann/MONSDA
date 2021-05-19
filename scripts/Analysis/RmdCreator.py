#! /usr/bin/env python3

import sys, argparse, os, inspect, gzip, glob, re, logging
import traceback as tb
from collections import defaultdict

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)
from Logger import *
scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Create Rmd snippets for Summary from downstream analyses results')
    parser.add_argument("--files", dest='files', required=True, type=str, nargs='*', help="Names of input files" )
    parser.add_argument("--output", dest='output', required=True, type=str, help="output Rmd File" )
    parser.add_argument("--env", dest='env', required=True, type=str, help="used conda environment" )
    parser.add_argument("--loglevel", default='INFO', help="Log verbosity" )

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

#Class
class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

def makeoutdir(outdir):
    if not os.path.isabs(outdir):
        outdir =  os.path.abspath(outdir)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    return outdir

def pretty(d, indent=0):
    for key, value in d.items():
        try:
            print('\t' * indent + str(key))
        except:
            print("not printable")
        if isinstance(value, dict):
           pretty(value, indent+1)
        else:
           try:
               print('\t' * (indent+1) + str(value))
           except:
               print("not printable")


# Analyses
compile = {
    "DE":  "Differential Gene Expression",
    "DEU": "Differential Exon Usage",
    "DAS": "Differential Alternative Splicing",
    "DTU": "Differential Transcript Usage"
    }

def check_workflow(files):
    logid = scriptname+'.check_workflow: '
    files = [os.path.basename(file) for file in files]
    works = set()
    for file in files:
        log.info(logid+str(file))
        w = file.split("_")[0]
        if w not in ["Sig", "SigDOWN", "SigUP"]:
            works.add(w)
    log.info(logid+str(works))
    if len(works) > 1:
        log.error(logid+"several Workflows found: > "+str(works)+" < check for wrong file naming.")
        sys.exit(1)
    try:
        compile[list(works)[0]]
    except:
        log.error(logid+"workflow > "+str(works)+" < not found")
        sys.exit(1)
    log.debug(logid+'workflow found: '+str(works))
    return list(works)[0]

def create_file_tree(files):
    logid = scriptname+'.create_file_tree: '
    tree = NestedDefaultDict()
    for file in files:
        setting = os.path.basename(file).split(".", 1)[0]
        setting = setting.replace("Sig_","").replace("SigUP_","").replace("SigDOWN_","")
        if "SESSION" in file:
            WF, TOOL, COMBI, NAME = setting.split('_')
            tree[WF][TOOL][COMBI]["DataSet"][NAME] = file
            continue
        WF, TOOL, COMBI, COMP, RES, NAME = setting.split('_')
        EXT = os.path.basename(file).split(".", 1)[1]
        tree[WF][TOOL][COMBI][COMP][RES][NAME] = file
    log.info(pretty(tree))
    return tree

def integrate_table(file):
    chunk = f"""
```{{r, table:{file}, include = TRUE, echo = FALSE, message = TRUE, warning = FALSE}}
library(DT)
#table <- read.table(gzfile(paste(params$root,'{file}',sep='/')), header = TRUE)
table <- read.table(gzfile('{file}'), header = TRUE)
DT::datatable(table)
```

"""
    return chunk

def integrate_figures(files):
    listlines = list()
    img = list()
    listlines.append(f"```{{r, echo=FALSE, out.width='{100/len(files)}%', out.height='{100/len(files)}%',fig.cap='caption',fig.show='hold',fig.align='center'}}")
    counter = 1
    for file in files:
        f = file.replace("\n","")
        #listlines.append(f"fig.{counter} <- paste(params$root,'{f}', sep='/')")
        listlines.append(f"fig.{counter} <- '{f}'")
        img.append(f'path.expand(fig.{counter})')
        counter += 1
    listlines.append(f"knitr::include_graphics(c({','.join(img)}))")
    listlines.append("```\n\n")
    return "\n".join(listlines)

def integrate_list(file):
    listlines = list()
    elements = NestedDefaultDict()
    with open(file, "r") as listfile:
        for line in listfile:
            if "geneID" in line: continue
            id = line.split("\t")[0]
            name = line.split("\t")[1]
            file = line.split("\t")[2]
            elements[id][name][file]
    for id in elements.keys():
        listlines.append(f"{id} : {list(elements[id].keys())[0]} \n")
        for name in elements[id].keys():
            listlines.append(integrate_figures(list(elements[id][name])))
    return "\n".join(listlines)

def create_Rmd(files, output, env):
    logid = scriptname+'.create_Rmd: '
    outdir = os.path.dirname(output)
    makeoutdir(outdir)
    workflow = check_workflow(files)
    tree = create_file_tree(files)
    lines = list()

    # Add title
    lines.append(f"# {compile[workflow]} {{.tabset}} \n\n")

    # Create file structure
    for tool in tree[workflow].keys():
        lines.append(f"## {tool}  {{.tabset}} \n\n")

        for combi in tree[workflow][tool].keys():
            lines.append(f"### {combi}  {{.tabset}} \n\n")

            for comparison in tree[workflow][tool][combi].keys():
                lines.append(f"#### {comparison}  {{.tabset}} \n\n")

                if "figure" in tree[workflow][tool][combi][comparison].keys():
                    lines.append(f"##### FIGURES  {{.tabset}} \n\n")
                    lines.append(f"###### OVERVIEW \n\n")
                    for name in tree[workflow][tool][combi][comparison]["figure"].keys():
                        lines.append(f"\n<br />\n")
                        lines.append(f"\n####### {name} \n")
                        lines.append(integrate_figures([tree[workflow][tool][combi][comparison]["figure"][name]]))

                if "list" in tree[workflow][tool][combi][comparison].keys():
                    for name in tree[workflow][tool][combi][comparison]["list"].keys():
                        lines.append(f"\n###### {name} \n")
                        lines.append(f"\n<br />\n")
                        lines.append(integrate_list(tree[workflow][tool][combi][comparison]["list"][name]))

                if "table" in tree[workflow][tool][combi][comparison].keys():
                    lines.append(f"##### TABLES  \n")
                    for name in tree[workflow][tool][combi][comparison]["table"].keys():
                        lines.append(f"\n<br />\n")
                        lines.append(f"\n{name}  \n")
                        # lines.append(integrate_table(tree[workflow][tool][combi][comparison]["table"][name]["tsv.gz"]))
                        lines.append(f"\n```\n{tree[workflow][tool][combi][comparison]['table'][name]}\n```\n")

                if "SESSION" in tree[workflow][tool][combi][comparison].keys():
                    lines.append(f"##### R-SESSION  {{.tabset}} \n\n")
                    lines.append(f"\n<br />\n")
                    lines.append(f"To access the {tool} R-session, follow these steps\n\n")
                    lines.append(f"1) Create the corresponding conda environment and activate it:  \n")
                    lines.append(f"```\nconda env create -f {env}\n```  \n")
                    lines.append(f"2) Start R and load the workingdir:  \n")
                    lines.append(f"```\nload('{tree[workflow][tool][combi][comparison]['SESSION']}')\n```\n")


    log.debug(logid+"lines: "+str(lines))

    if os.path.exists(output):
        os.rename(output, output+'.bak')
    with open(output, 'a') as writefile:
        for line in lines:
            writefile.write(line)
        writefile.write('\n\n')


####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        try:
            log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            #log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        log.debug(logid+str(log.handlers))
        print(os.getcwd())

        create_Rmd(args.files, args.output, args.env)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


# RmdCreator.py ends here
