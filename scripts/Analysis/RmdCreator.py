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
        works.add(file.split("_")[0])
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
        print(setting.split('_'))
        WF, TOOL, COMBI, COMP, RES, NAME = setting.split('_')
        EXT = os.path.basename(file).split(".", 1)[1]
        tree[WF][TOOL][COMBI][COMP][RES][NAME][EXT] = file
    log.info(logid+str(tree))
    return tree

def integrate_table(file):
    chunk = f"""
```{{r, table:{file}, include = TRUE, echo = FALSE, message = TRUE, warning = FALSE}}
library(DT)
table <- read.table(gzfile(paste(params$root,'{file}',sep='/')), header = TRUE)
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
        listlines.append(f"fig.{counter} <- paste(params$root,'{f}', sep='/')")
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

def create_Rmd(files, output):
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
                lines.append(f"##### FIGURES  {{.tabset}} \n\n")
                lines.append(f"###### OVERVIEW \n\n")

                if "figure" in tree[workflow][tool][combi][comparison].keys():
                    for name in tree[workflow][tool][combi][comparison]["figure"].keys():
                        lines.append(f"####### {name} \n")
                        lines.append(integrate_figures([tree[workflow][tool][combi][comparison]["figure"][name]["png"]]))

                if "list" in tree[workflow][tool][combi][comparison].keys():
                    for name in tree[workflow][tool][combi][comparison]["list"].keys():
                        lines.append(f"###### {name} \n")
                        print(str(tree[workflow][tool][combi][comparison]["list"][name][".tsv"]))
                        lines.append(integrate_list(tree[workflow][tool][combi][comparison]["list"][name]["tsv"]))

                lines.append(f"##### TABLES  {{.tabset}} \n\n")

                if "table" in tree[workflow][tool][combi][comparison].keys():
                     for name in tree[workflow][tool][combi][comparison]["table"].keys():
                         lines.append(f"###### {name} \n")
                         # lines.append(integrate_table(tree[workflow][tool][combi][comparison]["table"][name]["tsv.gz"]))
                         lines.append(tree[workflow][tool][combi][comparison]["table"][name]["tsv.gz"])
                         lines.append("\n\n")

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
        makelogdir('LOGS')
        try:
            log = setup_logger(name=scriptname, log_file='stderr', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            #log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        print(os.getcwd())
        create_Rmd(args.files, args.output)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))


# RmdCreator.py ends here
