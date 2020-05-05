#! /usr/bin/env python3

import sys, argparse, os, inspect, gzip, glob, re, logging
import traceback as tb

cmd_subfolder = os.path.join(os.path.dirname(os.path.realpath(os.path.abspath(inspect.getfile( inspect.currentframe() )) )),"../../lib")
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

scriptname=os.path.basename(__file__)

def parseargs():
    parser = argparse.ArgumentParser(description='Wrapper around snakemake to run config based jobs automatically')
    #parser.add_argument("-l", "--list", type=str, required=True, help="List of samples")
    parser.add_argument("-n", "--sample_name", action="store_true", help=" provide -n if sample names instead of group names should be used for header" )
    parser.add_argument("-o", "--order", action="store_true", help="if wanted the order of conditions can be given as comma separated list" )
    parser.add_argument("-c", "--conditions", required=True, type=str, help="Conditions to compare" )
    parser.add_argument("-t", "--types", required=False, type=str, help="Sequencing types to compare" )
    parser.add_argument("-r", "--replicates", required=True, type=str, help="Replicates belonging to conditions" )
    parser.add_argument("--cutoff", dest='cutoff', type=int, default=0 ,help="cutoff for minimum count" )
    parser.add_argument("-p", "--paired", required=False, type=str, default=None, help="Sequencing strategy for sample name processing" )
    parser.add_argument("--table", dest='table', required=True, type=str, default='counts.table' ,help="Name of table to write to" )
    parser.add_argument("--anno", dest='anno', required=True, type=str, default='counts.anno' ,help="Name of anno to write to" )
    parser.add_argument("--loglevel", default='INFO', help="Log verbosity" )

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()

class Sample_list(object):

    group_name=""
    replicate_names=[]
    replicate_paths=[]
    # the class constructor
    def __init__(self, group_name):
        self.group_name = group_name
        self.replicate_names = list()
        self.replicate_paths = list()
        self.replicate_types = list()

def prepare_table(conditions, replicates, types, paired, table, anno, sample_name=None, order=None, cutoff=None):
    try:
        #slist,
        logid = scriptname+'.prepare_table: '
        my_groups={}
        list_size=0

        # CLEANUP
        oldtab = os.path.abspath(table)
        oldanno = os.path.abspath(anno)
        for oldfile in glob.glob(oldtab):
            os.rename(oldfile,oldfile+'.bak')
            log.warning(logid+'Found old table file'+oldfile+', was moved to '+oldfile+'.bak')
        for oldfile in glob.glob(oldanno):
            os.rename(oldfile,oldfile+'.bak')
            log.warning(logid+'Found old anno file'+oldfile+', was moved to '+oldfile+'.bak')

        replist = str(replicates).strip().split(',')
        condlist = str(conditions).strip().split(',')
        typelist = str(types).strip().split(',') if types is not None else None
        pairedlist = str(paired).strip().split(',') if paired is not None else None

        log.debug(logid+'REPS: '+str(replist)+'\tLEN: '+str(len(replist)))
        log.debug(logid+'CONDS: '+str(condlist)+'\tLEN: '+str(len(condlist)))

        if types is not None:
            log.debug(logid+'TYPES: '+str(typelist)+'\tLEN: '+str(len(typelist)))

        if paired is not None:
            log.debug(logid+'PAIRED: '+str(pairedlist)+'\tLEN: '+str(len(pairedlist)))

        for i in range(len(replist)):
            rep = None
            cond = None
            typ = None

            rep = str(replist[i])
            cond = str(condlist[i])
            typ = str(typelist[i]) if types is not None else None

            if not rep or not cond:
                log.warning(logid+'No rep/cond found for sample '+str(replist[i]))

            log.debug(logid+'rep/cond: '+str([rep,cond]))

            list_size+=1

            if cond in my_groups:
                my_groups[cond].replicate_paths.append(rep)
                my_groups[cond].replicate_names.append(str.split(os.sep,rep)[-1])
                if typ is not None:
                    my_groups[cond].replicate_types.append(typ)
            else:
                my_groups[cond]=make_sample_list(cond)
                my_groups[cond].replicate_paths.append(rep)
                my_groups[cond].replicate_names.append(str.split(os.sep,rep)[-1])
                if typ is not None:
                    my_groups[cond].replicate_types.append(typ)

        log.debug(logid+'MyGroups: '+str(my_groups.keys()))

        myMatrix = []
        myMatrix.append([])
        myMatrix[0].append('')
        sample_counter=0

        conds = []
        if order:
            conds = order.split(',')
        else:
            conds = [x for x in my_groups.keys()]

        log.debug(logid+'CONDS: '+str(conds))
        typeanno = list()
        for gruppies in conds:
            condition_index=-1
            rep_nr=0
            for replicate in my_groups[gruppies].replicate_paths:
                log.info(logid+'Processing: '+str(replicate))
                condition_index +=1
                sample_counter+=1
                rep_nr+=1

                if (sample_name):
                    myMatrix[0].append(my_groups[gruppies].replicate_names[condition_index])
                    typeanno.append(my_groups[gruppies].replicate_types[condition_index])
                else:
                    myMatrix[0].append(str(my_groups[gruppies].group_name)+'_'+str(rep_nr))
                    typeanno.append(my_groups[gruppies].replicate_types[condition_index])
                if '.gz' in replicate:
                    myInput = gzip.open(replicate,'r')
                else:
                    myInput = open(replicate,'r')

                lineNumber=0
                for line in myInput:
                    if '#' in line[0:5] or '.bam' in line[-10:]:
                        log.debug(logid+'Rep: '+str(replicate)+'\t'+'Line: '+str(line))
                        continue
                    columns = line.strip().split('\t')
                    if columns[0] != "name" and columns[0] != "Geneid" and columns[1]!="count":
                        lineNumber+=1
                        if sample_counter==1:
                            newListi=[]
                            myMatrix.append(newListi)
                            myMatrix[lineNumber].append(str(columns[0]))
                        myMatrix[lineNumber].append(round(float(columns[-1])))

        line = "\t".join(myMatrix[0])
        annos = list()

        for i in range(1,len(myMatrix[0])):
        #for c in myMatrix[0][1:]:
            c = myMatrix[0][i]
            #a = ''.join([i for i in c if not i.isdigit()])
            a = str.join('_',str(c).split('_')[:-1])
            a += '\t'+str(typeanno[i-1]) if types is not None else None
            annos.append(str(c)+'\t'+str(a))

        with gzip.open(table, 'wb') as t:
            t.write(bytes(str(line)+'\n',encoding='UTF8'))
        with gzip.open(anno, 'wb') as a:
            a.write(bytes('\n'.join(annos)+'\n',encoding='UTF-8'))

        for z in range(1,len(myMatrix)):
            zeilen = myMatrix[z]
            willprint = False
            line = str(zeilen[0])+"\t"
            for x in range(1,len(zeilen)):
                line = line + str(zeilen[x])+"\t"
                if (int(zeilen[x]) >= cutoff):
                    willprint = True
            if willprint:
                with gzip.open(table, 'ab') as t:
                    t.write(bytes(str(line)+'\n',encoding='UTF8'))

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

def make_sample_list(group_name):
    sample_list=Sample_list(group_name)
    return sample_list


####################
####    MAIN    ####
####################

if __name__ == '__main__':

    logid = scriptname+'.main: '
    try:
        args=parseargs()
        makelogdir('LOGS')
        try:
            log = setup_logger(name=scriptname, log_file='LOGS/'+scriptname+'.log', logformat='%(asctime)s %(name)-12s %(levelname)-8s %(message)s', datefmt='%m-%d %H:%M', level=args.loglevel)
            log.addHandler(logging.StreamHandler(sys.stderr))  # streamlog
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        prepare_table(args.conditions, args.replicates, args.types, args.paired, args.table, args.anno, args.sample_name, args.order, args.cutoff)
    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        log.error(logid+''.join(tbe.format()))

#
# build_count_table_simple.py ends here
