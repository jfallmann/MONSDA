for conditions in samplecond(SAMPLES,config)):
    tool = subDict(config['QC'],conditions)['QCENV']
    src, treat, setup = conditions
    tempconf = defaultdict()
    for key in ['REFERENCE', 'BINS','MAXTHREADS']:
        tempconf[key] = config[key]
    for key in ['GENOME', 'NAME', 'SOURCE', 'SAMPLES', 'TRIMMING']:
        tempconf[key][src][treat][setup] = config[key][src][treat][setup]

    log.debug(tool, src, treat, setup, tempconf)
    with open('qcsubworkflow.json', 'w') as outfile:
        json.dump(tempconf, outfile)

    makeoutdir("subworkflow")
    shutil.copy(str(tool)+".smk", "subworkflow/Snakefile")


subworkflow sampleqc:
    snakefile: tool+".smk"
    configfile: "qcsubworkflow.json"

include: 'multiqc.smk'
