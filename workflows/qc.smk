for conditions in samplecond(SAMPLES,config)):
    tool = subDict(config['QC'],conditions)['QCENV']
    src, treat, setup = conditions
    tempconf = defaultdict()
    for key in ['REFERENCE', 'BINS','MAXTHREADS']:
        tempconf[key] = config[key]
    for key in ['GENOME', 'NAME', 'SOURCE', 'SAMPLES', 'QC']:
        tempconf[key][src][treat][setup] = config[key][src][treat][setup]

    with open('subworkflow.json', 'w') as outfile:
        json.dump(tempconf, outfile)

    subworkflow sampleqc:
        snakefile: "./"+str(tool)+".smk"
        configfile: 'subworkflow.json'
