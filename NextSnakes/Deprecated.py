########################################
############## DEPRECATED ##############
########################################

@check_run
def sources(config):
    logid = scriptname+'.Collection_sources: '
    ret = list()
    search =  [x[0] for x in keysets_from_dict(config["SOURCE"]) if x[0] != 'last']
    if len(getFromDict(config['SAMPLES'], search)) > 0:
        ret.extend(search)
    log.debug(logid+str(ret))
    return ret

@check_run
def genomepath(s, config):
    logid=scriptname+'.Collection_genomepath: '
    sa = os.path.basename(str(s))
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config["SAMPLES"])
    log.debug(logid+'GENOMEPATH: '+str([sa, cond, sk]))
    for skey in sk:
        klist = value_extract(skey, config["SOURCE"])
        for k in klist:
            for x, y in config["GENOME"].items():
                log.debug(logid+'GENOMEPATH: '+str([x, y, k]))
                if str(k) == str(y) or str(k) == str(x):
                    return os.path.join(str(x), str(y))

@check_run
def genome(s, config):
    logid=scriptname+'.Collection_genome: '
    sa = os.path.basename(str(s))
    sp = source_from_sample(str(s), config)
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config['SAMPLES'])
    for skey in sk:
        klist = value_extract(skey, config['SOURCE'])
        for k in klist:
            for x, y in config['GENOME'].items():
                log.debug(logid+str([k, x, y]))
                if str(k) == str(x):
                    return str(y)

@check_run
def fullgenomepath(sa, config):
    ret=list()
    for s in sa:
        l = config["GENOME"][s]
        ret.append(os.path.join(str(s), str(l)))
    return ret

@check_run
def genomename(s, config):
    s = os.path.basename(str(s))
    for k, v in config["SAMPLES"].items():
        for g, l in v.items():
            if s in l:
                for x, y in config["GENOME"].items():
                    if g == y:
                        return str(x)

@check_run
def transcriptomepath(s, config):
    logid=scriptname+'.Collection_transcriptomepath: '
    sa = os.path.basename(str(s))
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config["SAMPLES"])
    log.debug(logid+'TRANSCRIPTOMEPATH: '+str([sa, cond, sk]))
    for skey in sk:
        klist = value_extract(skey, config["SOURCE"])
        for k in klist:
            for x, y in config["TRANSCRIPTOME"].items():
                log.debug(logid+'TRANSCRIPTOMEPATH: '+str([x, y, k]))
                if str(k) == str(y) or str(k) == str(x):
                    return os.path.join(str(x), str(y))

@check_run
def transcriptome(s, config):
    logid=scriptname+'.Collection_transcriptome: '
    sa = os.path.basename(str(s))
    sp = source_from_sample(str(s), config)
    cond= s.split(os.sep)[-2]
    sk = find_key_for_value(sa, config['SAMPLES'])
    for skey in sk:
        klist = value_extract(skey, config['SOURCE'])
        for k in klist:
            for x, y in config['TRANSCRIPTOME'].items():
                log.debug(logid+str([k, x, y]))
                if str(k) == str(x):
                    return str(y)

@check_run
def fulltranscriptomepath(sa, config):
    ret=list()
    for s in sa:
        l = config["TRANSCRIPTOME"][s]
        ret.append(os.path.join(str(s), str(l)))
    return ret

@check_run
def transcriptomename(s, config):
    s = os.path.basename(str(s))
    for k, v in config["SAMPLES"].items():
        for g, l in v.items():
            if s in l:
                for x, y in config["TRANSCRIPTOME"].items():
                    if g == y:
                        return str(x)


@check_run
def namefromfile(s, config):
    if 'NAME' not in config:
        return ''
    else:
        sa = os.path.basename(str(s))
        cond= s.split(os.sep)[-2]
        sk = find_key_for_value(sa, config["SAMPLES"])
        for skey in sk:
            klist = value_extract(skey, config["NAME"])
            for k in klist:
                if str(skey) == str(cond):
                    return str(k)


@check_run
def namefrompath(p, config):
    p = os.path.dirname(p).split(os.sep)
    klist = getFromDict(config["NAME"], p) if 'NAME' in config else list('')
    for k in klist:
        return str(k)

@check_run
def pathstogenomes(samples, config):
    ret = list()
    for s in samples:
        s = os.path.basename(s)
        for k, v in config["SAMPLES"].items():
            for g, l in v.items():
                if s in l:
                    for x, y in config["GENOME"].items():
                        if g == y:
                            ret.append(os.path.join(str(x), str(y)))
    return sorted(list(set(ret)))


@check_run
def source_from_sample(sample, config):
    logid=scriptname+'.Collection_source_from_sample: '
    s = os.path.dirname(str(sample))
    cond= s.split(os.sep)
    log.debug(logid+str([s, cond]))
    ret = getFromDict(config["SOURCE"], cond)[0]
    return ret


@check_run
def anno_from_file(sample, config, step):
    logid = scriptname+'.Collection_anno_from_file: '
    p = os.path.dirname(genomepath(sample, config))
    s = source_from_sample(sample, config)
    ret = os.path.join(config["REFERENCE"], p, subDict(config["ANNOTATE"], s)[step])
    log.debug(logid+str(ret))
    return ret

@check_run
def anno_from_source(source, config, step):
    logid = scriptname+'.Collection_anno_from_source: '
    s = source.split(os.sep)[0:-1]
    p = s[0]
    samp = source.split(os.sep)[-1]
    log.debug(logid+str(s))
    runstate = runstate_from_sample([samp], config)[0]
    lst = list()
    lst.extend(s)
    lst.append(runstate)
    log.debug(logid+str(lst))
    ret = os.path.join(config["REFERENCE"], p, subDict(config["ANNOTATE"], lst)[step])
    log.debug(logid+str(ret))
    return ret
