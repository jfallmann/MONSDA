#! /usr/bin/env python3

import os
import json
import copy
import readline
import glob
from snakemake import load_configfile
from collections import defaultdict
import argparse
from NextSnakes.Logger import *
from functools import reduce
import operator
import datetime

parser = argparse.ArgumentParser(
    description='Helper to create initial config file used for workflow processing'
)
parser.add_argument(
    "-q", "--quickmode", action='store_true', default=False, help='choose quickmode to hide explanatory text'
)
parser.add_argument(
    "-d",
    "--debug",
    action='store_true',
    default=False,
    help='runnign in debugging mode for showing interim results',
)
args = parser.parse_args()

## PRINTERS
def prRed(skk):
    print("\033[91m{}\033[00m".format(skk))


def prGreen(skk):
    print("\033[92m{}\033[00m".format(skk))


def prYellow(skk):
    print("\033[93m{}\033[00m".format(skk))


def prLightPurple(skk):
    print("\033[94m{}\033[00m".format(skk))


def prPurple(skk):
    print("\033[95m{}\033[00m".format(skk))


def prCyan(skk):
    print("\033[96m{}\033[00m".format(skk))


def prLightGray(skk):
    print("\033[97m{}\033[00m".format(skk))


def prBlack(skk):
    print("\033[98m{}\033[00m".format(skk))


def print_rst(rst):
    if not args.quickmode:
        print('\n')
        source = "./docs/source/"
        path = os.path.join(source, rst)
        with open(path, 'r') as f:
            for l in f.readlines():
                prGreen(l.replace('\n', ''))
        print('\n')


def print_dict(dict, indent=6, col=None, gap=''):
    if col == 'Cyan':
        for line in json.dumps(dict, indent=indent).split('\n'):
            prCyan(gap + line)
    else:
        for line in json.dumps(dict, indent=indent).split('\n'):
            print(gap + line)


## DICT MANIPULATORS
def get_by_path(root, items):
    """Access a nested object in root by item sequence."""
    return reduce(operator.getitem, items, root)


def set_by_path(root, items, value):
    """Set a value in a nested object in root by item sequence."""
    if isinstance(get_by_path(root, items[:-1])[items[-1]], list):
        gp = get_by_path(root, items[:-1])[items[-1]]
        gp.append(value)
    else:
        get_by_path(root, items[:-1])[items[-1]] = value


def del_by_path(root, items):
    """Delete a key-value in a nested object in root by item sequence."""
    if isinstance(get_by_path(root, items[:-1])[items[-1]], str):
        get_by_path(root, items[:-1])[items[-1]] = ""
    else:
        del get_by_path(root, items[:-1])[items[-1]]


def setInDict(dataDict, maplist, value):
    first, rest = maplist[0], maplist[1:]
    if rest:
        try:
            if not isinstance(dataDict[first], dict):
                dataDict[first] = {}
        except KeyError:
            dataDict[first] = {}
        setInDict(dataDict[first], rest, value)
    else:
        try:
            if isinstance(dataDict[first], list):
                dataDict[first].append(value)
            else:
                dataDict[first] = value
        except:
            dataDict[first] = value


def get_conditions_from_dict(root, keylist=[]):
    for k, v in root.items():
        keylist.append(
            k,
        )
        if not v:
            yield ':'.join(keylist)
        else:
            yield from get_conditions_from_dict(v, keylist)
        keylist.pop()


def getPathesFromDict(d, value=None):
    def yield_func(d):
        q = [(d, [])]
        while q:
            n, p = q.pop(0)
            yield p
            if isinstance(n, dict):
                for k, v in n.items():
                    q.append((v, p + [k]))
            elif isinstance(n, list):
                for i, v in enumerate(n):
                    q.append((v, p + [i]))

    ret = [p for p in yield_func(d)]
    if value:
        sel = []
        for path in ret:
            if value in path:
                sel.append(path)
        return sel
    return ret


def decouple(d):
    string = json.dumps(d)
    return json.loads(string)


class NestedDefaultDict(defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)


class PROJECT:
    def __init__(self):
        self.name = "moin"
        self.path = ""
        self.cores = 1
        # self.configs = []
        self.baseDict = NestedDefaultDict()
        self.conditionDict = NestedDefaultDict()
        self.samplesDict = NestedDefaultDict()
        self.settingsDict = NestedDefaultDict()
        self.settingsList = []
        self.workflowsDict = NestedDefaultDict()
        self.commentsDict = NestedDefaultDict()
        # self.sampleDict = NestedDefaultDict()


class GUIDE:
    def __init__(self):
        self.mode = ''
        self.question = ""
        self.answer = ""
        self.toclear = 0
        self.debugging = False

    def clear(self, number=None):
        if not number:
            number = self.toclear
        os.system(f'echo -e "\e[{number}A\033[2K"')
        for i in range(number - 1):
            os.system(f'echo -e "\e[-1A\033[2K"')
        os.system(f'echo -e "\e[{number}A\03\c"')
        self.toclear = 0

    def proof_input(self, proof=None, spec=None):
        allowed_characters = [
            'a',
            'b',
            'c',
            'd',
            'e',
            'f',
            'g',
            'h',
            'i',
            'j',
            'k',
            'l',
            'm',
            'n',
            'o',
            'p',
            'q',
            'r',
            's',
            't',
            'u',
            'v',
            'w',
            'x',
            'y',
            'z',
            'A',
            'B',
            'C',
            'D',
            'E',
            'F',
            'G',
            'H',
            'I',
            'J',
            'K',
            'L',
            'M',
            'N',
            'O',
            'P',
            'Q',
            'R',
            'S',
            'T',
            'U',
            'V',
            'W',
            'X',
            'Y',
            'Z',
            '1',
            '2',
            '3',
            '4',
            '5',
            '6',
            '7',
            '8',
            '9',
            '0',
            '(',
            ')',
            '_',
            '-',
            ',',
            '.',
            ':',
            '/',
            ' ',
        ]
        while True:
            if spec:
                a = rlinput(">>> ", spec)
            else:
                a = input(">>> ").strip().replace(" ", "")
            if a == "outbreak":
                outbreak()
                return
            if any(x not in allowed_characters for x in a):
                safe = self.toclear
                self.clear(2)
                self.toclear = safe
                print("You used unallowed letters, try again")
                continue
            if proof is not None and proof != "only_numbers" and any(x not in proof for x in a.split(",")):
                safe = self.toclear
                self.clear(2)
                self.toclear = safe
                print(f"available are: {proof}")
                continue
            if proof == "only_numbers":
                try:
                    [float(x) for x in a.split(',')]
                    self.answer = a
                    break
                except:
                    self.clear(2)
                    print("please enter integer or float")
                    continue
            else:
                self.answer = a
                break

    def display(self, options=None, question=None, proof=None, spec=None):
        if options:
            for k, v in options.items():
                prYellow(f"   {k}  >  {v}")
            # print()
        if question:
            print('\n' + question)
        self.proof_input(proof, spec)


def complete(text, state):
    return (glob.glob(text + '*') + [None])[state]


def rlinput(prompt, prefill=''):
    readline.set_startup_hook(lambda: readline.insert_text(prefill))
    try:
        return input(prompt)
    finally:
        readline.set_startup_hook()


def print_dict_pointer(dict, path, copy, indent=6):
    text = json.dumps(dict, indent=indent)
    route = ['step'] + path.copy()
    out = f"\n{'='*60}\n\n"
    stepper = 1
    for line in text.split('\n'):
        level = int(((len(line) - len(line.lstrip(' '))) - indent) / indent)
        key = (
            line.replace('"', '').replace('{', '').replace('}', '').replace(':', '').replace(',', '').strip()
        )
        if level + 1 >= len(route):
            out += line + '\n'
        elif not key:
            out += line + '\n'
        elif route[level + 1] == key and route[level] == 'step':
            route[stepper] = 'step'
            stepper += 1
            if len(route) == level + 2:
                if route[level - 1] == 'step':
                    if copy and copy != ['']:
                        out += f"{line}    <-\n"
                        option = f"enter ID's on condition level comma separated \n\nor copy {copy} with 'cp'"
                        guide.toclear += 2
                    else:
                        out += f"{line}    <-\n"
                        option = "enter ID's on conditions comma separated "
                else:
                    out += f"{line}    <-\n"
                    option = "enter ID's on conditions comma separated "

            else:
                out += line + '\n'
        else:
            out += line + '\n'
    out += f"\n{'='*60}\n\n"
    return out, option


def rec_tree_builder(subtree, leafes, path=[], tree=None):
    if tree == None:
        tree = subtree
    if not leafes[0]:
        path.pop()
        return
    for leaf in leafes:
        if not leaf in subtree.keys():
            subtree[leaf] = NestedDefaultDict()
    copy = []
    for k, v in subtree.items():
        path.append(k)
        text, opt = print_dict_pointer(tree, path, copy)
        for line in text.split('\n'):
            if '<-' in line:
                prYellow("  " + line)
            else:
                prRed("  " + line)
        guide.toclear += len(text.split('\n')) + 3
        guide.display(question=opt)
        guide.clear()
        if guide.answer != '':
            guide.answer = guide.answer.rstrip(',')
            leafes = [x for x in guide.answer.split(',')]
        elif guide.answer == '' and isinstance(v, dict) and v != {}:
            leafes = list(v.keys())
        else:
            leafes = ['']
        if leafes == ["cp"]:
            leafes = copy
        rec_tree_builder(subtree[k], leafes, path, tree)
        copy = leafes
        leafes = ['']
    if len(path) > 0:
        path.pop()
    return


def depth(d):
    if isinstance(d, dict):
        return 1 + (max(map(depth, d.values())) if d else 0)
    return 0


def select_id_to_set(cdict, i, indent=6):
    text = json.dumps(cdict, indent=indent)
    project.settingsList = []
    d = depth(cdict)
    path = []
    reminder = ''
    counter = 0
    prRed(f"\n{'='*60}\n")
    for line in text.split('\n'):
        level = int(((len(line) - len(line.lstrip(' '))) - indent) / indent)
        key = (
            line.replace('"', '').replace('{', '').replace('}', '').replace(':', '').replace(',', '').strip()
        )
        if key:
            if len(path) > level:
                path = path[: -(len(path) - level)]
            path.append(key)
        if level == i and ':' in line:
            if reminder != key:
                counter += 1
                reminder = key
                project.settingsList.append([])
        elif level < i and '{}' in line:
            counter += 1
            project.settingsList.append([])
        if '{}' in line:
            if ',' in line:
                out = f"{line}{' '*(14-len(key) + indent*(d-2)-indent*level)} <-{' '*((counter+1)%2)*2}  {counter}"
                if counter % 2:
                    prPurple("  " + out)
                else:
                    prLightPurple("  " + out)
            else:
                out = f"{line}{' '*(15-len(key) + indent*(d-2)-indent*level)} <-{' '*((counter+1)%2)*2}  {counter}"
                if counter % 2:
                    prPurple("  " + out)
                else:
                    prLightPurple("  " + out)
            project.settingsList[-1].append(path)
        else:
            prRed("  " + line)
        guide.toclear += 1
    prRed(f"\n{'='*60}\n")
    guide.toclear += 6


def optionsDictToString(d):
    return ','.join(map(' '.join, d.items()))


def stringToOptionsDict(s):
    optsDict = NestedDefaultDict()
    s = ' '.join(s.split())
    pairs = s.split(',')
    pairs = [p.strip() for p in pairs]
    for pair in pairs:
        key = pair.split(' ')[0]
        try:
            value = pair.split(' ')[1]
            optsDict[key] = value
        except:
            optsDict[key] = ""
    return optsDict


def location(dictionary, setting, indent=6):
    prRed(f"\n{'='*60}\n")
    spots = copy.deepcopy(setting)
    d = depth(dictionary)
    text = json.dumps(dictionary, indent=indent)
    for line in text.split('\n'):
        switch = True
        level = int(((len(line) - len(line.lstrip(' '))) - indent) / indent)
        key = (
            line.replace('"', '').replace('{', '').replace('}', '').replace(':', '').replace(',', '').strip()
        )
        if key:
            for path in spots:
                if not path:
                    continue
                if path[0] == key:
                    path.pop(0)
                    if not path:
                        if ',' in line:
                            prYellow(f"  {line}{' '*(14-len(key) + indent*(d-2)-indent*level)} <-")
                        else:
                            prYellow(f"  {line}{' '*(15-len(key) + indent*(d-2)-indent*level)} <-")
                        switch = False
        if switch:
            prRed("  " + line)
        guide.toclear += 1
    prRed(f"\n{'='*60}\n")
    guide.toclear += 4


def interims(name, d):
    if guide.debugging == False:
        return
    print(f'\n***DEBUG INTERIMS:***\n')
    print(f'{name} =')
    print_dict(d)
    print("*********************")


readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: menu-complete")
readline.set_completer(complete)

project = PROJECT()
guide = GUIDE()


################################################################################################################
####                                            CONVERSATION                                                ####
################################################################################################################


def prepare_project(template):
    # add template dict to project object
    project.baseDict = template
    # add comments to commentsdict and remove them from baseDict
    comment_pathes = getPathesFromDict(project.baseDict, "comment")
    for path in comment_pathes:
        if path[-1] != "comment":
            continue
        opt = get_by_path(project.baseDict, path)
        set_by_path(project.commentsDict, path, opt)
        del_by_path(project.baseDict, path)
    interims("project.baseDict", project.baseDict)
    interims("project.commentsDict", project.commentsDict)
    return intro()


def intro():
    # print_rst("intro.rst"),
    prGreen("NEXTSNAKES GUIDE\n")
    opts = {'1': 'create new project', '2': 'modify existing config-file'}
    guide.display(
        question="choose an option",
        options=opts,
        proof=opts.keys(),
    )
    if guide.answer == '1':
        guide.mode = 'new'
        return new()
    if guide.answer == '2':
        guide.mode = 'modify'
        return modify()


def new():
    # print_rst("howitworks.rst")
    prGreen('\nCREATE NEW PROJECT')
    prCyan('\n   Directory:')
    er = 0
    while True:
        if er == 0:
            ques = "Enter the absolute path where your project-folder should be created"
        if er == 1:
            ques = "couldn't find this directory"
        if er == 2:
            ques = f"In the directory you entered, a folder with the name '{project.name}' already exist."
        guide.display(question=ques, spec=os.getcwd())
        project.name = os.path.basename(guide.answer)
        project.path = guide.answer
        if os.path.isdir(project.path):
            er = 2
            continue
        if os.path.isdir(os.path.dirname(project.path)):
            guide.clear(4)
            prCyan(f'  Directory: {project.path}')
            break
        else:
            er = 1
    # return explain_condition_tree()
    return create_condition_tree()


def explain_condition_tree():
    print_rst("ics.rst")
    guide.display(question="enter to continue")
    return create_condition_tree()


def create_condition_tree():
    prCyan('\n  Condition-Tree:\n')
    while True:
        rec_tree_builder(project.conditionDict, [project.name])
        # print("\nthis is your new condition tree:\n")
        # guide.clear(2)
        # print('\n')
        print_dict(project.conditionDict, gap="   ")
        guide.display(
            question="press enter to continue or type 'no' to create it again",
            proof=['', 'no'],
        )
        if guide.answer == 'no':
            project.conditionDict = NestedDefaultDict()
        else:
            # guide.clear(3)
            guide.clear(len(json.dumps(project.conditionDict, indent=4).split('\n')) + 3)
            print_dict(project.conditionDict, col='Cyan', gap="         ")
            project.settingsDict = decouple(project.conditionDict)
            break
    interims("project.settingsDict", project.settingsDict)
    return add_sample_dirs()


def add_sample_dirs():
    # print_rst("samples.rst")
    prCyan('\n  Add Sample-Files')
    path_to_samples_dict = NestedDefaultDict()
    er = 0
    while True:
        if er == 0:
            ques = "Enter an absolute path where samples are stored"
            sp = os.getcwd()
        if er == 3:
            ques = "Add another path or press enter to continue"
            sp = ''
        if er == 1:
            ques = f"Sorry, couldn't find '{dir}'. Enter an absolute path where samples are stored"
            sp = dir
        if er == 2:
            ques = f"Samples must be specified to continue. Enter an absolute path where samples are stored"
            sp = os.getcwd()
        guide.display(question=ques, options=path_to_samples_dict, spec=sp)
        dir = guide.answer
        if os.path.isdir(dir):
            for dirpath, dirnames, filenames in os.walk(dir):
                for filename in [f for f in filenames if f.endswith(".fastq.gz")]:
                    if filename.endswith('_R1.fastq.gz'):
                        filename = filename.replace('_R1.fastq.gz', '')
                    if filename.endswith('_R2.fastq.gz'):
                        filename = filename.replace('_R2.fastq.gz', '')
                    project.samplesDict[os.path.join(dirpath, filename)] = {}
            path_to_samples_dict[dir] = f"{len(project.samplesDict) } Files found"
            guide.clear(len(path_to_samples_dict) + 3)
            er = 3
            continue
        if guide.answer == '' and len(project.samplesDict) == 0:
            er = 2
            guide.clear(2)
            continue
        if guide.answer == '' and len(project.samplesDict) != 0:
            break
            switch = False
        if not os.path.isdir(dir):
            er = 1
            if path_to_samples_dict:
                guide.clear(len(path_to_samples_dict) + 2)
            guide.clear(2)
            continue
    counter = 1
    interims("project.samplesDict", project.samplesDict)
    return assign_samples()


def assign_samples():
    prCyan('\n  Assign Samples to Conditions')
    conditions = [pattern for pattern in get_conditions_from_dict(project.conditionDict)]
    samples = NestedDefaultDict()
    if project.settingsDict:
        safety_dict = decouple(project.settingsDict)
    else:
        safety_dict = decouple(project.conditionDict)
        project.settingsDict = decouple(project.conditionDict)
    while True:
        number = 1
        for samp in project.samplesDict.keys():
            samples[number] = samp
            number += 1
        for condition in conditions:
            if 'SAMPLES' in get_by_path(project.settingsDict, condition.split(':')).keys():
                continue
            cond_as_list = [x for x in condition.split(':')]
            location(project.settingsDict, [cond_as_list])
            guide.display(
                question=f"enter all sample-numbers according to the displayed condition comma separated",
                options=samples,
                proof=[str(i) for i in samples.keys()],
            )
            select = []
            for num in guide.answer.split(','):
                s = os.path.basename(samples[int(num)]).replace('.fastq.gz', '')
                set_by_path(project.samplesDict, [samples[int(num)]], condition)
                select.append(s)
                samples.pop(int(num))
            set_by_path(
                project.settingsDict,
                cond_as_list,
                {"SAMPLES": select},
            )
            guide.toclear += len(samples)
            guide.toclear += len(select)
            guide.toclear += 6
            guide.clear(guide.toclear)
        # print("\nthis is your sample assignment:\n")
        print_dict(project.settingsDict, gap="   ")
        guide.display(question="\npress enter to continue or type 'no' to assign samples again")
        if guide.answer == 'no':
            project.settingsDict = decouple(safety_dict)
        else:
            guide.clear(len(json.dumps(project.settingsDict, indent=4).split('\n')) + 3)
            print_dict(project.settingsDict, col='Cyan', gap="         ")
            break
    interims("project.settingsDict", project.settingsDict)
    return set_settings()


def set_settings():
    # print_rst("settings.rst")
    prGreen('\nGENERAL SETTINGS')
    safety_dict = decouple(project.settingsDict)
    conditions_list = [pattern.split(':') for pattern in get_conditions_from_dict(project.conditionDict)]
    while True:
        for maplist in conditions_list:
            if 'GROUPS' in get_by_path(project.settingsDict, maplist).keys():
                continue
            opt_dict = NestedDefaultDict()
            for key in ['SEQUENCING', 'REFERENCE', 'INDEX', 'PREFIX']:
                setInDict(project.settingsDict, maplist + [key], "")
                opt_dict[key] = ''
            setInDict(project.settingsDict, maplist + ['ANNOTATION', "GTF"], "")
            setInDict(project.settingsDict, maplist + ['ANNOTATION', "GFF"], "")
            opt_dict["GTF"] = ''
            opt_dict["GFF"] = ''

            for key in ['SEQUENCING', 'REFERENCE', 'INDEX', 'PREFIX', 'GTF', 'GFF']:
                guide.toclear = 0
                prCyan(f"\n   Settings for Condition: {':'.join(maplist)}")
                location(project.conditionDict, [maplist])
                if key == 'SEQUENCING':
                    p = ["single", "paired"]
                else:
                    p = None
                if key in ['REFERENCE', 'GTF', 'GFF']:
                    s = os.getcwd()
                else:
                    s = None
                guide.display(
                    question=f"comment: {project.commentsDict['SETTINGS']['comment'][key]}",
                    options=opt_dict,
                    proof=p,
                    spec=s,
                )
                # print_dict(project.settingsDict)
                # samps_number = len(get_by_path(project.settingsDict,maplist+["SAMPLES"]))
                if key in ['GTF', 'GFF']:
                    setInDict(project.settingsDict, maplist + ['ANNOTATION', key], guide.answer)
                    opt_dict[key] = guide.answer
                if key in ['SEQUENCING', 'REFERENCE', 'INDEX', 'PREFIX']:
                    setInDict(project.settingsDict, maplist + [key], guide.answer)
                    opt_dict[key] = guide.answer
                guide.toclear += 12
                guide.clear(guide.toclear)

            last_sample = ''
            opt_dict = NestedDefaultDict()
            for key in ['GROUPS', 'TYPES', 'BATCHES']:
                setInDict(project.settingsDict, maplist + [key], [])
                opt_dict[key] = ''
            for sample in get_by_path(project.settingsDict, maplist + ["SAMPLES"]):
                if last_sample:
                    guide.toclear = 0
                    location(project.conditionDict, [maplist])
                    print("Make Settings for Samples individual:\n")
                    for s in get_by_path(project.settingsDict, maplist + ["SAMPLES"]):
                        if s == sample:
                            prYellow(f"   >> {s} <<")
                        else:
                            prRed(f"      {s} ")
                    print('\n')
                    guide.display(
                        question=f"Press 'cp' to accept entries shown or press enter to make new settings",
                        options=opt_dict,
                    )
                    if guide.answer == 'cp':
                        for key in ['GROUPS', 'TYPES', "BATCHES"]:
                            to_add = get_by_path(project.settingsDict, maplist + [key])[-1]
                            get_by_path(project.settingsDict, maplist + [key]).append(to_add)
                        guide.toclear += 13 + len(get_by_path(project.settingsDict, maplist + ["SAMPLES"]))
                        guide.clear(guide.toclear)
                        continue
                    else:
                        for key in ['GROUPS', 'TYPES', 'BATCHES']:
                            opt_dict[key] = ''
                        guide.toclear += 13 + len(get_by_path(project.settingsDict, maplist + ["SAMPLES"]))
                        guide.clear(guide.toclear)

                last_sample = sample
                for key in ['GROUPS', 'TYPES', 'BATCHES']:
                    guide.toclear = 0
                    location(project.conditionDict, [maplist])
                    print("Make Settings for Samples individual:\n")
                    for s in get_by_path(project.settingsDict, maplist + ["SAMPLES"]):
                        if s == sample:
                            prYellow(f"   >> {s} <<")
                        else:
                            prRed(f"      {s} ")
                    print('\n')
                    guide.display(
                        question=f"comment: {project.commentsDict['SETTINGS']['comment'][key]}",
                        options=opt_dict,
                    )
                    get_by_path(project.settingsDict, maplist + [key]).append(guide.answer)
                    opt_dict[key] = guide.answer
                    guide.toclear += 13 + len(get_by_path(project.settingsDict, maplist + ["SAMPLES"]))
                    guide.clear(guide.toclear)
        prCyan("\n   SETTINGS Key:")
        print_dict(project.settingsDict, gap="   ")
        guide.display(
            question="\npress enter to continue or type 'no' to make settings again", proof=['', 'no']
        )
        if guide.answer == '':
            guide.clear(len(json.dumps(project.settingsDict, indent=4).split('\n')) + 4)
            print_dict(project.settingsDict, col='Cyan', gap="         ")
            interims("project.settingsDict", project.settingsDict)
            if guide.mode == 'new':
                return choose_workflows()
            if guide.mode == 'modify':
                return fillup_workflows()
        if guide.answer == 'no':
            project.settingsDict = decouple(safety_dict)
            continue


def modify():
    er = 0
    while True:
        if er == 0:
            ques = "Enter the absolut path of the config file to be modified"
        if er == 1:
            ques = "couldn't find the file"
        if er == 2:
            ques = "can't read the file, check if it's the right file or if it's corrupted"
        print("\n")
        guide.display(question=ques, spec=os.getcwd())
        if not os.path.isfile(guide.answer):
            er = 1
            continue
        try:
            modify_config = load_configfile(guide.answer)
            project.path = os.path.dirname(guide.answer)
            project.name = list(modify_config['SETTINGS'].keys())[0]
            project.cores = modify_config["MAXTHREADS"]
            project.settingsDict = modify_config['SETTINGS']
            project.conditionDict = NestedDefaultDict()
        except:
            er = 2

        condition_pathes = getPathesFromDict(modify_config, "SAMPLES")
        for path in condition_pathes:
            if path[-1] == "SAMPLES":
                path = path[1:-1]
                set_by_path(project.conditionDict, path, {})
        prRed("\nFollowing configuration was found :\n")
        prRed("Condition-Tree:")
        print_dict(project.conditionDict)
        prRed("\nactive WORKFLOWS:")
        active_workflows = modify_config['WORKFLOWS'].split(',')
        print(active_workflows)
        prRed("\ninactive WORKFLOWS:")
        inactive_workflows = list(modify_config.keys())
        for e in ["WORKFLOWS", "BINS", "MAXTHREADS", "SETTINGS"]:
            inactive_workflows.remove(e)
        for wf in inactive_workflows:
            project.workflowsDict[wf] = decouple(modify_config[wf])
        for e in modify_config['WORKFLOWS'].split(','):
            inactive_workflows.remove(e)
        print(inactive_workflows)

        print("\nWhat do you want to do...\n")
        opts = NestedDefaultDict()
        opts['1'] = "add workflows"
        opts['2'] = "remove workflows"
        opts['3'] = "add conditions (consider that all workflows has to be customized)"
        opts['4'] = "remove conditions"
        opts['5'] = "choose another file"
        guide.display(question="choose an option", options=opts, proof=opts.keys())
        if guide.answer == '1':
            return choose_workflows(inactive_workflows + active_workflows)
        if guide.answer == '2':
            return remove_workflows(active_workflows, inactive_workflows)
        if guide.answer == '3':
            return add_conditions()
        if guide.answer == '4':
            return remove_conditions()
        if guide.answer == '5':
            project.settingsDict = NestedDefaultDict()
            project.conditionDict = NestedDefaultDict()
            continue


def add_conditions():
    prGreen("\n\nAdd Conditions")
    while True:
        new_conditions = decouple(project.conditionDict)
        rec_tree_builder(new_conditions, [project.name])
        print("\nthis is your new condition tree:\n")
        print_dict(new_conditions)
        print('\n')
        guide.display(
            question="press enter to apply changes or type 'no' to make changes again",
            proof=['', 'no'],
        )
        if guide.answer == 'no':
            continue
        else:
            project.conditionDict = new_conditions
            break

    conditions = [pattern.split(':') for pattern in get_conditions_from_dict(project.conditionDict)]
    for condition in conditions:
        try:
            get_by_path(project.settingsDict, condition)
        except:
            setInDict(project.settingsDict, condition, NestedDefaultDict())
        for wf in project.workflowsDict.keys():
            try:
                get_by_path(project.workflowsDict, [wf] + condition)
            except:
                setInDict(project.workflowsDict, [wf] + condition, NestedDefaultDict())
    # print_dict(project.settingsDict)
    # print_dict(project.workflowsDict)
    interims("project.settingsDict", project.settingsDict)
    return add_sample_dirs()


def remove_workflows(actives, inactives):
    prGreen("\n\nRemove Workflows")
    workflows = project.workflowsDict.keys()
    opts = NestedDefaultDict()
    number = 1
    for wf in workflows:
        opts[str(number)] = wf
        number += 1
    while True:
        print("\nfollowing workflows exist:\n")
        guide.display(
            question="enter numbers of workflows to be removed comma separated",
            options=opts,
            proof=opts.keys(),
        )
        new_workflows = decouple(project.workflowsDict)
        new_actives = copy.deepcopy(actives)
        new_inactives = copy.deepcopy(inactives)
        for key, value in opts.items():
            if key in guide.answer.split(','):
                new_workflows.pop(value)
                if value in new_actives:
                    new_actives.remove(value)
                if value in new_inactives:
                    new_inactives.remove(value)
        prRed("\nactive WORKFLOWS:")
        print(new_actives)
        prRed("\ninactive WORKFLOWS:")
        print(new_inactives)
        guide.display(
            question="press enter to apply changes or type 'no' to make changes again", proof=['', 'no']
        )
        if guide.answer == '':
            project.workflowsDict = new_workflows
            return finalize()
        if guide.answer == 'no':
            continue


def remove_conditions():
    prGreen("\n\nRemove Conditions")
    conditions = get_conditions_from_dict(project.conditionDict)
    opts = NestedDefaultDict()
    number = 1
    for c in conditions:
        opts[str(number)] = c
        number += 1
    while True:
        print("\nfollowing conditions exist:\n")
        guide.display(
            question="enter numbers of conditions to be removed comma separated",
            options=opts,
            proof=opts.keys(),
        )
        new_settings = NestedDefaultDict()
        new_workflows = NestedDefaultDict()
        new_conditions = NestedDefaultDict()
        for key, value in opts.items():
            if key in guide.answer.split(','):
                continue
            path = value.split(':')
            setInDict(new_conditions, path, {})
            toadd = get_by_path(project.settingsDict, path)
            setInDict(new_settings, path, toadd)
            for wf in project.workflowsDict.keys():
                toadd = get_by_path(project.workflowsDict, [wf] + path)
                setInDict(new_workflows, [wf] + path, toadd)
        prRed("\nNew condition-Tree:")
        print_dict(new_conditions)
        guide.display(
            question="press enter to apply changes or type 'no' to make changes again", proof=['', 'no']
        )
        if guide.answer == '':
            project.settingsDict = new_settings
            project.workflowsDict = new_workflows
            project.conditionDict = new_conditions
            return finalize()
        if guide.answer == 'no':
            continue


def choose_workflows(existing_workflows=None):
    print_rst("workflows.rst")
    possible_workflows = list(project.baseDict.keys())
    for e in ["WORKFLOWS", "BINS", "MAXTHREADS", "SETTINGS"]:
        possible_workflows.remove(e)
    if existing_workflows:
        for e in existing_workflows:
            possible_workflows.remove(e)
    posWorkDict = NestedDefaultDict()
    counter = 1
    for w in possible_workflows:
        posWorkDict[counter] = w
        counter += 1
    guide.display(
        question="Enter which WORKFLOWS you would like to add",
        options=posWorkDict,
        proof=list(str(i) for i in posWorkDict.keys()),
    )
    for number in guide.answer.split(','):
        wf = posWorkDict[int(number)]
        project.workflowsDict[wf]
    print(f"\nSelected Workflows: {list(project.workflowsDict.keys())}")
    interims("project.workflowsDict", project.workflowsDict)
    return select_conditioning()


def select_conditioning():
    while True:
        d = depth(project.conditionDict)
        for i in range(d - 1):
            select_id_to_set(project.conditionDict, i)
            guide.display(
                question="In the following steps you will make different settings for each condition.\nTo avoid repetitions, specify which conditions should get the same settings\nYou will set all conditions with the same number at once afterwards\n\nTo loop through the possible selections press enter\n\nFinally enter 'ok' to make the settings",
                proof=["ok", ""],
            )
            guide.toclear += 8
            if guide.answer == 'ok':
                return set_workflows()
            else:
                guide.clear()


def fillup_workflows():
    conditions = [pattern.split(':') for pattern in get_conditions_from_dict(project.conditionDict)]
    fillupDict = NestedDefaultDict()
    for wf in project.workflowsDict.keys():
        fillupDict[wf]['on'] = []
        fillupDict[wf]['off'] = []
        for condition in conditions:
            if not get_by_path(project.workflowsDict, [wf] + condition):
                fillupDict[wf]['off'].append(condition)
            if get_by_path(project.workflowsDict, [wf] + condition):
                fillupDict[wf]['on'].append(condition)

    for wf in project.workflowsDict.keys():
        for off_con in fillupDict[wf]['off']:
            location(project.workflowsDict[wf], [off_con])
            prGreen(f'Fillup settings for {wf}\n')
            opts = NestedDefaultDict()
            number = 1
            for on_con in fillupDict[wf]['on']:
                opts[str(number)] = f"copy from {':'.join(on_con)}"
                number += 1
            opts[str(number)] = "make new settings"
            guide.display(question="enter option", options=opts, proof=opts.keys())
            if opts[guide.answer] == 'make new settings':
                project.settingsList = [[off_con]]
                set_workflows(wf)
            else:
                copy_con = opts[guide.answer].split(' ')[-1].split(':')
                toadd = get_by_path(project.workflowsDict, [wf] + copy_con)
                setInDict(project.workflowsDict, [wf] + off_con, toadd)
    return finalize()


def set_workflows(wf=None):
    for workflow in project.workflowsDict.keys():
        if wf:
            workflow = wf
        if not project.workflowsDict[workflow] or wf:
            prGreen(f"\nMake Settings for {workflow}\n")
            opt_dict = NestedDefaultDict()
            tools_to_use = NestedDefaultDict()
            if 'TOOLS' in project.baseDict[workflow].keys():
                number = 1
                for k in project.baseDict[workflow]['TOOLS'].keys():
                    opt_dict[number] = k
                    number += 1
                guide.display(
                    question='Select from these available Tools comma separated:',
                    options=opt_dict,
                    proof=[str(i) for i in opt_dict.keys()],
                )
                for number in guide.answer.split(','):
                    tools_to_use[opt_dict[int(number)]] = project.baseDict[workflow]["TOOLS"][
                        opt_dict[int(number)]
                    ]

            if 'CUTOFFS' in project.baseDict[workflow].keys() and not wf:
                project.workflowsDict[workflow].update({'CUTOFFS': {}})
                project.workflowsDict[workflow]["CUTOFFS"].update(project.baseDict[workflow]["CUTOFFS"])
                for key, value in project.baseDict[workflow]['CUTOFFS'].items():
                    print('\n')
                    guide.display(question=f'Set {key}', proof="only_numbers", spec=value)
                    set_by_path(project.workflowsDict, [workflow, "CUTOFFS", key], str(guide.answer))

            if 'COMPARABLE' in project.baseDict[workflow].keys():
                print_rst("comparable.rst")
                project.workflowsDict[workflow]['COMPARABLE']
                while True:
                    guide.display(
                        question="\nPress enter to add a differential contrast or type 'no'", proof=['', 'no']
                    )
                    if guide.answer == 'no':
                        break
                    groups = set()
                    conditions = get_conditions_from_dict(project.conditionDict)
                    for condition in conditions:
                        for g in get_by_path(project.settingsDict, condition.split(':') + ["GROUPS"]):
                            groups.add(g)
                    number = 1
                    opt_dict = NestedDefaultDict()
                    for g in sorted(groups):
                        opt_dict[number] = g
                        number += 1
                    contrast = [[], []]
                    for i in [0, 1]:
                        print('\n')
                        guide.display(
                            question=f"select {'base-condition' if i == 0 else 'contrast-condition'}",
                            options=opt_dict,
                            proof=[str(i) for i in opt_dict],
                        )
                        contrast[i].append(opt_dict[int(guide.answer)])
                    comp_name = '-VS-'.join(x[0] for x in contrast)
                    project.workflowsDict[workflow]['COMPARABLE'][comp_name] = contrast

            for setting in project.settingsList:
                for tool, bin in tools_to_use.items():
                    project.workflowsDict[workflow]['TOOLS'][tool] = bin
                    for maplist in setting:
                        setInDict(project.workflowsDict, [workflow] + maplist + [tool, "OPTIONS"], [])
                    for i in range(len(project.baseDict[workflow][tool]['OPTIONS'])):
                        if project.baseDict[workflow][tool]['OPTIONS'][i]:
                            call = optionsDictToString(project.baseDict[workflow][tool]['OPTIONS'][i])
                            guide.toclear = 0
                            location(project.conditionDict, setting)
                            print(f"Tool: {tool}\n")
                            guide.display(
                                question=f"comment: {project.commentsDict[workflow][tool]['comment'][i]}",
                                spec=call,
                            )
                            optsDict = stringToOptionsDict(guide.answer)
                            for maplist in setting:
                                setInDict(
                                    project.workflowsDict, [workflow] + maplist + [tool, "OPTIONS"], optsDict
                                )
                            guide.toclear += 6
                            guide.clear()
                        else:
                            for maplist in setting:
                                setInDict(project.workflowsDict, [workflow] + maplist + [tool, "OPTIONS"], {})
            if wf:
                return
    interims("project.workflowsDict", project.workflowsDict)
    if guide.mode == "new":
        return set_cores()
    if guide.mode == "modify":
        return finalize()


def set_cores():
    print('\n')
    guide.display(
        question="Several workflows use multithreading, set the maximum number of cores that can be used",
        proof='only_numbers',
    )
    project.cores = str(guide.answer)
    return finalize()


def finalize():
    # print(project.name)
    # print(project.path)
    # print(project.cores)
    # print_dict(project.conditionDict)
    # print(project.settingsList)
    # print_dict(project.settingsDict)
    # print_dict(project.workflowsDict)
    # print_dict(project.samplesDict)
    # print_dict(project.commentsDict)
    final_dict = NestedDefaultDict()

    final_dict["WORKFLOWS"] = ','.join(project.workflowsDict.keys())
    final_dict["BINS"] = project.baseDict["BINS"]
    final_dict["MAXTHREADS"] = project.cores
    final_dict["SETTINGS"] = project.settingsDict
    final_dict.update(project.workflowsDict)

    configfile = f"config_{project.name}.json"
    print('\n')
    prRed('Final Configuration:')
    print('\n')
    print_dict(final_dict)
    print('\n')

    if guide.mode == "new":
        space = len(configfile)
        print(
            "Above is your final configuration of NextSnakes. The Guide will create this directory as new project:\n"
        )
        prGreen(f"  {os.path.dirname(project.path)}")
        prGreen(f"  └─{os.path.basename(project.path)}")
        prGreen(f"     ├─NextSnakes{' '*(space-10)}   >  symlink to {os.getcwd()}")
        prGreen(f"     ├─FASTQ{' '*(space-5)}   >  contains symlinks of your samplefiles")
        prGreen(f"     ├─GENOMES{' '*(space-7)}   >  contains symlinks of your reference files")
        prGreen(f"     └─{configfile}   >  your brand new configuration file")

        guide.display(
            question="\npress enter to create your project or type 'abort' before it gets serious",
            proof=['', 'abort'],
        )
    if guide.mode == "modify":
        print("Above is your updated configuration of NextSnakes\n")
        print("The old config-file will be preserved (if will have a timestamp in it's name)\n")
        guide.display(
            question=f"\npress enter to update {configfile} or type 'abort' before it gets serious",
            proof=['', 'abort'],
        )

    if guide.answer == 'abort':
        quit()
    else:
        return create_project(final_dict)


def create_project(final_dict):
    if guide.mode == "new":

        # create Project Folder
        cwd = os.getcwd()
        os.mkdir(project.path)
        fastq = os.path.join(project.path, "FASTQ")
        gen = os.path.join(project.path, "GENOMES")
        os.mkdir(fastq)
        os.mkdir(gen)
        os.symlink(cwd, os.path.join(project.path, 'NextSnakes'))

        # LINK samples into FASTQ and insert samplenames in dict
        for sample, condition in project.samplesDict.items():
            if condition:
                cond_as_list = [x for x in condition.split(':')]
                os.chdir(fastq)
                for dir in cond_as_list:
                    if not os.path.exists(os.path.join(dir)):
                        os.mkdir(os.path.join(dir))
                    os.chdir(os.path.join(dir))
                path = '/'.join(cond_as_list)
                cond_dir = os.path.join(fastq, path)
                os.symlink(os.path.realpath(sample), os.path.join(cond_dir, os.path.basename(sample)))

        # link reference and annotation
        for setting in project.settingsList:
            for condition in setting:
                ref = get_by_path(project.settingsDict, condition + ['REFERENCE'])
                if os.path.isfile(ref):
                    if not os.path.exists(os.path.join(gen, os.path.basename(ref))):
                        os.symlink(os.path.realpath(ref), os.path.join(gen, os.path.basename(ref)))
                    f = os.path.join(gen, os.path.basename(ref))
                    rel = os.path.os.path.relpath(f, start=project.path)
                    setInDict(project.settingsDict, condition + ['REFERENCE'], rel)
                else:
                    prRed(
                        f"WARNING: reference path at {condition} is not correct, could not symlink, please do by hand"
                    )
                    setInDict(project.settingsDict, condition + ['REFERENCE'], "EMPTY")

                gtf = get_by_path(project.settingsDict, condition + ['ANNOTATION', 'GTF'])
                if os.path.isfile(gtf):
                    if not os.path.exists(os.path.join(gen, os.path.basename(gtf))):
                        os.symlink(os.path.realpath(gtf), os.path.join(gen, os.path.basename(gtf)))
                    f = os.path.join(gen, os.path.basename(gtf))
                    rel = os.path.os.path.relpath(f, start=project.path)
                    setInDict(project.settingsDict, condition + ['ANNOTATION', 'GTF'], rel)
                else:
                    prRed(
                        f"WARNING: GTF path at {condition} is not correct, could not symlink, please do by hand"
                    )
                    setInDict(project.settingsDict, condition + ['ANNOTATION', 'GTF'], "EMPTY")

                gff = get_by_path(project.settingsDict, condition + ['ANNOTATION', 'GFF'])
                if os.path.isfile(gff):
                    if not os.path.exists(os.path.join(gen, os.path.basename(gff))):
                        os.symlink(os.path.realpath(gff), os.path.join(gen, os.path.basename(gff)))
                    f = os.path.join(gen, os.path.basename(gff))
                    rel = os.path.os.path.relpath(f, start=project.path)
                    setInDict(project.settingsDict, condition + ['ANNOTATION', 'GFF'], rel)
                else:
                    prRed(
                        f"WARNING: GFF path at {condition} is not correct, could not symlink, please do by hand"
                    )
                    setInDict(project.settingsDict, condition + ['ANNOTATION', 'GFF'], "EMPTY")

    if guide.mode == "modify":
        file = os.path.join(project.path, f"config_{project.name}.json")
        bakfile = os.path.join(
            project.path, f"config_{project.name}_{datetime.datetime.now().strftime('%Y%m%d_%H_%M_%S')}.json"
        )
        os.system(f"cp {file} {bakfile}")
    with open(os.path.join(project.path, f"config_{project.name}.json"), 'w') as jsonout:
        print(json.dumps(final_dict, indent=4), file=jsonout)

    configfile = f"config_{project.name}.json"
    print(f"\nStart RunSnakemake with")
    prGreen(f"\n  python3 NextSnakes/RunSnakemake.py -c {configfile} --directory ${{PWD}}\n\n")


####################
####    TEST    ####
####################
# project.baseDict = template

'''
project.name = "TEST"
project.path = "/homes/brauerei/robin/Projects/NextSnakes/TEST"
project.cores = 2
project.conditionDict = {
      "TEST": {
            "a": {
                  "1": {}
            },
            "b": {}
      }
}

project.settingsList = []
project.settingsDict = {
      "TEST": {
            "a": {
                  "1": {
                        "SAMPLES": [
                              "hcc1395_tumor_rep2_R1",
                              "hcc1395_tumor_rep2_R2",
                              "hcc1395_tumor_rep1_R1",
                              "hcc1395_tumor_rep1_R2"
                        ]
                  },
            },
            "b": {
                  "SAMPLES": [
                        "hcc1395_tumor_rep1_R1",
                        "hcc1395_tumor_rep3_R2"
                  ]
            }
      }
}
project.workflowsDict = {"DEU":{}}
project.samplesDict = {
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep2_R1.fastq.gz": "TEST:a:1",
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep2_R2.fastq.gz": "TEST:a:1",
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep3_R1.fastq.gz": "TEST:a:2",
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep1_R2.fastq.gz": "TEST:a:2",
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep1_R1.fastq.gz": "TEST:b",
      "/homes/brauerei/robin/Projects/moin/FASTQ/practical/tumor/a/hcc1395_tumor_rep3_R2.fastq.gz": "TEST:b"
}
project.commentsDict = {
      "SETTINGS": {
            "comment": {
                  "GROUPS": "set groups",
                  "TYPES": "set types",
                  "BATCHES": "set batches",
                  "SEQUENCING": "paired or single",
                  "REFERENCE": "set path to reference",
                  "INDEX": "set index",
                  "PREFIX": "set prefix",
                  "GTF": "set path to gtf",
                  "GFF": "set path to gff"
            }
      },
      "TRIMMING": {
            "trimgalore": {
                  "comment": [
                        "??? options"
                  ]
            },
            "cutadapt": {
                  "comment": [
                        "??? options"
                  ]
            }
      },
      "DEDUP": {
            "umitools": {
                  "comment": [
                        "??? options",
                        "??? options"
                  ]
            }
      },
      "MAPPING": {
            "star": {
                  "comment": [
                        "??? options",
                        "??? options"
                  ]
            }
      },
      "COUNTING": {
            "countreads": {
                  "comment": [
                        "??? options"
                  ]
            }
      },
      "ANNOTATE": {
            "annotatebed": {
                  "comment": [
                        "??? options"
                  ]
            }
      },
      "UCSC": {
            "ucsc": {
                  "comment": [
                        "??? options"
                  ]
            }
      },
      "PEAKS": {
            "macs": {
                  "comment": [
                        "??? options"
                  ]
            },
            "peaks": {
                  "comment": [
                        "??? options"
                  ]
            },
            "piranha": {
                  "comment": [
                        "??? options"
                  ]
            }
      },
      "DE": {
            "deseq2": {
                  "comment": [
                        "Set featureCounts options",
                        "Set Cutoff ???"
                  ]
            },
            "edger": {
                  "comment": [
                        "Set featureCounts options",
                        "Set Cutoff ???"
                  ]
            }
      },
      "DEU": {
            "dexseq": {
                  "comment": [
                        "Set featureCounts options",
                        "Set Cutoff ???"
                  ]
            },
            "edger": {
                  "comment": [
                        "Set featureCounts options",
                        "Set Cutoff ???"
                  ]
            }
      },
      "DAS": {
            "diego": {
                  "comment": [
                        "set featureCounts options",
                        "set Diego options"
                  ]
            },
            "edger": {
                  "comment": [
                        "Set featureCounts options",
                        "Set Cutoff ???"
                  ]
            }
      },
      "DTU": {
            "drimseq": {
                  "comment": [
                        "set Salmon INDEXing options",
                        "set Salmon MAPPING options"
                  ]
            },
            "dexseq": {
                  "comment": [
                        "set Salmon INDEXing options",
                        "set Salmon MAPPING options"
                  ]
            }
      }
}

'''

####################
####    MAIN    ####
####################

# os.chdir(os.path.realpath(os.path.abspath(__file__)))
os.chdir(os.path.dirname(os.path.abspath(__file__)))
if __name__ == '__main__':

    template = load_configfile("configs/template_base_commented.json")

    if args.debug:
        guide.debugging = True

    if args.quickmode:
        prRed("\nrunning in quickmode\n")
    else:
        prRed("\nrunning in explanation mode\n")

    prepare_project(template)
    # set_settings()
