#! /usr/bin/env python3

import sys
import os
import inspect
import gzip
import re
import shutil
import pickle
import logging
import subprocess
import traceback as tb

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

cmd_subfolder = os.path.join(
    os.path.dirname(
        os.path.realpath(os.path.abspath(inspect.getfile(inspect.currentframe())))
    ),
    "../lib",
)
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from Logger import *

scriptname = os.path.basename(inspect.stack()[-1].filename).replace(".py", "")


def read_quant(path):
    logid = scriptname + ".read_quant: "
    if os.path.isdir(path):
        qf = os.path.join(path, "quant.sf.gz")
        if not os.path.exists(qf):
            qf = os.path.join(path, "quant.sf")
    else:
        qf = path
    log.debug(logid + "reading " + str(qf))
    comp = "gzip" if str(qf).endswith(".gz") else None
    df = pd.read_csv(qf, sep="\t", compression=comp)
    df["Name"] = [re.sub(r"\.[0-9]+$", "", str(x)) for x in df["Name"]]
    s = df.set_index("Name")["NumReads"]
    s = s.groupby(level=0).sum()
    return s


def parse_tx2gene(gtf):
    logid = scriptname + ".parse_tx2gene: "
    tx2gene = {}
    gene2name = {}
    opn = gzip.open if str(gtf).endswith(".gz") else open
    with opn(gtf, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            attr = cols[8]
            tid = re.search(r'transcript_id "?([^";]+)"?', attr)
            gid = re.search(r'gene_id "?([^";]+)"?', attr)
            gname = re.search(r'gene_name "?([^";]+)"?', attr)
            if tid and gid:
                t = re.sub(r"\.[0-9]+$", "", tid.group(1))
                g = re.sub(r"\.[0-9]+$", "", gid.group(1))
                tx2gene[t] = g
                if gname:
                    gene2name[g] = gname.group(1)
    log.debug(logid + "found " + str(len(tx2gene)) + " transcripts")
    return tx2gene, gene2name


def run_spit_cmd(subcmd, args_list):
    logid = scriptname + ".run_spit_cmd: "
    cmd = ["spit", subcmd] + [str(x) for x in args_list]
    log.info(logid + " ".join(cmd))
    proc = subprocess.run(cmd, capture_output=True, text=True)
    log.debug(logid + "STDOUT: " + proc.stdout)
    if proc.returncode != 0:
        log.error(logid + "STDERR: " + proc.stderr)
        raise RuntimeError("spit " + subcmd + " failed: " + proc.stderr)
    else:
        log.debug(logid + "STDERR: " + proc.stderr)
    return


def compute_lfc(ifs, ctrl_samples, case_samples):
    eps = 1e-6
    num = ifs.select_dtypes(include=[np.number])
    a = num[case_samples].mean(axis=1)
    b = num[ctrl_samples].mean(axis=1)
    lfc = np.log2((a + eps) / (b + eps))
    lfc[~np.isfinite(lfc)] = 0
    return lfc


def main(anname, gtf, outdir, cmp, combi, cores, spitopts):
    logid = scriptname + ".main: "

    pre_extra = []
    dtu_extra = []
    if spitopts and spitopts.strip():
        if "||" in spitopts:
            pre_part, dtu_part = spitopts.split("||", 1)
            pre_extra = pre_part.split()
            dtu_extra = dtu_part.split()
        else:
            dtu_extra = spitopts.split()

    tx2gene, gene2name = parse_tx2gene(gtf)

    sampleData_all = pd.read_csv(anname, sep="\t", compression="gzip", dtype=str)
    sampleData_all["sample_id"] = sampleData_all["sample_id"].astype(str)

    counts = {}
    for _, row in sampleData_all.iterrows():
        counts[row["sample_id"]] = read_quant(row["path"])
    cts_all = pd.DataFrame(counts)
    cts_all = cts_all.fillna(0)
    cts_all = cts_all[cts_all.index.isin(tx2gene.keys())]
    cts_all = cts_all[cts_all.sum(axis=1) > 0]
    cts_all = cts_all.round().astype(int)
    cts_all.index.name = "tx_id"

    if combi == "none":
        combi = ""

    os.makedirs(os.path.join(outdir, "Tables"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "Figures"), exist_ok=True)
    olddir = os.getcwd()
    os.chdir(outdir)

    session = {}

    for contrast in cmp.split(","):
        contrast_name = contrast.split(":")[0]
        groups = contrast.split(":")[1].split("-vs-")
        A = groups[0].split("+")
        B = groups[1].split("+")
        conds = A + B
        log.info(logid + "contrast " + contrast_name + " A=" + str(A) + " B=" + str(B))

        sub = sampleData_all[sampleData_all["condition"].isin(conds)]
        case_samples = sub[sub["condition"].isin(A)]["sample_id"].tolist()
        ctrl_samples = sub[sub["condition"].isin(B)]["sample_id"].tolist()
        keep = case_samples + ctrl_samples
        cts = cts_all[keep]
        cts = cts[cts.sum(axis=1) > 0]

        tmpdir = os.path.abspath("spit_tmp_" + contrast_name)
        if os.path.exists(tmpdir):
            shutil.rmtree(tmpdir)
        os.makedirs(tmpdir)

        cts_file = os.path.join(tmpdir, "tx_counts.txt")
        tx2gene_file = os.path.join(tmpdir, "tx2gene.txt")
        labels_file = os.path.join(tmpdir, "labels.txt")

        cts.to_csv(cts_file, sep="\t")
        t2g = pd.DataFrame(
            {"tx_id": cts.index, "gene_id": [tx2gene[t] for t in cts.index]}
        )
        t2g.to_csv(tx2gene_file, sep="\t", index=False)

        labels = pd.DataFrame(
            {
                "id": keep,
                "condition": [1] * len(case_samples) + [0] * len(ctrl_samples),
            }
        )
        batches = sub.set_index("sample_id").loc[keep, "batch"].tolist()
        if len(set(batches)) > 1:
            labels["batch"] = batches
        labels.to_csv(labels_file, sep="\t", index=False)

        shared = [
            "-i", cts_file,
            "-m", tx2gene_file,
            "-l", labels_file,
            "-O", tmpdir,
            "--no_clusters",
        ]
        run_spit_cmd("preprocess", shared + pre_extra)
        run_spit_cmd("dtu", shared + dtu_extra)

        spitout = os.path.join(tmpdir, "SPIT_analysis", "spit_out.txt")
        pvals = os.path.join(tmpdir, "SPIT_analysis", "all_p_values.txt")
        ifsf = os.path.join(tmpdir, "SPIT_analysis", "dominance_selected_ifs.txt")
        if not os.path.exists(ifsf):
            ifsf = os.path.join(tmpdir, "SPIT_analysis", "filtered_ifs.txt")

        spit_df = (
            pd.read_csv(spitout, sep="\t")
            if os.path.exists(spitout)
            else pd.DataFrame(columns=["gene_id", "flag", "wasserstein_distance"])
        )
        pv_df = (
            pd.read_csv(pvals, sep="\t", index_col=0)
            if os.path.exists(pvals)
            else pd.DataFrame(columns=["p_value"])
        )
        ifs_df = pd.read_csv(ifsf, sep="\t", index_col=0)

        lfc = compute_lfc(ifs_df, ctrl_samples, case_samples)

        tx_tab = pd.DataFrame(
            {
                "feature_id": ifs_df.index,
                "pvalue": [
                    pv_df["p_value"].get(t, np.nan) if "p_value" in pv_df else np.nan
                    for t in ifs_df.index
                ],
                "lfc": lfc.values,
                "gene_id": ifs_df["gene_id"].values,
            }
        )
        tx_tab["Gene"] = [gene2name.get(g, g) for g in tx_tab["gene_id"]]
        tx_tab = tx_tab.sort_values("pvalue", na_position="last")

        gene_p = tx_tab.groupby("gene_id")["pvalue"].min()
        maxidx = tx_tab.groupby("gene_id")["lfc"].apply(lambda s: s.abs().idxmax())
        gene_lfc = tx_tab.loc[maxidx.values].set_index("gene_id")["lfc"]
        gene_tab = pd.DataFrame({"gene_id": gene_p.index})
        gene_tab["pvalue"] = gene_p.values
        gene_tab["lfc"] = gene_lfc.reindex(gene_tab["gene_id"]).values
        wd = spit_df.set_index("gene_id")["wasserstein_distance"] if len(spit_df) else pd.Series(dtype=float)
        fl = spit_df.set_index("gene_id")["flag"] if len(spit_df) else pd.Series(dtype=str)
        gene_tab["wasserstein_distance"] = gene_tab["gene_id"].map(wd).fillna("NA").values
        gene_tab["flag"] = gene_tab["gene_id"].map(fl).fillna("").values
        gene_tab["Gene"] = [gene2name.get(g, g) for g in gene_tab["gene_id"]]
        gene_tab = gene_tab.sort_values("pvalue", na_position="last")

        prefix = "_".join(["DTU", "SPIT", combi, contrast_name, "table"])
        gene_tab.to_csv(
            os.path.join("Tables", prefix + "_genes.tsv.gz"),
            sep="\t", index=False, compression="gzip",
        )
        tx_tab.to_csv(
            os.path.join("Tables", prefix + "_transcripts.tsv.gz"),
            sep="\t", index=False, compression="gzip",
        )
        ifs_out = ifs_df.copy()
        ifs_out.insert(0, "feature_id", ifs_out.index)
        ifs_out.to_csv(
            os.path.join("Tables", prefix + "_ifs.tsv.gz"),
            sep="\t", index=False, compression="gzip",
        )

        figpre = "_".join(["DTU", "SPIT", combi, contrast_name, "figure"])
        pv_clean = tx_tab["pvalue"].dropna()
        plt.figure(figsize=(8, 6))
        if len(pv_clean) > 0:
            plt.hist(pv_clean, bins=20, color="steelblue", edgecolor="black")
        plt.xlabel("p-value")
        plt.ylabel("Frequency")
        plt.title("SPIT transcript p-values " + contrast_name)
        plt.tight_layout()
        plt.savefig(os.path.join("Figures", figpre + "_PValues.png"), dpi=150)
        plt.close()

        plt.figure(figsize=(8, 6))
        if len(spit_df) > 0:
            top = spit_df.head(20)
            plt.barh(
                [str(x) for x in top["gene_id"]][::-1],
                pd.to_numeric(top["wasserstein_distance"], errors="coerce").fillna(0).values[::-1],
                color="darkorange",
            )
        plt.xlabel("Wasserstein distance")
        plt.ylabel("Gene")
        plt.title("SPIT top DTU genes " + contrast_name)
        plt.tight_layout()
        plt.savefig(os.path.join("Figures", figpre + "_Wasserstein.png"), dpi=150)
        plt.close()

        session[contrast_name] = {
            "genes": gene_tab,
            "transcripts": tx_tab,
            "spit_out": spit_df,
        }
        shutil.rmtree(tmpdir)
        log.info(logid + "done " + contrast_name)

    with gzip.open("_".join(["DTU_SPIT", combi, "SESSION.gz"]), "wb") as fh:
        pickle.dump(session, fh)

    os.chdir(olddir)


if __name__ == "__main__":
    logid = scriptname + ".main: "
    try:
        args = sys.argv[1:]
        anname = args[0]
        gtf = args[1]
        outdir = args[2]
        cmp = args[3]
        combi = args[4]
        cores = int(args[5])
        spitopts = args[6] if len(args) > 6 else ""

        makelogdir("LOGS")
        try:
            log = setup_logger(
                name=scriptname,
                log_file="stderr",
                logformat="%(asctime)s %(name)-12s %(levelname)-8s %(message)s",
                datefmt="%m-%d %H:%M",
                level="INFO",
            )
        except:
            log = logging.getLogger(os.path.basename(inspect.stack()[-1].filename))

        main(anname, gtf, outdir, cmp, combi, cores, spitopts)
    except Exception:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(exc_type, exc_value, exc_tb)
        try:
            log.error(logid + "".join(tbe.format()))
        except Exception:
            print("".join(tbe.format()), file=sys.stderr)
        sys.exit(1)

#
# SPIT.py ends here
