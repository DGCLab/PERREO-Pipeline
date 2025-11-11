###############################################################################
####                                                                       ####
####                      PRIMER CONCENSOUS FINDER                         ####
####                                                                       ####
####             Authors:                                                  ####
####                                                                       ####
####                                                                       ####
###############################################################################

import os
import csv

from math import log2

def info_from_probs(p):
    H = 0
    for b in ('A','C','G','T'):
        if p[b] > 0:
            H -= p[b]*log2(p[b])
    return 2-H

def parse_hmm_probs(hmm_path):
    probs = []
    in_hmm = False
    with open(hmm_path) as fh:
        for line in fh:
            if line.startswith('HMM'):
                in_hmm = True
                next(fh)
                continue
            if not in_hmm:
                continue
            if line.strip() == '//' or not line.strip():
                break
            parts = line.split()
            if not parts[0].isdigit():
                continue
            pa, pc, pg, pt = map(float, parts[1:5])
            s = pa+pc+pg+pt
            if s == 0: p = {'A':.25,'C':.25,'G':.25,'T':.25}
            else: p = {b:v/s for b,v in zip('ACGT',(pa,pc,pg,pt))}
            probs.append(p)
    return probs

def best_window_within_consensus(hmm_path, fasta_path, w=80):

    probs = parse_hmm_probs(hmm_path)

    with open(fasta_path) as f:
        f.readline()
        consensus = f.read().replace('\n','').upper()
    L = min(len(probs), len(consensus))
    info = [info_from_probs(probs[i]) for i in range(L)]

    best = (-1,-1,-1)
    for i in range(L-w):
        mean_ic = sum(info[i:i+w])/w
        if mean_ic > best[2]:
            best = (i, i+w, mean_ic)
    start, end, mean_ic = best
    seq = consensus[start:end]
    return {'start': start+1, 'end': end, 'mean_IC_bits': round(mean_ic,3), 'sequence_nt': seq}

## Download Top10 DEGs from Dfam:
import requests
from pathlib import Path
import time

DFAM_BASE = "https://dfam.org/api"


def get_accession_from_name(family_name):
    """
    Busca una familia por nombre en Dfam y devuelve (accession, nombre_dfam).
    Usa el endpoint /families con el parámetro name_accession.
    """
    url = f"{DFAM_BASE}/families"
    params = {
        "name_accession": family_name,  
        "format": "summary",
        "limit": 1
    }
    r = requests.get(url, params=params)
    r.raise_for_status()
    data = r.json()

    results = data.get("results", [])
    if not results:
        print(f"❌ No se encontró familia con nombre/accesion '{family_name}'")
        return None, None

    fam = results[0]
    accession = fam["accession"]
    dfam_name = fam["name"]
    print(f"✅ {family_name} → {dfam_name} ({accession})")
    return accession, dfam_name




def download_dfam_file(family_id, filetype, outdir="."):
    """
    Descarga el archivo HMM o FASTA de una familia Dfam.
    """
    if filetype == "hmm":
        url = f"{DFAM_BASE}/families/{family_id}/hmm"
        params = {"format": "hmm"}
        ext = ".hmm"
    elif filetype == "fasta":
        url = f"{DFAM_BASE}/families/{family_id}/sequence"
        params = {"format": "fasta"}
        ext = ".fa"
    else:
        raise ValueError("filetype debe ser 'hmm' o 'fasta'.")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    outfile = outdir / f"{family_id}{ext}"

    if outfile.exists():
        print(f"🟡 Ya existe: {outfile}")
        return str(outfile)

    print(f"⬇️  Descargando {filetype.upper()} para {family_id} ...")
    r = requests.get(url, params=params)
    r.raise_for_status()
    outfile.write_text(r.text)
    return str(outfile)


def download_dfam_by_name(name_list, outdir="."):
    """
    Descarga los archivos HMM y FASTA a partir de nombres de familias.
    Además genera un CSV con el mapeo nombre_entrada ↔ nombre_dfam ↔ accession.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mapping = []  

    for name in name_list:
        acc, dfam_name = get_accession_from_name(name)
        if not acc:
            continue

        
        download_dfam_file(acc, "hmm", outdir)
        download_dfam_file(acc, "fasta", outdir)

        
        mapping.append((name, dfam_name, acc))

        time.sleep(1)  

    # escribe CSV con el mapeo
    map_path = outdir / "dfam_name_to_accession.csv"
    with open(map_path, "w") as out:
        out.write("input_name,dfam_name,accession\n")
        for input_name, dfam_name, acc in mapping:
            out.write(f"{input_name},{dfam_name},{acc}\n")

    print(f"\n📄 Mapeo guardado en: {map_path}")


if __name__ == "__main__":
    
    with open("names.txt") as f:
        names = [line.strip() for line in f if line.strip()]
    download_dfam_by_name(names, outdir="dfam_models")

## Busca las secuencias concenso


def batch_best_windows(folder="dfam_models", w=80, out_csv="dfam_best_windows.csv"):
    results = []
    hmm_files = [f for f in os.listdir(folder) if f.endswith(".hmm")]

    for hmm_file in hmm_files:
        base = os.path.splitext(hmm_file)[0]
        fasta_file = os.path.join(folder, base + ".fa")
        hmm_path = os.path.join(folder, hmm_file)

        if not os.path.exists(fasta_file):
            print(f"⚠️  No se encontró FASTA para {base}, se omite.")
            continue

        try:
            res = best_window_within_consensus(hmm_path, fasta_file, w=w)
            res["family"] = base
            results.append(res)
            print(f"✅ {base}: {res['start']}-{res['end']}  IC={res['mean_IC_bits']}")
        except Exception as e:
            print(f"❌ Error procesando {base}: {e}")

    # Guardar en CSV
    with open(out_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=["family", "start", "end", "mean_IC_bits", "sequence_nt"])
        writer.writeheader()
        writer.writerows(results)

    print(f"\n📄 Resultados guardados en: {out_csv}")
    return results


if __name__ == "__main__":
    batch_best_windows(folder="dfam_models", w=80)
