import requests
from bs4 import BeautifulSoup
import os
import re

BASE_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/therasabdab/search/?all=true"

def clean(text):
    return text.replace("\n", "").replace("\t", "").strip()

def extract_sequences(soup):
    sections = soup.find_all("div", {"class": "row"})
    sequences = []

    for sec in sections:
        headers = sec.find_all("h4")
        if not headers:
            continue

        header_text = headers[0].get_text()
        if "Target" not in header_text:
            continue

        lines = sec.get_text("\n").split("\n")
        lines = [l.strip() for l in lines if l.strip()]

        target = None
        heavy = None
        light = None

        for line in lines:
            if line.startswith("Target"):
                target = line.split("\t")[-1]

            elif line.startswith("Heavy Chain"):
                heavy = line.split("\t")[-1]

            elif line.startswith("Light Chain"):
                light = line.split("\t")[-1]

        if target and heavy and light:
            sequences.append((target, heavy, light))

    return sequences

def write_fasta_file(drug_name, seqs):
    fname = f"{drug_name.lower()}.py"

    fasta_lines = []

    for i, (target, heavy, light) in enumerate(seqs, start=1):
        fasta_lines.append(f"> {drug_name}|target_{i}|heavy")
        fasta_lines.append(heavy)
        fasta_lines.append("")
        fasta_lines.append(f"> {drug_name}|target_{i}|light")
        fasta_lines.append(light)
        fasta_lines.append("")

    fasta_block = "\n".join(fasta_lines).rstrip()

    with open(fname, "w", encoding="utf-8") as f:
        f.write('fasta = """\n')
        f.write(fasta_block)
        f.write('\n"""')

    print(f"[saved] {fname}")

def fetch_therapeutic(drug_name):
    url = BASE_URL.replace("?all=true", "") + f"?search_name={drug_name}&all=true"
    response = requests.get(url)

    if response.status_code != 200:
        print(f"[error] Could not fetch {drug_name}")
        return

    soup = BeautifulSoup(response.text, "html.parser")
    
    print("[debug] First 3000 chars of HTML:")
    print(soup.prettify()[:3000])
    print("\n" + "="*80 + "\n")
    
    seqs = extract_sequences(soup)

    if not seqs:
        print(f"[error] No sequences found for {drug_name}")
        return

    write_fasta_file(drug_name, seqs)

if __name__ == "__main__":
    print("Enter INN drug name (e.g., Emicizumab, Blinatumomab)")
    drug = input("> ").strip()

    fetch_therapeutic(drug)