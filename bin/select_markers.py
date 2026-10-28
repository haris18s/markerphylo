#!/usr/bin/env python3

"""
Select marker gene hits for one genome and write one FASTA per genome (marker | genome).

Inputs:
    - hmmsearch --tblout table (plain or .gz)
    - protein FASTA from gene calling

For each marker (HMM query):
    - exactly one hit --> keep it
    - more than one hit  --> mode: exclude: drop marker for this genome
                             mode: best hit: keep the highest scoring hit

Output headers use the format >marker|genome, so a later merge step can group sequences by marker across all genomes for alignment and concatenation.
An optional TSV summary records what happened to each marker.
"""

import argparse
import gzip
import sys
from collections import defaultdict


def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")

def parse_tbl_file(path):
    """
    Parse a hmmsearch --tblout file

    Input: a tbl out file
    Output: (dict), as marker --> list of (protein_id,score).
        Columns (white space -seperated): 1 target name, 3 query name, 6 full sequence score
   """

    hits = defaultdict(list)
    with open_maybe_gz(path) as file:
        for line in file:
            if line.startswith("#") or not line.rstrip():
                continue
            fields = line.split()
            protein_id = fields[0]
            marker = fields[2]
            score = float(fields[5])
            hits[marker].append((protein_id, score))
    return hits

def select_hits(hits, mode):
    """
    Choose one protein per marker.

    Returns:
        selected: dict marker -> protein_id
        status: dict marker -> (satus, number of hits)
    """
    selected = {}
    status = {}
    for marker, marker_hits in hits.items():
        n_hits = len(marker_hits)

        if n_hits == 1: # if one hit, keep it
            selected[marker] = marker_hits[0][0]
            status[marker] = ("single", n_hits)
        elif mode == "best_hit": # if several hits, keep the best
            best_protein , _ = max(marker_hits, key=lambda hit:hit[1])
            selected[marker] = best_protein
            status[marker] = ("best_hit", n_hits)
        else: # if several hits, exclude them frm analysis - cant tell ortholog from paralog, so skip this marker
            status[marker] = ("excluded_multicopy", n_hits)
    return selected, status


def read_fasta(path, selected_hits):
    """Read only the sequences whose ID is in the wanted ids
    Returns a dict: protein_id --> sequence
    """
    seqeunces = {}
    current_id = None
    chunks = []
    selected_hits = set(selected_hits.values()) # use a set becaue it is much faster the searching uses a hash to find.
    with open_maybe_gz(path) as file:
        for line in file:
            line = line.rstrip()
            if not line:
                continue
            if line.startswith(">"):
                #current id (previous protein's sequnece) is from whole proteome check againat what proteins are found in HMMER hits
                if current_id in selected_hits:
                    seqeunces[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                #if not a >, just append the sequence
                chunks.append(line)
        # Write the least header
        if current_id in selected_hits:
            seqeunces[current_id] = "".join(chunks)
    return seqeunces

def write_outputs(selected_hits, sequences, genome, outfile):
    """Write all selected marker proteins of one genome into one FASTA file
    Input:
        - selected_hits: dict marker -> protein id
        - sequence: dict protein id -> protein sequence
        - genome: genome name (str)
        - outfile: output FASTA pathname of outdir file

    Returns:
        FASTA with headers in the format >marker | genome
    """
    with open(outfile, "w") as out_file:
        for marker, protein_id in selected_hits.items():
            # pyrodigal uses * for stop codon, so remove it in case aligners complain downstream
            sequence = sequences[protein_id].rstrip("*")
            out_file.write(f">{marker}|{genome}\n{sequence}\n")

def write_summary(status, selected_hits, genome, outfile):
    """
    Write a per genome table summarising what happened to each marker.
    Input:
        -status: dict marker -> (status, number of hits)
        -selected_hits : dict marker -> protein_id (only kept markers)
        - genome: genome name (str)
        - outfile: output TSV file

    Output:
        TSV with columns genomme, marker, statu, n_hits, protein_id
        (protein_id is NA for excluded markers)

    """
    with open(outfile, "w") as summary_file:
        #write headers
        summary_file.write("genome\tmarker\tstatus\tn_hits\tprotein_id\n")
        for marker in sorted(status):
            marker_status, n_hits = status[marker]
            protein_id = selected_hits.get(marker, "NA")
            summary_file.write(f"{genome}\t{marker}\t{marker_status}\t{n_hits}\t{protein_id}\n")


def parse_args():
    parser =  argparse.ArgumentParser(description=__doc__, formatter_class=  argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--proteins", required=True, help = "Protein FASTA (plain or .gz)")
    parser.add_argument("--hits", required=True, help="hmmsearch --tblout file (plain or .gz)")
    parser.add_argument("--genome", required=True, help = "Genome ID used as sequence name")
    parser.add_argument("--outfile", required=True, help =" Output FASTA file with all selected markers of this genome")
    parser.add_argument("--summary", default=None, help = "Output TSV summarising marker selection")
    parser.add_argument(
        "--mode",choices= ["exclude", "best_hit"],
        default="exclude",
        help = "How to handle markers with multiple hits (default: exclude)",
    )
    return parser.parse_args()


def main():
    #args
    args= parse_args()


    hits_tbl = parse_tbl_file(args.hits)
    selected_hits, status = select_hits(hits_tbl, args.mode)

    protein_sequences = read_fasta(args.proteins, selected_hits)
    missing = [p for p in selected_hits.values() if p not in protein_sequences]
    if missing:
        sys.exit(f"ERROR: {len(missing)} selected protein not found in {args.proteins}, e.g. {missing[0]}")

    if args.summary:
        write_summary(status, selected_hits,genome=args.genome,outfile = args.summary)
    write_outputs(selected_hits, protein_sequences , genome = args.genome, outfile= args.outfile)



if __name__== "__main__":
    main()