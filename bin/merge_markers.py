#!/usr/bin/env python3

"""
A script that merges a per-genome file FASTA into file per marker with all genomes

input: all per-genome marker files
Output: One filer per marker with only the genome as header

e.g Output file
L YicC_N.faa
   >L_pneumophila
   MIHSMTAFSRIQ...
   >L_longbeachae
   MIHSLTAFSRVQ...
"""
import argparse
import gzip
from collections import defaultdict
import os

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_markerFilers(input_paths):
    """Read per genome marker FASTAs
    Returns: marker dict marker -> list of (genome, sequence)
    """

    markers = defaultdict(list)
    for path in input_paths:
        #re initiate after each file ends, because will have the info from previous file
        current_marker = None
        current_genome = None
        chunks = []

        with open_maybe_gz(path) as f:
            for line in f:
                print(line)
                line = line.strip()
                if not line:
                    continue

                elif line.startswith(">"):
                    #new header: save the previous seq first
                    if current_marker is not None:
                        markers[current_marker].append((current_genome, "".join(chunks)))

                    current_marker, current_genome = line[1:].split("|")
                    chunks = []
                else:
                    chunks.append(line)
            #end of file save the last sequence
            if current_marker is not None:
                markers[current_marker].append((current_genome, "".join(chunks)))

    return markers

def write_markers_file(markers, outdir):
    """Write one FASTA per marker, with one sequence per genome.
    Input:
        - markers: dict marker -> list of (genome, sequence)
        - outdir: output folder
    """

    os.makedirs(outdir, exist_ok=True)
    for marker_name, entries in markers.items():
        with open(os.path.join(outdir, f"{marker_name}.faa"), "w") as  output_file:
            for genome, sequence in entries:
                output_file.write(f">{genome}\n{sequence}\n")



def parse_args():
    parser = argparse.ArgumentParser(description=__doc__ , formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--inputs", nargs = "+", required=True, help="Per genome marker FASTA files") # + _> one or more inouts required
    parser.add_argument("--outdir", required=True, help="Output foler for per-marker FASTA files")

    return parser.parse_args()


def main():
    args = parse_args()

    input_paths = args.inputs

    markers = parse_markerFilers(input_paths)

    write_markers_file(markers, args.outdir)

if __name__ == "__main__":
    main()


