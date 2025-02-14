#!/usr/bin/env python3

import sys
from Bio import AlignIO
from collections import defaultdict

def justSpeciesNames(fasta):
    outputfile = fasta.split(".")[0] + "_speciesNames.fa"
    aln = AlignIO.read(fasta, "fasta")
    with open(outputfile, "w") as newFile:
        for sequence in aln:
            full_name = sequence.description.split("[", 1)[-1].split("]")[0]
            newFile.write(f">{full_name}\n{sequence.seq}\n")

def modifyAlignment(*fastas):
    species_set = set()
    sequences = defaultdict(dict)
    lengths = {}
    
    for fasta in fastas:
        spec_aln = fasta.split(".")[0] + "_speciesNames.fa"
        with open(spec_aln, "r") as f:
            lines = f.read().split("\n")
            for i in range(0, len(lines) - 1, 2):
                species, sequence = lines[i][1:], lines[i + 1]
                species_set.add(species)
                sequences[fasta][species] = sequence
                lengths[fasta] = len(sequence)
    
    for fasta in fastas:
        output_file = fasta.split(".")[0] + "_allSpecies.fa"
        with open(output_file, "w") as newFile:
            for species in species_set:
                newFile.write(f">{species}\n")
                seq = sequences[fasta].get(species, "?" * lengths[fasta])
                newFile.write(f"{seq}\n")

def concatAlignment(*fastas):
    spec_files = [f.split(".")[0] + "_allSpecies.fa" for f in fastas]
    concatenated = {}
    
    for spec_file in spec_files:
        with open(spec_file, "r") as f:
            lines = f.read().split("\n")
            for i in range(0, len(lines) - 1, 2):
                species, sequence = lines[i][1:], lines[i + 1]
                concatenated.setdefault(species, []).append(sequence)
    
    with open("concatenated_gene_alignments.fasta", "w") as outFile:
        for species, seq_list in concatenated.items():
            outFile.write(f">{species}\n{''.join(seq_list)}\n")

def main():
    fastas = sys.argv[1:]
    for fasta in fastas:
        justSpeciesNames(fasta)
    modifyAlignment(*fastas)
    concatAlignment(*fastas)

if __name__ == "__main__":
    main()
