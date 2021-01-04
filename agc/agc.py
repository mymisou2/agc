#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""
import argparse
import sys
import os
import statistics
from collections import Counter
from itertools import groupby
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"

def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """
    Genreateur de séquences
    Input:
    fastq_file: fichier fastq
    Return
    un générateur de séquences
    """
    with open(amplicon_file) as monf:
        seq = ''
        for line in monf:
            if (len(seq) >= minseqlen and '>' in line):
                yield seq
                seq = ''
            elif '>' not in line:
                seq = seq + line.replace(" ", "").replace("\n", "")
            else:
                seq = ''
            # Derniere sequence
        if len(seq) >= minseqlen:
            yield seq

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """
    input:
    - amplicon_file: fichier fasta
    - minseqlen: la longueur minimale des séquences
    - mincount: comptage minimum.
    output:
    - générateur des séquences uniques dans ordre décroissant d’occurrence
    - occurrence
    """
    generator_seq = read_fasta(amplicon_file, minseqlen)
    list_seq_occ =  list()
    for seq, gen in groupby(generator_seq):
        count_el = sum(1 for i in gen)
        if count_el >= mincount:
            list_seq_occ.append((seq, count_el))
            list_seq_occ.sort(key = lambda x: x[1], reverse=True)
    for seq in list_seq_occ:
        yield seq

def get_chunks(sequence, chunk_size):
    """
    prend une séquence et une longueur de segment l: chunk_size
    et retourne une liste de sous-séquences de taille l non chevauchantes.
    A minima 4 segments doivent être obtenus par séquence.
    """
    chunks = len(sequence)
    chunks_seq_list = [sequence[i:i+chunk_size] for i in range(0, chunks, chunk_size) if len(sequence[i:i+chunk_size])== chunk_size]
    if len(chunks_seq_list) >= 4:
        return chunks_seq_list

def get_unique(ids):
    """ tt """
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    """ tt """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    prend une séquence et une longueur de k-mer et retourne un générateur de
    tous les mots de longueur k présents dans cette séquence,  yield kmer
    """
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i : i + kmer_size]

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """ tt """
    for i in cut_kmer(sequence, kmer_size):
        if i not in kmer_dict:
            kmer_dict[i] = list()
        kmer_dict[i].append(id_seq)
    return kmer_dict

def get_identity(alignment_list):
    """ tt """
    len_seq = len(alignment_list[0])
    return sum([1 for i in range(len_seq) if alignment_list[0][i] ==  alignment_list[1][i]])/len_seq * 100

def search_mates(kmer_dict, sequence, kmer_size):
    """ tt """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size) if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]

def detect_chimera(perc_identity_matrix):
    """ tt """
    standard_deviation = 0
    similarity_chunk1 = set()
    similarity_chunk2 = set()

    for sim in perc_identity_matrix:
        standard_deviation += statistics.stdev(sim)
        similarity_chunk1.add(sim[0])
        similarity_chunk2.add(sim[1])

    if len(similarity_chunk2) >= 2 or len(similarity_chunk1) >= 2:
        standard_deviation_mean = standard_deviation/len(perc_identity_matrix)
        if standard_deviation_mean > 5:
            return True
    return False

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ tt """
    kmer_dict = {}

    for id_seq, sequence_ in enumerate(dereplication_fulllength(amplicon_file, minseqlen, mincount)):
        chimera = False
        chunk_mates = [search_mates(kmer_dict, sequence, kmer_size) for sequence in get_chunks(sequence_[0], chunk_size)]

        parentes = []
        for j in range(len(chunk_mates)):
            parentes = common(parentes, chunk_mates[j])

        if len(parentes) >= 2:
            for _ in parentes[0:2]:
                list_seq = get_chunks([], chunk_size)
                mat_id = [[]]
                mat_id = [get_identity(nw.global_align(chunk, list_seq[k])) for k, chunk in enumerate(chunks)]
            chimera = detect_chimera(mat_id)
        if not chimera:
            kmer_dict = get_unique_kmer(kmer_dict, sequence_[0], id_seq, kmer_size)
            yield sequence_

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """ tt """
    return [without_chimera for without_chimera in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)]

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(otu_list, output_file):
    """ tt """
    with open(output_file, "w") as file:
        for i in range(0,len(otu_list)):
            file.write(">OTU_" + str(i + 1) + " occurrence:"+ str(otu_list[i][1]) + "\n")
            file.write(fill(str(otu_list[i][0])))
            file.write("\n")

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()

    otu = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)

    write_OTU(otu, args.output_file)

if __name__ == '__main__':
    main()
