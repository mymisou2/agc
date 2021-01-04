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

__author__ = "Myriam Mendli"
__copyright__ = "EISTI"
__credits__ = ["MYRIAM MENDLI"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Myriam MENDLI"
__email__ = "mendlimyri@eisti.eu"
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
    input:
        - sequence
        - chunk_size: longueur du segment
    output:
        chunks_seq_list: liste de sous-séquences de taille chunk_size non chevauchantes.
    A minima 4 segments doivent être obtenus par séquence.
    """
    taille_seq = len(sequence)
    chunks_seq_list = [sequence[i:i+chunk_size] for i in range(0, taille_seq, chunk_size)
        if len(sequence[i:i+chunk_size]) == chunk_size]

    if len(chunks_seq_list) >= 4:
        return chunks_seq_list
    return None

def get_unique(ids):
    """
    input:
        ids: liste d'elements
    output:
        elements uniques d'une liste
    """
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    """
    input:
        - lst1: list
        - lst2: list
    ouput:
        element communs entre 2 listes
    """
    return list(set(lst1) & set(lst2))

def cut_kmer(sequence, kmer_size):
    """
    input:
        - sequence
        - kmer_size: longeueur de k-mer
    output:
        - un générateur de tous les mots de longueur k présents
        dans cette séquence
    """
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i : i + kmer_size]

def get_identity(alignment_list):
    """
    input:
        - alignment_list: prend un alignement (sous forme de liste)
    output:
        - calcule le pourcentage d’identité entre les deux séquences
    """
    len_seq = len(alignment_list[0])
    return sum([1 for i in range(len_seq)
        if alignment_list[0][i] ==  alignment_list[1][i]])/len_seq * 100

def search_mates(kmer_dict, sequence, kmer_size):
    """
    input:
        - dictionnaire de k-mer
        - sequence
        - kmer_size: taille du kmer
    output:
        8 séquences les plus communes du dictionnaire de k-mer
    """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
        if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]

def detect_chimera(perc_identity_matrix):
    """
    input:
        - perc_identity_matrix: matrice identité des séquences
    output:
        - True/False: est une chimère

    Si l’écart type moyen des pourcentages d’identité est supérieur à 5 et
    que 2 segments minimum de notre séquence montrent une similarité différente
    à un des deux parents, nous identifierons cette séquence comme chimérique
    """
    ecart_type_pourcentage = 0
    similarite_seq1, similarite_seq2 = set(), set()

    for pourcentage_identite in perc_identity_matrix:
        ecart_type_pourcentage += statistics.stdev(pourcentage_identite)
        similarite_seq1.add(pourcentage_identite[0])
        similarite_seq2.add(pourcentage_identite[1])

    if len(similarite_seq1) >= 2 or len(similarite_seq2) >= 2:
        ecart_type_moyen_p = ecart_type_pourcentage/len(perc_identity_matrix)
        if ecart_type_moyen_p > 5:
            return True
    return False

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """
    input:
        - sequence
        - id_seq: id de la séquence
        - kmer_dict: dictionnaire de kmer
        - kmer_size: taille du kmer
    output:
        - kmer_dict: retoure l'update du kmer
    Crée un dictionnaire de kmer avec les id des séquences coorespondantes
    """
    for i in cut_kmer(sequence, kmer_size):
        if i not in kmer_dict:
            kmer_dict[i] = list()
        kmer_dict[i].append(id_seq)
    return kmer_dict

def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    input:
        - amplicon_file:fichier fasta
        - minseqlen: la longueur minimale des séquences
        - mincount: comptage minimum.
        - chunk_size: longueur du segment
        - kmer_size: longueur du k-mer
    output:
        generateur de sequences non chimériques [sequence, count]:
        - sequence: une sequence non chimérique
        - count: son nombre d'occurence
    """
    generator_seq_unique = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    kmer_dict = {}
    list_non_chimere = list()
    perc_identity_matrix = []
    id_seq = 0
    for seq, compteur in generator_seq_unique:
        chunks_courant = get_chunks(seq, chunk_size)
        chunks_courant = chunks_courant[:4]
        chunk_mates = []
        for sub_seq in chunks_courant:
            mates = search_mates(kmer_dict, sub_seq, kmer_size)
            chunk_mates.append(mates)

        parents = []

        for j in range(len(chunk_mates)):
            parents = common(parents, chunk_mates[j])

        if len(parents) > 1:
            for parent in parents[0:2]:
                chunk_ref = get_chunks(list_non_chimere[parent], chunk_size)
                perc_identity_matrix = [[]]
                for element, chunk in enumerate(chunks_courant):
                    align = nw.global_align(chunk, chunk_ref[element])
                    identite =  get_identity(align)
                    perc_identity_matrix[element].append(identite)

        if not detect_chimera(perc_identity_matrix):
            kmer_dict = get_unique_kmer(kmer_dict, seq, id_seq, kmer_size)
            list_non_chimere.append(seq)
            id_seq += 1
            yield [seq, compteur]

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    input:
        - amplicon_file:fichier fasta
        - minseqlen: la longueur minimale des séquences
        - mincount: comptage minimum.
        - chunk_size: longueur du segment
        - kmer_size: longueur du k-mer
    ouput:
        - rm_chimer: une liste de séquences sans chimèreset son occurence
    """
    return [rm_chimer
        for rm_chimer in chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)]

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(otu_list, output_file):
    """
    input:
        - otu_list: liste d’OTU
        - output_file: chemin vers un fichier de sortie
    ouput:
        création d'un fichier avec les OTU et leurs occurences
    """
    with open(output_file, "w") as file:
        for i, _ in enumerate(otu_list):
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

if __name__ == '__main__':
    main()
