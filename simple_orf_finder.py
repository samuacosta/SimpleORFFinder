#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''SimpleORFFinder
Finds open reading frames in given nucleotide sequences.

The script searches for all the possible open reading frames (ORFs) from DNA
nucleotide sequences and writes them as a FASTA formatted file, to the standard
output or to a file if specified. It supports custom genetic code, custom start
and stop codons, multiple files and multiple sequences per file. Minimum ORF
length and limitation to specific subset of reading frames are also supported.
'''
import numpy
import sys
import random
import subprocess
import re
import decimal
import math
import os
import shutil
import time
import types
import uuid
import argparse
import datetime

__author__ = "Samuel Acosta"
__email__ = "samuel.acostamelgarejo@postgrad.manchester.ac.uk"
__license__ = "GPL"
__version__ = "1.0.0"


def main():
    '''Main method'''
    try:
        # Read command line parameters
        args = read_parameters()
        # Parameter warnings
        print_warnings(args)
        # Overwrite previously generated file if exists
        if args.output is not None and "" != args.output and os.path.exists(args.output):
            os.remove(args.output)
        genetic_code = load_genetic_code(args.gencode)
        start_codons = args.start
        stop_codons = args.stop
        # For each input file received as parameter
        for input_file in args.input:
            # For each different organism in each provided input file
            for organism in read_input(input_file):
                name_org = organism[0]
                nucleotides_seq = organism[1]
                # Search for ORFs
                orfs = get_orfs(args, nucleotides_seq, args.frames,
                                genetic_code, start_codons, stop_codons)
                # Print the results
                print_orfs(name_org, orfs, args)
        print("-----")
        print("ORF finder: Finished succesfully")
    except Exception as e:
        print("ORF finder: An error occurred during execution:")
        print(e)
        # raise


def read_parameters():
    '''Method for retrieving command line parameters

    Returns:
        args: The parsed command line arguments
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-i",
                        "--input",
                        nargs='*',
                        help="The name of the input file(s) in FASTA format. [Default=test_input.fasta]",
                        default=["test_input.fasta"]
                        )
    parser.add_argument("-o",
                        "--output",
                        help="The name of the optional output file. [Default=standard output]"
                        )
    parser.add_argument("-min", "--minimum",
                        help="The minimum size of ORFs to search, in nucleotides. [Default=150]",
                        type=int,
                        default=150
                        )
    parser.add_argument("-g", "--gencode",
                        help="The name of the file with the genetic code to use. [Default=standard_genetic_code.csv]",
                        default="standard_genetic_code.csv"
                        )
    parser.add_argument("-start", "--start",
                        nargs='*',
                        help="List of possible start codons to use. [Default=ATG]",
                        default=["ATG"]
                        )
    parser.add_argument("-stop", "--stop",
                        nargs='*',
                        help="List of possible stop codons to use. [Default=TAG,TAA,TGA]",
                        default=["TAG", "TAA", "TGA"]
                        )
    parser.add_argument("-f",
                        "--frames",
                        help="The frames that should be considered (1 to 6). [Default=All]",
                        type=int,
                        nargs='*',
                        default=[1, 2, 3, 4, 5, 6]
                        )
    return(parser.parse_args())


def print_warnings(args):
    print("ORF FINDER PARAMETERS:")
    print("Input file(s): " + ", ".join(args.input))
    if args.output is not None and "" != args.output:
        print("Output file: " + args.output)
    else:
        print("Using standard output")
    print("Minimum amount of nucleotides for ORFs: " + str(args.minimum))
    print("Frames to be considered for search: " +
          ", ".join(list(map(str, args.frames))))
    print("-----")


def load_genetic_code(genetic_code_filename):
    genetic_code = {}
    file = open(genetic_code_filename, 'r')
    for line in file.readlines():
        codon = line.split(",")[0]
        amino_acid = line.split(",")[1].strip()
        genetic_code[codon] = amino_acid
    return genetic_code


def read_input(input_filename):
    '''Method for reading input files

    Args:
        input_filename (str): The name of the file to be processed

    Returns:
        sequences[]: A list of sequences and their names
    '''
    sequences = []
    file = open(input_filename, 'r')
    is_first_line = True
    name_org = nucleotides_seq = ""
    # For each line of the file
    for line in file.readlines():
        if line.startswith(">"):
            if not is_first_line:
                sequences.append([name_org, nucleotides_seq.upper()])
            # Create a short name for the organism to be processed (for FASTA formatting)
            name_org = re.split(" | \|", line)[0].replace(
                ">", "").replace("\n", "").upper()
            nucleotides_seq = ""
            is_first_line = False
        else:
            nucleotides_seq += line.replace("\n", "")  # Removing newlines
    sequences.append([name_org, nucleotides_seq.upper()])
    return sequences


def get_orfs(args, nucleotide_seq, frames_considered, genetic_code, start_codons, stop_codons):
    '''Method for calculating all the ORFs from a particular nucleotide sequence

    Args:
        args: The object containing the command line arguments
        nucleotide_seq(str): The nucleotide sequence to be processed
        frames_considered([int]): The list of frames to be considered for ORF search

    Returns:
        orfs[]: A list of all found ORFs
    '''
    orfs = {}
    reverse_complement = get_reverse_complement(nucleotide_seq)
    # Only append ORFs in frames that were specified
    if 1 in frames_considered:
        orfs[1] = get_orfs_frame(args, nucleotide_seq,
                                 genetic_code, start_codons, stop_codons)
    if 2 in frames_considered:
        orfs[2] = get_orfs_frame(
            args, nucleotide_seq[1:], genetic_code, start_codons, stop_codons)
    if 3 in frames_considered:
        orfs[3] = get_orfs_frame(
            args, nucleotide_seq[2:], genetic_code, start_codons, stop_codons)
    if 4 in frames_considered:
        orfs[4] = get_orfs_frame(
            args, reverse_complement, genetic_code, start_codons, stop_codons)
    if 5 in frames_considered:
        orfs[5] = get_orfs_frame(
            args, reverse_complement[1:], genetic_code, start_codons, stop_codons)
    if 6 in frames_considered:
        orfs[6] = get_orfs_frame(
            args, reverse_complement[2:], genetic_code, start_codons, stop_codons)
    return orfs


def get_orfs_frame(args, nucleotide_frame, genetic_code, start_codons, stop_codons):
    '''Method for calculating the ORFs from a particular nucleotide subsequence

    Args:
        args: The object containing the command line arguments
        nucleotide_frame(str): The nucleotide subsequence to be processed

    Returns:
        orfs_frame[]: A list of the found ORFs in the specified frame
    '''
    orfs_frame = []
    orf_nucleotide = ""
    orf_start_protein = 0
    number_nucleotides_considered = len(
        nucleotide_frame)-len(nucleotide_frame) % 3
    for i in range(0, number_nucleotides_considered, 3):
        codon = nucleotide_frame[i:i+3]
        # If the codon contains an unspecified nucleotide, we ignore that ORF
        if "N" in codon:
            orf_nucleotide = ""
        else:
            if "" != orf_nucleotide:
                if codon not in stop_codons:
                    # We keep appending new codons
                    orf_nucleotide += codon
                else:
                    # We stop reading and append the translated sequence to the results
                    if args.minimum <= len(orf_nucleotide):
                        orfs_frame.append(
                            [get_protein_seq(orf_nucleotide, genetic_code), orf_start_protein])
                    orf_nucleotide = ""
            else:
                if codon in start_codons:
                    # We start reading
                    orf_nucleotide = codon
                    orf_start_protein = (i+3)//3
        if (i+3 == number_nucleotides_considered and "" != orf_nucleotide):
            # The sequence finishes without a stop codon
            if args.minimum <= len(orf_nucleotide):
                orfs_frame.append(
                    [get_protein_seq(orf_nucleotide, genetic_code), orf_start_protein])
    return orfs_frame


def get_protein_seq(nucleotides_seq, genetic_code):
    '''Method for calculating the resulting amino acid sequence from a nucleotide sequence

    Args:
        nucleotides_seq(str): The sequence of nucleotides to be translated

    Returns:
        protein_seq(str): The resulting amino acid sequence
    '''
    protein_seq = ""
    for i in range(0, len(nucleotides_seq)-len(nucleotides_seq) % 3, 3):
        protein_seq += get_amino_acid(nucleotides_seq[i:i+3], genetic_code)
    return protein_seq


def get_amino_acid(codon, genetic_code):
    '''Method for translating each codon to the corresponding amino acid

    Args:
        codon(str): The codon to be translated

    Returns:
        amino_acid: The translated amino acid
    '''
    amino_acid = genetic_code.get(codon)
    if amino_acid is not None and "" != amino_acid:
        return amino_acid
    else:
        # It's an unknown amino acid
        print("Warning: Unknown amino acid found for codon: " + codon)
        return "X"


def get_reverse_complement(nucleotide_seq):
    '''Method for calculating the reverse complement of a nucleotide sequence

    Args:
        nucleotide_seq(str): The nucleotide sequence

    Returns:
        complement[::-1]: The reverse complement of the nucleotide sequence
    '''
    trans_tab = str.maketrans('ATGC', 'TACG')
    # We calculate the complement
    complement = nucleotide_seq.translate(trans_tab)
    # We reverse the complement
    return complement[::-1]


def print_orfs(name_org, orfs, args):
    '''Method for printing the found ORFs to file or standard output

    Args:
        name_org(str): The short name of the organism (for FASTA formatting)
        orfs[]: The list of found ORFs
        args: The object containing the command line arguments
    '''
    # If an output file was specified, we open it
    if args.output is not None and "" != args.output:
        outputfile = open(args.output, 'a')

    for frame, orfs_frame in orfs.items():
        for i in range(0, len(orfs_frame)):

            # We calculate the header for FASTA formatting
            header = ""
            name = "_".join([">" + name_org, "F"+str(frame), '%04d' % (i+1)])
            framenr = "FRAME" + str(frame)
            length = str(len(orfs_frame[i][0]))
            start = str(orfs_frame[i][1])
            header = " ".join([name, framenr, length, start])

            # We append the amino acid sequence in fixed length subsequences
            orf = "\n".join(
                [header] + list(fixed_length_subsequences(orfs_frame[i][0], 70)))

            if args.output is not None and "" != args.output:
                # Print to file
                print(orf, file=outputfile)
            else:
                # Standard output
                print(orf)


def fixed_length_subsequences(string, length):
    '''Method for splitting sequences in fixed lenght substrings (for FASTA formatting)

    Args:
        string(str): The string to be splitted
        length(int): The maximum length of the resulting substrings
    '''
    return (string[0+i:length+i] for i in range(0, len(string), length))


# Entry
if __name__ == "__main__":
    main()
