#!/usr/bin/python

import h5py
import argparse
from Bio import SeqIO


def sequence_to_array(fasta_file, output_path):
    """
    Convert a genomic sequence to a dictionary of 5xn arrays of
    nucleotides.

    Parameters:
    -----------
    fasta_file: string
        The name of the fasta file that we wish to extract the sequence
        from

    Returns
    -------
    sequence_arrays: {'sequence_name': [[0,1,0,0,0]...]...}
        A dictionary indexed by the sequence names containing 5xn arrays
        where the position of the letter is 1 and the other values are 0
    """
    f = h5py.File(output_path, "w")  # root group
    record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    # aggregated_dict = {}

    for key in record_dict.keys():
        sequence = record_dict[key]
        nucleotide_array = f.create_dataset(str(key), (len(str(sequence.seq)), 6), "f", compression='gzip')  # is value in dict
        temporary_array = []
        # order: [A, T, G, C, N, Other]
        for i in range(0, len(str(sequence.seq))):
            letter = str(sequence.seq)[i]
            if (letter == 'A') | (letter == 'a'):
                temporary_array.append([1, 0, 0, 0, 0, 0])
            elif (letter == "T") | (letter == "t"):
                temporary_array.append([0, 1, 0, 0, 0, 0])
            elif (letter == "G") | (letter == "g"):
                temporary_array.append([0, 0, 1, 0, 0, 0])
            elif (letter == "C") | (letter == "c"):
                temporary_array.append([0, 0, 0, 1, 0, 0])
            elif (letter == "N") | (letter == "n"):
                temporary_array.append([0, 0, 0, 0, 1, 0])
            else:
                temporary_array.append([0, 0, 0, 0, 0, 1])
        nucleotide_array[:] = temporary_array

    f.close()


def main():
    parser = argparse.ArgumentParser(description="""description""")
    parser.add_argument('fasta_file')
    parser.add_argument('output_path')
    args = parser.parse_args()
    sequence_to_array(args.fasta_file, args.output_path)


if __name__ == '__main__':  # if inside terminal this gets called. if in python interpreter call main directly
    main()

