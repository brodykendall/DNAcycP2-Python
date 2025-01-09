"""DNAcycP2 - DNA Sequence Cyclizability Prediction

Usage:
    dnacycp2-cli -f <inputfile> <basename> [-L <chunk_length>] [-n <num_threads>]
    dnacycp2-cli -f -s <inputfile> <basename> [-L <chunk_length>] [-n <num_threads>]
    dnacycp2-cli -t <inputfile> <basename>
    dnacycp2-cli -t -s <inputfile> <basename>
    dnacycp2-cli (-h | --help)
    
Arguments:
    <inputfile> Input file name.
    <basename>  Output file name base.
    
Options:
    -h --help           Show this screen.
    -f                  FASTA mode.
    -t                  TXT mode.
    -s                  SMOOTH mode.
    -L <chunk_length>   Chunk length [default: 100000].
    -n <num_threads>    Number of threads [default: 1].
"""
from docopt import docopt
from dnacycp2 import cycle_fasta, cycle_txt
import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO

def main():
    arguments = docopt(__doc__)
    print("Input file: "+arguments['<inputfile>'])

    smooth = True if arguments['-s'] else False

    chunk_size = int(arguments['-L']) if arguments['-L'] else 100000
    num_threads = int(arguments['-n']) if arguments['-n'] else 1

    if arguments['-f']:
        cycle_fasta(arguments['<inputfile>'],
                    arguments['<basename>'],
                    smooth,
                    chunk_size,
                    num_threads)
    elif arguments['-t']:
        cycle_txt(arguments['<inputfile>'],
                  arguments['<basename>'],
                  smooth)
            
if __name__ == "__main__":
    main()
