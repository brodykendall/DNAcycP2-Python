"""DNAcycP - DNA Sequence Cyclizability Prediction

Usage:
    dnacycpv2-cli -f <inputfile> <basename>
    dnacycpv2-cli -t <inputfile> <basename>
    dnacycpv2-cli (-h | --help)
    
Arguments:
    <inputfile> Input file name.
    <basename>  Output file name base.
    
Options:
    -h --help   Show this screen.
    -f          FASTA mode.
    -t          TXT mode.
    
"""
from docopt import docopt
from dnacycpv2 import cycle_fasta, cycle_txt
import keras
import pandas as pd
import numpy as np
from numpy import array
from Bio import SeqIO

def main():
    arguments = docopt(__doc__)
    print("Input file: "+arguments['<inputfile>'])

    if arguments['-f']:
        cycle_fasta(arguments['<inputfile>'],
            arguments['<basename>'])
    elif arguments['-t']:
        cycle_txt(arguments['<inputfile>'],
            arguments['<basename>'])
            
if __name__ == "__main__":
    main()
