#!/usr/bin/env python

##############################################################################################
#### SAMI CHAABAN v. 2017-12-12 ##############################################################
##############################################################################################

##############################################################################################
#Input: two fasta files with multiple sequences
#Output: txt file formated to be read by chimera's attributes
#
#1. All sequences must be the same length. Generate the gaps with mafft:
#   > mafft-linsi  --retree 2 --inputorder "all-beta.fasta" > "all-beta-gaps.fasta"
#2. Separate sequences into two groups of interest
#3. Run this script
#4. Load attributes on atomic model in chimera using the commandline.
#   > Make sure alpha and beta tubulin are separate models
#   > cd /path/to/similarity/files
#   > defattr name-of-alpha-similarity.txt spec #1
#   > defattr name-of-beta-similarity.txt spec #0
##############################################################################################

import sys
import os
from Bio import motifs
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import optparse

def setupParserOptions():
    parser = optparse.OptionParser()
    parser.add_option("--i1",dest="input1",type="string",metavar="FILE",
                help="sequences1.fasta")
    parser.add_option("--i2",dest="input2",type="string",metavar="FILE",
        help="sequences2.fasta")

    options,args = parser.parse_args()

    if len(args) > 1:
            parser.error("Unknown commandline options: " +str(args))

    if len(sys.argv) < 2:
            parser.print_help()
            sys.exit()

    params={}

    for i in parser.option_list:
            if isinstance(i.dest,str):
                    params[i.dest] = getattr(options,i.dest)
    return params

def mainloop(params):

    ################
    #Remove dashes

    input1filename = params["input1"]
    input2filename = params["input2"]

    input1filename_nodash = input1filename + '.tmp'
    input2filename_nodash = input2filename + '.tmp'

    file1 = open(input1filename, 'r')
    file1tmp = open(input1filename_nodash, 'w')
    file2 = open(input2filename, 'r')
    file2tmp = open(input2filename_nodash, 'w')

    for line in file1:
        file1tmp.write(line.replace('-', 'X'))
    for line in file2:
        file2tmp.write(line.replace('-', 'X'))

    file1.close()
    file2.close()
    file1tmp.close()
    file2tmp.close()

    ################
    #Create matrices

    ###INPUT1

    input1seqs = []

    for seq_record in SeqIO.parse(input1filename_nodash, "fasta"):
        input1seqs.append(Seq(str(seq_record.seq),IUPAC.extended_protein))
        
    input1_motifs = motifs.create(input1seqs)
    input1_pwm = input1_motifs.pwm
    input1_consensus = input1_motifs.consensus

    print("consensus1: " + input1_consensus + "\n")

    ###INPUT2

    input2seqs = []

    for seq_record in SeqIO.parse(input2filename_nodash, "fasta"):
        input2seqs.append(Seq(str(seq_record.seq),IUPAC.extended_protein))
        
    input2_motifs = motifs.create(input2seqs)
    input2_pwm = input2_motifs.pwm
    input2_consensus = input2_motifs.consensus

    print("consensus2: " + input2_consensus + "\n")

    ################
    #Calculate similarity

    total = len(input2_motifs)

    similarity = []

    for p in range(0,total-1):
        sim = 0
        if (input2_consensus[p] == 'X' or input1_consensus[p] == 'X'):
            sim = 1
        else:
            for a in range(0,25):
                sim = sim + input2_pwm[a][p]*input1_pwm[a][p]
        similarity.append(sim*100)

    ################
    #Account for gaps and output each file and output attribute

    for filename in [input1filename_nodash, input2filename_nodash]:
      
        for seq_record in SeqIO.parse(filename, "fasta"):
            
            print("Working on " + seq_record.id)
            similarity_fix = []
            
            for p in range(0,total-1):
                if seq_record[p] != "X":
                    similarity_fix.append(similarity[p])
                
            print(" - " + str(len(similarity_fix)) + " residues")
            
            outputname = "similarity-maped-on-" + str(seq_record.id) + ".txt"
            outputname = outputname.replace("|", "-")
            
            output = open(outputname,"w")

            output.write("attribute: similarity\n")
            output.write("match mode: 1-to-1\n")
            output.write("recipient: residues\n")

            for p in range(0,len(similarity_fix)-1):
                output.write("\t:" + str(p+1) + "\t" + str(similarity_fix[p]) + "\n")

            output.close()

    os.remove(input1filename_nodash)
    os.remove(input2filename_nodash)

if __name__ == "__main__":
    params = setupParserOptions()
    mainloop(params)
