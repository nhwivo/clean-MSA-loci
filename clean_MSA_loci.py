#!/usr/bin/env python

"""
This program takes in a multiple sequence aignment (MSA) file and a reference genome name as input and cleans the region between the loci based on the specified reference genome. 
Ensure that the version of Python used is 2+

Example command: clean_MSA_loci.py FcC_supermatrix.fas 'BMORI'
The command above will result in a cleaned file named 'cleaned_loci_FcC_supermatrix.fas'

To specify name of output file: clean_MSA_loci.py FcC_supermatrix.fas 'BMORI' --out desired_name.fas

Nhi Vo    01 August 2022
"""

import argparse
import os

# COMMAND LINE ARGUMENTS:
parser = argparse.ArgumentParser()
# required parameters/inputs:
parser.add_argument("file", help="path to fasta file containing the multiple sequence alignment to be cleaned.")
parser.add_argument("ref", help="name of the reference genome in the MSA fasta file (does not have to be the whole name, just needs to be unique enough). Other sequences will be trimmed based on the gap of this sequence. ex: 'BMORI'")
# optional parameters/inputs:
parser.add_argument("--out", help="path/name of file output. False = output would be in current directory.") 
args = parser.parse_args()

class Sequence:
    def __init__(self, filename=''):
        """
        Initializes the object
        
        param mode: string, path of file to be edited
        """
        self.filename = filename  # name of file to be read 
        self.seqdict = {}  # dictionary of sequence names (dict keys) and genomic sequences (dict values) 
        
        if filename:
            self.open_file('r') 
            
    def open_file(self, mode):
        """
        Opens file contained in self.filename and sets the filehandle as self.file 
        If the file can't be opened, exit status will equal 1

        param mode: string, mode for opening file, usually 'r' or 'w'
        return: filehandle of opened file
        """
        try:
            self.file = open(self.filename, mode)
        except OSError:
            print('Error opening file named' + self.filename)
            # exit with status = 1
            exit(1)
            
    def next_line(self):
        """
        Returns the next line from the file
        
        return: logical, True if line is read, False when end of file 
        """
        line = self.file.readline().strip() 
        if not line: 
            # no more line, file is finished
            return False 
        
        self.line = line 
        return True  #return True when there is still line to be read 
    
    def obtain_reference(self, ref_name):
        """
        Obtain the reference genome from the MSA file
        
        param mode: string, name of reference genome.
        """
        basename = os.path.basename(self.filename)
        self.ref_filename = "reference_genome_" + basename  # name of the new reference genome file 
        bash_command = "grep -A1 '" + ref_name + "' " + self.filename + " > " + self.ref_filename
        os.system(bash_command)  # obtain the reference genome
        
    
    def read_file(self, seq_type, flank_ranges=''): 
        """
        Process the sequence data by determining flanking regions from the reference genome, and using those ranges to
        clean aligned sequences. 
        
        param mode: string, type of sequence, either 'reference_genome' (used to determine flanking regions)
        or 'aligned_sequences' (sequences to be edited). 
        """
        while self.next_line():  # iterate through every lines in the file
            if self.line.startswith('>'):  # create new dictonary key using sequence name/description 
                self.seqdict[self.line] = ''  # set value of dictionary key to be empty 
            else:
                last_key = list(self.seqdict.keys())[-1]  # access the last key created in the dictionary
                self.seqdict[last_key] += self.line  # append line (sequence) to the value of key 
        
        if seq_type == 'reference_genome':
            self.process_ref_genome()
        
        if seq_type == 'aligned_sequences':
            self.process_aligned_seq(flank_ranges)
    
    def process_ref_genome(self):
        """
        Create ranges of flanking regions based on the reference genome. The resulting ranges (coordinates) will then 
        determine the regions in other sequences to be cleaned. 
        """
        ref_seq_num = len(self.seqdict.keys())
        if ref_seq_num > 1:
            raise Exception("Argument for name of reference genome is not unique enough, resulting in multiple genomes being used. Please try another argument that is more specific.")
            
        if ref_seq_num == 0:
            raise Exception("Could not find the specified reference genome. Please check to make sure the argument matches name of reference genome.") 
        
        seq = list(self.seqdict.items())[0][1]  # access the reference genome sequence
        self.flank_range_list=[]  # list of dictionaries of flank ranges [{'id' : 1, 'start' : 20, 'end' : 25, 'total' : 6}]
        pos = 0  # keep track of position of characters in sequence string (first letter is 0)
        prev_char = False  # keep track of the previous character (true if '-', false if not '-')
        range_num = 1  # keep track of number of ranges 
        # iterate through each character in the reference sequence - keep count of character's position
        for char in seq:
            # check to see if character is '-'
            if char == '-':
                current_char = True  # set the current character = '-'
                if prev_char == False:  # if previous character isn't '-', then start a new range 
                    self.flank_range_list.append({'id' : range_num, 'start' : pos})
                    range_num+=1  

            else:  # if not '-'
                current_char = False  # set the current character as not '-'
                if prev_char == True:  # if previous character was '-', then prev position was the end 
                    temp_dict = self.flank_range_list[-1]
                    temp_dict['end'] = pos-1  # record the end of the range 
                    temp_dict['total'] = temp_dict['end'] - temp_dict['start'] + 1
                    
            pos+=1  # increase character position counter 
            prev_char = current_char  # set condition for previous character 
    
    def process_aligned_seq(self, flank_ranges):
        """
        The rest of the alignment is cleaned using the ranges/coordinates from the reference genome. 
        
        param mode: list of dictionaries of flanking regions/ranges (from the reference genome)
        """
        for key, value in self.seqdict.items():
            for each_range in flank_ranges:
                start = each_range['start']
                end = each_range['end'] 
                total = each_range['total']
                add = '-'*total
                self.seqdict[key] = value[:start]+add+value[end+1:]
                value = self.seqdict[key]
                
    def save_output(self, output_name):
        """
        Save the cleaned loci result
        
        param mode: logical, False if user did not specify name for output
        """
        print("Cleaning completed, saving result")
        if not output_name:  # user did not specify name for output 
            basename = os.path.basename(self.filename)
            output_name = 'cleaned_loci_'+basename  # name for output file 
        
        with open(output_name, 'w') as f:
            for key, value in self.seqdict.items():
                f.write(key)
                f.write('\n')
                f.write(value)
                f.write('\n')
                
        os.system("rm " + self.ref_filename)  # delete reference genome file 
        
        print("Resulting file named '" + output_name + "' saved - program is finished.")

if __name__ == '__main__':
    print("Processing <" + args.file + "> using <" + args.ref + "> as the reference genome.")
    MSA = Sequence(args.file)  # create object and opens the file
    MSA.obtain_reference(str(args.ref))  # obtain the reference genome from the file
    ref = Sequence(MSA.ref_filename)  # create object and open file
    ref.read_file('reference_genome')  # process the reference genome file 
    MSA.read_file('aligned_sequences', flank_ranges=ref.flank_range_list)  # process the MSA file 
    
    # SAVE FILE: 
    output_name = False  # no user specified name
    if args.out:  # user specified output file name 
        output_name=args.out
    MSA.save_output(output_name)
                
