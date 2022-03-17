#!/usr/bin/env python
import argparse

#Constants:
global DNA_BASES 
global RNA_BASES 

DNA_BASES ="ABDSD"
RNA_BASES ="GTSDSG"

# def get_args():

#         parser = argparse.ArgumentParser(description ="A program to parse protein fasta files and filter only the longest protein per gene")
#         parser.add_argument("-f", required=True)
#         #parser.add_argument("-o", required=True)
#         #parser.add_argument("-b", required=True)
#         #parser.add_argument("-z", help="specifies the path to the Zebra file" , required=True)

#         return parser.parse_args()

# args=get_args()
# input_file = args.f
#output_file = args.o
#biomart_table = args.b




#Functions:
def convert_phred(letter):
    """Converts a single character into a phred score"""
    value = ord (letter)
    score = value - 33
    return score

def qual_score(phred_score : str) ->int:
    # YOUR CODE HERE
    res = []
    for letter in phred_score:
         res.append(convert_phred(letter))
         avg_score = sum(res)/len(phred_score)
    return avg_score

def multi_to_single_line (f):
        First_Line=True
        seq = ""
        head_list = [] #list of all the headers
        seq_list = [] # list of all the sequences converted to single line
        with open(f, 'r') as f:
            print("BIOINFO MODULE")
            while True:
                    line = f.readline().strip() #Stripping all the lines
                    if line == '':
                            break
                    if line[0]=='>':
                            if First_Line:
                                    header = line
                                    First_Line=False
                            else:
                                    seq_list.append(header)
                                    seq_list.append(seq)
                                    header = line
                                    seq = ""

                    else:
                            seq+=line

        seq_list.append(header)
        seq_list.append(seq)
        return seq_list
        #writing the list to a file
        # with open("output_file", 'w') as fo:
        #     for item in seq_list:
        #         fo.write("%s\n" %item)
        # return output_file

def validate_base_seq(seq, RNAflag=False):
    DNA_bases = set('ATGCatcgNn')
    RNA_bases = set('AUGCaugc')
    return set(seq)<=(RNA_bases if RNAflag else DNA_bases)

def gc_content(input_file):
    sequence =""
    with open(input_file, 'r') as tf:
        for line in tf:
            if line.startswith('>'):
                seq_id = line.rstrip()[0:]
                
            else:
                sequence += line.rstrip()
    GC_content = float((sequence.count('G') + sequence.count('C'))) / len(sequence) * 100
    with open("output_file", 'w') as file_out:
        file_out.write("The GC content of file is\t %.2f%%" % (GC_content))
    return output_file 


if __name__ == "__main__":
    #Test case to check convert_phred
    assert convert_phred("A") == 32, "Phred score incorrect"
    assert convert_phred("@") == 31, "Phred score incorrect"
    assert convert_phred("#") == 2, "Phred score incorrect"
    print("Phred scores converted correctly")
    
    #Test case to check qual_score
    phred_score = "FFHHHHHJJJJIJIJJJIJJJJJJIIIJJJEHJJJJJJJIJIDGEHIJJFIGGGHFGHGFFF@EEDE@C??DDDDDDD@CDDDDBBDDDBDBDD@"
    assert qual_score(phred_score) == 37.62105263157895, "wrong average phred score"
    print("You calcluated the correct average phred score")

    #Test case to check qual_score
    #input_file = ~/bioinfo/Bi621/PS/ps4-sneharchalla
    #multi_to_single_line(input_file)

    #Test case to check validate_base-seq:
    print(validate_base_seq("ATTGCG"))
    print(validate_base_seq("ATTGpyCG"))

    #Test case to check gc_content:
    gc_content(input_file)