#Imports
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import os

# Functions
def parsing_fasta(input_file):
    """A function that parses a fasta file.
    Args:
    input_file: [str]
    The path of the parsed fasta file.
    Returns:
    sequences_str: [list]
    A list of strings, the elements are the sequences of the fasta file.
    sequences: [list]
    A list of Bio.Seq.Seq objects (Sequence Records).
    """
    sequences=[]
    sequences_str = []
    for Seq_record in SeqIO.parse(input_file, 'fasta'):
        format_string = "%s" % Seq_record.seq
        sequences_str.append(str(format_string))
        sequences.append(Seq_record)

    return(sequences_str,sequences)

def consensus_sequence(sequences_list):
    """ A function that gets the consensus sequence of a set of DNA sequences.
    Args:
    sequences_list: [list]
    The list of DNA sequences.
    n: [integer]
    The length of the sequences.
    Returns:
    consensus_sequence: [string]
    The required consensus sequence.
    """
    n= len(sequences_list[0])

    profile_matrix = { 'A': [0]*n , 'C': [0]*n,
    'G' : [0]*n , 'T': [0]*n, 'N': [0]*n,'-':[0]*n}

    #Calculating the profile matrix
    for dna in sequences_list:
        for position, nucleotide in enumerate(dna):
                profile_matrix[nucleotide][position] += 1
    output = []

    #Getting the consensus sequence from the profile matrix.
    for position in range(n):
        max_count = 0
        max_nucleotide = ''
        for nucleotide in ['A','C','G','T']:
            count = profile_matrix[nucleotide][position]
            if count > max_count:
                max_count = count
                max_nucleotide = nucleotide
        output.append(max_nucleotide)
    consensus_sequence = "".join(output)
    return(consensus_sequence)

def sequence_content(sequence, flag):
    """ A function that calculates the A, C, G, T, CG content of a sequence.
    Args:
    sequence: [str]
    The sequence for which the content is calculated.
    flag: [str]
    Which content is required, the flag is either A,C,G,T, CG or ALL.
    Returns:
    content: [float or str]
    The percentage of the required content.
    If the flag is ALL, a string shows all contents.
    all_list: [list] (optional)
    If the flag is ALL, a list of floats has the contents in the order A, C, G, T and CG.
    """
    sequence_length = len(sequence)
    #Counting all the A, C, G, T and CG
    A_count = sequence.count('A')
    C_count = sequence.count('C')
    G_count = sequence.count('G')
    T_count = sequence.count('T')
    CG_count = C_count + G_count
    #Calculating the contents --> 100*(count/sequence_length)
    A_content= round(100*(A_count/sequence_length),2)
    C_content= round(100*(C_count/sequence_length),2)
    G_content= round(100*(G_count/sequence_length),2)
    T_content= round(100*(T_count/sequence_length),2)
    CG_content= round(100*(CG_count/sequence_length),2)
    #Check for the flag, and return the asked content
    if flag == 'A':
        return(A_content)
    elif flag == 'C':
        return(C_content)
    elif flag == 'G':
        return(G_content)
    elif flag == 'T':
        return(T_content)
    elif flag == 'CG':
        return(CG_content)
    elif flag == 'ALL':
        ALL_content = f"A content = {A_content}%.\nC content = {C_content}%.\nG content = {G_content}%.\nT content = {T_content}%.\nCG content = {CG_content}%."
        ALL_list = [A_content, C_content, G_content, T_content, CG_content]
        return(ALL_content, ALL_list)
    else:
        print("The flag is incorrect")

def nucleotide_check(nucleotide_column,percentage):
    """A function that checks for the occurence of any nucleotide
    more than a certain threshold.
    Args: 
    nucleotide_column [list]:
    A list of nucleotides, which represents a column of aligned sequences.
    percentage [float]:
    A float that should vary from 0.1 to 1, which is the accepted percentage 
    for the nucleotide to be assigned to the return value.
    Returns:
    A string which is either: 'A', 'C', 'G', 'T', 'N', '-' or 'D'
    D is assigned when there is no nucleotide that is present more than the asked percentage.
    """
    threshold = percentage*len(nucleotide_column)
    count_A = nucleotide_column.count('A')
    count_C = nucleotide_column.count('C')
    count_G = nucleotide_column.count('G')
    count_T = nucleotide_column.count('T')
    count_N = nucleotide_column.count('N')
    count_gap = nucleotide_column.count('-')
    #Check for the dominant nucleotide, if exists.
    if count_A >= threshold: 
        return('A')
    elif count_C >= threshold:
        return('C')
    elif count_G >= threshold:
        return('G')
    elif count_T >= threshold:
        return('T')
    elif count_N >= threshold:
        return('N')
    elif count_gap >= threshold:
        return('-')
    else:
        return('D')

def conserved_regions_seq(sequences_list,percentage):
    """A function that takes a list of sequences and returns a sequence 
    that contains the conserved regions between all the sequences.
    This function is used as an intermediate step for calculating the dissimilar
    regions.
    It uses the function "nucleotide_check" above. 
    Example: ['ACTG', 'ACCA' , 'ACTC','ATTG'] ---> result = 'ACTD'
    Args: 
    sequences_list: [list]
    A list of aligned sequences.
    percentage: [float]
    A float that should vary from 0.1 to 1, which is the accepted percentage 
    for the nucleotide to be assigned to the return value.
    It is an argument for "nucleotide check" function.
    Returns:
    result: [str]
    The resulting sequence.
    """
    mycolumn = []
    result = []
    #Loop over all the columns of the sequences and call nucleotide check,
    #to get the nucleotide at each column.
    for i in range(len(sequences_list[0])):
        for j in range(len(sequences_list)):
            mycolumn.append(sequences_list[j][i])
        result.append(nucleotide_check(mycolumn,percentage))
        mycolumn=[]
    #Get the resulting representative sequence.
    result = "".join(result)
    return(result)

def dissimilarRegions(seq1, seq2):
    """A function that extracts the dissimilar regions between two sequences.
    Args:
    seq1,seq2: [str]
    The two sequences.
    Returns:
    allRegion:
    A list of lists, the lists contain the two dissimilar regions in the two 
    sequences, and the indecies of starting and ending of these regions.
    """
    i,k = 0,0
    flag = False
    region, allRegion = [],[]
    for count in range(len(seq1)):
        if seq1[i]  == seq2[k] or seq1[i] == 'D' or seq2[k] == 'D':
            if flag: 
                region.append(i-1)
                region.append(seq1[region[0]:i])
                region.append(seq2[region[0]:i])
                allRegion.append(region)
                region = []
                flag = False
        else:
            if region == []:
                flag = True
                region.append(i)
            if i == len(seq1) - 1:
                region.append(i)
                region.append(seq1[region[0]:i+1])
                region.append(seq2[region[0]:i+1])
                allRegion.append(region)
                break
        i += 1
        k += 1
    return allRegion


# Code for calculations on delta and omicron variants

#For the file to run properly, please run the file in the folder its sent in (for the pathes to be correct).
#We use the first 10 sequences in fasta files for our analysis.

current_directory = os.getcwd()

omicron_str, omicron_seq = parsing_fasta(f"{current_directory}\Data\SARS-Cov2-Data\Ghana\Delta\gisaid_hcov-19_2021_12_31_12.fasta")
delta_str, delta_seq = parsing_fasta(f"{current_directory}\Data\SARS-Cov2-Data\Ghana\Omicron\gisaid_hcov-19_2021_12_31_12.fasta")

#Fasta file of multiple alignment of omicron sequences.
omicron_alignment,_ = parsing_fasta(f"{current_directory}\Data\omicron_aligment.fas")

#Fasta file of multiple alignment of delta sequences.
delta_alignment,_ = parsing_fasta(f"{current_directory}\Data\delta_alignment.fas")

#First: Consensus Sequence:

delta_consensus = consensus_sequence(delta_alignment)
record = SeqRecord(Seq(delta_consensus), id ="DeltaCons123", description ="Consensus sequence")
SeqIO.write(record, "CodeResults/Consensus.fasta", "fasta")

#After getting the result, we took the "Consensus.fasta" file and applied multiple alignment 
#with the omicron sequences and the delta consensus sequence.
#Parsing that file:
casevsconsensus,_ = parsing_fasta(f"{current_directory}\Data\CaseVsConsensus.fas")

#Second: A, C, G, T and CG content:

omicron_dictionary = {}
for i in range(10):
    _, omicron_results = sequence_content(omicron_str[i],"ALL")
    omicron_dictionary[f'{omicron_seq[i].id}'] = omicron_results

delta_dictionary={}
_, delta_results = sequence_content(delta_consensus,"ALL")
delta_dictionary['Consensus Sequence'] = delta_results

omicron_dataframe = pd.DataFrame.from_dict(omicron_dictionary, orient='index',
                       columns=['A content', 'C Content', 'G content', 'T content', 'CG content'])
delta_dataframe = pd.DataFrame.from_dict(delta_dictionary, orient='index',
                       columns=['A content', 'C Content', 'G content', 'T content', 'CG content'])

#Extracting an excel sheet for analysis.
omicron_dataframe.to_excel("CodeResults/Omicron_Content.xlsx")
delta_dataframe.to_excel("CodeResults/Delta_Content.xlsx")

#Third : Extracting dissimilar regions:

#Getting the omicron sequences aligned with the consensus sequence in a list.
omicron_for_comparison = []
for i in range(0,len(casevsconsensus)-1):
    omicron_for_comparison.append(casevsconsensus[i])

#Getting the omicron representitive sequence.
conserved_omicron = conserved_regions_seq(omicron_for_comparison,0.7)

#Getting the dissimilar regions between the omicron representitive sequence and the consensus sequence.
dissimilar_regions = dissimilarRegions(conserved_omicron, casevsconsensus[10])
dissimilar_dataframe = pd.DataFrame(dissimilar_regions[0:], columns=['Staring Index','Ending Index','Region in Omicron','Region in Consensus'])
dissimilar_excel = dissimilar_dataframe.to_excel("CodeResults/DissimilarRegions.xlsx")