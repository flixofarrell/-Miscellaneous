"""

This is a simplified smith-waterman algorithm which 
inputs two sys input fasta files and outputs an alignment of the 
two.

Author: Felix O'Farrell

Date: April 2020

"""


#Import numpy for matrix generation.
import numpy as np
import sys 


def fasta_reader(inp):  
   
    """
    Fasta File Reader.
    
    Returns string of nucleotide seqeunce present in the fasta file.
    This needs to be called once for each FASTA file. Inp is the FASTA filename 
    which is to be inputted.
    
    """
    #inp is hard coded as "Sequence1/2.fasta in this script".
    with open(inp) as in_file:  
        for line in in_file.readlines():
            #Guarantees sequence is pulled from the FASTA file not the title 
            if line[0].isalpha():
                seq = line.rstrip()
                return (seq)
                
                
def scr_calc(m,i,j):  
   
    """
    Matrix Score Calculator.
    
    Updates the zeros matrix with scores based on matching and non-matching
    indices of strA and strB. Matches are given a score of +1, mismatches and
    gaps are given a score of -1. 
    

    """
    #Indices contributing value at any given position in matrix.   
    diag = m[i-1][j-1] 
    up = m[i-1][j]
    left = m[i][j-1]
    #Takes highest out of all 3 options above.
    highest = max(diag,up,left)
                    
    #Matches.
    if strA[i] == strB[j]: 
        m[i,j] = diag + 1
           
    #Non matches.
    if strA[i] != strB[j]:
        m[i,j] = highest - 1
        #Guarantees no value in matrix falls below 0.
        if m[i,j] < 0:
            m[i,j] = 0


def backtracer(i, j, i_align, j_align):
    
    """
    Backtracer.
    
    Finds highest scoring path in the initialised matrix and relates that to 
    the indices in strA and strB. In each back-step, the highest scoring move 
    is identified and the new values of i and j are updated accordingly.
    
    """
    #Directions to move during backtracing.
    diag = matrix[i-1][j-1]
    up = matrix[i-1][j]
    left = matrix[i][j-1]
    lowest = min(diag,up,left)

    #Guarantees the loop continues when atleast 1 of 3 possible moves > 0.
    while lowest > 0:
        #Updates new directions to move in matrix for each loop.
        diag = matrix[i-1][j-1]
        up = matrix[i-1][j]
        left = matrix[i][j-1]
        
        #Prioritse the diagonal move over up and left.
        if diag >= up and diag >= left:
            i_align.append(strA[i])
            j_align.append(strB[j])
            #Update new positions in matrix of i and j.
            i = i-1
            j = j-1
        
        #When the up move is the highest scoring append a '-' to j sequence.                            
        if up > diag and up > left:
            i_align.append(strA[i])
            j_align.append("-")
            i = i-1
        
        #When the left move is the highest scoring append a '-' to i sequence.                             
        if left > diag and left > up:           
            i_align.append("-")
            j_align.append(strB[j])
            j = j-1

        #Break loop when all next moves = 0.
        if diag == up == left == 0:
            i_align.append(strA[i])
            j_align.append(strB[j])
            break



#Read in the two sequences as inputs to fasta_reader function.           
inp1 = sys.argv[1]
inp2 = sys.argv[2]


#Create strA and strB for each input sequence.
strA = ''
strB = ''       


#Call fasta_reader function to read in each fasta input.
strA = fasta_reader(inp1)
strB = fasta_reader(inp2)

       
#Add space to sequences. This is done to create the correct number of rows 
#and columns in the scoring matrix. 
strA = ' '+ strA
strB = ' '+ strB


#Reverse the order of the sequence. This is done to simplify the backtracer 
#function.
strA = strA[::-1]
strB = strB[::-1]


#Create matrix with zeros of correct number of rows and columns.
matrix = np.zeros((len(strA),len(strB)))

 
#Initialise matrix using scr_calc function.   
for i in range(1, len(strA)):
    for j in range(1, len(strB)):
        scrs = scr_calc(matrix, i, j)


#Locate maximum scoring coordinates in matrix.
max_scrs = np.amax(matrix)
start = np.where(matrix == max_scrs)
 

#Take integer values of coordinates with max score. Numpy functions return 
#array.
i = int(start[0])
j = int(start[1])    


#Create empty alignment lists of i and j coordinates to fill using backtracer 
#function.
i_align = []
j_align = []


#Fill above lists with coordinates of each step in backtracer.
backtracer(i, j, i_align, j_align)
 

#Create strings from alignment lists.
ali_1 = ''.join(i_align)
ali_2 = ''.join(j_align)


#Print alignment to screen.
print ("Sequence1", ali_1)
print ("Sequence2", ali_2)
