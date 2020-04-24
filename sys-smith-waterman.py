"""

This is a simplified smith-waterman algorithm which 
inputs two sys input fasta files and outputs an alignment of the 
two.

Author: Felix O'Farrell

Date: April 2020

"""


#PELICAN and COELACANTH 


#import numpy for matrix generation and handling
import numpy as np
import sys


def fasta_reader(inp):  
   
    """
    Fasta file reader.
    
    Returns string of nucleotide sequnce present in the fasta file.
    This needs to be called once for each fasta file. 
    
    """

    with open(inp) as in_file:  
        for line in in_file.readlines():
            if line.startswith("A" or "T" or "C" or "G"):
                seq = line.rstrip()
                return(seq)
                
                
def scr_calc(m,i,j):  
   
    """
    Matrix Score Calculator.
    
    Updates matrix with coordinates of alignments between
    the two strings
    """
    
    diag = m[i-1][j-1] 
    up = m[i-1][j]
    left = m[i][j-1]
    highest = max(diag,up,left)
                    
    #matches
    if strA[i] == strB[j]:
            m[i,j] = diag + 1
           
    #non matches
    if strA[i] != strB[j]:
            m[i,j] = highest - 1
            if m[i,j] < 0:
                m[i,j] = 0


def backtracer(i, j, i_align, j_align):
    
    """
    Backtracer.
    
    Finds highest scoring path in the intitialised 
    matrix and relates that to the indexes 
    in strA and strB
    """
    
    diag = matrix[i-1][j-1]
    up = matrix[i-1][j]
    left = matrix[i][j-1]


    while diag > 0 or up > 0 or left > 0:
        
        
        diag = matrix[i-1][j-1]
        up = matrix[i-1][j]
        left = matrix[i][j-1]
        
      
        if diag > up and diag > left:
            i_align.append(strA[i])
            j_align.append(strB[j])
            i = i-1
            j = j-1
            
        if diag > up and diag == left:
            i_align.append(strA[i])
            j_align.append(strB[j])
            i = i-1
            j = j-1        
           
        if diag == up and diag > left:
            i_align.append(strA[i])
            j_align.append(strB[j])
            i = i-1
            j = j-1
           
        if diag == up and diag == left:
            i_align.append(strA[i])
            j_align.append(strB[j])
            i = i-1
            j = j-1
            
        if up > diag and up > left:
            i_align.append(strA[i])
            j_align.append("-")
            i = i-1
            j = j
        
        if left > diag and left > up:           
            i = i
            j = j-1
            i_align.append("-")
            j_align.append(strB[j])
        
        if diag == up  == left == 0:
            i_align.append(strA[i])
            j_align.append(strB[j])
            break


#Read in the two sequences
            
inp1 = sys.argv[1]
inp2 = sys.argv[2]
 
strA = ''
strB = ''       

strA = fasta_reader(inp1)
strB = fasta_reader(inp2)


#Add space to sequences. This is done to 
#create correct number of rows and collumns
#in the scoring matrix


strB = ' '+strB
strA = ' '+strA


#Reverse the order of the sequence. This is done
#simpify the Backtracer function.


strA = strA[::-1]
strB = strB[::-1]


#Create matrix with zeros 
matrix = np.zeros((len(strA),len(strB)))
 
#Initialise matrix using scr_calc funciton   

   
for i in range(1, len(strA)):
    for j in range(1, len(strB)):
        scrs = scr_calc(matrix, i, j)

       
#print the matrix to screen
        
        
print ("SW alignment matrix = ")
print(matrix)


#locate maximum scoring coordinates 
#in matrix


max_scrs = np.amax(matrix)
start = np.where(matrix == max_scrs)
 

#take integer values of coordinates with max score

i = int(start[0])
j = int(start[1])    


#create empyty alignment lists of i and 
#j coordinates to fill using Backtracer 
#function

i_align = []
j_align = []

#fill above lists with coordinates of each
#step in Backtracer
backtracer(i, j, i_align, j_align) 

#create string from alignment lists
ali_1 =''.join(i_align)
ali_2 =''.join(j_align)


#print alignments to screen

print ("Sequence1",ali_1)
print ("Sequence2",ali_2)


quit()
