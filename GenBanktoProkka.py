#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri May 10 10:59:02 2019

@author: leslieannesmith
"""
''' Written by Leslie Smith. Function to take in a fasta file (downloaded from GenBank) and parse it into 
files containing 100 sequences each (files that will be small enough to run Prokka on).
If you want to change the amount of sequences in each parsed file change the variable countMax. Further
descriptions below.'''
#prokka takes in a .fasta file (would be downloaded from genBank)
import re
import sys
import os

from Bio import SeqIO

'''variable declartation and opening first files:'''
#i is the count used to change file name(when countMax is reached in each file and a new file is created)
i = 0 
#countMax is the maximum number of sequences user wants in each file 
countMax = 100
#counter is counting the number of sequences in each file, when this is equal to countMax a new file will be created 
counter =0
#inFile is name of file containing sequences downloaded from genBank
inFile = 'GenBankPlasmids2019.fasta'
#outName is the name of the file the sequences will be sent to after being parsed, they will be numbered sequentially
#as countMax is reached and a new file is created 
outName = 'GenBankInputForProkka'
#open up handle to file and read it 
inHandle = open(inFile,"rU")
#output from PROKKA/name of PROKKA annotated files
prokkaOutDir = 'GenBankAnnotatedPlamids'
#open up handle of outfile to write to (str(i) is what is counting each file sequentially)
#a seperate variable was made here for outFile to easily increase each folder name sequentially(at the end of the for loop)
outFile = outName+str(i)+'.fasta'
outHandle = open(outFile,"w")

'''INPUT: Fasta file downloaded from GenBank holding all sequences that will be run through PROKKA
   OUTPUT: Multiple files, first set of files will be the parsed files that hold the sequences from 
   GenBank(the outName variable holding sequences now parsed into several different files 
   easier for PROKKA to annotate). The second set of files are the files annotated from PROKKA. 
   The outName files and PROKKA annotated files will be coordinated by the numbers in the file title.
   
    Funtion: For each sequence in the inFile, write the sequence to an outfile given by user.
    OutFiles will have the same basename but will have a number (i in the code)
    that increases with each new file, the user will define how many sequences they 
    want in each file (labeled above as countMax).
    A new file will be created when the current file reaches the 
    maximum number of sequences. '''
for seq in SeqIO.parse(inHandle,"fasta"):
    #write the current sequence to new file 
    SeqIO.write(seq,outHandle,"fasta")
    #increase counter to keep track of # sequences in file
    counter = counter+1
    #if max number of files is reached call prokka on current file, send .tbl files from prokka
    #output to a new  file and close the current file and open a new one 
    if counter == countMax:
        print("Max number of sequences for folder reached, calling Prokka:")
       
        #command to call prokka 
        os.system('prokka --outdir '+prokkaOutDir+str(i)+' --proteins Summer19CompleteDatabase.faa --evalue 0.001 --addgenes '+outFile)
        print('Sending Prokka .tbl files to a new folder')
        #change folder after double arrows to change the folder where Prokka .tbl files will go
        '''os.system('cat prokkaFileOut'+str(i)+'/*.tbl >> ProkkaTBLFiles'+str(i)+'.tbl')#puts all prokka files together
        ^^^line above does not do anyting but i am not going to erase it quite yet'''
        #make sure that the file made by prokka is executable 
        os.system('chmod 755 prokkaFileOut'+str(i))
        print ('Complete')
        #close previous outhandle 
        outHandle.close()
        #increase i to keep track of new file names 
        i= i+1
        #open a new file (same name with new i)
        outFile = outName+str(i)+'.fasta'
        outHandle = open(outFile,"w")
        #set counter back to zero
        counter = 0        
        
        




        
       
        
        
        
    
    
    