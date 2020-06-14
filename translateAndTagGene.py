#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 30 10:16:37 2019

@author: leslieannesmith
"""
'''Written by Leslie Smith, summer2019.  Program has several defined functions.
***FURTHER EXPLANATIONS OF EACH FUNCTION THROUGHOUT*** 
Functions listed in order of appearance:
    
-1-the function translate() will translate DNA nucleotide sequence into 
into a protein sequence. Some the sequences given in the online 
databases are given in nucleotide sequences and must be translated into proteins before
being added to our databse.

-2-the function resistanceGeneTable() will take in a given directory that contains 
the resistance .fsa files and will write the Resistance type, Gene 
name, and Accession number to a csv file.

-3-the function tagResGeneName() will take in a given directory 
that contains the resistance .fsa files and alter the gene title to contain 
the accession number, gene name, resistance type, and which database the gene came from.
EX:>JX440350 ~~~blaCMY~~~blaCMY-70_1_:_beta-lactam-resFinderdb~~~

-4-the function tagToxinGeneName() will take in a given directory that contains the Toxin/Antitoxin .fsa files
and reformat the gene title to include the accession number, type of Toxin/Antitoxin, gene description, and 
which database the gene came from.
EX:>NP_718721.1 ~~~TypeIIToxin~~~ toxin-antitoxin system antidote Mnt family _:_TypeIIToxin-TAfinderdb~~~

-5-the function tagIntGeneName() will take in a given directory that contains the Integrase .fsa files. 
NOTE: we chose only to take the INTI1 genes from the I-VIP database(an online integrase database).  Funciton will
reformat the gene title of the Integrases to include the accession number, class of integrase, and integrase label.
All of the integrase genes here are from the I-VIP database.
EX:>YP_003675754.1 ~~~Int1~~~class1_:_Integrase~~~'''


from Bio.Alphabet import generic_dna
import os
import csv 
import re
import Bio
from Bio import SeqIO
from Bio import Seq

'''-1-INPUT: directory that holds files that hold nucleotide sequences that need to be translated into protein sequences
   OUTPUT: a directory that holds the files(same names as before), now translated into protein sequences. 
   
Function will translate nucleotide sequences that will be read in from a file directory given at function call.
It will create a new file that contains the sequences translated as proteins
this file will later be fed into the function(depending on which database you are working on) that will
reformat the database and write it to our database 
EX of execution of both functions at once: tagResGeneName(translate('resfinder_db'))'''


#OutFolder is folder that will be made to contain all of the protein sequences
OutFolder = 'ResistanceProteins624'

'''INPUT(directory1) HERE SHOULD BE A DIRECTORY THAT HOLDS THE FILES OF WHERE THE SEQUENCES FROM DOWNLOADED DATABASE
    CURRENTLY ARE'''

def translate(directory1):
    #counter to count how many sequences were not able to be translated(incorrect length for codon translation)
    counter = 0
    #create a folder that all newly translated sequences will be sent to
    #(use complete file path to make sure the DNA sequences are where you want them to be)
    os.system('mkdir /Users/leslieannesmith/SummerResearch19/'+OutFolder)
    #cycle through files in directory given by user that contains sequences to be translated
    #get the ones ending in .fsa to be translated into proteins and written to OutFolder
    
    for filename in os.listdir(directory1):
        if filename.endswith('.fsa'):
            
            #OutFileName is the path to sequences files which are now translated into proteins
            #they have the same name as they did before but will now be proteins in OutFolder
            OutFileName = OutFolder+'/'+filename
            OutFile = open(OutFileName, 'w')
            
            #filePath is a path to where original DNA sequences currently are 
            filePath = '/Users/leslieannesmith/SummerResearch19/'+directory1+'/'+filename
            with open(filePath) as handle:
                #cycling through each sequence one at a time using biopython
                for record in SeqIO.parse(handle, "fasta"):
                        
                    #make sure sequence can be translated
                    if len(record.seq)%3 ==0:
                        #translate the sequence and write it to OutFile in fasta format
                        OutFile.write('>'+record.id+'\n')
                        seqRec = record.seq.translate()
                        OutFile.write(str(seqRec)+'\n')
                    else:
                        counter = counter+1#counting the sequences not translated

            OutFile.close()  
    print('There were '+counter+' sequences unable to be translated')
    print('Finished')
    return OutFolder


##################################################################################################################

'''-2-INPUT: directory that holds the files with the resistance protein sequences BEFORE THE SEQUENCES HAVE BEEN REFORMATTED 
TO FIT OUR DATABASE
    OUTPUT:csv file containing all the resistance types, gene names, and accession number from the resFinder database

Function will make a csv file with the resistance gene accession number and name.
file will take in the output from tagGene (function below) and will find the accession 
and geneName/type of gene and will put them into a new csv file table 
table gives list of Resistance Type, GeneName, and Accesstion Number '''

#ResTableOutFile is the name of the csv file the resistance type, gene name, and accession number will be written to
ResTableOutFile = 'resistanceGeneTable1.csv'
def resistanceGeneTable(directory):
    #find all files in directory that end in .fsa
    for filename in os.listdir(directory):
        if filename.endswith('.fsa'):
            #create an absolute path to the current .fsa file found
            filePath = '/Users/leslieannesmith/SummerResearch19/'+filename
            #strip the .fsa off the file name so filename can be used as the resistance type
            geneName = filename.strip('.fsa')
            
            #label the columns in new csv file
            with open(ResTableOutFile, 'a') as writeFile:
                 writer = csv.writer(writeFile)
                 name1 = 'Resistance Type'
                 name2 = 'Gene Name'
                 name3 = 'Accession Number'
                 writer.writerow([name1, name2, name3])
                 
                 #open the current .fsa file
            with open(filePath) as f:
                #read the file line by line, find the title line of each sequence and use
                #regular expressions to match the resistance type, gene name, and accession number
                data = f.readlines()
                for line in data:
                    if(line.find('>') != -1):
                        #regExTitle matches the entire title line, then renamed as completeTitle
                        regExTitle = re.search("^>.+(\d)$|^>.+(\d\s)$", line)
                        completeTitle = regExTitle.group()
                        
                        #if there is a whitespace at the end of the line, remove it
                        if completeTitle[-1] == ' ':
                            completeTitle = completeTitle[:-1]
                        
                        #regExAccessNum matches the accession number within completeTitle, then renamed as acceessionNum
                        regExAccessNum = re.search("[A-Z]+\d+\d{4}$",completeTitle)
                        if regExAccessNum:
                            accessionNum = regExAccessNum.group()
                        
                        #RegEx matches the name within completeTitle, then renamed as name
                        regExName = re.search(">(.+_\d{1,2}_)",completeTitle)
                        if regExName:
                            name = regExName.group()
                            #get rid of the carat at the beginning of name
                            name = name.strip('>')
                            #get rid of the extra digits in the name for simplicity
                            excess = re.search("_\d{1,2}_", name)
                            if excess:
                               excess = excess.group()
                               name = name.strip(excess)
                        #create empty list so it can be appended and add the [geneName,name,accessionNum] to it
                        outList = []
                        outList.append([geneName,name,accessionNum])
                        writer.writerows(outList)##removed some stuff, test this
                        writeFile.close()

##################################################################################################################

'''-3-INPUT: directory that contains files with the resistance protein sequences you want to reformat and add to database. 
       NOTE: this function uses the names of the files to tag the type of resistance of each gene when writing
       to the database.
   OUTPUT: function will append a file, adding all genes with a reformatted title. User should input the file/database
   they are creating (should be the same file/database for all database writing functions in this file) 

Function reads in files from a directory that contain genomic sequences (in this case the sequences should be 
resistance genes) and reformats their title to fit into out own database which is compatibe with PROKKA(annotation tool).
Funciton will also tag each gene with the type of resistance as well as which database it came from.
EX: >JX440350 ~~~blaCMY~~~blaCMY-70_1_:_beta-lactam-resFinderdb~~~'''

#directory is the name of the directory containing files of resistance protein sequences you want to label 
def tagResGeneName(directory):
    #open a file to write to, all files ending in .fsa will be written to this file
    #(here I edited our current database)
    OutFile = open('Summer19CompleteDatabase.faa','a')#name of database to write to
    #cycle through the files in directory given at method call, find files that end in .fsa
    for filename in os.listdir(directory):
        if filename.endswith('.fsa'):
            
            #split the filename by . and get the first part of that split (gets the resistance type from the filename)
            filename1 = filename.split('.')
            filename2 = filename1[0]
            
            #create an absolute file path to where the sequences currently here
            #NOTE:earlier in the script we assigned OutFolder to the place where we wanted our newly translated 
            #protein sequences to go
            filePath = '/Users/leslieannesmith/SummerResearch19/'+OutFolder+'/'+filename 
            
            #open file to read from
            with open(filePath) as f:
                data = f.readlines()
                
                #loop through lines to find lines containing title name 
                for line in data:
                    if(line.find('>') != -1):
                        #RE will match entire title line and assigns to variable 'completeTitle'
                        regExTitle = re.search("^>.+(\d)$|^>.+(\d\s)$", line)
                        completeTitle = regExTitle.group()
                        #if the final character in the regExTitle is a white space, remove white space
                        if completeTitle[-1] == ' ':
                            completeTitle = completeTitle[:-1]
                        
                        #RE matches just the accession number and assigns to variable 'accessionNum'
                        regExAccessNum = re.search("[A-Z]+\d+\d{4}$",completeTitle)
                        if regExAccessNum:
                            accessionNum = regExAccessNum.group()
                        
                        #Re matches just the name and assigns to variable 'name'
                        regExName = re.search(">(.+_\d{1,2}_)",completeTitle)
                        
                        if regExName:
                            name = regExName.group()
                            #get rid of the '>' at the beginning (caret is grouped in the wrong place here)
                            if name[0]=='>':name = name[1:]
                            #get short version of gene name to fill in a section of the title PROKK needs
                            #assign to variable 'simpleName'
                            simpleName = re.search('^[A-Za-z]+', name)
                            if simpleName:
                                simpleName = simpleName.group()
                        #create that has new appropriate title for gene
                        newTitleName = '>'+accessionNum+' ~~~'+simpleName+'~~~'+name+':_'+filename2+'-resFinderdb~~~\n'
                        line = newTitleName
                    OutFile.write(line)
    print('Finished')
            
##################################################################################################################
'''-4-INPUT: directory that holds the files of the Toxin/Antitoxin sequences to reformat and add to database.
   OUTPUT:function will append a file, adding all genes with a reformatted title. User should input the file/database
   they are creating (should be the same file/database for all database writing functions in this file)

Function reads in files from a directory that contain genomic sequences (in this case the sequences should be 
toxin/antitoxin genes) and reformats their title to fit into out own database which is compatibe with PROKKA(annotation tool).
Function will also tag each gene with the type of resistance as well as which database it came from.
EX: >NP_718721.1 ~~~TypeIIToxin~~~ toxin-antitoxin system antidote Mnt family _:_TypeIIToxin-TAfinderdb~~~'''

#enter directory containing files of DNA sequences you want to label (put in directory where 
#the files from the ToxinAntitoxin database were downloaded to)
def tagToxinGeneName(directory):

    #open a file to write to, all files ending in .fsa will be written to this file(here I edited our currentdatabase),
    #again if you are creating a new database the file here should be the same name as used in the other functions
    OutFile = open('Summer19CompleteDatabase.faa','a')
    #NOTE: some of the sequences in the TAfinder database did not match any uniform format
    #or they were missing some information, any sequences that were this way I had written to a seperate file
    #then I went back and manually fixed them/got needed information. There should only be around 5-6 sequences that
    #this happens with
    notReadable = open('Summer19TAunmatched.faa','a')
    #find all files ending with .fsa in directory
    for filename in os.listdir(directory):
        if filename.endswith('.fsa'):
            
            #create an absolute path to the current filename where the TA protein sequences are
            filePath = '/Users/leslieannesmith/SummerResearch19/'+directory+'/'+filename 
            filename = filename.strip('.fsa')
            
            #open file to read as a handle so sequences may be parsed out 
            with open(filePath) as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    #get the record ID from the current sequence and cast it into a string
                    recordID = str(record.description)
                    
                    #Note for following code:
                    #There are several different formats in the TAfinder database. Based on specific formats/things found in the 
                    #description of each sequence the information needed to create a new title and write the sequence
                    #to our database is found different ways
                    
                    #find the number of '|' in the ID 
                    num = recordID.count('|')
                    if num>3:
                        
                    #split the record ID at each | and assign the correct name to indexes of the list 
                        Splits = recordID.split('|')
                        accessNum = Splits[4]
                        description = Splits[5]
                        
                    #cast description into a string to get rid of spaces in it
                        description = str(description)
                        description.strip(' ')
                        
                        #get only the first part of the gene name (other information is not needed)
                        if ',' in description:
                            detail = description.split(',')
                            nameDescript = detail[0]
                        else:
                            detail = description.split('[')
                            nameDescript = detail[0]
                            
                        #if the sequences has a random series of numbers in it beginning with :, get rid of the numbers
                        if nameDescript.startswith(':'):
                            newNameDescript = nameDescript.split(' ',2)
                            nameDescript = newNameDescript[2]
                            
                        #if the name of the gene has hypothetical in it, get rid of the hypothetical and replace it with
                        #the type of toxin/antitoxin it is
                        if 'hypothetical protein' in nameDescript:
                            nameDescript = filename
                        if 'Hypothetical protein' in nameDescript:
                            nameDescript = filename
                            
                        #build the new gene title to add to the database
                        Sequence = '>'+accessNum+' ~~~'+filename+'~~~'+nameDescript+'_:_'+filename+'-TAfinderdb~~~\n'+str(record.seq)+'\n'
                        OutFile.write(Sequence)
                        
                        
                    #if the sequence did not match any of the formats above/did not have enough information the sequence
                    #will be written to a seperate file(written as originally presented, not reformatted)
                    else:
                        Sequence = '>'+str(record.description)+'\n'+str(record.seq)+'\n'
                        notReadable.write(Sequence)
    
    print('Finished')
##################################################################################################################
'''-5-INPUT: directory that contains files with the integrase protein sequences you want to reformat and add to database.
    OUTPUT:function will append a file, adding all genes with a reformatted title. User should input the file/database
   they are creating (should be the same file/database for all database writing functions in this file) 
   
Function reads in files from a directory that contain genomic sequences (in this case the sequences should be 
integrase genes) and reformats their title to fit into our own database which is compatible with PROKKA(annotation tool).
Funciton will reformat the gene title of the Integrases to include the accession number, class of integrase, 
and integrase label.
EX:>YP_003675754.1 ~~~Int1~~~class1_:_Integrase~~~'''
def tagIntGeneName(directory):
    #open a file to write to, all files ending in .fsa will be written to this file(here I edited our currentdatabase),
    #again if you are creating a new database the file here should be the same name as used in the other functions
    OutFile = open('Summer19CompleteDatabase.faa','a')
    #NOTE: some of the sequences in the IVIP database did not match any uniform format
    #or they were missing some information, any sequences that were this way I had written to a seperate file
    #then I went back and manually fixed them/got needed information. There should only be around 5-6 sequences that
    #this happens with
    OutFileMisFits = open('Summer19INTunmatched.faa','a')
    #find all files ending with .fsa in directory
    for filename in os.listdir(directory):
        #there should only be one file found in this (the integrase) case
        if filename.endswith('.fsa'):
            #create an absolute path to the current filename
            filePath = '/Users/leslieannesmith/SummerResearch19/'+directory+'/'+filename 
            filename = filename.strip('.fsa')
            
            #open file to read as handle so sequences may be parsed out 
            with open(filePath) as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    
                    #Note for following code:
                    #There are several different formats in the IVIP database. Based on specific formats/things found in the 
                    #description of each sequence, the information needed to create a new title and write the sequence
                    #to our database is found different ways:
                    
                    if '|' in record.id:
                        #cast record ID into a string so can split it up by '|' and access needed information 
                        #at different indexes in the list
                        RecStr = str(record.id)
                        Splits = RecStr.split('|')
                        accessNum = Splits[3]
                        Sequence = '>'+str(accessNum)+' ~~~IntI1~~~class1_:_Integrase-IVIPdb~~~\n'+str(record.seq)+'\n'
                        #write sequence to new file (OutFile)
                        OutFile.write(Sequence)
                        
                    elif 'Gene_ID' in record.id:
                        #cast record into a string and write it to the seperate misfit file. 
                        #these genes give Gene_ID and accession num had to be manually looked up
                        OutFileMisFits.write(str(record.id)+'\n'+str(record.seq)+'\n')
                        
                    elif 'GCA' in record.id:
                        #cast record into a string to split it by '_' and access needed information 
                        #at different indexes in the list
                        RecStr = str(record.id)
                        Splits = RecStr.split('_')
                        accessNum = Splits[2]
                        #split the index again to get solely the necessary accession number
                        accessSplit = accessNum.split('.')
                        newAccessNum = accessSplit[0]
                        #write the new record title and sequence to the OutFile
                        Sequence = '>'+str(newAccessNum)+' ~~~IntI1~~~class1_:_Integrase~~~\n'+str(record.seq)+'\n'
                        OutFile.write(Sequence)
                     
                    else:
                        #cast record into a string to split it by '_' and access needed information 
                        #at different indexed in the list 
                        RecStr = str(record.id)
                        Splits = RecStr.split('_')
                        accessNum = Splits[4]
                        #write the new record title and sequence to the OutFile
                        Sequence = '>'+str(accessNum)+' ~~~IntI1~~~class1_:_Integrase~~~\n'+str(record.seq)+'\n'
                        OutFile.write(Sequence)

    print('Finished')
##################################################################################################################
tagResGeneName(translate('resfinder_db'))#functions 1&3

tagToxinGeneName('ToxinAntitoxin')#function 4

tagIntGeneName('dataBaseTester')#funciton 5