#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 13:10:48 2019

@author: ghovis
"""

##############################ProkkaReformat.py##############################
#  Finds genes in Prokka output annotated from a user defined database, noted by annotation name_db.faa
#  Input csv output from Prokka (.tbl file)
#  Output csv with the start index, stop index, ARG name, ARG type, Backbone Gene Name, Incompatibility Group, plasmid name
 
import sys
import re
import csv
import copy
import os
from collections import defaultdict
from Bio import SeqIO


def rowsFromProkka(inFile, dbDelim = ".faa"):
	### Note that infile is a .csv file of output from Prokka
    	# Input 
	#	inFile = csv file for reading Prokka output
	# 	dbDelim = the key indicating user defined database in the Prokka output
	# Output
	#	outlist = a list of each CDS's data formatted as
	# 	(start, stop, arg, arg type, bg, inc croup, plasmidName, keep (Boolean indicating if it was from a passed in database)
	with open(inFile, "r") as csvfile:

		csvreader = csv.reader(csvfile, delimiter= '\t')

		#create an empty dictionary to temporarily hold the information for one gene	
		tempout = {"start" : "NA", "stop" : "NA", "geneName" : "NA", "plasmidName" : "NA", 'keep' : False}
		outlist = []

		for line in csvreader: 
			# identify lines specifying the start of a new plasmid
			if len(line) == 1 and line[0].find(">")==0:
				# identify when the plasmids are switching and write out the previous record
				outlist.append([tempout[I] for I in tempout])
				# clear tempout
				tempout = {"start" : "NA", "stop" : "NA", "resName" : "NA", "resType" : "NA", "geneName" : "NA","incGroup" : "NA",\
				 "plasmidName" : line[0][line[0].index(" "):], 'keep' : False}

			# identify rows starting a new CDS
			elif len(line) == 3 and line[2] == 'CDS':
				outlist.append([tempout[I] for I in tempout])
				# replace values in tempout, note that we only need to 
				# replace the plasmid name when we get to a new plasmid
				tempout["start"] = line[0]
				tempout["stop"] = line[1]
				tempout["geneName"] = "NA" 
				tempout["resName"] = "NA" 
				tempout["resType"] = "NA"
				tempout["incGroup"] = "NA" 
				tempout['keep'] = False
		
			# identify if the CDS is from database and tag for keeping
			elif len(line)>4 and line[3] == 'inference' and ((line[4].find(dbDelim)>-1) or (line[4].find('ISfinder'))):
				tempout['keep'] = True
		
			# identify the gene name if it has one
			elif len(line)>4 and line[3] == 'product':
				if "group" not in line[4] and "RES" not in line[4]:
					tempout['geneName'] = line[4]
					tempout['incGroup'] = "NA"
					tempout['resName'] = "NA"
					tempout["resType"] = "NA"                    
				if "group" in line[4]:
					tempout['incGroup'] = line[4]
					tempout['geneName'] = "NA"
					tempout['resName'] = "NA"
					tempout["resType"] = "NA"                    
				if "RES" in line[4]:
					tempout['resName'] = line[4]
					tempout['geneName'] = "NA"
					tempout['incGroup'] = "NA"
					tempout["resType"] = "NA"
				if "_:_" in line[4]:
					resToList = line[4].split("_:_")
					tempout['resName'] = resToList[0]
					tempout['geneName'] = "NA"
					tempout['incGroup'] = "NA"
					tempout['resType'] = resToList[1]					                    
	return outlist
def cleanList(CDSlist):
	### function removes rows that were not from the current database
	# Input
	#	CDSlist = a list of the information contained for each CDS in one row. 
	#		Last element determines if it was from user defined database
	# Out
	#	returns a list of only the rows from CDS's from our database, 
	#	removes column identifying keeper rows
	return [row[:-1] for row in CDSlist if row[-1]]

#########################findNeighbors2.py#######################
#takes gene table and finds the nearest backbone genes up and dwnstream of every ARG
#input file:file created by prokkaReformat.py
#output file is .csv of ARG (with ARG type), upstream BG, downstream BG, inc group, and plasmid name 
def findNeighbors2(filePath):
    with open(filePath,'r') as  csvreadfile:
        csvreader = csv.reader(csvreadfile)
        tempN1 = ""
        N1 = ""
        N2 = ""
        outList = []
        resList = []
        afterRES = False 
        for line in csvreader:
        	#set gene in BG column to upstream neighbor temporarily
            if not line[4] == "NA" and afterRES == False:
                tempN1 = line[4]
            #if theres an ARG, set temp upstream neighbor to neighbor1 and assign ARG
            if not line[2] == "NA":
                N1 = tempN1
                resName = line[2]
                resType = line[3]#inserted by Leslie
                afterRES = True 
                #resLIst accounts for the possibility of multiple ARGS in row
                resList.append([resName, resType])#resType inserted by Leslie 
                
            #first BG that occurs after the res gene is neighbor 2 (downstream)
            if not line[4] == "NA" and afterRES == True:
                N2 = line[4]
                afterRES = False
                #every ARG in a cluster will share the same up/downstream BG neighbors
                for gene in resList:
                    outList.append([gene,N1,N2,line[5],line[6]])#resType inserted by Leslie 
                tempN1 = line[4]
                resList = []
        #print(outList)
    return outList

######################geneTableEditor.py#############################
#function takes in a parameter, the parameter will be an absolute file pathway
# to the gene table file (should be .csv file from prokkaReformat.py)
#(ex:/Users/marielelensink/Documents/H2/GeneTable1Edit.csv)
#output is  another gene table that is compatible with R
#Gabrielle altered to add incompatibility group 
#leslie just changes all the numbers up one basically
def geneTableEditor(filePath):
    with open(filePath, 'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)
        incList= []
        #Create default dictionary for incompatibility (inc) groups 
        dDict= defaultdict(list)
        outList = []
        plasmid = ""
        for line in csvreader:
            plasmid = line[6]
            plasmid = plasmid.strip()
            if not line[5] == 'NA':
                incOld= line[5]
                #Splice inc group name to make it easier to read
                incSave= incOld[5:10]
                #Add plasmid and inc group pairs to list
                incList.append([plasmid, incSave])
        #Add all the inc groups to each plasmid key in the default dictionary
        for plasmid, incSave in incList:
            dDict[plasmid].append(incSave)
    #Write the data for each plasmid to a new file
    with open(filePath, 'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)  
        resName1 = ""
        for line in csvreader:
            start= line[0]
            stop = line[1]
            if "RES" in line[2]:
                resName1 = line[2]
                resName = resName1[3:]
                resType = line[3]
            else:
                resName = line[2]
                resType = line[3]
            geneName = line[4]
            plasmid = line[6]
            parseinc = line[5]
            parseinc = parseinc[5:10]
            countplasmid = 0
            #Ensure that all the inc groups will be listed for each plasmid 
            for i in dDict.keys():
                key= str(plasmid)
                key= key.strip()
                if (key == str(i)):
                    countplasmid += 1
                    curincs= dDict[key]
                    lencurincs= len(curincs)
                    count= 0
                    if (lencurincs > 1):
                        for j in range(0,lencurincs):
                            oneinc= str(curincs[j])
                            if (oneinc == parseinc):
                                count += 1
                        if (count == 0):
                            incgroup= ",".join(dDict[key])
                            break
                    else:
                        if (str(curincs) != str(parseinc)):
                            incgroup= ",".join(dDict[key])
                            break
            if (countplasmid == 0):
                incgroup= "other"
            outList.append([start,stop,resName,resType,geneName,incgroup,plasmid])
    return outList    

############################ResGeneEdit.py##########################
#takes output.csv file from GeneTableEdit.py and edits to make more compatible with R
#changes instances of "NA" to "NoGene", cleans ARG names 
#output .csv file of the same format 
#Leslie added ResType and moved most index numbers up 1 
def resGeneEdit(inputFile):
    with open(inputFile,'r') as csvreadfile:
        csvreader = csv.reader(csvreadfile)
        start = ""
        stop = ""
        ResGene1 = ""
        ResGene = ""
        ResType = ""
        BackboneGene = ""
        IncGroup = ""
        PlasmidName = ""
        outlist = []
        for line in csvreader:
            start = line[0]
            stop = line[1]
            #cleans/standardizes ARGS by removing subgroups 
            ResGene1 = line[2]
            ResType = line[3]
            #abrvName = /^[A-Za-z]#currently working on this line
            ResGene = ResGene1[:4]
            #change "NA's" to NoGene in ARG column 
            if ResGene == "NA":
                ResGene = "NoGene"
            if ResType == "NA":
                ResType == "NoGene"
            BackboneGene = line[4]
            #change "NA"s to NoGene in BG column 
            if BackboneGene == "NA":
                BackboneGene = "NoGene"
            IncGroup = line[5]
            PlasmidName = line[6]
            outlist.append([start,stop,ResGene,ResType,BackboneGene,IncGroup,PlasmidName])
    return outlist



###################step1(ProkkaReformat.py)#########################
#INPUT: Prokka will output 10 files, this function takes in the file from Prokka ending in .tbl.
#OUTPUT: Function will output a .csv file with the start index, stop index, ARG name, Backbone Gene Name, 

# Incompatibility Group, and plasmid name
cleanedUp = (cleanList(rowsFromProkka("PROKKA_06062019.tbl")))#insert name of input file (from prokka) here 
fName = "6-10resTypeTestNine"

with open(fName+'TableOutput.csv','w') as csvfile: #name of output file here
	
	csvwriter = csv.writer(csvfile)
	csvwriter.writerows(cleanedUp)		
#################step2(findNeighbors2.py)##########################
newNeighbors = findNeighbors2('/Users/ghovis/Documents/Research2019/Databases/'+fName+'TableOutput.csv')
with open(fName+'Neighbors.csv','w') as csvwritefile:
    csvwriter = csv.writer(csvwritefile)
    csvwriter.writerows(newNeighbors)
#################step3(geneTableEditor.py)#########################

editedTable = geneTableEditor("/Users/ghovis/Documents/Research2019/Databases/"+fName+"TableOutput.csv")#name of input file from prokkareformat      
        
with open(fName+'EditedGeneTable.csv','w') as csvwritefile:#need to change the input here 
        csvwriter = csv.writer(csvwritefile)
        csvwriter.writerows(editedTable)
######################step4(ResGeneEdit.py)########################
geneOutputList = resGeneEdit(fName+'EditedGeneTable.csv')
with open(fName+'RCompatibleTable.csv','w') as csvwritefile: #write in output file name 
    csvwriter = csv.writer(csvwritefile)
    csvwriter.writerows(geneOutputList)