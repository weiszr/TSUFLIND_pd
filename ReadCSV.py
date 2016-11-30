'''
Created on 2013-5-24;
===============================================Read data From CSV document==================================
#=================================================readCSV function one======================================
#=================================================readCSV function two======================================
@author: Hui Tang tanghui@vt.edu
'''
from pylab import *
from types import *
import csv
#============================================================================================================
#============================================================================================================
#===========================================Read Data From CSV document Script===============================
#============================================================================================================
#============================================================================================================
#This Script is used to read data from csv document 
#Writen by Hui Tang 05/24/2013
#==================================================readCSV one===============================================
#This is a function to read data from CSV document,two items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction)
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu
def readCSV(name,separator=';'):
    data=zeros([1,2]); start=True # initialize data array
    with open(name,'rU') as f:
        reader = csv.reader(f) # read CSV file
        for row in reader: # read each line and store them into a n by 2 array
            try: rs=row[0].split(separator); d1=float(rs[0]); d2=float(rs[1]);
            except IndexError: continue # if line is not delimited by separator, continue
            except ValueError: continue # if not numbers, continue
            if (d1!="\x00"):
                if not start: data=append(data,[[d1,d2]],axis=0); 
                else: start=False; data[0]=[d1,d2]
    return data
#===============================================readCSV two===================================================
#This is a function to read data from CSV document,five items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction),sediment thickness, water
#                        depth and location 
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu
def readCSV1(name,separator=';'):
    data=zeros([1,5]); start=True # initialize data array
    with open(name,'rU') as f:
        reader = csv.reader(f) # read CSV file
        for row in reader: # read each line and store them into a n by 2 array
            try: rs=row[0].split(separator); d1=float(rs[0]); d2=float(rs[1]);d3=float(rs[2]);d4=float(rs[3]);d5=float(rs[4])
            except IndexError: continue # if line is not delimited by separator, continue
            except ValueError: continue # if not numbers, continue
            if (d1!="\x00"):
                if not start: data=append(data,[[d1,d2,d3,d4,d5]],axis=0); 
                else: start=False; data[0]=[d1,d2,d3,d4,d5]
    return data
#===============================================readCSV three===================================================
#This is a function to read data from CSV document,five items
#Input:             filename-The data file's name
#                   separator-The separator between different types of data
#Output:            data-grain size class (phi) and weight (actual, %, or fraction),sediment thickness, water
#                        depth and location 
#writed by Hui Tang,Virginia Tech, Jan 31, 2013 tanghui@vt.edu
def readCSV2(name,separator=';'):
    data=zeros([1,3]); start=True # initialize data array
    with open(name,'rU') as f:
        reader = csv.reader(f) # read CSV file
        for row in reader: # read each line and store them into a n by 2 array
            try: rs=row[0].split(separator); d1=float(rs[0]); d2=float(rs[1]);d3=float(rs[2]);
            except IndexError: continue # if line is not delimited by separator, continue
            except ValueError: continue # if not numbers, continue
            if (d1!="\x00"):
                if not start: data=append(data,[[d1,d2,d3]],axis=0); 
                else: start=False; data[0]=[d1,d2,d3]
    return data
