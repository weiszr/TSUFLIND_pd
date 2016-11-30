'''
Created on 2013-5-26;
=====================Joint Tsunami Sediment Speed model pre-processing Script=====================
#==================================Grainsize Distribution Part====================================
#==========================================Water Density Part=====================================
#==========================================Water Depth Part=======================================
@author: Hui Tang, tanghui@vt.edu
'''
#=================================================================================================
#=================================================================================================
#====================Joint Tsunami Sediment Speed Model Pre-processing Script=====================
#=================================================================================================
#=================================================================================================
#This function is designed to preprocess the field data to input the joint model.
#Inputs:            filename-The name of field data document, must be csv documnet contain grainsize
#                            and percentage of each grainsize,separated by ";" and in a row.
#Function used:     sedstat.py,sw_dens0.py
#Output:            Dl: Largest grainsize in phi
#                   Ds: Smallest grainsize in phi
#                   Dm: Medean grainsize in phi
#                   Nc: The number of grainsize class
#                   Rho:Density of sea water
#                   H:Max water depth
from pylab import *
from function import *
from sedstats import *
from ReadCSV import *
from output2CSV import *
import csv
import os
def preproc(filename):
#===================================================================================================
#===========================================Input Parameter=========================================
#===================================================================================================
#This part is designed for inputing all parameter for all model and major function
    separator=':'
    with open('parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
	        if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='water_temperature': wtemp=float(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_run_up':Rz=float(v)
                elif n=='slope': m=float(v)
                elif n=='Largest_distance': l=float(v)
                elif n=='Depth': h=float(v)
                elif n=='Sediment_thickness':th=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
#=========================================Grainsize Distribution Part==============================
    filename1='data_process.csv'
    data=readCSV(filename,separator=';')
    phi=data[:,0]
    fr=data[:,1]
    Dl=0
    Ds=data[-1,0]
    Nc=len(data[:,0])
    fr=fr.reshape(len(fr),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    if(len(fr)==1):
            Dm=phi[0]
    else:
            sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
            Dm=sedstat[0]
#==========================================Water Density Part=====================================
    fr=fr.reshape(len(fr))
    rho=sw_dens0(sal,wtemp)
#==========================================Water Depth Part=======================================
    H=h*m*Rz/(m*Rz-l)
    fp=output2CSV1(filename1,phi,100*fr,th,h,l,'phi','fraction','thickness','depth','location')
    return(Dl,Ds,Dm,Nc,rho,H)

