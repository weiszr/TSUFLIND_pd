'''
Created on 2013-5-24
========================The Root Error Caculate Part for sediments thickness=============================
#===============================================Input Data===============================================
#==========================================The calculate part============================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=====================================The Root Error Calculate part model================================
#========================================================================================================
#========================================================================================================
#This part is designed to caculate the sediment difference between field data and model result, if the 
# difference between two data is larger than 5%, change the average velocity, if less than 5%, stop go to 
# next step directly.
#Hui Tang 05/24/2013
#Input:         name-Document name of field data document
#               th2-Model result of suspended sediemnts thickness m
#               N-The number of sample location 
#Output:        averthRSE-Average thickness difference between field data and model result 
#=======================================================================================================
#==========================================Input Data===================================================
#This part is designed to input data to calculate thickness difference
from pylab import *
import os
from ReadCSV import *
from function import *
def thcalrse(name,th2,N):
    rse=zeros(N)#Average thickness difference between field data and model result for each sample location 
    data1=readCSV(name,separator=';')
    th1=data1[:,1]
    k=0
    n=len(th1)
    m=len(th2)
    if(n==m):
        for j in range(N):
            if (th1[j]!=0):
                if(th2[j]!=0):
                    rse[j]=(th1[j]-th2[j])*100/th1[j]
                    k=k+1
    else:
        print("These Dimension of two document are not same, Please check the data structure")
    averthRSE=sum(rse)/N
    #print th1
    #print th2
    #print rse
    return(averthRSE)
