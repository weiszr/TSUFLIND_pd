'''
Created on 2013-5-24
===========================================The Root Error Caculate Part==================================
#===============================================Input Data===============================================
#==========================================The calculate part============================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=====================================The Root Error Calculate part model================================
#========================================================================================================
#========================================================================================================
#This part is designed to calculate the root square error between field data and model result, if the RSE
# is larger than 10, change the grainsize distribution of sediment source, if less than 10, stop go to nextt
# step directly.
#This part is based on the paper published on Sedimentary Geology 282(2012) 90-109
# Flow Speed estimated by inverse modeling of sandy tsunami deposits: results from the 11 March tsunami on 
# the coastal plain near the Sendai Airport, Honshu, Japan
#               Bruce E. Jaffe and Kazuhisa Goto
#Hui Tang 05/24/2013
#Input:		phi-grain size
#           Fr2-Model result Distribution 
#           Num-The number of Fraction of sub-intervals in each sample location
#Output:	RSE-Average Root Sqaure Root
#=======================================================================================================
#==========================================Input Data===================================================
#This part is designed to input data to calculate RSE

from pylab import *
import os
from function import *
from ReadCSV import *
def calrse(name,phi,Fr2,Num):
    sumerr=zeros(Num)#RSE for each sub-intervals 
    for i in range(Num):
        data1=readCSV(name,separator=';')
    	phi=data1[:,0]
        Fr1=data1[:,1]
#=============================================The calculate part========================================
        n=len(Fr1)
        m=len(Fr2)
        err=zeros(n)
        sumerr=zeros(n)
        if(m!=n):
            print('These Dimension of two document are not same, Please check the data structure')
        else:
            for j in range(n):
                err[j]=((Fr1[j]-Fr2[j]))**2
                sumerr[i]=sum(err)
                sumerr[i]=sqrt(sumerr[i]/n)
    RSE=sum(sumerr)
    return(RSE)

