'''
Created on 2013-5-21
=====================================The source grain size distribution===================================
=====================================The calculate model part=============================================
#=============================================Input parameters============================================
#=====================================Coculate the pecerntage of each grain size==========================
#=====================plot original grain size distibution================================================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#=======The original tsunami sediment grain size distibution model based normalize distribution==========
#========================================================================================================
#========================================================================================================
#This model is based on idealize grainsize distribution which will be used for an input for joint model,
#it will be set as normalize distribution and can be other dstribution based on the problem.
#Input:             Dl-Largest grain diameter in phi
#                   Dm-Medean grain diameter in phi
#                   Ds-Smallest grain diameter in phi
#                   nc-Number of grain size class
#Output:            Se-Grain size in phi
#                   fr-Percentage of each grain size %
# Hui Tang 2013-5-21
#========================================================================================================
#=============================================Input parameters===========================================
#========================================================================================================
from pylab import *
import os
from output2CSV import *
from function import *
from ReadCSV import *
def source_distribution(Dl,Ds,Dm,nc):
#========================================================================================================
#=====================================Coculate the pecerntage of each grain size=========================
#========================================================================================================
#Conculating the pencertage for different grain size and set different grian size for plot figure in the 
#phi scale
    se=linspace(Dl,Ds,nc)
    #Dm=2.8
    filename1='sample_P14abc.csv'
    Field_data=readCSV(filename1,separator=';')
    Field_se=Field_data[:,0] 
    fr=exp(-(Field_se-Dm)**2/2*1)/sqrt(2*(3.14))
    fr=fr/sum(fr)
#========================================================================================================
#=====================plot original grain size distibution===============================================
#========================================================================================================
    #fname1='original grain size distibution.csv'
    #fp2=output2CSV(fname1,Field_se,fr,'phi','fraction')
    #f=figure()
    #fname='original grain size distibution'
    #ax=subplot(111)
    #ax.set_xlim(0,7)
    #ax.set_ylim(0,6)
    #ax.set_title(fname)
    #ax.plot(Field_se,100*fr,'-')
    #ax.set_xlabel('Grain size')
    #ax.set_ylabel('Pecentage (%)')
    #grid()
    #show()
    return(Field_se,fr)
def source_distribution1(Dl,Ds,Dm,nc):
#========================================================================================================
#=====================================Coculate the pecerntage of each grain size=========================
#========================================================================================================
#Conculating the pencertage for different grain size and set different grian size for plot figure in the 
#phi scale
    se=linspace(Dl,Ds,nc)
    #Dm=2.8
    filename1='sample_P14abc.csv'
    Field_data=readCSV(filename1,separator=';')
    Field_se=Field_data[:,0] 
    fr=exp(-(Field_se-Dm)**2/2*1)/sqrt(2*(3.14))
    fr=fr/sum(fr)
#========================================================================================================
#=====================plot original grain size distibution===============================================
#========================================================================================================
    #fname1='original grain size distibution.csv'
    #fp2=output2CSV(fname1,Field_se,fr,'phi','fraction')
    #f=figure()
    #fname='original grain size distibution'
    #ax=subplot(111)
    #ax.set_xlim(0,7)
    #ax.set_ylim(0,6)
    #ax.set_title(fname)
    #ax.plot(Field_se,100*fr,'-')
    #ax.set_xlabel('Grain size')
    #ax.set_ylabel('Pecentage (%)')
    #grid()
    #show()
    return(se,fr)
