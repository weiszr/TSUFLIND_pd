'''
Created on 2013-5-21
=============================The original tsunami sediment concentration model===========================
#=========================================Input parameters===============================================
#=====================================The calculate model part===========================================
##==============================Caculate critical velocity and critical shear stress=====================
##=================================Caculate the concentration of each grain size=========================
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#==============================The original tsunami sediment concentration model=========================
#========================================================================================================
#========================================================================================================
#This model is based on the paper published Proceedings of the 24th International Conferrence on Coastal 
#Engeering, 1993
#SUSPENDED SEDIMENT SUSPENSION ON THE INNER SHELF DURING EXTREME STROMS 
#              Madsen, O.S, Chisholm, T.A and Wright, L.D
#Hui Tang 05/21/2013
#Input:             Se-Grain size in phi
#                   fr-Percentage of each grain size %
#                   Vi-Smallest possible average velocity m/s
#                   H-water depth at still water level m
#                   Cb-Bed concentration M^3/m^3
#                   rho-Density of sea water
#Used constant:     k-Karman's constant 0.41
#                   rho-Density of fluid: 1050 kg/m^3
#                   rhos-Density of sediments: 2770 kg/m^3
#                   g0-Resuspension coefficient 0.0004
#Function used:     UstarCrit
#Output:            Se-Grain size in phi
#                   C0-initial concentration for sediments source m^3/m^3
#========================================================================================================
#==========================================Input parameters==============================================
#This part is designed to input the parameter for this model 
from pylab import *
import os
from output2CSV import *
from function import *
def inconcentration(se1,fr,Vi,H,Cb,rho):
    with open('parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                if n=='resuspension_coefficient': g0=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
		    	continue
#========================================================================================================
#======================Deal with the data from original grainsize distribution===========================
#This part is designed to deal with the data from orignal grainsize distribution file, in order to use in 
# caculating concentration.
    se=2**(-se1)*0.001
    nc=len(se)
#========================================================================================================
#=====================================The calculate model part===========================================
#This part is designed for calculate the original sediment concentration base on the grain size
##==============================Caculate critical velocity and critical shear stress=====================
    Ucrit=zeros(nc)#Critical shear velocity for initiation of sediment transport
    tcrit=zeros(nc)#Critical shear stress (N/m^2)
    C=zeros(nc)# Original concentration for each grain size
    for i in range(nc):
        Ucrit[i]=UstarCrit(se[i],KinVisc(wtemp),Rhos/rho)# caculate critical shear velocity for initiation of sediment transport
        tcrit[i]=rho*Ucrit[i]**2# Caculate critical shear stress (N/m^2)
    taub=(Vi*vk/(log(30*H/(2**(-1)-2**(-9))-(1-(H*(2**(-1)-2**(-9)))))))**2*rho # bottom stress(?)
    taustar=taub/tcrit
    S=taustar-1     #normalized excess shear stress
    for i in range(nc):#returns zero when taub<=tcrit
        if(S[i]<0):
            S[i]=0
    for i in range(len(S)):
        C[i]=fr[i]*g0*Cb*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        C0=Rhos*C
    return(C0)
