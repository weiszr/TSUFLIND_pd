'''
Created on 2012-12-10
=========================The tsunami water depth and average velocity model==========================
#==================================The calculaet model part==========================================
##===================================handle input parameter==========================================
##================================preprocess the data for equation===================================
#=================================The solve equation part============================================
@author:Hui Tang-tanghui@vt.edu
'''
#====================================================================================================
#====================================================================================================
#============The tsunami water depth and average velocity model based on Andrew L. Moore=============
#====================================================================================================
#====================================================================================================
#This model ia bsed on a paper from Sedimentary Geology 200(2007)336_346
#LANDWARD FINING FROM MULTIPY SOURCES IN A SAND SHEET DEPOSITED BY THE 1929 GRAND BANKS TSUNAMI
#NEWFOUNDLAND
#Andrew L. Moore, Brian G. McAdoo,Alan Ruffman
#Input:         Dl-Largest grain's diameter m
#               Dm-medean grain's diameter m
#               rho-water density
#Used constant: tcr1-Dimensionless critical shear stress 0.045-0.06
#               g-gravity acceleration: 9.8m/s^2
#               rho-Density of fluid: 1050 kg/m^3
#               rhos-Density of sediments: 2770 kg/m^3
#               k-Karmon's constant: 0.41
#               R-Rouse number: 0.55<R<0,63 for fine peak; 1.8<R<2.1 for coaser peak
#Output:        Vi-Smallest possible average velocity m/s
#Hui Tang 05/22/2013
#====================================================================================================
#=================================The calculate model part===========================================
#This part is designed for calculatethe critical shear stress,shear velocity,particle velocity and set
#these equations for this model
from pylab import *
import os
def Moore_function(Dl,Dm,rho):
#========================================= Handle input parameters===================================
#This part is used to read parameters for this model frome document
    separator=':'
    with open('parameter_P16.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
	        if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='Dimensionless_shear_velovity': Tstarcr=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Rouse_number': P=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='Largest_distance': l=float(v)
                elif n=='Gravity_acceleration': g=float(v)
            except IndexError as e:
	        continue
	    except ValueError as e:
	        continue
#======================================Preprocess the data for equations===========================
#This function accumulate the critical shear stress,shear velocity and particle settling velocity
#Take the result as these equation's parameters
    def processdata1(tcr1,D,rhos1,rhos2,g):#get the critical shear stress and shear velocity
        if tcr1>0.045:
            if tcr1<0.06:
                tcr2=tcr1*D*g*(rhos1-rhos2)# get the critical shear stress
                U1=sqrt(tcr2/rhos2)# get the shear velocity
            else: 
                    print "please make dimensionless critical shear stress less than 0.06"
        else: print "Please make dimensionless critical shear stress larger than0.045"
        return U1
    def processdata2(R,Vc,u1):#get the particle settling velocity
        if R>0.05:
            if R<0.63:
                ws=R*Vc*u1#get the particle settling velocity
            else:
                    print "Please make Rouse number smaller than 0.63"
        else: print "Please make Rouse number larger than 0.05"
        return ws
    u1=processdata1(Tstarcr,Dl,Rhos,rho,100*g)
    Ws=processdata2(P,vk,u1)
    h=arange(1,10000,1) #Depth of water
    z0=Dm/30 #Roughtness of land surface
    u2=l*Ws/h
    u3=u1*(log(h/z0-(1-z0/h)))/P 
#=================================Solve the equation system==================================
#This part is designed for solving equation and get the tsunami's average depth and velocity
    def fn(h,Ws,P,u1,l,z0):
        s=l*Ws/h-u1*(log(h/z0-(1-z0/h)))/P
        return(s)
    def solve_function(s,h,delta=0.0001):
        s1=s
        if s1>0:
            while (s1>0):
                h=h+delta;
                s1=fn(h,Ws,P,u1,l,z0);
            H=h;
        else:
            while(s1<0):
                h=h-delta;
                s1=fn(h,Ws,P,u1,l,z0);
            H=h;
        return H
    h=5
    s=fn(h,Ws,P,u1,l,z0)
    H=solve_function(s,h,delta=0.001) #The average depth of water
    U1=u1*(log(H/z0-(1-z0/H)))/P #The average velocity of flow
    U2=l*Ws/H
    return(U2)
