'''
Created on 2013-1-13
======================================The tsunami sediment model=========================================
#=====================================The calculate model part===========================================
#=============================================Input parameters===========================================
##==============================Conculate the Deposit Velosity===========================================
##================Conculate some different parameters for tsunami sediment model ========================
##=====================================Calculate the thickness at x=0====================================
##==================Calculate the horizontal and vertical sediment run-up along the slope================
##==Calculate different grain size's thickness,the number of grains and area everywhere along the slope==
##=============================calculate percentage of different grian size==============================
#=========================================The plot figure part===========================================
##=================================plot the location of sediments end====================================
##===================plot the different grain size's thickness on the still water level==================
##========================plot different grains thickness on different location==========================
##===================================plot pecentation of different grian size============================
##===================================plot D50 vs different location======================================
##=======================plot thickness vs different location in different grain size====================
##=======================plot total thickness vs different location in different grain size==============
@ author: Hui Tang- tanghui@vt.edu
'''
#========================================================================================================
#========================================================================================================
#===========================The tsunami sediment model based on R L. Soulsby=============================
#========================================================================================================
#========================================================================================================
#This model is based on the paper published in sixth International Symposium on Coastal Sediment Processes
#RECONSTRUCTING TSUNAMI RUN-UP FROM SEDIMENTARY CHARACTERISTICS- A SIMPLE MATHMATICAL MODEL
#             Richard L. Soulsby, David E. Smith, Alan Ruffman
#Input:             V-Average velocity m/s
#                   C0-Initial sediment concentration kg/m^3
#                   se-Sediment grainsize in phi
#                   H-Water depth at still water level
#                   nc-Number of grain size classes 
#Constant used:     Rz-Vertical water run-up m
#                   m-slope
#                   Rhos-Sediments density
#                   wtemp-Water tempertature C
#                   sal-Salinity of sea water (psu)
#                   q-sediment deposit rate
#                   nc-Number of grainsize classes
#                   N-Number of sample location 
#Fuction used:      Tubesetlevel
#                   inconcerntration
#Output:            th2-Sediment load for each sample location 
#                   se1-grain size in phi
#                   fr-Grain size distribution for each sample point
#========================================================================================================
#===========================================The calculate model part=====================================
#This part is designed for calculate the sediment grains deposit velosity, horizontal run-up limits of water,
#horizontal and vertical run-up limits of tsunami sediments, sediment thickness, the number of sediment grains
#and percentage of different grain size in different location along the slope
#========================================================================================================
#=============================================Input parameters===========================================
#========================================================================================================
#This part is designed to input the parameter for this model
from pylab import *
import os
from output2CSV import *
from function import *
from ReadCSV import *
def sousby(se,C0,V,H,nc):
    with open('parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                if n=='sediment_deposit_rate': q=float(v)
                elif n=='water_run_up':Rz=float(v)
                elif n=='slope': m=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='Number_sample': N=int(v)
                elif n=='Filename': filename=str(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='salinity': sal=float(v)
                elif n=='pi': pi=float(v)
                elif n=='Filename_Sousby': filename=str(v)
                elif n=='settling_velocity_type': setType=float(v)
            except IndexError as e:
                continue
            except ValueError as e:
			continue
    C=sum(C0)
    #print C
    se2=phi2mm(se)*0.001
#===================================Calculate the Deposit Velosity=======================================
#Calculating the deposit velosity for different grain size and set different grian size for plot figure in
#the phi scale
    rhow=sw_dens0(sal,wtemp)
    if(setType==1.0):
        ws=tubesetvel(se,Rhos,wtemp,sal)#Settling velocity of each grain size
    elif(setType==2.0):
        ws=100*dietrichWs(wtemp,se2,Rhos,rhow,C/2650.0,0.7,3.5)
    elif(setType==3.0):
        ws[i]=Setvel1(se2,kinvisc,s[i])
    #ws=tubesetvel(se,Rhos,wtemp,sal)#Settling velocity of each grain size 
#====================Conculate some different parameters for tsunami sediment model =====================
    Rw=Rz*m #The horizontal run-up limit of tsunami water
    T=Rw/V#run up time
    a=ws*T/H #constant a in model 
    Rs=Rw/(1+a*q) #The horizontal run-up limit of tsunami sediment
    #print "Rs"
    #print Rs
    x=linspace(0,350,int(N))
#======================================Calculate the sediment load at x=0=================================
#This part is designed for calculate different grain size's thickness at still water level
    th=zeros(nc)# The sediment load of different grain size sediments in still water level
    th=a*(1+a*q)*C0*H/((1+a)*Rhos)
#=============================Calculate different grain size's sediment load along the slope==================
#This part is designed to calculate sediment load, the area of the sample location, the area of sediment
#will be used to calculate the percentage of different grain size
    i=0 #index for grain size
    j=0 #index for sample location 
    th2=zeros(shape=(N,nc)) #The sediment thickness in different sample location
    A=zeros(shape=(N,nc)) #The area of sediment grain in different sample location
    sumA=zeros(N) # The total area of sediment grain in different sample location  
    totalth2=zeros(N)# The total sediment thickness in different sample location
    Dt=0.01 #The width of area
    while(j<N): #calculate the thickness, number and area
        i=0
        while(i<nc):
            th2[j,i]=th[i]*(1-x[j]/Rs[i])
            if(th2[j,i]>0):
                A[j,i]=th2[j,i]*Dt
                sumA[j]=sumA[j]+A[j,i]
                totalth2[j]=sumA[j]/Dt
            i=i+1
        j=j+1
#================================calculate percentage of different grian size=================================
#This part is design to use the area to calculate the percentage of different grain size in different sample
#location
    #print th2[9]
    p=0 #index for sample location
    pecentA=zeros(shape=(N,nc)) #The percentage of different grain size
    med=zeros(N)# The D50 of sample location
    while(p<N): # calculate the percentage
        n=0
        sumpecent=0
        while(n<nc):
            pecentA[p,n]=100*A[p,n]/sumA[p]
            sumpecent=sumpecent+pecentA[p,n]
            if(sumpecent<=55):
                med[p]=se[n]
            n=n+1
        p=p+1
#===========================================================================================================
#==================================The plot figure part=====================================================
#This part is designed to plot figure for this model include:total sediment load, mean grain size,flow depth
##==================================Output Data and parameters==========================================
#This part is designed to output and plot Data and parameter
    h=x*(Rz-H)/(Rz*m)-x/m+H# Water depth
    for i in range(N):
        pe=pecentA[i,:]
        filename2="sample%02d"%i+".csv"
        fp2=output2CSV1(filename2,se,pecentA[i,:],totalth2[i],h[i],x[i],'phi','fraction','thickness','depth','location')
        f4=figure()
        filename="Grainsize distribution at different location x=%f"%x[i]+"m"
        ax4=subplot(111)
        ax4.plot(se,pecentA[i,:],'-')
        ax4.set_xlim(0,10)
        ax4.set_ylim(0,40)
        ax4.set_xlabel('grain size($\phi$)')
        ax4.set_ylabel('The percentage of different sediment grains size(%)')
        #savefig(filename+'.png',dpi=100)
        close()
    return(totalth2)
