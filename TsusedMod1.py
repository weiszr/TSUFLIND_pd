'''
Created on 2013-1-25
========================================The TsuSedMod tsunami sediment transport model==============================
#=============================================The calculate model part==============================================
##==========================================Read data from file for model===========================================
##============================================Calculate water parameter=============================================
###=================================================calculate s=====================================================
###===============================calculate critical shear stress (N/m2) and settling velocity for each=============
###=================================================calculate zo====================================================
###==================================calculate shear velocity needed to suspend deposit=============================
#==============================================subfunction part=====================================================
##============================================Subfunction grading===================================================
###============================================Calculate part=======================================================
###===========================================plot figure part======================================================
####============================plot figure for Grain Statistics, 1st Moment for layers=============================
####===================================plot figure for Standard Deviation for layers================================
####=====================================plot figure for FW mean for layers=========================================
####=====================================plot figure for FW sorting for layers======================================
####=====================================plot figure for FW Skewness for layers=====================================
####======================================plot figure for FW Kurtosis for layers====================================
####==================================plot figure for Thickness vs. Time sediment hits bottom=======================
##=============================================subfunction sedstats=================================================
##=============================================subfunction readModelFile============================================
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
#===================================================================================================================
#===================================================================================================================
#=========================================The tsunami sediment model based on Bruce Jaffe===========================
#===================================================================================================================
#===================================================================================================================
#This model is based on the paper published in Sedimentary Geology
#A SIMPLE MODEL FOR CALCULATING TSUNAMI FLOW SPEED FROM TSUNAMI DEPOSITS
#               Bruce E. Jaffe, Guy Gelfenbuam
#===================================================================================================================
#Fuction:estimates water velocity from tsunami sediment deposit
#Assumes steady uniform flow
#Functions needed: phi2mm.py, sw_dens0.py, sw_smow.py, KinVisc.py, UstarCrit.py, zoWS.py,dietrichWs.py, Setvel1.py,
#                  tubesetvel.py, linearK.py, parabolicK.py, gelfK.py
#Inputs for this model:         th:deposit thickness (m)
#               	        	h:estimate of water depth (m)
#                       		sal:salinity (psu)
#               	        	wtemp:water temperature (deg. C)
#               	        	nclass:number of sed sizes to consider
#                       		Cb:total bed concentration (1-porosity), usually 0.65
#                       		fr:fractional portion of each grain size in bed
#                               phi:grain size in phi units for each size class
#                       		RhoS:sediment density (kg/m^3) for each size class        
#                         		Dmn:mean grain diamter(mm);
#                               ustrc:shear velocity, initial guess input and then allowed to change
#                               ustrcAdjFac:Factor that sets how quickly ustrc is adjusted during iterations
#                               nit:number of iterations to run
#                               diststep:how often to iterate size classes versus concentration
#                               nz:number of vertical bins for calculations
#                        		vk:Van Karmen's constant
#                        		g0:the resuspension coefficient  *** standard value 1.4 x 10^-4 from Hill et al., 1989
#                         		concconvfactor:difference between modeled and observed concentration difference allowed
#                               sizeconvfactor:set the difference between modeled and observed size distribution allowed
#                               var:sets eddy viscosity shape; 1=linear, 2=parbolic,3=[Gelfenbaum]
#Used:                          ssfrwanted:suspended size distribution wanted
#                               ssconverge:keeps track of whether amount of sediment in each size class is acceptable
#                               ssoffa:array to track of how far off the suspended sediment concentration is from desired concentration 
#                               ssfra:array to track of suspended sediment concentration in each size class
#                               fra:array to track of bed sediment concentration in each size class
#                               offa:array to track how far total suspended loasource_distribution.pyd is off from load required to create deposit
#                               ustrca:array to track how  u*c changes with iteration
#Parameter calculated by import function:
#                               d:grain diameter (set to m)
#                         		rhow:density of water (kg/m^3)sedstat
#                       		kinvisc:kinematic viscosity of water (m^2/s), usually about 10^-6
#                               UcritDmn:critical shear velocity for mean grain size
#                               zotot:Z naut total, from Nickaradse bed roughness and saltation bed roughness
#                         		ws:settling velocity (m/s)
#                        		ucrit:critical shear velocity for initiation of sediment transport
#                               z:log space elevation values where concentrations and velocities are calculated
#                          		K:eddie viscosity (m/s^2)
#                       		zo:Z naut, zero velocity intercept
#parameter calculated in this model:
#                               zoN- Nickaradse grain roughness
#      	                        ustrc:current shear velocity, u*c [m/s], initial value 0.1 m/s
#		                        s:ratio of density of sediment to density of fluid density for each size class
#   	                        thload:sediment load needed to get observed tsunami deposit thickness
#		                        tcrit:critical shear stress (N/m^2)
#		                        taub:bottom shear stress due to currents (N/m^2)
#		                        taustar:bottom shear stress normalized by critical shear stress
#		                        S:excess shear stress
#                               Ca:reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters 
#		                        Ki:integral of 1/eddie viscosity
#                               C:Suspended-sediment profile (calculated for each size class)
#                               G:Suspended load (calculated for each size class)
#                               off:normalized difference between suspended load and load required to create deposit
#		                        spd:mean current speed calculated using current shear velocity (m/s)
#                               froude:Froude number calculated using maximum speed
#Options:
#                               'infile':specifies filename of model input file. 
#                               'outfile':writes results to comma-separated (.csv) file. If the specified output file already exists,
#                               results will be appended to current file.  
#                               'grading':plots results of grading function . The default output depth of this function is 0.01 m.  To
#                                modify the defaut, the user can specify an additional arguement of a vector of depths or a single value.
#Usage examples:
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5;       
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('infile','C:\exampleFile.xls');
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',0.01);
#                                   where 0.01 specifies set interval (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('grading',[0.01,0.03,0.04,0.8]); 
#                                   where [ ] specifies variable intervals (m) for grain statistics for the synthetic deposit
#                               [modelP,siteInfo,datain,details,results]=Tsunami_InvVelModel_V3p5('outfile','modelResults.csv');
#Version information:
#written in March 12, 2001 by Bruce Jaffe, based on code written by Chris Sherwood, documented in July 3,2001
#Significant code clean-ups and user friendly input added by Andrew Stevens from March to October, 2007
#modfied on 11/29/07 by Bruce Jaffe, modification of Andrew Stevens' cleaned-up code including: allowing single size class runs
#modified 4/29/08 to output weight percent (Bruce Jaffe) and to fix grading function interval error (Mark Buckley)
#made in python 1/25/2013 by Hui Tang (tanghui@vt.edu)
#==============================================================================================================================================
#==============================================================The caculate model part=========================================================
from pylab import *
from function import *
from sedstats import *
from numpy import *
from grading import *
from ReadCSV import *
import csv
import os
def Tsusedmod1(name):
##============================================================Read data from file for model====================================================
    print("\n ====================================running==================================== \n")
    #read in model inputs from excel file
    separator=':'
    with open('parameter_P14abc.txt','r') as f:
        for line in f:
            try:
                s=line.split('=');n=s[0];v=s[1];#remove '\n'
                if os.name!='nt':v=v[0:len(v)-1] # remove '\r' for Unix-system
                #re-name input variables to match input list
                #model parameters
                if n=='number_of_iterations': nit=float(v)
                elif n=='Von_Karmen_constant': vk=float(v)
                elif n=='number_of_vertical_bins': nz=float(v)
                elif n=='resuspension_coefficient': g0=float(v)
                elif n=='concentration_convergence_factor': concconvfactor=float(v)
                elif n=='size_convergence_factor': sizeconvfactor=float(v)
                elif n=='shear_velocity': ustrc=float(v)
                elif n=='bed_roughness': zotot=float(v)
                elif n=='settling_velocity_type': setType=float(v)
                elif n=='eddy_viscosity_shape': var=float(v)
                #site info
                elif n=='site_description': t=str(v)
                elif n=='salinity': sal=float(v)
                elif n=='water_temperature': wtemp=float(v)
                elif n=='sediment_density': Rhos=float(v)
                elif n=='bed_concentration': Cb=float(v)
                #data
                elif n=='Number_sample':num=int(v)
                elif n=='Number_grainsize':nclass=int(v)
                elif n=='data_dimension':l=int(v)
            except IndexError as e:
                continue
            except ValueError as e:
                continue
    filename=name
    data=readCSV1(filename,separator=';')
    phi=data[:,0]
    fr1=data[:,1]
    fr=zeros(shape=(len(fr1),1))
    th=data[0,2]
    h=data[0,3]
    x=data[0,4]
    #end data input
    nclass=len(fr1) #number of sed sizes to consider
    fr=fr1.reshape(len(fr1),1) #make fr into a column vector
    fr=fr/sum(fr) #normalize size fractions to avoid rounding errors in sed. analysis causing problems
    ssfrwanted=fr #susp sed size distribution wanted, keep this for comparison with model results
    ssconverge=ones(shape=(nclass,1))#initiate ss converge to 1 for all size classes (does  converge)
    for i in range(len(fr)):
        if(fr[i]!=0):
            ssconverge[i]=0 #reset ss converge values to 0 (does not converge) for all size classes with sediments
    sizeconverge=sizeconvfactor*ones(shape=(nclass,1))#set percent difference in size distribution allowed, same value for all sizes
    D=phi2mm(phi) #convert grain size to mm for settling velocity calculation using Dietrich's formulae
    d=0.001*D #convert classes to m for calculating settling velocity using other forumlae and for UstarC calculation
    if(nclass==1):
        dMean_phi=phi[0]
    else:
        sedstat=sedstats(phi,transpose(fr))#calculate mean grain size
        dMean_phi=sedstat[0]
        m3a=sedstat[1]
        m4a=sedstat[2]
        stdeva=sedstat[3]
        fwmedian1=sedstat[4]
        fwsort1=sedstat[5]
        fwmean1=sedstat[6]
        fwskew1=sedstat[7]
        fwkurt1=sedstat[8]
    dMean_mm=phi2mm(dMean_phi) # convertmean grain size in phi to mm
    Dmn=0.001*dMean_mm  # convert mean grain size to m for calculating bed roughess
    rhos=Rhos*ones(shape=(nclass,1)) #set sediment density, same value for all sizes
#==========================================================Calculate water parameter=============================================================
    rhow=sw_dens0(sal,wtemp)  #calculate the density of seawater [kg/m^3]
    kinvisc=KinVisc(wtemp)    #calculate the kinematic viscosity of the water
###===============================================================calculate s=====================================================================
    s=rhos/rhow
###===============================calculate critical shear stress (N/m2) and settling velocity for each===========================================
    #size class
    ws=zeros(nclass)
    Ucrit=zeros(nclass)
    tcrit=zeros(nclass)
    if(setType==1.0):
        for i in range(nclass):
            ws[i]=-tubesetvel(phi[i],rhos,wtemp,sal)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==2.0):
        for i in range(nclass):
            ws[i]=-dietrichWs(wtemp,d[i],rhos[i],rhow,0.005,0.7,3.5)
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    elif(setType==3.0):
        for i in range(nclass):
            ws[i]=Setvel1(d[i],kinvisc,s[i])
            Ucrit[i]=UstarCrit(d[i],kinvisc,s[i])
            tcrit[i]=rhow*Ucrit[i]**2
    thload=th*Cb #suspended sediment load to make deposit of observed thickness, thickness times Cb = 1-porosity
###============================================================calculate zo=======================================================================
    zoN=Dmn/30.0  #Nickaradse bed roughness
    #guess ustrc to make deposit, need this to calculate zos, bed roughness from moving flat bed, saltation
    #Note:Make the initial ustrc an input, also make ustrcAdjFac an input and delelet lines 228 and 229
    ustrcAdjFac=0.01
    UcritDmn=UstarCrit(Dmn,kinvisc,2650/rhow) #need critical shear velocity for mean grain size to calculate bed roughness
    if(zotot!=1):
        zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #zo total (Nickaradse bed roughness and Wiberg/Smith bed roughness from saltation,
###===============================================calculate shear velocity needed to suspend deposit===============================================
    #Begin loop to determine shear velocity needed to suspend deposit
    #This loops used to assume that the bed sediment disribution is the same as the suspended sediment distribution
    #This gives reasonable results when the size of the bed material is number of iterations to run
    diststep=1 #adjusts sediment grain size distribution every 2nd iteration
    ssoffa=zeros(shape=(nit/diststep,nclass)) # initiate array tracking ss deviation from desired concentration for iteration=1:nit
    ssfra=ones(shape=(nit/diststep,nclass)) #initiate array tracking ss concentrations for iteration=1:nit
    fra=zeros(shape=(nit/diststep,nclass)) #initiate array tracking bed sediment concentrations for iteration=1:nit
    offa=zeros(nit)#inititate array to track how far total suspended load is off from desired suspended load
    ustrca=zeros(nit)#inititate array to track how  u*c changes with iteration                                            #ref. Wiberg and Rubin, 1989 JGR
    cconverge=0
    for iteration in range(int(nit)):
        fr=fr/sum(fr) #make sure fractions add up to 1
        if(all(ssconverge) and cconverge):
            break
        taub=rhow*ustrc**2 # bottom stress
        taustar=taub/tcrit
        S=taustar-1     #normalized excess shear stress
        for i in range(len(S)):#returns zero when taub<=tcrit
            if S[i]<0:
                S[i]=0
            else:s[i]=s[i]
        S=S.reshape(len(S),1)
        Ca=zeros(shape=(len(s),1))
        for i in range(len(S)):
            Ca[i]=g0*Cb*fr[i]*S[i]/(1+g0*S[i])# reference concentration (volume conc; m3/m3),ref. elevation is zo, not 2 grain diameters
        # resuspension coefficient(g0) for flat moveable bed Madsen 94a "Susp. Sed. Transport on the inner shelf waters during extreme storms"
        # ICCE pg 1849 Madsen, Chisholm, Wright based on Wikramanayake 92 PhD thesis first shot, might want to change this using mean current litterature
        # log-spaced z grid
        z=logspace(log10(zotot),log10(h),num=nz,endpoint=True,base=10.0)
        z=z.reshape(len(z),1)
#*************************************************elevation of suspended load, SL******************************************
        def diff(z):
            diff=zeros(shape=(len(z)-1,1))
            for i in range(0,len(z)-1):
                diff[i]=z[i+1]-z[i]
            return (diff)
        zsl=z[0:nz-1]+diff(z)
        diff1=diff(z)
#**************************************elevation of speed calculated from eddy viscosity***********************************
        zmid=(z[0:nz-1]+z[1:nz])/2
        zmid1=zmid.reshape(len(zmid))
#*******************************************calculate eddy viscosity profile***********************************************
        #var=1 is linear eddy viscosity profile, var=2 is parabolic profile, var=3 is gelf profile
        if var==1.0:
            K=linearK(z,h,ustrc)
        elif var==2.0:
            K=parabolicK(z,h,ustrc)
        elif var==3.0:
            K=gelfK(z,h,ustrc)
        #K=(var==1)*linearK(z,h,ustrc)+(var==2)*parabolicK(z,h,ustrc)+(var==3)*gelfK(z,h,ustrc)
#*******************************************Calculate current and eddie viscosity profiles*********************************
        #Integral of 1/K
        K=K.reshape(len(K),1)
        Ki=cumsum(2.0*diff(z)/(K[0:nz-1]+K[1:nz]))
#*********************************calculate the suspended-sediment profile and sediment load for each class****************
        C=zeros(shape=(nz,nclass))
        sl=zeros(shape=(nz-1,nclass))
        sl1=zeros(shape=(nz,nclass))
        G=zeros(nclass)
        for i in range(nclass):
            #Suspended-sediment profile
            C[:len(Ca[0]),i]=Ca[i][0]
            C[len(Ca[0]):,i]=Ca[i][0]*exp(ws[i]*K[i])
            for j in range(int(nz-2)):
                sl[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
                sl1[j][i]=0.5*((C[j][i]+C[j+1][i])*diff1[j])#used to create deposit by sediment settling
#******************************************Depth-intergral of sediment profile*********************************************
            G[i] = sum(sl[:,i])
#************************************************** set cconverge to 0*****************************************************
        cconverge=0
        off=sum(G)/thload-1
        offa[iteration]=off
        if(abs(off)*100<concconvfactor):
            cconverge=1
        if(all(ssconverge) and cconverge): #if both conc and size are good enough, get outta here
            break
        ustrca[iteration]=ustrc
        ustrc=ustrc*(1-ustrcAdjFac*off) #change ustrc based on how far off suspended load is from desired value
        if(zotot!=1):
            zotot=zoN+zoWS(Dmn,ustrc,UcritDmn,rhow) #recalculate bed roughness
#*******check for convergence for size of suspended sediment and adjust bed sediment distribution every diststep***********
        def fix(n):
            if(n>0):
                if(round(n)>n):
                    fix=round(n)-1
                else:fix=round(n)
            if(n<0):
                if(round(n)<n):
                    fix=round(n)+1
            return(n)    
        doss=diststep*fix(iteration/diststep)
        if(doss==iteration): #redo distribution every diststep steps
            ssfr=G/sum(G) #calculate susp sed size distribution for model
            ssfra[doss/diststep,:]=ssfr
            fra[doss/diststep,:]=fr.reshape(1,len(fr))
            ssfr=ssfr.reshape(len(ssfr),1)
            #keep track of whether model is converging ** array indexes
            #different from iteration
            for k in range(nclass):
                if(ssfrwanted[k]!=0): # only adjust for classes with material in them
                    ssoff=ssfr[k][0]/ssfrwanted[k]-1
                    ssoffa[doss/diststep,k]=ssoff #change bed grain size distribution based on suspended load distribution
                    fr[k]=fr[k]*(1-0.01*ssoff)
                    if (abs(ssoff)*100<sizeconverge[k]): #check for convergence for each size class
                        ssconverge[k]=1.0 #set to 1 if OK
                if(all(ssconverge) and cconverge):
                    break #if both conc and size are good enough, get outta here
    spd=zeros(len(z)-1)
    diff1=diff(z)
    diff1=diff1.reshape(len(diff1),)
    spd[:]=cumsum(diff1*ustrc*ustrc/gelfK(zmid,h,ustrc))#speed calculated using eddy viscosity profile integration
    if(iteration==nit):
        print('Model may not have converged, use extreme caution with these resuls')
    #stop clock
    predictedssload=sum(G[:])
    maximumspeed=max(spd)
    spd1=0.5*(spd[0:len(spd)-1]+spd[1:len(spd)])
    spd1=spd1.reshape(len(spd1),1)
    #print zsl
    #print spd1
    #print shape(spd)
    #print shape(z)
    #fig3=figure()
    #plot(spd,zsl)
    #title("Velocity Ditribution in vertical direction at x=%d"%x+"m")
    #ylabel("Distance from bed (m)")
    #xlabel("Veloctiy (m/s) ")
    #xlim(0,8)
    #savefig(filename+' .png',dpi=100)
    diff2=zeros(shape=(len(diff1)-1,1))
    diff2=diff1[1:len(diff1)]
    diff2=diff2.reshape(len(diff2),1)
    avgspeed=sum((spd1) * diff2)/h
    MaxFroude=maximumspeed/sqrt(9.8*h)
    AvgFroude=avgspeed/sqrt(9.8*h)
    if (var==1): ed='Linear'
    elif (var==2): ed='Parabolic'
    elif (var==3): ed='Gelfenbaum'
    print('Deposit thickness:th=%4.3f'%th+'m \n')
    print('Water depth:h=%4.1f'%h+'m \n')
    print('Size classes: %i'%nclass+'\n')
    print('D mean: %.3e'%Dmn+'m \n')
    print('Viscosity Profile:%s'%ed+'\n')
    print('\n Model Results \n')
    print('Shear velocity:%4.2f'%ustrc+'m/s \n')
    print('Bed Roughness:%.3e'%zotot+'\n')
    print('Sediment load needed:%4.3f'%thload+'m \n')
    print('Sediment load predicted:%4.3f'%predictedssload+'m \n')
    print('Max. speed:%5.2f'%maximumspeed+'m/s \n')
    print('Avg. speed:%5.2f'%avgspeed+'m/s \n')
    print('Max. Froude:%4.2f'%MaxFroude+'\n')
    print('Avg. Froude:%4.2f'%AvgFroude+'\n')
    fp=open('result.txt','wb')
    fp.write("Model stats \n")
    fp.write("datestr(now) \n")
    fp.write("iterations:")
    fp.write("%i \n"%iteration)
    fp.write("Input Data \n")
    fp.write("Deposit thickness (m):")
    fp.write("%4.3f \n"%th)
    fp.write("Depth (m):")
    fp.write("%4.1f \n"%h)
    fp.write("size classes:")
    fp.write("%i \n"%nclass)
    fp.write("Mean grain diameter (m):")
    fp.write("%.3e \n"%Dmn)
    fp.write("Eddy viscosity profile:")
    fp.write("%s \n"%ed)
    fp.write("Model Results: \n")
    fp.write("ustarc (shear velocity) (m/s):")
    fp.write("%4.2f \n"%ustrc)
    fp.write("ztot (bed roughness):")
    fp.write("%.3e \n"%zotot)
    fp.write("thload (sediment load need) (m):")
    fp.write("%4.3f \n"%thload)
    fp.write("predictload (sediment load predicted) (m): ")
    fp.write("%4.3f \n"%predictedssload)
    fp.write("max. speed (m/s):")
    fp.write("%5.2f \n"%maximumspeed)
    fp.write("avg. speed (m/s):")
    fp.write("%5.2f \n"%avgspeed)
    fp.write("max. Froude:")
    fp.write("%4.2f \n"%MaxFroude)
    fp.write("avg. Froude:")
    fp.write('%4.2f \n'%AvgFroude)
    fp.write("Detail of structure: \n")
    fp.write("ssfrwanted (suspended size distribution wanted): \n")
    fp.write("%s \n"%ssfrwanted)
    fp.write("ssconverge (keeps track of whether amount of sediment in each size class is acceptable): \n")
    fp.write("%s \n"%ssconverge)
    fp.write("cconverge: \n")
    fp.write("%s \n"%cconverge)
    fp.write("ssfra(array to track of suspended sediment concentration in each size class): \n")
    fp.write("%s \n"%ssfra)
    fp.write("offa(array to track how far total suspended load is off from load required to create deposit): \n")
    fp.write("%s \n"%offa)
    fp.write("C(Suspended-sediment profile (calculated for each size class)): \n")
    fp.write("%s \n"%C)
    fp.write("G(Suspended load (calculated for each size class)): \n")
    fp.write("%s \n"%G)
    fp.write("ws(settling velocity (m/s)): \n")
    fp.write("%s \n"%ws)
    fp.write("SL: \n")
    fp.write("%s \n"%sl)
    fp.write("z(log space elevation values where concentrations and velocities are calculated): \n")
    fp.write("%s \n"%s)
    fp.write("zsl: \n")
    fp.write("%s \n"%zsl)
    fp.write("zmid: \n")
    fp.write("%s \n"%zmid)
    fp.write("ustrca(array to track how  u*c changes with iteration): \n")
    fp.write("%s \n"%ustrca)
    fp.write("S(excess shear stress): \n")
    fp.write("%s \n"%S)
    fp.write("spd(mean current speed calculated using current shear velocity (m/s)): \n")
    fp.write("%s \n"%spd)
    fp.write("results of the structure \n")
    fp.write("zoN(Nickaradse grain roughness): \n")
    fp.write("%s \n"%zoN)
    fp.write("tcrit(critical shear stress (N/m^2)): \n")
    fp.write("%s \n"%tcrit)
    fp.write("taub(bottom shear stress due to currents (N/m^2)): \n")
    fp.write("%s \n"%taub)
    fp.write("taustar(bottom shear stress normalized by critical shear stress): \n")
    fp.write("%s \n"%taustar)
    fp.write("phi(grain size in phi units for each size class): \n")
    fp.write("%s \n"%phi)
    fp.write("Ca(reference concentration (volume conc; m3/m3) *** ref. elevation is zo, not 2 grain diameters): \n")
    fp.write("%s \n"%Ca)
    fp.write("Ki(integral of 1/eddie viscosity): \n")
    fp.write("%s \n"%Ki)
    fp.write("off(normalized difference between suspended load and load required to create deposit): \n")
    fp.write("%s \n"%off)
    fp.write("\n End")
    fp.write('General results \n')
    fp.write('Predictload (kg/m3) \n')
    fp.write('%f \n'%predictedssload)
    fp.write('maximum speed \n')
    fp.write('%f \n'%maximumspeed)
    fp.write('average speed \n')
    fp.write('%f \n'%avgspeed)
    fp.write('Max froude \n')
    fp.write('%f \n'%MaxFroude)
    fp.write('average froude \n')
    fp.write('%f \n'%AvgFroude)
    fp.close()
    return(maximumspeed,avgspeed,MaxFroude,AvgFroude)
