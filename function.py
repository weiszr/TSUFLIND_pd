'''
Created on 2013-1-25
=====================================The tsunami sediment model3 (function part)================================
#===============================================function dietrichWs=============================================
#Input: 	T- temperature of fluid in deg. C[Array2]
#               D- grain size in mm [Num]
#               rhos- density of sediment in gm/cm**3 [Num]
#               rhow- density of water in gm/cm**3 [Num]
#               Cs- volume concentration of sediment [Num]
#               csf- Corey shape factor, usually taken as 0.7 [Num]
#               P- roundness of grains, usaully taken as 3.5 [Num]
#Output         ws- settling velocity (cm/s) [Array2]
#function       KinVisc(T)
#=================================================Function gelfK================================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array1]
#   		h- water depth (m) [Num]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K- parabolic eddy viscosity [Array1]
#===============================================Function kinvisc================================================
#Input: 	T- temperature of fresh water in degrees celsius [Array2]
#Output:        kinvisc- kinematic viscosity in m**2/s[Num]
#==================================================function linearK=============================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array1]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K-linear eddy viscosity profile [Array1]
#================================================Function parabolick============================================
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated [Array 1]
#   		h- water depth (m) [Num]
#		ustrc- current shear velocity (m/s) [Num]
#Output:        K-parabolic eddy viscosity profile [Array1]
#====================================================Function phi2mm============================================
#Input: 	phi- Grain size in phi [Num]
#Output:        mm-Grain size in mm [Num]
#====================================================function Setvel1===========================================
#Input: 	d- grain diamter (m) [Num]
#		kinvisc- kinematic viscosity (m**2/s) [Array2]
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#Output:        ws- particle settling velocity (m/s) [Array2]
#=================================================Function sw_dens0=============================================
#Input: 	S = salinity    [psu      (PSS-78)] [Array2]
#               T = temperature [degree C (IPTS-68)][Array2]
#Output:        dens0 = density  [kg/m**3] of salt water with properties S,T, [Array2]
#               P=0 (0 db gauge pressure)  [Num]
#===================================================Function sw_smow============================================
Input: 	        T = temperature [degree C (IPTS-68)][Array2]
#Output:        dens = density  [kg/m**3] [Array2]
#================================================function tubesetvel============================================
#Input: 	phi- grain size in phi [Num]
#               rhos- density of sediment (or sphere) in gm/cm**3 [Num]
#               T- temperature of water in settling tube, deg. C [Array2]
#               S- salinity of water in settling tube (added to explore real world) [Array2]
#Output:        ws- settling velocity of grain (cm/s) [Array2]
#Functions:     KinVisc(T) and sw_dens0(S,T)
#================================================function UstarCrit=============================================
#Input: 	d- grain diamter (m) [Num]
#		kinvisc- kinematic velocity (m**2/s) [Array2]
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#               shieldscrit- Shield's parameter for initiation of sediment transport [Num]
#		Shield's parameter= ustarcrit^2/((s-1)*g*d) 
#Output:        ustarcrit- critical shear velocity for initiation of sediment transport (m/s) [Array2]
#Function:      KinVisc(T)
#===============================================function zoD====================================================
#Inputs:    D- nominal grain diameter (m) [Num]
#	    ustarc- current shear velocity (m/s) (may be a vector) [Num]
#	    rhow- density of water (kg/m3) [Num]
#           s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#           T- temperature of water in settling tube, deg. C [Array2]
#Output:    zosD- bed roughness for saltating layer (m)  *** must be added to zoN (Nikuradse
#	    bed roughness) to get total roughness  (zo=zoN+zosD) [Array2]
#Function:  KinVisc(T),UstarCrit(D,T,s)
#===================================================function zoWS===============================================
#Inputs:    D- nominal grain diameter (m) [Num]
#	    ustarc- current shear velocity (m/s) (may be a vector) [Num]
#	    rhow- density of water (kg/m3) [Num]
#           s- ratio of density of sediment to density of fluid (2.65 for quartz and water) [Num]
#           T- temperature of water in settling tube, deg. C [Array2]
#Output:        zos- bed roughness for mobile flat bed (m)  *** must be added to zoN (Nikuradse
#		bed roughness) to get total roughness  (zo=zoN+zos)[Array]
#Function: UstarCrit(D,T,s)
First version written March 12, 2001by Bruce Jaffe in matlab
originate from Bruce Jaffe,USGS, Made into python by Hui Tang, Virginia Tech, tanghui@vt.edu
'''
from pylab import *
from types import *
import csv
#===============================================function dietrichWs===========================================
#This is a fuction that calculates settling velocity from Dietrich'82
#Input: 	T- temperature of fluid in deg. C
#               D- grain size in mm
#               rhos- density of sediment in gm/cm**3
#               rhow- density of water in gm/cm**3
#               Cs- volume concentration of sediment
#               csf- Corey shape factor, usually taken as 0.7
#               P- roundness of grains, usaully takefunction.pyn as 3.5
#Used
#               g- acceleration due to gravity, 980 cm/sec**2
#               k- 1/von karmen's constant, 2.5
#Calculated
#               mu- dynamic viscosity of water (poise)
#               musubs- dynamic viscosity of fluid (poise)
#               rhof- density of fluid in gm/cm^3
#               nu- kinematic viscosity of fluid (poise)
#               dstar- dimensionless particle size
#               wstar- dimensionless settling velocity
#Output         ws- settling velocity (cm/s)
#need function  KinVisc
#originate from Bruce Jaffe, USGS, Agu 6,2003
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def dietrichWs(T,D,rhos,rhow,Cs,csf,P):
    g=9.8  # acceleration due to gravity, cm/sec**2
    k=2.5
    nargin=5
    if(nargin==5):# if csf and P are not set, use defaults
        csf=0.7    # set Corey shape factor of grains to default value
        P=3.5   # set roundness of grains to default value
    lu=KinVisc(T)
    mu=KinVisc(T)*rhow*10000         # convert kinematic viscosity in m**2/s to dynamic viscosity in poise 
    #calculate dynamic viscosity for water and sediment
    musubs=mu*(1+k*Cs)
    #calculate fluid density
    rhof=rhow+(rhos-rhow)*Cs
    #calculate nu for the sediment conentration
    nu=musubs/rhof
    #calculate settling velocity, converts grain size from mm to cm
    dstar=(rhos-rhof)*g*(D/10)**3/(rhof*nu**2)             # eq. 6, Dietrich '82 (D82), converting size from mm to cm
    dstarl=log10(dstar)
    r1=-3.76715+1.92944*dstarl-0.09815*dstarl**2-0.00575*dstarl**3+0.00056*dstarl**4 #fitted equation for size and density effects, eq. 9 D82
    omc=1.0-csf
    r2=log10(1-omc/0.85)-omc**2.3*tanh(dstarl-4.6)+0.3*(0.5-csf)*omc**2*(dstarl-4.6) #fitted equation for shape effects. eq. 16 D82
    r3=(0.65-(csf/2.83)*tanh(dstarl-4.6))**(1+(3.5-P)/2.5) #fitted equation for roundness effects, eq. 18 D82
    wstar= r3*10**(r1+r2)                                       # eq. 19, D82
    ws=((rhos-rhof)*g*nu*wstar/rhof)**(1.0/3) # rearranged eq. 5 from D82
    return(ws)
#=================================================Function gelfK=============================================
#This is a function that calculates the parabolic eddy viscosity profile using
#Input: 	z- vector with elevations (m) whsedstats.pyere the eddy viscosity is calculated
#   		h- water depth (m)
#		ustrc- current shear velocity (m/s)
#Output:        K- eddy viscosity
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def gelfK(z,h,ustrc):
    vk=0.41
    K=zeros(len(z))
    for i in range(len(z)):
        K[i]=vk*ustrc*z[i]*exp((-z[i]/h)-3.2*(z[i]/h)**2+2.13*(z[i]/h)**3)
    return K
#===============================================Function kinvisc=============================================
#This is a function that calculates kinematic viscosity of fresh water given temperature
#Uses eq. 3.1.4 from Van Rijn, 1993, "Principles of Sediment Transport in Rivers, Estuaries, and Coastal Seas."
#Input: 	T- temperature of fresh water in degrees celsius
#Output:        kinvisc- kinematic viscosity in m^2/s   		
#originate from Bruce Jaffe, USGS, Jan 30,2001
#Made into python by Hui Tang, Virginia Tech, Jan 30, 2013 tanghui@vt.edu
def KinVisc(T):
    kinvisc=(1.14-0.031*(T-15)+0.00068*(T-15)**2)*10**-6
    return kinvisc
#==================================================function linearK==========================================
#This is a fuction that calculates the linear eddy viscosity profile
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated
#		ustrc- current shear velocity (m/s)
#originate from Bruce Jaffe, USGS, April 3,2sedstats.py001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def linearK(z,h,ustrc):
    K=zeros(len(z))
    vk=0.41
    for i in range(len(z)):
        K[i]=vk*z[i]*ustrc
    return K
#================================================Function parabolick=========================================
#This is a function that calculates the parabolic eddy viscosity profile
#Input: 	z- vector with elevations (m) where the eddy viscosity is calculated
#   		h- water depth (m)
#		ustrc- current shear velocity (m/s)
#Output:        K-parabolic eddy viscosity profile
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def parabolicK(z,h,ustrc):
    K=zeros(len(z))
    vk=0.41
    for i in range(len(z)):
        K[i]=vk*z[i]*(1-z[i]/h)*ustrc
    return K
#====================================================Function phi2mm==========================================
#This is a function that converts phi to mm
#Input: 	phi- Grain size in phi
#Output:        mm-Grain size in mm
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def phi2mm(phi):
    mm=2**(-phi)
    return mm
#====================================================function Setvel1=========================================
#This is a fuction calculate settling velocity of sediment.  Formula from 
#Ole Madsen's short course at Coastal Seds '99.  These are fits to data from 
#Dietrich '82, van Rijn '90, Cheng '97.
#Input: 	d- grain diamter (m)
#		kinvisc- kinematic viscosity (m**2/s)
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water)
#Output:        ws- particle settling velocity (m/s)
#Variation of Setvel allowing diameters smaller than 0.00006 m
#originate from Bruce Jaffe, USGS, April 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def Setvel1(d,kinvisc,s):
    g=9.80 #gravitational acceleration
    # calculate sstar
    sstar=(d/(4*kinvisc))*sqrt((s-1)*g*d)
    if(d>= 0.0001 and d <= 0.001 ): #diametesedstats.pyr between 0.1 and 1mm
        wstar=1/(5.40/sstar + 0.92)
    elif(d>= 0.0001 and d <= 0.005): #diameter between 0.1 and 5mm
        wstar=1/(5.21/sstar + 0.88)
    ws=wstar*sqrt((s-1)*g*d) #settling velocity in m/s (?)
    return ws
#=================================================Function sw_dens0==========================================
#This is a function that calculates Density of Sea Water at atmospheric pressure using UNESCO 1983 (EOS 1980) polynomial.
#DISCLAIMER:
#This software is provided "as is" without warranty of any kind.  
#REFERENCES:
#Unesco 1983. Algorithms for computation of fundamental properties of seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#UNESCO 1983 p17  Eqn(14) Millero, F.J & Poisson, A. INternational one-atmosphere equation of state for seawater. Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
#Input: 	S = salinity    [psu      (PSS-78)]
#               T = temperature [degree C (IPTS-68)]
#Output:        dens0 = density  [kg/m**3] of salt water with properties S,T,
#               P=0 (0 db gauge pressure) 
#originate from Phil Morgan, Sep 05,1992 morgan@ml.csiro.au
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def sw_dens0(S,T):
    nargin=2
    if(nargin!=2):
        print('error:sw_dens0.m: Must pass 2 parameters')
    else:
        #if (type(S)=='float' or type(T)=='float'):
        b0=8.24493e-1
        b1=-4.0899e-3
        b2=7.6438e-5
        b3=-8.2467e-7
        b4=5.3875e-9
        c0=-5.72466e-3
        c1=+1.0227e-4
        c2=-1.6546e-6
        d0=4.8314e-4
        dens1=sw_smow(T)
        dens=dens1+(b0+(b1+(b2+(b3+b4*T)*T)*T)*T)*S+(c0+(c1+c2*T)*T)*S*sqrt(S)+d0*S**2
    return(dens)
#===================================================Function sw_smow==========================================
#This is a function that calculates Denisty of Standard Mean Ocean Water (Pure Water) using EOS 1980.
#DISCLAIMER:
#This software is provided "as is" without warranty of any kind.  
#REFERENCES:
#Unesco 1983. Algorithms for computation of fundamental properties of seawater, 1983. _Unesco Tech. Pap. in Mar. Sci._, No. 44, 53 pp.
#UNESCO 1983 p17  Eqn(14) Millero, F.J & Poisson, A. INternational one-atmosphere equation of state for seawater. Deep-Sea Research Vol28A No.6. 1981 625-629.    Eqn (6)
#Input: 	T = temperature [degree C (IPTS-68)]
#Output:        dens = density  [kg/m^3] 
#originate from Phil Morgan, Sep 05,1992 morgan@ml.csiro.au
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def sw_smow(T):
    #check input arguments
    #test inputs
    nargin=1
    if(nargin!=1):
        print('error:sw_smow.m: Only one input argument allowed')
    else:
        a0=999.842594
        a1=6.793952e-2
        a2=-9.095290e-3
        a3=1.001685e-4
        a4=-1.120083e-6
        a5=6.536332e-9
        dens=a0+(a1+(a2+(a3+(a4+a5*T)*T)*T)*T)*T
    return(dens)
#================================================function tubesetvel=========================================
#This is a function to determine settling velocity in USGS settling
#tubes from grain size output of sedsize (phi).  Uses Gibbs'
#settling formula (from excel spreadsheet provided by John
#Penscil, a contrator who wrote the programs for Mike
#Torresan).  This spreadsheet is GibbsEquations.xls.  I obtained
#a copy in 2003.
#Input: 	phi- grain size in phi
#               rhos- density of sediment (or sphere) in gm/cm**3
#               T- temperature of water in settling tube, deg. C
#               S- salinity of water in settling tube (added to explore real world)
#Used:
#               g- acceleration due to gravity (cm/s**2)
#               rhof- density of fluid in settling tube (gm/cm**3)
#               kv- kinematic viscosity of fluid in settling tube(m/s**2)
#               dynvisc- dynamic viscosity of fluid in settling tube (poise)
#               r- radius of grain (or sphere) in mm
#Output:        ws- settling velocity of grain (m/s)
#only need to input grain size in phi, will use defaults of
#rhos=2.65, T=21, S=0
#Need functions KinVisc and sw_dens0sedstats.py
#originate from Bruce Jaffe, USGS, Aug 1,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def tubesetvel(phi,rhos,T, S):
    g=979.9 #gravitational acceleration in cm/s^2
    nargin=1
    if(nargin==1): #use defaults if only phi size is provided
        rhos=2.65 #quartz density
        T=11.6  #21 deg. C water temperature
        S=0#fresh water in settling tubes
    #calculate dynamic viscosity of water in settling tube (poise)
    kv=KinVisc(T)      # m/s^2
    #kv=1.004e-6
    #print kv
    rhof=sw_dens0(S,T) # kg/m^3, 
    #print rhof
    dynvisc=10*kv*rhof #dynamic viscosity, converted to poise (dyne-sec/cm^2)
    #print dynvisc
    rhof=rhof/1000.0#convert fluid density to gm/cm^3 for use in formula
    #reproduce spreadsheet columns and calculate settling velocity
    r=(1.0/(2.0**phi))/20.0
    term= sqrt((9.0*dynvisc**2)+(g*r**2*rhof*(rhos-rhof)*(0.015476+(0.19841*r))))
    ws=(term-(3*dynvisc))/((rhof)*(0.011607+(0.14881*r)))/100.0
    return ws
#================================================function UstarCrit==========================================
#This is a fuction that calculates the critical shear velocity for initiation of sediment transport. Formula from Ole Madsen's short course
#at Coastal Seds '99.
#Input: 	d- grain diamter (m)
#		kinvisc- kinematic velocity (m**2/s)
#		s- ratio of density of sediment to density of fluid (2.65 for quartz and water)
#               shieldscrit- Shield's parameter for initiation of sediment transport
#		Shield's parameter= ustarcrit^2/((s-1)*g*d)
#Output:        ustarcrit- critical shear velocity for initiation of sediment transport (m/s)
#Function:      KinVisc(T)
#originate from Bruce Jaffe, USGS, Fed 8,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def UstarCrit(d,kinvisc,s):
    #kinvisc=KinVisc(T)    
    if(type(kinvisc)=='float'):
        mK=len(kinvisc)
        nK=len(kinvisc[0][:])
        sstar=zeros(shape=(mK,nK))
        ustarcrit=zeros(shape=(mK,nK))
        for i in range(mK) :
            for j in range(nK):
                g=9.80
                sstar[i][j]=d/(4*kinvisc[i][j])*sqrt((s-1)*g*d)
                if (sstar[i][j]>= 0.8): 
                    shieldscrit=0.095*sstar[i][j]**(2.0/3.0) + 0.056*(1-exp((-sstar[i][j]**(3.0/4.0))/20.0))
                else:
                    shieldscrit=0.1*sstar[i][j]**(-2.0/7.0)
                ustarcrit[i][j]=sqrt((s-1)*g*d)*sqrt(shieldscrit) # critical shear velocity in m/s
        
    else:
        g=9.80 #gravitational acceleration
        #calculate sstar
        sstar=d/(4*kinvisc)*sqrt((s-1)*g*d)
        #print kinvisc
        #calculate Shield parameter for initiation of sediment transport
        if (sstar >= 0.8): 
            shieldscrit=0.095*sstar**(-2.0/3.0) + 0.056*(1-exp((-sstar**(3.0/4.0))/20.0))
        else:
            shieldscrit=0.1*sstar**(-2.0/7.0)
        ustarcrit=sqrt((s-1)*g*d)*sqrt(shieldscrit) # critical shear velocity in m/s
    return (ustarcrit)
#===============================================function zoD=================================================
#This is a fuction that calculates zosD, Dietrich's 1982 expression for bottom
#roughness created by saltating layer.
#Inputs:    D- nominal grain diameter (m)
#	    ustarc- current shear velocity (m/s) (may be a vector)
#	    rhow- density of water (kg/m3)
#Used:      alpha- constant, 0.077, from Muddy Creek, Wyoming data
#           (averages about 0.10, can be a max of 0.13 according
#           to Wiberg and Rubin 1989)
#Output:    zosD- bed roughness for saltating layer (m)  *** must be added to zoN (Nikuradse
#	    bed roughness) to get total roughness  (zo=zoN+zosD)
#References:Dietrich, W.D., Flow, boundary shear stress, and sediment transport in a river meander,
#           Ph.D. dissertation, 261 p., Univ. of Wash., Seattle, 1982        
#           Wieberg, P. L. and Rubin, D. M., 1989, Bed roughness produced by
#	    saltating sediment. J. Geophys. Res., V. 94, no. C4, pp. 5011-5016.
#originate from Bruce Jaffe, USGS, Fed 8,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def zosD(T,D,ustarc,rhow,s):
    alpha=0.077         # constant, set using Muddy Creek, Wyoming data
    s=2.65
    kinvisc=KinVisc(T)
    ustarcrit=UstarCrit(D,T,s)
    alpha=0.077
    taub=rhow*ustarc**2
    mU=len(ustarcrit)
    nU=len(ustarcrit[0][:])
    taucrit=zeros(shape=(mU,nU))
    Tstar=zeros(shape=(mU,nU))
    zosD=zeros(shape=(mU,nU))
    for i in range(mU):
        for j in range(nU):
            taucrit[i][j]=rhow*ustarcrit[i][j]**2
            Tstar[i][j]=taub/taucrit[i][j]
            zosD[i][j]=alpha*D*(0.6*Tstar[i][j]/(1+0.2*Tstar[i][j])) #Dietrich 1982 from Wiberg and Rubin, 1989 equation 3
    return zosD
#===================================================function zoWS===========================================
#This is a fuction that calculates zos, the bed roughness parameter for mobile
#flat bed, using the formula from Wieberg and Rubin, 1989.
#ws stands for "Wiberg and Smith"
#Input: 	D- nominal grain diameter (m)
#		ustarc- current shear velocity (m/s) (may be a vector)
#		ustarcrit- critical shear velocity for initiation of motion (m/s) (may be a vector)
#		rhow- density of water (kg/m3))
#Used:          gammawWS- a constant determined by Wieberg and Rubin to fit the data
#		taub- bottom shear stress (N/m2)
#		taucrit- critical shear stress for initiation of motion (N/m2)
#		Tstar- transport stage, ratio of taub to taucrit
#		delb- average saltation height
#		a1- constant used in fomulation of average saltation height
#		a2- constant used in fomulation of average saltation height
#Output:        zos- bed roughness for mobile flat bed (m)  *** must be added to zoN (Nikuradse
#		bed roughness) to get total roughness  (zo=zoN+zos)
#Reference:     Wieberg, P. L. and Rubin, D. M., 1989, Bed roughness produced by
#		saltating sediment. J. Geophys. Res., V. 94, no. C4, pp. 5011-5016.
#originate from Bruce Jaffe, USGS, Mar 3,2001
#Made into python by Hui Tang, Virginia Tech, Jan 24, 2013 tanghui@vt.edu
def zoWS(D,ustarc,ustarcrit,rhow):
    gammaWS = 0.056
    a1=0.68
    a2=0.0204*(log(D*100))**2+0.0220*log(D*100)+0.0709 	#equation 6, D converted to cm, Wiberg and Rubin, 1989
    taub=rhow*ustarc**2
    taucrit=rhow*ustarcrit**2
    Tstar=taub/taucrit
    delb=D*a1*Tstar/(1+a2*Tstar)
    zos=gammaWS*delb   # equation 7 Wiberg and Rubin, 1989
    return zos

