'''
Created on 2013-5-26
==========================The TsuSedForm Tsunami Sediment formation Model==========================
#========================================Calculate Part============================================
##=========================Calculate SL and z for sediment Formation===============================
##================================Calculate Settling velocity======================================
##=================================Calculate Settling Time=========================================
##=========================Calculate hittingsize and concertration=================================
##=================================Determine Boundary==============================================
##=================================Calculate grainsize distribution================================
##==============================Calculate stats for each interval==================================
#========================================Plot figure part==========================================
##========================Plot grainsize distribution for each interval============================
##===================plot figure for Grain Statistics, 1st Moment for layers=======================
##========================plot figure for Standard Deviation for layers============================
##============================plot figure for FW mean for layers===================================
##==========================plot figure for FW sorting for layers==================================
##==========================plot figure for FW Skewness for layers=================================
##==========================plot figure for FW Kurtosis for layers=================================
##====================plot figure for Thickness vs. Time sediment hits bottom======================
#=========================================OUTPUT Data Part=========================================
@ author Hui Tang- tanghui@vt.edu
'''
#==================================================================================================
#==================================================================================================
#============================The TsuSedForm Tsunami Sediments Formation Model======================
#==================================================================================================
#==================================================================================================
#This model is based on Jaff's TsuSedMod model's grading and static psrts, grading and static parts
#orignate from Jaffe 9/29/2003, Made into python by Hui Tang, Virginia Tech 1/25/2013
#This function is desinged  to form suspended sediments grainsize distribution based on whole grainsize
#distribution in tsunami sediments.
#Input:         name1-Filename of data document from Sousby model
#               h-Water depth m
#               th1-Sediment load m
#               intv-sub-interval bounadry 
#               x-Sample location
#Function Used: Tsusemod.py,tubesetvel.py and sedstatics.py
#Output:        se-Suspened sediments grainsize in phi
#               Fr-Suspended seidments percentage for each grain size
from pylab import *
from numpy import *
from ReadCSV import *
from function import *
from sedstats import *
from TsusedMod import *
from output2CSV import *
from ReadCSV import *
def Tsusedform(name1,intv,h,x):
    with open('parameter_P14abc.txt','r') as f:
            for line in f:
                try:
                    s=line.split('=');n=s[0];v=s[1];#remove '\n'
                    if os.name!='nt':v=v[0:len(v)] # remove '\r' for Unix-system
                    if n=='porosity': porosity=float(v)
                    elif n=='sediment_density': Rhos=float(v)
                    elif n=='Number_sample': N=int(v)
                    elif n=='Filename': filename=str(v)
                    elif n=='water_temperature': wtemp=float(v)
                    elif n=='salinity': sal=float(v)
                    elif n=='pi': pi=float(v)
                    elif n=='Filename_Sousby': filename=str(v)
                except IndexError as e:
                    continue
                except ValueError as e:
			        continue
#========================================Calculate Part============================================
##=========================Calculate SL and z for sediment Formation===============================
    data=Tsusedmod(name1)
    sl=data[0]# Sediment concentration profile
    z=data[1]# Elevation
    phi=data[2]# Grain size
    numrow_sl=len(sl)
    numcol_sl=len(sl[0]) #get matrix sl size
    if((len(z)!=numrow_sl) or (len(phi)!=numcol_sl)):
        print('Number of size classes or number of elevations do not match profiles')
    else:
        z=z.reshape(len(z),1) #z must be a column ector
##================================Calculate Settling velocity======================================
        setvel=tubesetvel(phi,Rhos,wtemp,sal)#calculate settling velocity (m/s)
        setvel=setvel.reshape(1,len(setvel))
        numrow_setvel=len(setvel)
        numcol_setvel=len(setvel[0])
        if(numrow_setvel>numcol_setvel):
            setvel=setvel #make sure setvel is a row vector
##=================================Calculate Settling Time=========================================
        settlingtime=zeros(shape=(numrow_sl,numcol_setvel))#Settling time for each grain size
        stime=zeros(shape=(numrow_sl*numcol_setvel,3))#store for sorting based on settling time
        for i in range(numrow_sl):
            for j in range(numcol_setvel):
                settlingtime[i][j]=z[i]/setvel[0][j]
##=========================Calculate hittingsize and concertration=================================
###==========================Sort by time it takes to hit bottom===================================
        stime[:,0]=settlingtime.reshape(int(numrow_sl*numcol_setvel),)
        ind=argsort(stime[:,0],axis=0)
        k=0
        for i in range(len(ind)):#keep track of grain size and elevation
            stime[k][1]=int(int(ind[i])/int(numcol_setvel))#grain size
            stime[k][2]=int(int(ind[i])%int(numcol_setvel))#location
            k=k+1
        stime1=sorted(stime[:,0])
###===========================find size of sediment hitting the bed================================
        phiind=stime[:,2] # find size of sediment
        hitsize=zeros(len(phiind))#initiate array with size of particles hitting bed
        for i in range(len(phiind)):
            hitsize[i]=phi[phiind[i]]    # keep track of order sediment size hitting the bed 
###=========================find elevations that sediment started from=============================
        zind=stime[:,1]# find elevation of sediment
        hitheight=zeros(len(zind)) # initiate array with height of particles hitting bed
        for i in range(len(zind)):
            hitheight[i]=z[zind[i]] #keep track of order of height sediment came from
        hitload=zeros(len(zind)) # initial array with concentration of particles hitting bed
        for i in range(len(zind)):
            hitload[i]=sl[zind[i],phiind[i]]
##=========================Calculate sediment thickness for each hitting grain size=================
        thickness=cumsum(hitload/(1-porosity))#cumulative thickness of deposit
        thickness2=zeros(len(thickness))
##=================================Determine Boundary===============================================
        intbind=zeros(len(thickness))
        indthick=zeros(len(intv))
        if(len(intv)==1):#calculate interval boundaries for set intervals
            for intbound in range(intv,thickness[len(thickness)]-intv/2,intv):
                for i in range(len(thickness)):
                    if(thickness[i]<intv[intbound]):
                        thickness2[j]=thickness[i]
                        j=j+1
                intbind[round(intbound/intv)]=max(thickness2)
            intbind2=zeros(len(intbind)+2)
            intbind2[0]=0
            i=1
            while(i<len(intbind)+1):
                intbind2[i]=intbind[i-1]
                i=i+1
            intbind[i]=len(thickness)#add break at first and last thickness value
        else:#if interval boundaries are specified
            for intbound in range(len(intv)):
                for i in range(len(thickness)):
                    j=0
                    if(thickness[i]<intv[intbound]):
                        thickness2[j]=thickness[i]
                        j=j+1
                        indthick[intbound]=i
                intbind[intbound]=max(thickness2)
            indthick2=zeros(len(intv)+1)
            intbind2=zeros(len(intbind)+1)
            intbind2[0]=1
            i=1
            while(i<len(intbind)+1):
                intbind2[i]=intbind[i-1]
                i=i+1
            j=0
            while(intbind2[j]!=0):
                j=j+1
            intbind3=zeros(int(j))
            for i in range(int(j)):
                intbind3[i]=intbind2[i]
            intbind3[0]=0
            k=1
            while(k<len(indthick)+1):
                indthick2[k]=indthick[k-1]
                k=k+1
##=================================Calculate grainsize distribution================================
        st=zeros(len(intbind3)-1)# start point 
        last=zeros(len(intbind3)-1) # end point 
        sload=zeros(shape=(len(indthick),len(phi)))# sediment load
        m1=zeros(len(st))# mean grain size
        m3=zeros(len(st))
        m4=zeros(len(st))
        stdev=zeros(len(st))# Standard deviation 
        fwmedian=zeros(len(st))# Median grain size 50%
        fwsort=zeros(len(st))# sorting phi
        fwmean=zeros(len(st))# Mean grain size
        fwskew=zeros(len(st))# Skewness
        fwkurt=zeros(len(st))# kurtosis
        p=zeros(shape=(len(st),len(phi)))# grain size
        s=zeros(shape=(len(st),len(phi))) # sediment load
        wpc=zeros(shape=(len(st),len(phi)))
        midint=zeros(len(indthick2)-1)
        for i in range(0,len(intbind3)-1):
            st[i]=intbind3[i]
        j=0
        for i in range(1,len(intbind3)):
            last[j]=intbind3[i]
            j=j+1
        thickness1=zeros(shape=(len(indthick2)-1,len(phi)))
        Fr=zeros(shape=(len(indthick2)-1,len(phi)))#Pecentage for each grain size
        totalth=zeros(len(indthick2)-1)
        for j in range(0,len(indthick2)-1):
            for k in range(len(phi)): #get cumulative sediment loads for each phi class
                #sload[k][j]=sum(hitload[int(indthick[j]):int(indthick[j+1])])
                for l in range(int(indthick2[j]),int(indthick2[j+1])):
                    if(phi[stime[l][2]]==phi[k]):
                        sload[j][k]=hitload[l]+sload[j][k]
                        thickness1[j][k]=sload[j][k]/(1-porosity)
##==============================Calculate stats for each interval==================================
            sedstats1=sedstats(phi,sload[j])
            m1[j]=sedstats1[0]
            m3[j]=sedstats1[1]
            m4[j]=sedstats1[2]
            stdev[j]=sedstats1[3]
            fwmedian[j]=sedstats1[4]
            fwsort[j]=sedstats1[5]
            fwmean[j]=sedstats1[6]
            fwskew[j]=sedstats1[7]
            fwkurt[j]=sedstats1[8]
            p[j,:]=phi
            s[j,:]=sload[j]
            wpc[j,:]=100*sload[j]/(sum(sload[j])) #weight percent in each phi interval
            midint[j]=(intbind3[j]+intbind3[j+1])/2.0
        for j in range(0,len(indthick2)-1):
            for k in range(len(phi)): #get cumulative sediment loads for each phi class
                Fr[j][k]=thickness1[j][k]*100/sum(thickness1[j])
                totalth[j]=sum(thickness1[j])
#========================================Plot figure part==========================================
##========================Plot grainsize distribution for each interval============================
        #This part is designed for plot grainszie distribution and ouput data for each interval 
        for i in range(0,int(len(intv))):
            pe=Fr[i,:]
            name=name1.split(".")
            filename=name[0]+"_suspended_sample%02d"%i+".csv"
            fp=output2CSV1(filename,phi,pe,totalth[i],h,x,'phi','fraction','thickness','depth','location')
            f1=figure()
            ax1=subplot(111)
            ax1.plot(phi,pe,'-')
            ax1.set_xlim(0,10)
            ax1.set_ylim(0,40)
            ax1.set_xlabel('grain size($\phi$)')
            ax1.set_ylabel('The percentage of different sediment grains size(%)')
            #savefig(filename+'.png',dpi=100)
            close()
##============================plot figure for Grain Statistics, 1st Moment for layers================
        ##This figure is desinged to show the relationship between location and first hit grain size
        fig1=figure()
        fname1='Grain Statistics, 1st Moment for layers'
        ax1=subplot(111)
        ax1.set_xlabel('1st Moment (phi)')
        ax1.set_ylabel('Midpoint of layer (m)')
        ax1.set_title(fname1)
        ax1.plot(m1,midint,'-')
        close()
        #savefig(fname1+' Sample%i.png'%ni,dpi=100)
##============================plot figure for Standard Deviation for layers=============================
        ##This figure is designed to show the relationship between location and standard deviation 
        fig2=figure()
        fname2='Grain Statistics, Standard Deviation for layers'
        ax2=subplot(111)
        ax2.set_xlabel('Standard Deviation (phi)')
        ax2.set_ylabel('Midpoint of layer (m)')
        ax2.set_title(fname2)
        ax2.plot(stdev,midint,'-')
        close()
        #savefig(fname2+' Sample%i.png'%ni,dpi=100)
##=================================plot figure for FW mean for layers===================================
        ##This figure is designed to show the relationship between location and mean
        fig3=figure()
        fname3='Grain Statistics, FW mean for layers'
        ax3=subplot(111)
        ax3.set_xlabel('FW Mean (phi)')
        ax3.set_ylabel('Midpoint of layer (m)')
        ax3.set_title(fname3)
        ax3.plot(fwmedian,midint,'-')
        close()
        #savefig(fname3+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW sorting for layers==================================
        ##This figure is designed to show the relationship between location and sorting
        fig4=figure()
        fname4='Grain Statistics, FW sorting for layers'
        ax4=subplot(111)
        ax4.set_xlabel('Sorting (phi)')
        ax4.set_ylabel('Midpoint of layer (m)')
        ax4.set_title(fname4)
        ax4.plot(fwsort,midint,'-')
        close()
        #savefig(fname4+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Skewness for layers=================================
        ##This figure is designed to show the relationship between location and skewness
        fig5=figure()
        fname5='Grain Statistics, FW Skewness for layers'
        ax5=subplot(111)
        ax5.set_xlabel('Skewness')
        ax5.set_ylabel('Midpoint of layer (m)')
        ax5.set_title(fname5)
        ax5.plot(fwskew,midint,'-')
        close()
        #savefig(fname5+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Kurtosis for layers=================================
        ##This figure is designed to show the reltaionship between location and kurtosis
        fig6=figure()
        fname6='Grain Statistics, FW Kurtosis for layers'
        ax6=subplot(111)
        ax6.set_xlabel('Kurtosis')
        ax6.set_ylabel('Midpoint of layer (m)')
        ax6.set_title(fname6)
        ax6.plot(fwskew,midint,'-')
        close()
        #savefig(fname6+' Sample%i.png'%ni,dpi=100)
##==========================plot figure for Thickness vs. Time sediment hits bottom=====================
        ##This figures is designed to show the relationship between thickness and time sediment hits bottom
        fig7=figure()
        fname7='Thickness vs. Time sediment hits bottom'
        ax7=subplot(111)
        ax7.set_xscale('log',basex=10)
        ax7.set_xlabel('Time(s)')
        ax7.set_ylabel('Thickness (m)')
        ax7.set_title(fname7)
        ax7.plot(stime1,thickness,'-')
        print shape(stime1)
        close()
        #savefig(fname7+' Sample%i.png'%ni,dpi=100)
        return (phi,Fr)
