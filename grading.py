from pylab import *
from numpy import *
from function import *
from sedstats import *
#============================================Subfunction grading==============================================
#This is a function to calculate grading of sediment based on concentration profiles for different grain size
#classes.
#Input:             sl:Sediment load for each size class (columns) for given z(rows)
#                   phi:size classes in phi(classes=columns of C),row vector
#                   z:elevation (m) above bottom of sediment loads (elevations=row of sl), column vector
#                   intv:interval (m) of sediment texture stats, maybe a set interval or specified(e.g thickness
#                        1cm,2cm,1cm,4cm intv[0.01 0.03 0.04 0.08]),default is 0.01 if not specified
#                   porosity:porosity of deposit,default is 0.35 if not specified
#Function used:     tubesetvel.py, sedstats.py
#Note:set defaults for intervals and porosity if not specified
#originate from Bruce Jaffe 9/29/2003
#Made into python by Hui Tang, Virginia Tech, Jan 25, 2013 tanghui@vt.edu
def grading(sl,phi,z,intv,porosity,ni):
##============================================Calculate part============================================
    numrow_sl=len(sl)
    numcol_sl=len(sl[0]) #get matrix sl size
    hitsize=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with size of particles hitting bed
    hitload=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with concentration of particles hitting bed
    hitheight=zeros(shape=(numrow_sl,numcol_sl)) #initiate array with height of particles hitting bed
    if((len(z)!=numrow_sl) or (len(phi)!=numcol_sl)):
       print('Number of size classes or number of elevations do not match profiles')
    else:
       z=z.reshape(len(z),1) #z must be a column vector
       #calculate time for sediment to reach reach bottom
       setvel=0.01*tubesetvel(phi,rhos=2.65,T=21.6,S=0)#calculate settling velocity (m/s)
       setvel=setvel.reshape(1,len(setvel))
       numrow_setvel=len(setvel)
       numcol_setvel=len(setvel[0])
       if(numrow_setvel>numcol_setvel):
           setvel=setvel #make sure setvel is a row vector
       settlingtime=zeros(shape=(numrow_sl,numcol_setvel))
       stime=zeros(shape=(numrow_sl*numcol_setvel,3))
       k=0
       #print int(numcol_setvel)
       for i in range(numrow_sl):
           for j in range(numcol_setvel):
               settlingtime[i][j]=z[i]/setvel[0][j]
                      #stime[k][0]=settlingtime[i][j]
       stime[:,0]=settlingtime.reshape(int(numrow_sl*numcol_setvel),)
       #print shape(stime)
#***********************************sort by time it takes to hit bottom*****************************************
       stime1=sorted(stime[:,0])
       ind=argsort(stime[:,0],axis=0)
       for i in range(len(ind)):
           stime[k][1]=int(int(ind[i])/int(numcol_setvel))
           stime[k][2]=int(int(ind[i])%int(numcol_setvel))
           k=k+1
#**********************************find size of sediment hitting the bed****************************************
       #phiind=floor(((stime[:,2]-0.001)/numrow_sl))+1 # find size of sediment
       phiind=stime[:,2] # find size of sediment
       hitsize=zeros(len(phiind))
       for i in range(len(phiind)):
           hitsize[i]=phi[phiind[i]]    # keep track of order sediment size hitting the bed 
#***********************************find elevations that sediment started from**********************************
       zind=stime[:,1]
       #print ind
       hitheight=zeros(len(zind))
       for i in range(len(zind)):
           hitheight[i]=z[zind[i]] #keep track of order of height sediment came from
       hitload=zeros(len(zind))
       for i in range(len(zind)):
           hitload[i]=sl[zind[i],phiind[i]]
       thickness=cumsum(hitload/(1-porosity))#cumulative thickness of deposit
       #print sum(thickness)
#********************************************determine interval breaks*******************************************
       thickness2=zeros(len(thickness))
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
           #print indthick
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
#**************************************calculate stats for each interval*****************************************
       st=zeros(len(intbind3)-1)
       last=zeros(len(intbind3)-1)
       sload=zeros(shape=(len(indthick),len(phi)))
       m1=zeros(len(st))
       m3=zeros(len(st))
       m4=zeros(len(st))
       stdev=zeros(len(st))
       fwmedian=zeros(len(st))
       fwsort=zeros(len(st))
       fwmean=zeros(len(st))
       fwskew=zeros(len(st))
       fwkurt=zeros(len(st))
       p=zeros(shape=(len(st),len(phi)))
       s=zeros(shape=(len(st),len(phi)))
       wpc=zeros(shape=(len(st),len(phi)))
       midint=zeros(len(indthick2)-1)
       for i in range(0,len(intbind3)-1):
           st[i]=intbind3[i]
       j=0
       for i in range(1,len(intbind3)):
           last[j]=intbind3[i]
           j=j+1
       #print last,st
       for j in range(0,len(indthick2)-1):
           for k in range(len(phi)): #get cumulative sediment loads for each phi class
               #sload[k][j]=sum(hitload[int(indthick[j]):int(indthick[j+1])])
               for l in range(int(indthick2[j]),int(indthick2[j+1])):
                    if(phi[stime[l][2]]==phi[k]):
                        sload[j][k]=hitload[l]+sload[j][k]
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
       #print phi,midint,m1,stdev,m3,m4,fwmedian,fwmean,fwsort,fwskew,fwkurt
#=================================================plot figure part======================================
##============================plot figure for Grain Statistics, 1st Moment for layers===================
##This figure is desinged to show the relationship between location and first hit grain size
       fig1=figure()
       fname1='Grain Statistics, 1st Moment for layers'
       ax1=subplot(111)
       ax1.set_xlabel('1st Moment (phi)')
       ax1.set_ylabel('Midpoint of layer (m)')
       ax1.set_title(fname1)
       ax1.plot(m1,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/m1/'+fname1+' Sample%i.png'%ni,dpi=100)
##============================plot figure for Standard Deviation for layers=============================
##This figure is designed to show the relationship between location and standard deviation 
       fig2=figure()
       fname2='Grain Statistics, Standard Deviation for layers'
       ax2=subplot(111)
       ax2.set_xlabel('Standard Deviation (phi)')
       ax2.set_ylabel('Midpoint of layer (m)')
       ax2.set_title(fname2)
       ax2.plot(stdev,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/stdev/'+fname2+' Sample%i.png'%ni,dpi=100)
##=================================plot figure for FW mean for layers===================================
##This figure is designed to show the relationship between location and mean
       fig3=figure()
       fname3='Grain Statistics, FW mean for layers'
       ax3=subplot(111)
       ax3.set_xlabel('FW Mean (phi)')
       ax3.set_ylabel('Midpoint of layer (m)')
       ax3.set_title(fname3)
       ax3.plot(fwmedian,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwmedian/'+fname3+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW sorting for layers==================================
##This figure is designed to show the relationship between location and sorting
       fig4=figure()
       fname4='Grain Statistics, FW sorting for layers'
       ax4=subplot(111)
       ax4.set_xlabel('Sorting (phi)')
       ax4.set_ylabel('Midpoint of layer (m)')
       ax4.set_title(fname4)
       ax4.plot(fwsort,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwsorting/'+fname4+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Skewness for layers=================================
##This figure is designed to show the relationship between location and skewness
       fig5=figure()
       fname5='Grain Statistics, FW Skewness for layers'
       ax5=subplot(111)
       ax5.set_xlabel('Skewness')
       ax5.set_ylabel('Midpoint of layer (m)')
       ax5.set_title(fname5)
       ax5.plot(fwskew,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwskew/'+fname5+' Sample%i.png'%ni,dpi=100)
##===============================plot figure for FW Kurtosis for layers=================================
##This figure is designed to show the reltaionship between location and kurtosis
       fig6=figure()
       fname6='Grain Statistics, FW Kurtosis for layers'
       ax6=subplot(111)
       ax6.set_xlabel('Kurtosis')
       ax6.set_ylabel('Midpoint of layer (m)')
       ax6.set_title(fname6)
       ax6.plot(fwskew,midint,'-')
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/fwkurt/'+fname6+' Sample%i.png'%ni,dpi=100)
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
       savefig('/Users/tanghui/Dropbox/project1/project1/model01/model1_2/figure/hittime/'+fname7+' Sample%i.png'%ni,dpi=100)

