from pylab import *
#=============================================function sedstats===========================================
#This function is designed to compute moment measures and Folk and Ward statistics from sediment grain size
#weight % and phi size classes.
#Inputs: phi:grain size (phi) for each size class
#        weight:weight (or weight percent, vol concentration) for size classes
#Calculated: phimdpt:phi midpoints used for moment measures for tsunami samples, [-1.125:.25:10 12]
#Outputs:    m1:mean grain size (phi), 1st moment
#            stdev:standard deviation (phi) similar to sorting, is equal to square root of 2nd moment
#            m3:third moment
#            m4:fourth moment
#            fwmedian:Folk and Ward medium grain size (phi), 50th percentile
#            fwmean:Folk and Ward mean grain size (phi)
#            fwsort:Folk and Ward sorting (phi)
#            fwskew:Folk and Ward skewness
#            fwkurt:Folk and Ward Kurtosis
def sedstats(phi,weight):
    if(len(phi)==1):
        print('Moment measures not valid because there is only one size class')
        #phi=
        #create a cumulative weight with 0% weight in it to allow calculation
    else:
        mpt1=phi[0]+0.5*(phi[0]-phi[1])
        phimdpt=zeros(len(phi))
        phimdpt[0]=mpt1
        phimdpt[1:len(phi)]=phi[1:len(phi)]+0.5*(phi[0:len(phi)-1]-phi[1:len(phi)])
        phimdpt=transpose(phimdpt)
        #p=phimdpt*weight
        #print p
        #print phimdpt
        #print weight
        #calculate the 1st moment (mean grain size)
        m1=sum(phimdpt*weight)/sum(weight)
        dev=phimdpt-m1
        #calculate the standard deviation, similar to sorting, square root of the 2nd moment
        var=sum(weight*(dev**2))/sum(weight)
        stdev=sqrt(var)
        #calculate the third moment
        m3=sum(weight*((dev/stdev)**3))/sum(weight)
        #calculate the fourth moment
        m4=sum(weight*((dev/stdev)**4))/sum(weight)
        #Use phi intervals, not midpoints to calculate Folk and Ward stats- for tsunami samples use
        cw=100*cumsum(weight)/sum(weight)#calculate normalized cumulative weight percent
        cp=[5,16,25,50,75,84,95]
        if(cw[0]>=5):
            print('Folk and Ward Statistics suspect because first size class has >=5% cumulative %')
            cw1=zeros(len(cw)+1)
            cw1[1:len(cw)+1]=cw[:] #create a cumulative weight with 0% weight in it to allow calculation
        #calculate cumulative percents used in Folk and Ward statistics using
        #linear interpolation (could use different interpolation scheme)
        cumphi=zeros(7)
        cw1=zeros(int(len(phi)))
        for i in range(7):
            k=0
            while(cw[k]<=cp[i]):
                cw1[k]=cw[k]
                k=k+1
            lp=argmax(cw1)
            slp=(cp[i]-cw[lp])/(cw[lp+1]-cw[lp])
            cumphi[i]=phi[lp]+slp*(phi[lp+1]-phi[lp])
        #make some names that have meaning
        phi5=cumphi[0]
        phi16=cumphi[1]
        phi25=cumphi[2]
        phi50=cumphi[3]
        phi75=cumphi[4]
        phi84=cumphi[5]
        phi95=cumphi[6]
        #calculate Folk and Ward stats
        fwmedian=phi50
        fwmean=(phi16+phi50+phi84)/3
        fwsort=(phi84-phi16)/4.0+(phi95-phi5)/6.6
        fwskew=(phi16+phi84-2*phi50)/(2*(phi50-phi16))+(phi5+phi95-2*phi50)/(2*(phi95-phi5))
        fwkurt=(phi95-phi5)/(2.44*(phi75-phi25))
        wp=weight*phi
        #print len(weight[0,:])
        #percen=zeros(len(weight[0,:]))
        #print weight[0,39]
        #for j in range(len(weight[0,:])):
            #percen[j]=weight[0,j]
        #percen=percen.reshape(len(percen[:,0]))
        #print percen
        #sumw=sum(weight)
        #wp=zeros(len(phi))
        #percen=zeros(len(phi))
        #for j in range(len(phi)):
            #wp[j]=percen[j]*phi[j]
        fwmean=sum(wp)/100
        sedstat=[m1,m3,m4,stdev,fwmedian,fwsort,fwmean,fwskew,fwkurt]
        return(sedstat)
        
