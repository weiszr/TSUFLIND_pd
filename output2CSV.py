'''
Created on 2012-11-24;
===============================================Output data to CSV document==================================
#======================================================progressBar==========================================
#=================================================output2CSV function one===================================
#=================================================output2CSV function two===================================
@author: Hui Tang tanghui@vt.edu
'''
#===========================================================================================================
#===========================================================================================================
#==================================================Output Data to CSV Script================================
#===========================================================================================================
#===========================================================================================================
#This Script is used to output data to CSV document
from pylab import *
import csv
import sys
# ================================================= progressBar ============================================
# CLASS NAME: DLLInterface
# Author: Larry Bates (lbates@syscononline.com)
# Written: 12/09/2002
# Released under: GNU GENERAL PUBLIC LICENSE
class progressBar: 
    def __init__(self, finalcount, title='Progress' ,progresschar='>'):
        import sys
        self.finalcount=finalcount; self.blockcount=0
        if not progresschar: self.block=chr(178)
        else:                self.block=progresschar
        # Get pointer to sys.stdout so I can use the write/flush
        self.f=sys.stdout
        # If the final count is zero, don't start the progress gauge
        if not self.finalcount : return
        self.f.write('\n'+title+'\n[')
        return
    def progress(self, count):
        # Make sure I don't try to go off the end (e.g. >100%)
        count=min(count, self.finalcount)
        # If finalcount is zero, I'm done
        if self.finalcount:
            percentcomplete=int(round(100*count/self.finalcount))
            if percentcomplete < 1: percentcomplete=1
        else:            percentcomplete=100
        #print "percentcomplete=",percentcomplete
        blockcount=int(percentcomplete/2)
        #print "blockcount=",blockcount
        if blockcount > self.blockcount:
            for i in range(self.blockcount,blockcount):
                self.f.write(self.block)
                self.f.flush()
        if percentcomplete >= 100: self.f.write("]")
        self.blockcount=blockcount
        return
# ============================================== output2CSV =============================================
# this function output given data to csv file
def output2CSV(name,data1,data2,text1,text2):
    data=zeros([len(data1),2]);
    data[:,0]=data1; data[:,1]=data2;
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(name,'wt') as f:
        writer=csv.writer(f,delimiter=';') 
        writer.writerow((text1,text2)) # write description
        for row in data: writer.writerow((row[0],row[1]));pb.progress(n+1);n+=1;
    return data
# ============================================= output2CSV1 =============================================
# this function output given data to csv file
def output2CSV1(name,data1,data2,data3,data4,data5,text1,text2,text3,text4,text5):
    data=zeros([len(data1),5]);
    data[:,0]=data1; data[:,1]=data2;data[:,2]=data3;data[:,3]=data4;data[:,4]=data5
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(name,'wt') as f:
        writer=csv.writer(f,delimiter=';') 
        writer.writerow((text1,text2,text3,text4,text5)) # write description
        for row in data: writer.writerow((row[0],row[1],row[2],row[3],row[4]));pb.progress(n+1);n+=1;
    return data
# ============================================= output2CSV2 =============================================
# this function output given data to csv file
def output2CSV2(name,data1,data2,data3,text1,text2,text3):
    data=zeros([len(data1),3]);
    data[:,0]=data1; data[:,1]=data2;data[:,2]=data3;
    pb=progressBar(len(data1),'# Outputting to File['+name+']:',">"); n=0;
    with open(name,'wt') as f:
        writer=csv.writer(f,delimiter=';') 
        writer.writerow((text1,text2,text3)) # write description
        for row in data: writer.writerow((row[0],row[1],row[2]));pb.progress(n+1);n+=1;
    return data
