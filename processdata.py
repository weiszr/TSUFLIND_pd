from pylab import *
from function import *
from ReadCSV import *
from output2CSV import *
import os
import csv

filename='data_draw.csv'
filename1='data_draw1.csv'
Data=readCSV2(filename,separator=';')
se=Data[:,0]
data_x=Data[:,1]
data_y=Data[:,2]
for j in range(2):
    n=0
    for i in range(len(se)-1):
        if(se[i]==nan):
            if(se[i+1]!=nan):
                if(i!=len(se)-1):
                    se[i]=(se[i-1]+se[i+1])/2
                else:
                    se[i]=se[i-1]
            else:
                se[i]=(se[i-1]+se[i+2])/2
        if(se[i]!=0):
            n=n+1
data1_se=zeros(n)
data1_x=zeros(n)
data1_y=zeros(n)
for i in range(n):
    data1_se[i]=se[i]
    data1_x[i]=data_x[i]
    data1_y[i]=data_y[i]
            #data1_se=append(data1_se,se[i])
            #data1_x=append(data1_x,data_x[i])
            #data1_y=append(data1_y,data_y[i])
fp=output2CSV2(filename1,data1_se,data1_x,data1_y,'Mean_grainsize','x_location','y_location')
