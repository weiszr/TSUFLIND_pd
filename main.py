from pylab import *
from src.function import *
from src.output2CSV import *
from src.concentration import *
from src.Moore_model import *
from src.RSEcalcu import *
from src.thick_RSE import *
from src.ReadCSV import *
#from thickness import *
from src.source_distribution import *
from src.sousby_model import *
from src.TsusedForm import *
from src.TsusedMod import *
from src.preproc import *
from src.TsusedMod1 import *
from src.TsusedMod2 import *
import csv
import os
import matplotlib.colors as colors
import matplotlib.cm as cmx
import matplotlib.ticker as ticker
import shutil
#TODO:
#   [ ] save output data in data directory
#   [ ] create a new way to read in data
#   [ ] make input files and data directories a commandline input
def main():
    separator = ':'
    separator = ':'
    path = os.getcwd()
    print ("The current working directory is %s" % path)
    global path1
    path1 = '{t1}/data'.format(t1=path)
    if os.path.isdir(path1) !=False:
        shutil.rmtree(path1)
    try:
        os.mkdir(path1)
    except OSError:
        print ("Creation of the directory {t1} already exists".format(t1=path1))
    else:
        print ("Successfully created the directory {t1}".format(t1=path1))
    with open('parameter_P14abc.txt', 'r') as f:
        for line in f:
            try:
                s = line.split('=')
                n = s[0]
                v = s[1]  # remove '\n'
                if os.name != 'nt':
                    v = v[0:len(v)-1]  # remove '\r' for Unix-system
                # re-name input variables to match input list
                # model parameters
                if n == 'number_of_iterations':
                    nit = float(v)
                elif n == 'Von_Karmen_constant':
                    vk = float(v)
                elif n == 'number_of_vertical_bins':
                    nz = float(v)
                elif n == 'resuspension_coefficient':
                    g0 = float(v)
                elif n == 'concentration_convergence_factor':
                    concconvfactor = float(v)
                elif n == 'size_convergence_factor':
                    sizeconvfactor = float(v)
                elif n == 'shear_velocity':
                    ustrc = float(v)
                elif n == 'bed_roughness':
                    zotot = float(v)
                elif n == 'settling_velocity_type':
                    setType = float(v)
                elif n == 'eddy_viscosity_shape':
                    var = float(v)
                elif n == 'salinity':
                    sal = float(v)
                elif n == 'water_temperature':
                    wtemp = float(v)
                elif n == 'sediment_density':
                    Rhos = float(v)
                elif n == 'bed_concentration':
                    Cb = float(v)
                # data
                elif n == 'Dimensionless_shear_velovity':
                    Tstarcr = float(v)
                elif n == 'Gravity_acceleration':
                    g = float(v)
                elif n == 'Rouse_number':
                    P = float(v)
                elif n == 'water_run_up':
                    Rz = float(v)
                elif n == 'slope':
                    m = float(v)
                elif n == 'sediment_deposit_rate':
                    q = float(v)
                elif n == 'Number_sample':
                    N = int(v)
                # Input
                elif n == 'Largest_distance':
                    l = float(v)
                elif n == 'Filename':
                    filename = str(v)
                elif n == 'Depth':
                    h = float(v)
                elif n == 'Sediment_thickness':
                    th = float(v)
                elif n == 'Thickness_file':
                    th_file = str(v)
                elif n == 'Filename_Sousby':
                    filename_sousby = str(v)
                elif n == 'layer_thickness':
                    thb = float(v)

            except IndexError as e:
                continue
            except ValueError as e:
                continue
    para = preproc(filename)  # preprocess data from field data document
    Dl = 0.0
    Dm = 2.6860
    Ds = 4.0
    V = 7.743
    Nc = para[3]
    rho = para[4]  # Density of sea water
    H = 6.0  # Max water depth
    Dl1 = phi2mm(-Dl)*0.001  # Change to m in order to use in moore's model
    Dm1 = phi2mm(-Dm)*0.001  # Change to m in order to use in moore's model
    D = Dm
    Data_source = source_distribution(Dl, Ds, D, Nc)
    Se = Data_source[0]  # Grain size of original sediment source
    Fr = Data_source[1]  # Fraction for each grain size of sediment source
    Data_concentration = inconcentration(Se, Fr, V, H, Cb, rho)
    con = Data_concentration
    data2 = sousby(Se, con, V, H, Nc)
    thD = thcalrse(th_file, data2, N)
    x = linspace(0, 350, 16)
    h = x*(Rz-H)/(Rz*m)-x/m+H
    result = zeros(shape=(N, 3, 4))
    result1 = zeros(shape=(N, 3, 2, 100))
    for i in range(int(len(data2))):
        intv = linspace(data2[i]/3.0, data2[i], 3)
        name = 'sample%02d' % i+'.csv'
        name1 = 'sample%02d' % i
        data5 = Tsusedform(name, intv, h[i], x[i])
        for j in range(len(intv)):
            name2 = path1+'/'+name1+"_suspended_sample%02d" % j+".csv"
            data6 = readCSV1(name2, separator=';')
            result[i, j, :] = Tsusedmod1(name2)
            speedfile = Tsusedmod2(name2)
            result1[i, j, 0, :] = speedfile[:, 0]
            result1[i, j, 1, :] = speedfile[:, 1]
            name3 = path1+ '/result_'+'sample%02d' % i+'.txt'
            fp = open(name3, 'w')
            fp.write("%s \n" % result[i, j, :])
        close()
    fg1 = figure(1, figsize=(8, 30))
    p1, = plot(x, result[:, 0, 0])
    p2, = plot(x, result[:, 0, 1])
    title("Velocity Ditribution")
    fname1 = 'Velocity Ditribution'
    xlabel("Location")
    ylabel("Veloctiy")
    legend([p1, p2], ["Max Velocity", "Average Velocity"], loc=3)
    savefig(fname1+' Sample%i.png' % i, dpi=100)
    fg2 = figure(2, figsize=(8, 30))
    p3, = plot(x, result[:, 0, 2])
    p4, = plot(x, result[:, 0, 3])
    fname2 = "Froude Number Ditribution"
    title("Froude Number Ditribution")
    xlabel("Location")
    ylabel("Froude Number")
    legend([p3, p4], ["Max Froude Number", "Average Froude Number"], loc=4)
    savefig(fname2+' Sample%i.png' % i, dpi=100)




if __name__ == '__main__':
    main()
#    main('conv')
#    main_ind()
