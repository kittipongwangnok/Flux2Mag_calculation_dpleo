#PhD reseach 2020
#Kittipong Wangnok, D6010218
#School of Physics, Institute of Science, Suranaree University of Technology

#Import all module
#import sys
#import os
#from scipy import stats
import matplotlib.pyplot as plt
import numpy as np
#from statistics import stdev
#from statistics import mean
np.seterr(divide='ignore', invalid='ignore')
import barycorr
import itertools

#Latex font
import matplotlib as mpl
from matplotlib import rc
plt.rc('text', usetex=True)
plt.rc('font', family='serif',size=16)
#################################################################################
'''
1. Input file: Aperture_size_I
'''
#################################################################################
#Please change the input file
DP_Leo = open("2013_12_27_dpleo_run040_normal_gaussian_II.log",'r').readlines()
N_dpleo = len(DP_Leo)

#################################################################################
'''
2. MJD2BJD conversion
'''
#################################################################################
MJD_time=[]
counts_1 = []
countse_1 = []
counts_2 = []
countse_2 = []
counts_3 = []
countse_3 = []
counts_4 = []
countse_4 = []
counts_5 = []
countse_5 = []
for line in open("2013_12_27_dpleo_run040_normal_gaussian_II.log"):
    li=line.strip()
    #li=line.splt(" ")
    if not li.startswith("#"):
        MJD_time.append(float(li.split(" ")[2]))
        counts_1.append(float(li.split(" ")[15]))
        countse_1.append(float(li.split(" ")[16]))
        counts_2.append(float(li.split(" ")[31]))
        countse_2.append(float(li.split(" ")[32]))
        counts_3.append(float(li.split(" ")[47]))
        countse_3.append(float(li.split(" ")[48]))
        counts_4.append(float(li.split(" ")[63]))
        countse_4.append(float(li.split(" ")[64]))
        counts_5.append(float(li.split(" ")[79]))
        countse_5.append(float(li.split(" ")[80]))

#How to print JD, counts1, countse1, counts2, countse2, counts3, countse3
#print (MJD_time)
n_MJD=len(MJD_time)
#print (n_MJD)

JD = []
for i in range (n_MJD):
    JD_time = np.array(MJD_time)+2400000.5
N_Max = 300
N_JD=len(JD_time)
N_iteration = int(N_JD/N_Max)
N_extra = round((N_JD/N_Max-N_iteration)*N_Max)

BJD_time=[]
#RA and Dec of DP Leo (in deg)
for i in range (N_iteration):
    params = {
        'jd_utc': JD_time[i*N_Max:i*N_Max+N_Max].tolist(),
        'ra': 169.3163492,
        'dec': 17.96157789,
        'lat': 18.573725,
        'lon': 98.482194,
        'elevation': 2457,
    }
    BJD_time.append(barycorr.bvc(**params))
    
for j in range (1):
    params = {
        'jd_utc': JD_time[N_iteration*N_Max:N_iteration*N_Max+N_extra+1].tolist(),
        'ra': 169.3163492,
        'dec': 17.96157789,
        'lat': 18.573725,
        'lon': 98.482194,
        'elevation': 2457,
    }
    BJD_time.append(barycorr.bvc(**params))
    
BJD_time=(list(itertools.chain.from_iterable(BJD_time)))

#################################################################################
'''
3. BJD2Phase calculation
'''
#################################################################################
#We are using the ephemeris of TO and Period obtained by Beuermann et al. (2011)
#BJD(TT) = 2454914.8322920(20) + 0.06236285648(90) E
BJD_0 = 2454914.8322920
P_orb = 0.06236285648
#Arrays
BJD0 = [i for i in range(N_dpleo)]
P_orb1 = [i for i in range(N_dpleo)]
BJD1 = [i for i in range(N_dpleo)]
Cycle1 = [i for i in range(N_dpleo)]
Epoch1 = [i for i in range(N_dpleo)]
Phase1 = [i for i in range(N_dpleo)]

print ("No.\t\t BJD\t\t Phase\t Tar\tTar_err\t Ref\t\tRef_err\tChkI\tChkI_err\tDelta_Mag12\tDelta_Mag12_err")
#Number of BJD_time
n_BJD=len(BJD_time)
#print (n_MJD)
BJD = []
BJD2Phase = []
mean_std = []
phase_int = 0.9
phase_fnl = 1.1
Flux2Mag = []
Flux2Mag_in = []
for i in range (n_BJD):
    BJD = np.array(BJD_time)
    BJD_0 = np.array(BJD_0)
    BJD0[i] = BJD_0
    P_orb = np.array(P_orb)
    P_orb1[i] = P_orb
    BJD = np.array(BJD_time)
    BJD1[i] = BJD
    Cycle = (BJD_time[i] - BJD_0)/P_orb
    Cycle1[i] = Cycle
    Epoch = int(Cycle)
    Epoch1[i] = Epoch
    Phase = Cycle1[i] - Epoch1[0]
    Phase1[i] = Phase
    counts_1 = np.array(counts_1)
    countse_1 = np.array(countse_1)
    counts_2 = np.array(counts_2)
    countse_2 = np.array(countse_2)
    counts_3 = np.array(counts_3)
    countse_3 = np.array(countse_3)
    counts_4 = np.array(counts_4)
    countse_4 = np.array(countse_4)
    counts_5 = np.array(counts_5)
    countse_5 = np.array(countse_5)
    Delta_Mag12 = 2.5*np.log(counts_2/counts_1)
    Delta_Mag12_err = ((counts_2*countse_1) - (counts_1*countse_2))/(counts_1*counts_2)
#    if Phase1[i] >= phase_int and Phase1[i] <= phase_fnl:
#        mean_1 = np.mean(counts_1)
#        std_1 = np.std(counts_1)
#        BJD2Phase.append('%0.0f\t%0.10f\t%0.5f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f' %(i, BJD[i], Phase1[i], counts_1[i], countse_1[i], counts_2[i], countse_2[i], counts_3[i], countse_3[i], counts_4[i], countse_4[i], counts_5[i], countse_5[i]))
 #   print('%0.0f\t%0.10f\t%0.5f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f' %(i, BJD[i], Phase1[i], counts_1[i], countse_1[i], counts_2[i], countse_2[i], counts_3[i], countse_3[i], counts_4[i], countse_4[i], counts_5[i], countse_5[i]))
    print('%0.0f\t%0.10f\t%0.5f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t\t%0.2f\t%0.2f' %(i, BJD[i], Phase1[i], counts_1[i], countse_1[i], counts_2[i], countse_2[i], counts_3[i], countse_3[i], Delta_Mag12[i], Delta_Mag12_err[i]))
    Flux2Mag.append('%0.0f\t%0.10f\t%0.5f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t%0.1f\t\t%0.2f\t%0.2f' %(i, BJD[i], Phase1[i], counts_1[i], countse_1[i], counts_2[i], countse_2[i], counts_3[i], countse_3[i], Delta_Mag12[i], Delta_Mag12_err[i]))
    Flux2Mag_in.append('%0.10f\t%0.2f\t%0.2f' %(BJD[i], Delta_Mag12[i], Delta_Mag12_err[i]))
#print (mean_1, std_1)
#mean_std.append('%0.1f\t%0.10f' %(mean_1, std_1))

'''
4. Save the output file
'''
#Result_1 = BJD2Phase
Result_1 = Flux2Mag
f1 = open("2013_12_27_dpleo_run040_normal_gaussian_II_flux2mag.out", 'w')
for i in range(len(Result_1)):
    f1.write(str(Result_1[i])+ '\n')
f1.close()

#Result_2 = mean_std
Result_2 = Flux2Mag_in
#f2 = open("2013_12_27_dpleo_run040_normal_gaussian_II.std", 'w')
f2 = open("2013_12_27_dpleo_run040_normal_gaussian_II_flux2mag.dat", 'w')
for i in range(len(Result_2)):
    f2.write(str(Result_2[i])+ '\n')
f2.close()

##################################################################################
#'''
#5. Plot Relative flux vs. Orbital phase
#'''
##################################################################################
###Input_1
#InputFileInput_1  = "2013_12_27_dpleo_run040_normal_gaussian_II.bjd"
InputFileInput_1  = "2013_12_27_dpleo_run040_normal_gaussian_II_flux2mag.out"
Data_1   = np.genfromtxt(InputFileInput_1)
Number = Data_1[:,0]
BJD = Data_1[:,1]
Phase1 = Data_1[:,2]
counts_1 = Data_1[:,3]
countse_1 = Data_1[:,4]
counts_2 = Data_1[:,5]
countse_2 = Data_1[:,6]
counts_3 = Data_1[:,7]
countse_3 = Data_1[:,8]
Delta_Mag12 = Data_1[:,9]
Delta_Mag12_err = Data_1[:,10]
#counts_4 = Data_1[:,9]
#countse_4 = Data_1[:,10]
#counts_5 = Data_1[:,11]
#countse_5 = Data_1[:,12]
#
##Flux: Tar
fig, (ax0) = plt.subplots(nrows=1, sharex=True, sharey=True, figsize=(10, 5),tight_layout=True)
ax0.set_xlim(phase_int,phase_fnl)
ax0.set_ylim(22,8)
ax0.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')
ax0.errorbar(Phase1, Delta_Mag12, yerr= Delta_Mag12_err, color='black', fmt='o', markersize=3.0)
#ax0.hlines(y=mean_1, xmin=phase_int, xmax=phase_fnl, colors='red', linestyles='--', #lw=1)
ax0.set_xlabel('Orbital phase')
ax0.set_ylabel('Delta\_Mag')
ax0.grid(linestyle='dotted')
#ax0.legend(loc='lower left')
#plt.savefig('2013_12_27_dpleo_run040_normal_gaussian_II_flux_tar.png', fmt='PNG', #dpi=1000 )
#
#Relative flux: Tar/Chek1
#fig, (ax1) = plt.subplots(nrows=1, sharex=True, sharey=True, figsize=(10, 5), tight_layout=True)
##ax1.set_xlim(0.95,1.05)
#ax1.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')
#ax1.errorbar(Phase1, Relative_flux_13, yerr= Relative_fluxe_13, color='red', fmt='o', markersize=3.0)
#ax1.set_xlabel('Orbital phase')
#ax1.set_ylabel('Relative Flux')
#ax1.grid(linestyle='dotted')
###ax0.legend(loc='lower left')
#plt.savefig('2013_12_27_dpleo_run040_normal_gaussian_Relative_flux_13.png', fmt='PNG', dpi=1000 )
#
#Relative flux: Tar/Chek2
#fig, (ax2) = plt.subplots(nrows=1, sharex=True, sharey=True, figsize=(10, 5), tight_layout=True)
##ax2.set_xlim(0.95,1.05)
#ax2.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')
#ax2.errorbar(Phase1, Relative_flux_14, yerr= Relative_fluxe_14, color='blue', fmt='o', markersize=3.0)
#ax2.set_xlabel('Orbital phase')
#ax2.set_ylabel('Relative Flux')
#ax2.grid(linestyle='dotted')
###ax0.legend(loc='lower left')
#plt.savefig('2013_12_27_dpleo_run040_normal_gaussian_Relative_flux_14.png', fmt='PNG', dpi=1000 )
#
#Relative flux: Tar/Chek3
#fig, (ax3) = plt.subplots(nrows=1, sharex=True, sharey=True, figsize=(10, 5), tight_layout=True)
##ax3.set_xlim(0.95,1.05)
#ax3.tick_params(direction='in', which='both', bottom='on',top='on', right = 'on')
#ax3.errorbar(Phase1, Relative_flux_15, yerr= Relative_fluxe_15, color='orange', fmt='o', markersize=3.0)
#ax3.set_xlabel('Orbital phase')
#ax3.set_ylabel('Relative Flux')
#ax3.grid(linestyle='dotted')
###ax0.legend(loc='lower left')
#plt.savefig('2013_12_27_dpleo_run040_normal_gaussian_Relative_flux_15.png', fmt='PNG', dpi=100 )

plt.show()

