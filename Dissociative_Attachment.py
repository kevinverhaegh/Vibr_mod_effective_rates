
#read vibrational distribution

import csv
# Open the CSV file for reading
with open('VibrDistr_1eV_Holm_Digi_Yac.csv', 'r') as file:
    # Create a CSV reader object
    csv_reader = csv.reader(file)

    # Skip the header row
    next(csv_reader)

    # Initialize empty lists for X and Y values
    x_values = []
    y_values = []

    # Iterate over each row in the CSV file
    for row in csv_reader:
        # Assuming the X and Y columns are in the first and second positions
        x_values.append(float(row[0]))
        y_values.append(float(row[1]))

from scipy.interpolate import interp1d
import numpy as np
f = interp1d(x_values,np.log(y_values),bounds_error=False,fill_value='extrapolate')
import numpy as np
v = np.linspace(0,20,21)

vdist = np.exp(f(v))
vdist = vdist/sum(vdist)

#fit with vibrational temperature

import CRUMPET
X = CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr_custom.tex'})
tev = np.linspace(0.1,2,100)

def eval_1D(coeff, T):
    o = np.zeros(np.shape(T))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(T) ** i)
    return 1e-6 * np.exp(o)

da_rate = np.zeros(np.shape(tev))

for i in range(0,len(tev)):
    for j in range(0,21):
        da_rate[i] = da_rate[i] + eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(j)+'L4'], tev[i])*vdist[j]

cx_rate = np.zeros(np.shape(tev))

for i in range(0,len(tev)):
    for j in range(0,15):
        cx_rate[i] = cx_rate[i] + eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(j)+'Q6'], tev[i]/2)*vdist[j]


def eval_2D(coeff, T, n):
    '''
    Function to evaluate 2D polynomial logarithmic fits of reaction rate coefficients of the type in Eirene data files (e.g. AMJUEL).
    In this case, the rate coefficients have both a temperature and density dependence.

        Parameters
        ----------
        coeff: numpy array
            Fit constants of the rate coefficients. (For Eirene data, this is an 8x8 array)
        T: float
            Evaluation temperature in eV (may also be an array of temperatures).
        n: float
            Evaluation density in m^-3 (may also be an array of densities).

        Output
        ------
        Reaction rate coefficient(s) in m^3s^-1.
    '''

    o = np.zeros(np.shape(T))
    for i in range(0, np.shape(coeff)[0]):
        for j in range(0, np.shape(coeff)[1]):
            o = o + coeff[i, j] * (np.log(T) ** i) * (np.log(n * 1e-14) ** j)
    return 1e-6 * np.exp(o)

AMJ = CRUMPET.ratedata.RateData(rates={'AMJUEL' : '/rates/amjuel.tex'})

MAR_H2p = (cx_rate/(eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.14'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.12'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.11'],tev,np.ones(np.shape(tev))*1e19)))*eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.14'],tev,np.ones(np.shape(tev))*1e19)
MAR_Hm = (da_rate/(eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3A'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3B'],tev,np.ones(np.shape(tev))*1e19))) * eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3A'],tev,np.ones(np.shape(tev))*1e19)

import matplotlib.pyplot as plt
plt.figure()
plt.plot(tev,MAR_H2p,'m')
plt.plot(tev,MAR_Hm,'c')
print('...')

import csv
# Open the CSV file for reading
with open('VibrDistr_2eV_Holm_Digi_Yac.csv', 'r') as file:
    # Create a CSV reader object
    csv_reader = csv.reader(file)

    # Skip the header row
    next(csv_reader)

    # Initialize empty lists for X and Y values
    x_values = []
    y_values = []

    # Iterate over each row in the CSV file
    for row in csv_reader:
        # Assuming the X and Y columns are in the first and second positions
        x_values.append(float(row[0]))
        y_values.append(float(row[1]))

from scipy.interpolate import interp1d
import numpy as np
f = interp1d(x_values,np.log(y_values),bounds_error=False,fill_value='extrapolate')
import numpy as np
v = np.linspace(0,20,21)

vdist = np.exp(f(v))
vdist = vdist/sum(vdist)

#fit with vibrational temperature

import CRUMPET
X = CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr_custom.tex'})
tev = np.linspace(0.1,2,100)

def eval_1D(coeff, T):
    o = np.zeros(np.shape(T))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(T) ** i)
    return 1e-6 * np.exp(o)

da_rate = np.zeros(np.shape(tev))

for i in range(0,len(tev)):
    for j in range(0,21):
        da_rate[i] = da_rate[i] + eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(j)+'L4'], tev[i])*vdist[j]

cx_rate = np.zeros(np.shape(tev))

for i in range(0,len(tev)):
    for j in range(0,15):
        cx_rate[i] = cx_rate[i] + eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(j)+'Q6'], tev[i]/2)*vdist[j]


def eval_2D(coeff, T, n):
    '''
    Function to evaluate 2D polynomial logarithmic fits of reaction rate coefficients of the type in Eirene data files (e.g. AMJUEL).
    In this case, the rate coefficients have both a temperature and density dependence.

        Parameters
        ----------
        coeff: numpy array
            Fit constants of the rate coefficients. (For Eirene data, this is an 8x8 array)
        T: float
            Evaluation temperature in eV (may also be an array of temperatures).
        n: float
            Evaluation density in m^-3 (may also be an array of densities).

        Output
        ------
        Reaction rate coefficient(s) in m^3s^-1.
    '''

    o = np.zeros(np.shape(T))
    for i in range(0, np.shape(coeff)[0]):
        for j in range(0, np.shape(coeff)[1]):
            o = o + coeff[i, j] * (np.log(T) ** i) * (np.log(n * 1e-14) ** j)
    return 1e-6 * np.exp(o)

AMJ = CRUMPET.ratedata.RateData(rates={'AMJUEL' : '/rates/amjuel.tex'})

MAR_H2p = (cx_rate/(eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.14'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.12'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.11'],tev,np.ones(np.shape(tev))*1e19)))*eval_2D(AMJ.reactions['AMJUEL']['H.4']['2.2.14'],tev,np.ones(np.shape(tev))*1e19)
MAR_Hm = (da_rate/(eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3A'],tev,np.ones(np.shape(tev))*1e19) + eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3B'],tev,np.ones(np.shape(tev))*1e19))) * eval_2D(AMJ.reactions['AMJUEL']['H.4']['7.2.3A'],tev,np.ones(np.shape(tev))*1e19)

plt.plot(tev,MAR_H2p,'m--')
plt.plot(tev,MAR_Hm,'c--')

print('...')

plt.xlabel('T (eV)')
plt.ylabel('MAR rate (part. m^3 s^-1)')
plt.axis([0.2,2,0,1.5e-16])

plt.savefig('MAR_HmH2p.eps')

print('...')
