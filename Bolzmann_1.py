def vibr_dist(input_crm,iso_mass=1, Te_max=100,Te_reso=int(1e3),Te_min=0.1):
    import CRUMPET
    import numpy as np

    print('Running DEMO with input:')
    print(input_crm)

    indx_H2v = np.append(2,np.arange(3,17))

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    Tev = np.linspace(Te_min,Te_max,Te_reso)
    Tiv = Tev/iso_mass # Assume Ti=Te
    ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)
    crm.source[2] = 1e-100 #add small source for numerical stability (only needed if reactions that dissociate are included)

    #compute vibrational distribution H2
    fv_H2 = np.zeros([15,len(Tev)])

    #calculate vibrational distribution using Tiv = Tev/iso_mass
    for i in range(0,len(Tev)):
        fv_H2[:,i]=crm.steady_state(Tev[i],ne[i],Ti=Tiv[i],plot=False,dt=True)[indx_H2v]
        print(i)

    #normalise vibrational distribution by dividing the distribution values to the sum of the distribution
    fv_H2 = fv_H2/(np.sum(fv_H2,axis=0)[None,:])

    return fv_H2, Tev

def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E

def fit_eval(coeffs, x):
    o = np.zeros(np.shape(x))
    for i in range(0,len(coeffs)):
        o = o + coeffs[i]*x**i
    return np.exp(o)

Te_min = 0.1
Te_max = 100
Te_reso = int(1e2)

input_crm = 'input_false.dat'
fv, Tev = vibr_dist(input_crm,iso_mass=2, Te_reso=Te_reso)

# input_crm = 'input_ichi.dat'
# fv_ichi, Tev = vibr_dist(input_crm, iso_mass=2, Te_reso=Te_reso)

import numpy as np
import matplotlib.pyplot as plt

def err_bolzmann(fv_H2, T):
    x = 0.5+np.linspace(0,14,15)

    err = np.zeros(len(T))
    coeff = np.zeros(len(T))

    for i in range(len(T)):
        fit = np.flip(np.polyfit(x,np.log(fv_H2[:,i]), 1))
        err[i] = np.sqrt(np.sum(((fit_eval(fit,x)-fv_H2[:,i])/fit_eval(fit,x))**2/len(T)))
        coeff[i] = fit[1]
    return err, coeff


# err_ichi, coeff_ichi = err_bolzmann(fv_ichi, Tev)
err, coeff = err_bolzmann(fv, Tev)
E = D2_energies('Fantz/Table 1 Vib Eigenvalues/X1_EV.txt')


# Plot of Ichihara fit coefficients
# plt.plot(Tev, coeff_ichi, label='Ichihara')
plt.plot(Tev, coeff, label='H2VIBR')
plt.ylabel('Decay constant')
plt.xlabel('Temperature (eV)')


# Plot of the rms deviation from exponential fit
plt.figure()
plt.yscale('log')
plt.plot(Tev,err, label='H2VIBR')
# plt.plot(Tev,err_ichi, label='Ichihara')
plt.ylabel('RMS error')
plt.xlabel('Temperature (eV)')
plt.legend()
plt.title("Deviation from Bolxmann distribution as a function of temperature")


# Plot of distribution at given temperature
i = 90
plt.figure()
plt.yscale('log')
plt.plot(E,fv[:, i], 'o', label='Calculated data')

y = np.log(fv[:,i])
c = np.flip(np.polyfit(E,y, 1))

fit = fit_eval(c,E)
plt.plot(E,fit, '--', label='Exponential fit')
plt.title('T = %1.2f eV' %Tev[i])
plt.xlabel('Vibrational level')
plt.ylabel('Fractional abundance')
plt.legend()
plt.show()

print('hello')

