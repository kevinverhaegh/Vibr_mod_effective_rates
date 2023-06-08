import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import CRUMPET

e = 1.60217663e-19  #Elementary charge
me = 9.1093837e-31 #electron mass
mp = 1.67262192e-27 #proton mass
kb = 1.380649e-23 #Bolzmann constant

def eval_1D(coeff,T):
    o = np.zeros(np.shape(T))
    for i in range(0,len(coeff)):
        o = o + coeff[i]*(np.log(T)**i)
    return 1e-6*np.exp(o)

# Get cross sections from H2VIBR
E = 10**np.linspace(np.log10(10.2),np.log10(5000),10000)
X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr.tex', 'HYDHEL' : '/rates/HYDHEL.tex'})
# sigma = eval_1D(X.reactions['H2VIBR']['H.1']['2.0L2'],E) 
sigma = 1e2*eval_1D(X.reactions['HYDHEL']['H.1']['2.1.1'][:9],E) # Calculate in m^2
sigma[E<10.8] = 1e-4*1e-19
sigma[sigma>(1e-4*7.5e-17)] = 1e-4*7.5e-17

table = []
table.append(E)
table.append(sigma)
table = np.transpose(np.array(table))

plt.loglog(table[:,0], table[:,1])
plt.show()

print('hello')

def calc_rates(E, sigma, Tev):
    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return rate

def fit_temp(table, Tev):

    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1], fill_value = 'extrapolate')
    arr = np.transpose([np.linspace(max(logtable[:,0]),np.log(500),100), o(np.linspace(max(logtable[:,0]),np.log(500),100))])
    table = np.exp(np.concatenate((logtable,arr)))

    # Calculating rates as function of temperature
    rates = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates[i] = 1e-6*calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!

    fit = np.flip(np.polyfit(np.log(Tev),np.log(rates/1e-6),8))

    return fit, rates

Tev = 10**np.linspace(np.log10(1.26), np.log10(100), 100)
rates = np.zeros(len(Tev))
for i in range(len(Tev)):
    rates[i] = calc_rates(table[:,0],table[:,1],Tev[i]) 

# Make comparison
# y = eval_1D(X.reactions['H2VIBR']['H.2']['2.0L2'],Tev)
y = eval_1D(X.reactions['HYDHEL']['H.2']['2.1.1'],Tev)

plt.loglog(Tev,rates,'b', label='Calculated using cross section')
plt.loglog(Tev,y,'r', label = 'From rate data')
plt.legend()
plt.show()

# I have absolutely no idea how the units work in H2VIBR