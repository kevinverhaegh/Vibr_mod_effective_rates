import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate

e = 1.6e-19  #Elementary charge
me = 9.11e-31 #electron mass
kb = 1.38e-23 #Bolzmann constant


def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E

def eval_1D(coeff,T):
    o = np.zeros(np.shape(T))
    for i in range(0,len(coeff)):
        o = o + coeff[i]*(np.log(T)**i)
    return 1e-6*np.exp(o)

def test_integration(input_crm):
    print('Running DEMO with input:')
    print(input_crm)
    import CRUMPET
    

    indx_H2v = np.append(2,np.arange(3,17))

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    crm.source[2] = 1e10 #add small source for numerical stability (only needed if reactions that dissociate are included)

def calc_rates(E,sigma, Tev):
    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def extract_table(lines, initial_state, final_state):

    for i,line in enumerate(lines): 
        if f'PROCESS: E + D2(X,v={initial_state}) -> E + D2(X,v={final_state}), Excitation' in line:
            start_index = i+5

    table = []
    for line in lines[start_index:]:
        try: 
            number = [float(line.split()[0]), float(line.split()[1])]
            table.append(number)
        except ValueError:
            break
    
    table = np.array(table)
    return table


def fit_temp(Tev, E, i_state, f_state):

    with open('rates/Laporta/vibr_trans/Cross section.txt', 'r') as file:
        lines = file.readlines()
        table = extract_table(lines,0,1)

    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1], fill_value = 'extrapolate')
    arr = np.transpose([np.linspace(max(logtable[:,0]),np.log(500),100), o(np.linspace(max(logtable[:,0]),np.log(500),100))])
    table = np.exp(np.concatenate((logtable,arr)))

    # Calculating rates as function of temperature
    rates_ex = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates_ex[i] = calc_rates(table[:,0],table[:,1],Tev[i])
    
    
    # Calculating deexcitation rates using detailed balance 
    rates_deex = rates_ex*np.exp((E[f_state]-E[i_state])/Tev)

    fit_ex = np.flip(np.polyfit(np.log(Tev),np.log(rates_ex/1e-6),8))
    fit_deex = np.flip(np.polyfit(np.log(Tev),np.log(rates_deex/1e-6),8))

    return rates_ex, rates_deex, fit_ex, fit_deex

Tev = np.linspace(0.2,100, 100)
Tiv = Tev # Assume Ti=Te
ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)

E = D2_energies('Fantz/Table 1 Vib Eigenvalues/X1_EV.txt')
r_ex,r_deex,fit_ex,fit_deex = fit_temp(Tev,E,3,4)
r = eval_1D(fit_deex,Tev)
r_ = eval_1D(fit_deex,Tev)

plt.plot(Tev,r_ex)

plt.plot(Tev,r_deex)

plt.show()


print('Done')

