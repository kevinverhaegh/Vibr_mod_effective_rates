#%%
import numpy as np
import matplotlib.pyplot as plt
import CRUMPET
import pickle
#%%


def get_decay_table(crm, init_state, final_state, indx_i_state, Te, ne):
    fv = crm.steady_state(Te,ne,Ti=Te,plot=False,dt=True)[indx_i_state]
    # fv = fv/np.sum(fv)


    coeffs = np.zeros((15,15))

    for i in range(15):
        for j in range(15):
            coeffs[i,j] = crm.reactions['USER']['COEFFICIENT'][init_state+'_'+final_state+f'{i}to{j}'].coeffs

    
    E_i=np.zeros(15)
    E_f=np.zeros(15)
    for i in range(15):
        coeffs[i,:] = coeffs[i,:]*fv[i]
        E_i[i] = crm.species['H2(n='+init_state[:-2]+f',v={i})']['V']
        E_f[i] = crm.species['H2(n='+final_state[:-2]+f',v={i})']['V']

    table = np.zeros((2,15**2))
    k=0
    for i in range(15):
        for j in range(15):
            table[0,k] = E_i[i]-E_f[j]
            table[1,k] = coeffs[i,j]
            k+=1
    return table

def spectrum(input_crm, i_state, f_state, indx_i_state, Te, ne):

    crm = CRUMPET.Crumpet(input_crm)
    crm.source[2] = 1e-10
    # pickle.dump(crm.getM(10,1e19)[0],open('m.pkl','wb'))
    

    indx_d3 = np.append(92,np.arange(93,107))
    indx_X1 = np.append(2,np.arange(3,17))
    indx_mol = np.append(2,np.arange(3,107))
    
    
    fv = crm.steady_state(Te,ne,Ti=Te,plot=False,dt=True)
    n0 = np.sum(fv[indx_mol])

    fv_X1 = fv[indx_X1]/n0
    fv_d3 = fv[indx_d3]/np.sum(fv[indx_d3])

    R_0 = fv/(n0*ne)  # Population coefficients


    import fc_mapping
    fv_d3_mapped = fc_mapping.get_upper(ground=fv_X1)
    fv_d3_mapped = fv_d3_mapped/np.sum(fv_d3_mapped)

    table = get_decay_table(crm,i_state,f_state,indx_i_state, Te, ne)
    
    # hbar = 1.0546e-34
    # e = 1.602e-19
    # table[0,:] = 1e9*2*np.pi*hbar/(table[0,:]*e)

    return fv_X1, fv_d3, fv_d3_mapped, np.array(sorted(np.transpose(table), key=lambda row: row[0]))


def rad_power(crm, i_state, f_state, i_indx, Te, ne):
    indx_X1 = np.append(2,np.arange(3,17))
    indx_B1 = np.append(17,np.arange(18,32))
    indx_C1 = np.append(32,np.arange(33,47))

    table = get_decay_table(crm,i_state,f_state,i_indx,Te,1e19)
    
    power_loss = sum(table[0,:]*table[1,:])
    return power_loss
#%%

indx_X1 = np.append(2,np.arange(3,17))
indx_B1 = np.append(17,np.arange(18,32))
indx_C1 = np.append(32,np.arange(33,47))
indx_EF1 = np.append(47,np.arange(48,62))
indx_a3 = np.append(62,np.arange(63,77))
indx_c3 = np.append(77,np.arange(78,92))
indx_d3 = np.append(92,np.arange(93,107))

input_crm = 'input_fin.dat'
fv_X1, fv_d3, fv_d3_mapped, spect = spectrum(input_crm, 'd3Pu', 'a3Sg',indx_d3, 10,1e19)

#%%

e=1.602e-19
h=6.626e-34
hbar = h/(2*np.pi)
c=3e8
kb = 1.38e-23


# plt.plot(fv_X1, label='H2[v] in ground state')
# plt.plot(fv_d3_mapped, label = 'H2[v] in Fulcher (mapped from X1)')
# plt.plot(fv_d3, label='H2[v] in Fulcher (to be compared against experiment)')
# plt.xlabel('Vibrational state')
# plt.ylabel('H2[v]/H2[0]')
# plt.yscale('log')
# plt.legend()
# plt.show()

# wavelength = 1e9*2*np.pi*c*hbar/(e*spect[:,0])

# plt.plot(spect[:,0], spect[:,1])
# plt.xlabel('Energy (eV)')
# plt.ylabel('Intensity (a.u.)')
# # plt.xlim(400,900)
# plt.show()

Tev = np.linspace(1,100, 10)
ne = 1e19
p = np.zeros(len(Tev))

crm = CRUMPET.Crumpet(input_crm)
crm.source[2] = 1e10



for i in range(len(Tev)):
    p[i] = rad_power(crm,'B1Su', 'X1Sg', indx_B1, Tev[i], ne)
    p[i] +=rad_power(crm,'C1Pu', 'X1Sg', indx_C1, Tev[i], ne)
    print(i)

plt.plot(Tev,p)
plt.show()
# %%
print('hello')