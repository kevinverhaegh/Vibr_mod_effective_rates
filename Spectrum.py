import numpy as np
import matplotlib.pyplot as plt
import CRUMPET

def get_decay_table(crm, init_state, final_state, indx_f_state, Te, ne):
    fv = crm.steady_state(Te,ne,Ti=Te,plot=False,dt=True)[indx_f_state]
    # fv_B1 = fv_X1/np.sum(fv_X1)


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

def spectrum(input_crm, Te, ne):
    crm = CRUMPET.Crumpet(input_crm)
    crm.source[2] = 1e-10

    indx_d3 = np.append(92,np.arange(93,107))
    indx_X1 = np.append(2,np.arange(3,17))
    
    n=20
    fv_X1 = crm.steady_state(Te,ne,Ti=Te,plot=False,dt=True, gl=True)[indx_X1]
    fv_X1 = fv_X1/np.sum(fv_X1)

    fv_d3=crm.steady_state(Te,ne,Ti=Te,plot=False,dt=True)[indx_d3]
    fv_d3 = fv_d3/np.sum(fv_d3) 


    import fc_mapping
    fv_d3_mapped = fc_mapping.get_upper(ground=fv_X1)
    fv_d3_mapped = fv_d3_mapped/np.sum(fv_d3_mapped)



    coeffs = np.zeros((15,15))
    for i in range(15):
        for j in range(15):
            coeffs[i,j] = crm.reactions['USER']['COEFFICIENT'][f'd3Pu_a3Sg{i}to{j}'].coeffs
    

    E_d3=np.zeros(15)
    E_a3=np.zeros(15)
    for i in range(15):
        coeffs[i,:] = coeffs[i,:]*fv_d3[i]
        E_a3[i] = crm.species[f'H2(n=a3,v={i})']['V']
        E_d3[i] = crm.species[f'H2(n=d3,v={i})']['V']
    
    table = np.zeros((2,15**2))
    k=0
    for i in range(15):
        for j in range(15):
            table[0,k] = E_d3[i]-E_a3[j]
            table[1,k] = coeffs[i,j]
            k+=1
    
    # hbar = 1.0546e-34
    # e = 1.602e-19
    # table[0,:] = 1e9*2*np.pi*hbar/(table[0,:]*e)

    return fv_X1, fv_d3, fv_d3_mapped, np.array(sorted(np.transpose(table), key=lambda row: row[0]))

def rad_power(crm, Te, ne):
    indx_X1 = np.append(2,np.arange(3,17))
    indx_B1 = np.append(17,np.arange(18,32))
    indx_C1 = np.append(32,np.arange(33,47))

    table = get_decay_table(crm,'C1Pu','X1Sg',indx_C1,Te,1e19)
    
    power_loss = sum(table[0,:]*table[1,:])

    return power_loss



input_crm = 'input_res.dat'
fv_X1, fv_d3, fv_d3_mapped, spect = spectrum(input_crm, 100,1e19)
plt.plot(fv_X1, label='H2[v] in ground state')
plt.plot(fv_d3_mapped, label = 'H2[v] in Fulcher (mapped from X1)')
plt.plot(fv_d3, label='H2[v] in Fulcher (to be compared against experiment)')
plt.xlabel('Vibrational state')
plt.ylabel('H2[v]/H2[0]')
plt.yscale('log')
plt.legend()
plt.show()


plt.plot(spect[:,0], spect[:,1])
plt.xlabel('Energy (eV)')
plt.ylabel('Intensity (a.u.)')
plt.show()

# Tev = np.linspace(1,100, 100)
# p = np.zeros(len(Tev))

# crm = CRUMPET.Crumpet(input_crm)
# crm.source[2] = 1e-10

# for i in range(len(Tev)):
#     p[i] = rad_power(crm,Tev[i],1e19)
#     print(i)

# plt.plot(Tev,p)
# plt.show()