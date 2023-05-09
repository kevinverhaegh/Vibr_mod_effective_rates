def eval_1D(coeff, T):
    import numpy as np
    o = np.zeros(np.shape(T))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(T) ** i)
    return 1e-6 * np.exp(o)


def run_demo(input_crm,iso_mass=1,T_request=1, ichi=False,Te_max=100,Te_reso=int(1e3),Te_min=0.1):
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

    #Get vibrationally resolved molecular CX rates from H2VIBR

    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr_ichi.tex', 'AMJUEL' : '/rates/amjuel.tex'})

    #Get rates as function of Te


    vibr_resolved_CX = np.zeros([15,len(Tev)])
    vibr_resolved_Diss = np.zeros([15,len(Tev)])
    vibr_resolved_Ion = np.zeros([15,len(Tev)])

    states = ['a3Sg', 'B1Su', 'Bp1Su', 'C1Pu', 'c3Pu', 'D1Pu', 'd3Pu', 'e3Su', 'EF1Sg', 'g3Sg', 'GK1Sg', 'H1Sg', 'h3Sg', 'I1Pg', 'i3Pg', 'J1Dg', 'j3Dg']
    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L2'],Tiv)

        # for j in range(len(Tev)):
        #     vibr_resolved_Ion[i,j] = crm.reactions['MCCCDB']['ion'][f'{i}'].integrated(Tev[j],0)

            # for state in states:
            #     vibr_resolved_Diss[i,j] += crm.reactions['MCCCDB']['diss']['X_'+state+f'{i}'].integrated(Tev[j],0)

        vibr_resolved_Diss[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L1'],Tev)
        vibr_resolved_Ion[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L4'],Tev)

    
    def ichihara_rates(Ti):
        #calculate effective molecular CX rate using Ichihara rates
        import mat73, os
        ichi_tables = os.path.join(os.path.dirname(__file__), 'MolCX_Ichi.mat')

        ichi_tab = mat73.loadmat(ichi_tables)
        from scipy.interpolate import interp1d

        #make interpolation object
        f_ichi = interp1d(ichi_tab['T'], 1e-6*ichi_tab['MolCX_Ichi_01'], bounds_error=False,
                            fill_value=(1e-6*ichi_tab['MolCX_Ichi_01'][:, 0], 1e-6*ichi_tab['MolCX_Ichi_01'][:, -1])) #interpolate Ichihara tables at EH2=0.1 eV. Nearest neighbour extrapolation

        return f_ichi(Ti)
    
    #in the case ichihara data is used for the vibrational distribution, this also needs to be used for the calculation of effective rates. 
    if ichi: 
        vibr_resolved_CX = ichihara_rates(Tiv)
        #for i in range(0,15):
        #    vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'Q6'],Tiv)


    #Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
    eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2,axis=0)
    eff_mol_diss = np.sum(vibr_resolved_Diss*fv_H2, axis=0)
    eff_mol_ion = np.sum(vibr_resolved_Ion*fv_H2, axis=0)

    # Calculate logarithmic fit coefficients for the effective rates
    p_cx = np.flip(np.polyfit(np.log(Tev),np.log(eff_mol_cx/1e-6),8))
    p_diss = np.flip(np.polyfit(np.log(Tev),np.log(eff_mol_diss/1e-6),8))
    p_ion = np.flip(np.polyfit(np.log(Tev),np.log(eff_mol_ion/1e-6),8))

    # fit_cx = eval_1D(p_cx,Tev)
    # fit_diss = eval_1D(p_cx,Tev)
    # fit_ion = eval_1D(p_cx,Tev)

    return p_cx, p_diss, p_ion, eff_mol_cx, eff_mol_diss, eff_mol_ion, fv_H2 #, fit_cx, fit_diss, fit_ion

def replace_line(file_name, line_num, text, overwrite=False, new_file_name='bla.tex'):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()


def output_string(A):
    exponent = str(f"{A:E}").split('E')[1]
    output = f"{A / 10**int(exponent):.12f}D{exponent}"
    if A>0:
        output = ' ' + output
    return output

def replace_block_1d(file_name_i, start_line, X,overwrite=False,new_file_name = 'bla.tex'):
    #overwrite protection
    if not overwrite:
        import shutil
        shutil.copyfile(file_name_i, new_file_name)
        file_name = new_file_name
    else:
        file_name = file_name_i
    #overwrite lines
    replace_line(file_name, start_line, '  b0 ' + output_string(X[0]) + '  b1 ' + output_string(X[1]) + '  b2 ' + output_string(X[2])+'\n')
    replace_line(file_name, start_line+1, '  b3 ' + output_string(X[3]) + '  b4 ' + output_string(X[4]) + '  b5 ' + output_string(X[5])+'\n')
    replace_line(file_name, start_line+2, '  b6 ' + output_string(X[6]) + '  b7 ' + output_string(X[7]) + '  b8 ' + output_string(X[8])+'\n')

Te_min = 0.1
Te_max = 100
Te_reso = int(1e3)

input_crm = 'input_ichi.dat'
co_cx_ichi, co_diss, co_ion, eff_ichi_cx, eff_mol_diss, eff_mol_ion, vibr_ichi = run_demo(input_crm, iso_mass=2, ichi=True,Te_min=Te_min,Te_max=Te_max,Te_reso=Te_reso)

import numpy as np
Tev = np.linspace(Te_min,Te_max,Te_reso)

import matplotlib.pyplot as plt
plt.figure()
plt.loglog(Tev,eff_ichi_cx,label='Effective rate (NEW)')
# plt.plot(Tev, fit_ichi,  '--', label='Fit')
plt.xlabel('Temperature (eV)')
plt.title('CX effective rate Ichihara')
plt.legend()

input_crm = 'input.dat'
co_cx, co_diss, co_ion, eff_mol_cx, eff_mol_diss, eff_mol_ion, vibr_h2vibr = run_demo(input_crm, iso_mass=2,Te_min=Te_min,Te_max=Te_max,Te_reso=Te_reso)

replace_block_1d('rates/amjuel.tex', 3131, co_cx, new_file_name='rates/amjuel_altered.tex')
replace_block_1d('rates/amjuel_altered.tex', 3157, co_diss, overwrite=True)
replace_block_1d('rates/amjuel_altered.tex', 3183, co_ion, overwrite=True)
replace_block_1d('rates/amjuel_altered.tex', 3209, co_cx_ichi, overwrite=True)

import CRUMPET

plt.figure()
plt.loglog(Tev,eff_mol_cx,'b',label='Charge Exchange')
plt.plot(Tev, eff_mol_diss,'g', label='Dissociation')
plt.plot(Tev, eff_mol_ion,'r', label='Ionization')
# plt.plot(Tev, fit_h2vibr,  '--', label='Fit')
X = CRUMPET.ratedata.RateData(rates={'HYDHEL': '/rates/HYDHEL.tex'})
plt.plot(Tev, eval_1D(X.reactions['HYDHEL']['H.2']['2.2.9'],Tev),'r--',label='Ionization HYDHEL')
plt.plot(Tev, eval_1D(X.reactions['HYDHEL']['H.2']['2.2.5'],Tev),'g--',label='Dissociation HYDHEL')
plt.plot(Tev, eval_1D(X.reactions['HYDHEL']['H.2']['3.2.3'],Tev/2),'b--',label='Charge exchange HYDHEL')
plt.xlabel('Temperature (eV)')
plt.ylabel('Effective rate')
plt.title('Effective rates for different reactions')
plt.legend()

plt.figure()
plt.loglog(Tev,eff_ichi_cx, label='Ichihara')
plt.loglog(Tev,eff_mol_cx, label='H2VIBR')
plt.xlabel('Temperature (eV)')
plt.ylabel('Effective rate')
plt.title('CX effective rates')
plt.legend()

plt.figure()
plt.loglog(Tev,np.transpose(vibr_ichi))
plt.xlabel('Temperature (eV)')
plt.ylabel('Fractional abundance of vibrational distribution')
plt.title('Vibr. distribution as function of T, Ichihara')

plt.figure()
plt.loglog(Tev,np.transpose(vibr_h2vibr))
plt.xlabel('Temperature (eV)')
plt.ylabel('Fractional abundance of vibrational distribution')
plt.title('Vibr. distribution as function of T, H2VIBR')



plt.show()

print('Done')
# replace_block_1d('C:/Users/Gebruiker/OneDrive - TU Eindhoven/Documenten/Masters/Internship/Vibr_mod_effective_rates/rates/amjuel.tex',  3131, co)
