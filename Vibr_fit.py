def run_demo(input_crm,iso_mass=1,T_request=1):
    print('Running DEMO with input:')
    print(input_crm)
    import CRUMPET
    import numpy as np 
    from scipy.optimize import curve_fit

    indx_H2v = np.append(2,np.arange(3,17))

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    Tev = np.linspace(0.2,10,100)
    Tiv = Tev # Assume Ti=Te
    ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)
    crm.source[2] = 1e-100 #add small source for numerical stability (only needed if reactions that dissociate are included)

    #compute vibrational distribution H2
    fv_H2 = np.zeros([15,len(Tev)])

    #calculate vibrational distribution using Tiv = Tev/iso_mmass
    for i in range(0,len(Tev)):
        fv_H2[:,i]=crm.steady_state(Tev[i]/iso_mass,ne[i],plot=False,dt=True)[indx_H2v]

    #normalise vibrational distribution by dividing the distribution values to the sum of the distribution
    fv_H2 = fv_H2/(np.sum(fv_H2,axis=0)[None,:])

    #Get vibrationally resolved molecular CX rates from H2VIBR

    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr.tex', 'AMJUEL' : '/rates/amjuel.tex'})

    #Get rates as function of Te

    def eval_1D(coeff,T):
        o = np.zeros(np.shape(T))
        for i in range(0,len(coeff)):
            o = o + coeff[i]*(np.log(T)**i)
        return 1e-6*np.exp(o)

    vibr_resolved_CX = np.zeros([15,len(Tev)])
    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L2'],Tiv/iso_mass)

    #Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
    eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2,axis=0)

    p = np.flip(np.polyfit(np.log(Tev),np.log(eff_mol_cx/1e-6),8))
    fit = eval_1D(p,Tev)

    # show a plot comparing different effective rates
    import matplotlib.pyplot as plt
    plt.figure()
    plt.loglog(Tev,eff_mol_cx,label='Effective rate (NEW)')
    plt.plot(Tev, fit,  '--', label='Fit')
    plt.xlabel('Temperature (eV)')
    plt.legend()
    print('Generating Plot')
    plt.show()
    return p 

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

input_crm = 'input_ichi.dat'
#simple input file with additional dissociation & molecular depletion channels
co = run_demo(input_crm, iso_mass=2)

# replace_block_1d('C:/Users/Gebruiker/OneDrive - TU Eindhoven/Documenten/Masters/Internship/Vibr_mod_effective_rates/rates/amjuel.tex',  3131, co)
