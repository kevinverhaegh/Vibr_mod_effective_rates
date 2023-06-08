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
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def extract_table(lines, initial_state, final_state, el_state):

    for i,line in enumerate(lines): 
        if 'PROCESS: E + D2('+el_state[0]+f',v={initial_state}) -> E + D2('+el_state[0]+f',v={final_state}), Excitation' in line:
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


def fit_temp_laporta(table, Tev, E, i_state, f_state):

    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1], fill_value = 'extrapolate')
    table = np.exp(np.transpose([np.linspace(min(logtable[:,0]),np.log(1000),1000), o(np.linspace(min(logtable[:,0]),np.log(1000),1000))]))

    # Calculating rates as function of temperature
    rates_ex = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates_ex[i] = calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!
    
    
    # Calculating deexcitation rates using detailed balance 
    rates_deex = rates_ex*np.exp((E[f_state]-E[i_state])/Tev)

    fit_ex = np.flip(np.polyfit(np.log(Tev),np.log(rates_ex/1e-6),8))
    fit_deex = np.flip(np.polyfit(np.log(Tev),np.log(rates_deex/1e-6),8))

    return fit_ex, fit_deex, rates_ex

def numpy_to_string(state,i,j): 
    part_1 = "\subsection{\n" +\
                f"Reaction 3.{i}t{j}\n"+\
                f"$ e + H_2(n="+state+f",v={i}) \\rightarrow e + H_2(n="+state+f",v={j}) (direct transition)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived from Laporta cross sections.\n"+\
                "Taken at $E(H_2) = 0.1 \\approx 0.0$ eV,  and fit is for temperature $T_p=T$ with $H_2$ at rest.\n"+\
                "\n"+\
                "\\begin{small}\\begin{verbatim}\n"+\
                "\n"
                
    part_2 =   "\n"+\
                "\\end{verbatim}\end{small}\n"+\
                "\n"+\
                "\\newpage\n"
    return part_1, part_2

def output_string(A):
    exponent = str(f"{A:E}").split('E')[1]
    output = f"{A / 10**int(exponent):.12f}D{exponent}"
    if A>0:
        output = ' ' + output
    return output

def block_string(X):
    block = '  b0 ' + output_string(X[0]) + '  b1 ' + output_string(X[1]) + '  b2 ' + output_string(X[2])+'\n'+\
            '  b3 ' + output_string(X[3]) + '  b4 ' + output_string(X[4]) + '  b5 ' + output_string(X[5])+'\n'+\
            '  b6 ' + output_string(X[6]) + '  b7 ' + output_string(X[7]) + '  b8 ' + output_string(X[8])+'\n'
    return block

def gen_full_string(coeffs):
    full_string = ''
    for i in range(15):
        str1, str2 = numpy_to_string(i)
        full_string += str1+block_string(coeffs[i,:])+str2
    return full_string

def insert_string(file_name, file_name_new, string, line_num):
    
    # Copy file into new file
    import shutil
    shutil.copyfile(file_name, file_name_new)

    # Create lines variable with all the line numbers of the file
    with open(file_name_new,'r') as f:
        lines = f.readlines()
    
    with open(file_name_new, 'r+') as f: 
        for i, line in enumerate(lines):
            if i == line_num:
                f.write(string + '\n')
            f.write(line)

Tev = np.linspace(0.2,100, 100)
Tiv = Tev # Assume Ti=Te
ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)

def create_tex(state):

    E = D2_energies('Fantz/Table 1 Vib Eigenvalues/'+state[:2]+'_EV.txt')

    full_string = ''
    with open('rates/Laporta/vibr_trans_'+state+'/Cross section.txt', 'r') as file:
        lines = file.readlines()
    for i in range(15):
        print(i)
        for j in range(i,15):
            table = extract_table(lines,i,j,state)
            fit_ex,fit_deex, rates_ex = fit_temp_laporta(table,Tev,E,i,j)

            # Excitation part
            str1,str2=numpy_to_string(state,i,j)
            full_string += str1+block_string(fit_ex)+str2
            
            # De-excitaion part
            if i is not j:
                str1,str2=numpy_to_string(state,j,i)
                full_string += str1+block_string(fit_deex)+str2
    return full_string

full_string = create_tex('B1Su')

# insert_string('rates/h2vibr_ichi.tex', 'rates/h2vibr_ichi_lap.tex',full_string, 2805)

# full_string = '\section{H.2 :  Fits for $<\sigma v> (T)$}\n\n'+full_string

with open('rates/h2vibr_custom.tex', 'a') as file:
    file.write(full_string)

print('Done')

