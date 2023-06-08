import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate

e = 1.6e-19  #Elementary charge
me = 9.11e-31 #electron mass
kb = 1.38e-23 #Bolzmann constant

def get_coeffs(l_state, h_state):

    with open('Fantz/Table 2 Franck-Condon Factors/D2_'+h_state[:-2]+'-'+l_state[:-2]+'_FCF.dat', 'r') as f:
        lines  = f.read().splitlines()
        data = []
        for line in lines:
            columns = line.split()
            data.append(columns)

    arr = np.array(data)[1:, 1:]

    table = np.zeros(np.shape(arr))
    for i in range(np.shape(arr)[0]):
        for j in range(np.shape(arr)[1]):
            table[i,j] = float(arr[i,j])

    return table

def extract_table(i_state, f_state, vibr_num):
    a0 = 5.29177210903e-11
    with open('rates/MCCC/'+i_state+f'-excitation/vi={vibr_num}/MCCC-el-D2-'+f_state+f'_bound.'+i_state+f'_vi={vibr_num}.txt','r') as file:
        table = []
        lines = file.readlines()
        for line in lines:
            if '#' not in line:
                number = [float(line.split()[0]), float(line.split()[1])]
                table.append(number)
    table = np.array(table)
    table[:,1] = table[:,1]*a0**2
    table[0,1] = 1e-99
    return table

def extract_table_res(i_state,f_state,i_vibr,f_vibr):
    a0 = 5.29177210903e-11
    with open('rates/MCCC/'+i_state+f'-excitation/vi={i_vibr}/MCCC-el-D2-'+f_state+f'_vf={f_vibr}.'+i_state+f'_vi={i_vibr}.txt','r') as file:
        table = []
        lines = file.readlines()
        for line in lines:
            if '#' not in line:
                number = [float(line.split()[0]), float(line.split()[1])]
                table.append(number)
    table = np.array(table)
    table[:,1] = table[:,1]*a0**2
    table[0,1] = 1e-99
    return table


def calc_rates(E, sigma, Tev):
    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def fit_temp(table, Tev, E_init, E_final, i_vibr, f_vibr):

    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1], fill_value = 'extrapolate')
    table = np.exp(np.transpose([np.linspace(min(logtable[:,0]),np.log(1000),100), o(np.linspace(min(logtable[:,0]),np.log(1000),100))]))

    # Calculating rates as function of temperature
    rates = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates[i] = calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!

    # Deexcitation using detailed balance
    rates_deex = rates*np.exp((E_final[f_vibr]-E_init[i_vibr])/Tev)

    

    fit = np.flip(np.polyfit(np.log(Tev),np.log(rates/1e-6),8)) 
    fit_deex = np.flip(np.polyfit(np.log(Tev),np.log(rates_deex/1e-6),8))

    return fit, fit_deex, rates

def numpy_to_string(i_state, f_state,i): 
    part_1 = "\subsection{\n" +\
                f"Reaction 2.{i}c1\n"+\
                f"$ e + H_2(v={i}) \\rightarrow e + H2(n="+f_state[0]+")+ $ (electronic excitation)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived by integrating MCCC cross sections over Bolzmann distribution\n"+\
                "Excitation rates from "+i_state+" to "+f_state+" \n"+\
                "\n"+\
                "\\begin{small}\\begin{verbatim}\n"+\
                "\n"
                
    part_2 =   "\n"+\
                "\\end{verbatim}\end{small}\n"+\
                "\n"+\
                "\\newpage\n"
    return part_1, part_2

def numpy_to_string_res(i_state, f_state,i,j, direction): 
    part_1 = "\subsection{\n" +\
                f"Reaction 12.{i}"+direction+f"{j}\n"+\
                f"$ e + H_2(n="+i_state+f",v={i}) \\rightarrow e + H2(n="+f_state+f",v={j}) $ (electronic excitation)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived by integrating MCCC cross sections over Bolzmann distribution\n"+\
                "Excitation rates from "+i_state+f", v={i} to "+f_state+f", v={j}\n"+\
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



def eval_1D(coeff,T):
    o = np.zeros(np.shape(T))
    for i in range(0,len(coeff)):
        o = o + coeff[i]*(np.log(T)**i)
    return 1e-6*np.exp(o)




Tev = 10**np.linspace(np.log10(0.2),np.log10(100),100)



def unresolved(i_state,f_state):
    full_string = ''
    for i in range(15):
        table = extract_table(i_state,f_state,i)
        fit, rates = fit_temp(table,Tev)

        str1,str2=numpy_to_string_res(i_state,f_state,i,j)
        full_string += str1+block_string(fit)+str2
        
        plt.plot(Tev, rates)
    plt.show()
    return full_string

def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2].replace('*','0'))
    return E

def resolved(i_state, f_state):
    E_init = D2_energies('Fantz/Table 1 Vib Eigenvalues/'+i_state[:-2]+'_EV.txt')
    E_final = D2_energies('Fantz/Table 1 Vib Eigenvalues/'+f_state[:-2]+'_EV.txt')
    
    full_string = ''
    for i in range(15):
        for j in range(15):
            table = extract_table_res(i_state,f_state,i,j)
            fit, fit_deex, rates = fit_temp(table,Tev, E_init,E_final,i,j)

            str1,str2=numpy_to_string_res(i_state,f_state,i,j,'u')
            full_string += str1+block_string(fit)+str2

            str1,str2=numpy_to_string_res(f_state,i_state,j,i,'d')
            full_string += str1+block_string(fit_deex)+str2
    return full_string

def resolved_exc(i_state,f_state):
    FCFs = get_coeffs(i_state, f_state)
    factor = np.sum(FCFs[:15,0])
    a0 = 5.29177210903e-11
    with open('rates/MCCC/MCCC-el-H2/MCCC-el-H2-'+i_state+'/excitation/vi=0/MCCC-el-H2-'+f_state+'_total.'+i_state+'_vi=0.txt','r') as file:
        table_0 = []
        lines = file.readlines()
        for line in lines:
            if '#' not in line:
                number = [float(line.split()[0]), float(line.split()[1])]
                table_0.append(number)
    table_0 = np.array(table_0)
    table_0[:,1] = table_0[:,1]*a0**2
    table_0[0,1] = 1e-99

    E_init = D2_energies('Fantz/Table 1 Vib Eigenvalues/'+i_state[:-2]+'_EV.txt')
    E_final = D2_energies('Fantz/Table 1 Vib Eigenvalues/'+f_state[:-2]+'_EV.txt')

    full_string = ''
    for i in range(15):
        for j in range(15):
            table = np.transpose([table_0[:,0],table_0[:,1]*FCFs[j,i]/factor])  # i and j indexes switched because of orientation of original table
            fit, fit_deex, rates = fit_temp(table,Tev, E_init,E_final,i,j)

            str1,str2=numpy_to_string_res(i_state,f_state,i,j,'u')
            full_string += str1+block_string(fit)+str2

            str1,str2=numpy_to_string_res(f_state,i_state,j,i,'d')
            full_string += str1+block_string(fit_deex)+str2
    
    return full_string





i_state = 'c3Pu'
f_state = 'd3Pu'
full_string = resolved_exc(i_state, f_state)      




# for i in range(4):
#     table = extract_table_res(i_state,f_state,0,i)
#     fit, rates = fit_temp(table,Tev)
#     plt.plot(Tev,rates, label = f'{i}')

# plt.legend()
# plt.show()

# x = eval_1D(fit,Tev)
# plt.plot(Tev,rates)
# plt.plot(Tev,x,'--')
# plt.show()


with open('rates/h2vibr_custom.tex','a') as file:
    file.write(full_string)


print('Done')