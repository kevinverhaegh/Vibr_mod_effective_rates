import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate

e = 1.6e-19  #Elementary charge
me = 9.11e-31 #electron mass
kb = 1.38e-23 #Bolzmann constant

def extract_table(i_state, i_num, f_num):
    a0 = 5.29177210903e-11
    with open('rates/MCCC/'+i_state+f'-ERD/vi={i_num}/MCCC-el-D2-'+i_state+f'_vf={f_num}.ERD.'+i_state+f'_vi={i_num}.txt','r') as file:
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

def fit_temp(table, Tev):

    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1], fill_value = 'extrapolate')
    table = np.exp(np.transpose([np.linspace(min(logtable[:,0]),np.log(1000),100), o(np.linspace(min(logtable[:,0]),np.log(1000),100))]))

    # Calculating rates as function of temperature
    rates = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates[i] = calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!

    fit = np.flip(np.polyfit(np.log(Tev),np.log(rates/1e-6),8))

    return fit, rates

def numpy_to_string(i_state, i,j): 
    part_1 = "\subsection{\n" +\
                f"Reaction 2.{i}p{j}\n"+\
                f"$ e + H_2(v={i}) \\rightarrow e + H_2(v={j}) $ (electronic excitation + radiative decay)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived by integrating MCCC cross sections over Bolzmann distribution\n"+\
                "Excitation rates from "+i_state+" to "+i_state+" \n"+\
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

i_state = 'X1Sg'
f_state = 'B1Su'

full_string = ''
for i in range(15):
    for j in range(15):
        table = extract_table(i_state,i,j)
        fit, rates = fit_temp(table,Tev)

        str1,str2=numpy_to_string(i_state,i,j)
        full_string += str1+block_string(fit)+str2
    


with open('rates/h2vibr_custom.tex','a') as file:
    file.write(full_string)


print('Done')