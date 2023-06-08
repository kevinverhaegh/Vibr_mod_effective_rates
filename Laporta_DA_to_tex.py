import numpy as np 
import matplotlib.pyplot as plt
import scipy.integrate as integrate
import scipy.interpolate as interpolate

e = 1.6e-19  #Elementary charge
me = 9.11e-31 #electron mass
kb = 1.38e-23 #Bolzmann constant


def eval_1D(coeff,T):
    o = np.zeros(np.shape(T))
    for i in range(0,len(coeff)):
        o = o + coeff[i]*(np.log(T)**i)
    return 1e-6*np.exp(o)

def calc_rates(E,sigma, Tev):
    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def extract_table(state, lines, i_num):

    for i,line in enumerate(lines): 
        if 'PROCESS: E + D2('+state[0]+f',v={i_num}) -> ' in line:
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


def fit_temp_laporta(table, Tev):
    test = table
    # Add extrapolation in loglog space 
    logtable = np.log(table)
    o = interpolate.interp1d(logtable[:,0], logtable[:,1],bounds_error=False, fill_value = (-1e99,-1e99))
    table = np.exp(np.transpose([np.linspace(min(logtable[:,0]),np.log(1000),1000), o(np.linspace(min(logtable[:,0]),np.log(1000),1000))]))

    # Calculating rates as function of temperature
    rates = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates[i] = calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!
    

    fit = np.flip(np.polyfit(np.log(Tev),np.log(rates/1e-6),8))

    return fit,rates

def fit_temp_laporta_2(table1,table2, Tev):

    # Add extrapolation in loglog space 
    logtable1 = np.log(table1)
    logtable2 = np.log(table2)

    o = interpolate.interp1d(logtable1[:,0], logtable1[:,1], fill_value = 'extrapolate')
    p = interpolate.interp1d(logtable2[:,0], logtable2[:,1], fill_value = (-1e99,-1e99), bounds_error=False)
    

    table = np.transpose([np.exp(np.linspace(min(logtable1[:,0]),max(logtable1[:,0]),10000)), np.exp(o(np.linspace(min(logtable1[:,0]),max(logtable1[:,0]),10000)))+np.exp(p(np.linspace(min(logtable1[:,0]),max(logtable1[:,0]),10000)))])

    # Calculating rates as function of temperature
    rates = np.zeros(len(Tev))
    for i in range(len(Tev)):
        rates[i] = calc_rates(table[:,0],table[:,1],Tev[i]) # 1e-6 because units!!!!!
    

    fit = np.flip(np.polyfit(np.log(Tev),np.log(rates/1e-6),8))
    return fit, rates

def numpy_to_string(i, state): 
    part_1 = "\subsection{\n" +\
                f"Reaction 2.{i}z4\n"+\
                "$ e + H_2(n="+state[:-2]+f"v={i}) \\rightarrow H + H^- (dissociative attachment from "+state+")\n" +\
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



full_string = ''
with open('rates/Laporta/diss_attachment_B1Su/Cross section.txt', 'r') as file:
    lines1 = file.readlines()

# with open('rates/Laporta/diss_attachment_X1Sg/Cross section exc.txt', 'r') as file:
#     lines2 = file.readlines()


state = 'B1Su'
for i in range(15):
    # table1 = extract_table(lines1,i)
    # table2 = extract_table(lines2,i)

    # fit,rates = fit_temp_laporta_2(table1,table2,Tev)

    table = extract_table(state,lines1,i)
    fit,rates = fit_temp_laporta(table,Tev)

    # Excitation part
    str1,str2=numpy_to_string(i,state)
    full_string += str1+block_string(fit)+str2
    plt.loglog(Tev,rates)

plt.show()

with open('rates/h2vibr_custom.tex', 'a') as file:
    file.write(full_string)


print('Done')
