import numpy as np

def insert_string(file_name, file_name_new, string, line_num, overwrite=False):
    
    # Copy file into new file
    if not overwrite: 
        import shutil
        shutil.copyfile(file_name, file_name_new)
    else: 
        file_name_new = file_name

    # Create lines variable with all the line numbers of the file
    with open(file_name_new,'r') as f:
        lines = f.readlines()
    
    with open(file_name_new, 'r+') as f: 
        for i, line in enumerate(lines):
            if i == line_num:
                f.write(string + '\n')
            f.write(line)

def rea_X1_B1Su(arr):
    t = ''
    for i in range(15):
        t +=f'* MCCCDB B {i} \n'+\
            f'e + H2(v={i}) > e + H2(n=B)\n'+\
            f'rates/MCCC/X1Sg-excitation/vi={i}/MCCC-el-D2-B1Su_bound.X1Sg_vi={i}.txt\n\n'+\
            f'* USER COEFFICIENT B_X{i}\n'+\
            f'H2(n=B) > H2(v={i})\n'+\
            f'{arr[i]}\n\n'
    return t

def rea_X1_C1Pu(arr):
    t = ''
    for i in range(15):
        t +=f'* MCCCDB C {i} \n'+\
            f'e + H2(v={i}) > e + H2(n=C)\n'+\
            f'rates/MCCC/X1Sg-excitation/vi={i}/MCCC-el-D2-C1Pu_bound.X1Sg_vi={i}.txt\n\n'+\
            f'* USER COEFFICIENT C_X{i}\n'+\
            f'H2(n=C) > H2(v={i})\n'+\
            f'{arr[i]}\n\n'
    return t

def rea_X1_ion():
    t = ''
    for i in range(15):
        t +=f'* MCCCDB ion {i} \n'+\
            f'e + H2(v={i}) > 2*e + H2+\n'+\
            f'rates/MCCC/X1Sg-ionization/vi={i}/MCCC-el-D2-TICS.X1Sg_vi={i}.txt\n\n'
    return t

def rea_X1_diss():
    states = ['a3Sg', 'B1Su', 'Bp1Su', 'C1Pu', 'c3Pu', 'D1Pu', 'd3Pu', 'e3Su', 'EF1Sg', 'g3Sg', 'GK1Sg', 'H1Sg', 'h3Sg', 'I1Pg', 'i3Pg', 'J1Dg', 'j3Dg']
    t = ''
    for i in range(15):
        for j in states:
            t +=f'* MCCCDB diss X_'+j+f'{i} \n'+\
                f'e + H2(v={i}) > e + 2*H(n=1)\n'+\
                f'rates/MCCC/X1Sg-excitation/vi={i}/MCCC-el-D2-'+j+f'_DE.X1Sg_vi={i}.txt\n\n'
    return t

def rea_vibr_trans():
    t = ''
    for i in range(15):
        for j in range(15):
            t +=f'* MCCCDB trans {i}to{j}\n'+\
                f'e + H2(v={i}) > e + H2(v={j})\n'+\
                f'rates/Laporta/vibr_trans/vi={i}_vf={j}.txt\n\n'
    return t


            
def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E


def FC_factors(file_name, nu_eff):

    with open(file_name, 'r') as f:
        lines  = f.read().splitlines()
        data = []
        for line in lines:
            columns = line.split()
            data.append(columns)

    arr = np.array(data)[1:16, 1:]

    table = np.zeros(np.shape(arr))
    for i in range(np.shape(arr)[0]):
        for j in range(np.shape(arr)[1]):
            table[i,j] = float(arr[i,j])

    coeffs = nu_eff*np.sum(table, axis = 1)/15

    return coeffs
    

def gen_input(new_file_name, B_X = True, C_X = True, ion = True, diss = True, vibr = True, cx_ichi = True, diss_ion = True):

    ## TITLE 
    string = '# This is an input file for the UEDGE Python CRM\n'+\
                '# Created from scratch by Stijn Kobussen\n'+\
                '# May 2023\n\n'

    ## SPECIES
    # Input species atomic hydrogen and ions
    string+='** SPECIES\n'+\
            '* H(n=1)\n'+\
            '    V 2.375\n'+\
            '* H2+\n'+\
            '    V 15.56\n'

    #Input species ground state molecular hydrogen (D2 energies)
    E = D2_energies('Fantz/Table 1 Vib Eigenvalues/X1_EV.txt')
    for i in range(15):
        string+=f'* H2(v={i})\n'+\
	            f'  V {E[i]}\n'
        
    # Input excited state hydrogen molecules
    if B_X:
        string+='* H2(n=B)\n'+\
	            '   V 11.36832\n'
    if C_X:
        string+='* H2(n=C)\n'+\
	            '   V 12.41104\n'
    
    string +='\n\n'

    ## BACKGROUND SPECIES
    string +='** BACKGROUND\n'+\
                '* e\n'+\
                '    V 0\n'+\
                '* p\n'+\
                '    V 15.975\n'

    string+='\n\n'

    ## REACTIONS
    string+='** REACTIONS\n\n'

    if vibr:
        string+=rea_vibr_trans()
    
    if cx_ichi:
        string +='* H2VIBR H.2 2.$q6\n'+\
                    'p + H2(v=$) > H2+ + H(n=1)\n\n'

    # Excitation and emission from B1Su to ground
    if B_X:
        B_X_rate = FC_factors('Fantz/Table 2 Franck-Condon Factors/D2_B1-X1_FCF.dat', 7.7771e+08)
        string += rea_X1_B1Su(B_X_rate)

    # Excitation and emission from C1Pu to ground
    if C_X: 
        C_X_rate = FC_factors('Fantz/Table 2 Franck-Condon Factors/D2_C1-X1_FCF.dat', 1.0532e+09)
        string += rea_X1_C1Pu(C_X_rate)

    # Ionization from ground
    if ion:
        string +=rea_X1_ion()

    # Dissociation from ground
    if diss: 
        string += rea_X1_diss()
    
    if diss_ion: 
        string += '* HYDHEL H.2 2.2.10\n'+\
                    'e + H2(v=0) > 2*e + p + H(n=1)\n\n'
    
    string +='\n\n'

    ## RATES
    string += '** RATES\n'+\
                '# Define the files for the standard inputs\n'+\
                'H2VIBR  rates/h2vibr_ichi.tex\n'+\
                'HYDHEL rates/HYDHEL.tex\n\n'

    ## SETTINGS
    string +='** SETTINGS\n'+\
            '* vmax      14\n'+\
            '* n0\n'+\
            'H(n=1)      0e12\n'+\
            'H2(v=0)     1e10\n'+\
            '* verbose   0   # Show verbose output\n'
    
    with open(new_file_name, 'w') as f: 
        f.write(string)

    
gen_input('input_new.dat')

        


