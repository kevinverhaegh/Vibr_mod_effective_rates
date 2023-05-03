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
    

def gen_input(new_file_name, B_X = True, C_X = True, ion = True, diss = True):
    
    n=0
    string = ''

    # Excitation and emission from B1Su to ground
    if B_X:
        B_X_rate = FC_factors('Table 2 Franck-Condon Factors/D2_B1-X1_FCF.dat', 7.7771e+08)
        string += rea_X1_B1Su(B_X_rate)

        H2B = '* H2(n=B)\n'+\
	            '   V 11.36832\n'
        
    # Excitation and emission from C1Pu to ground
    if C_X: 
        C_X_rate = FC_factors('Table 2 Franck-Condon Factors/D2_C1-X1_FCF.dat', 1.0532e+09)
        string += rea_X1_C1Pu(C_X_rate)

        H2C = '* H2(n=C)\n'+\
	            '   V 12.41104\n'

    # Ionization from ground
    if ion:
        string +=rea_X1_ion()

    # Dissociation from ground
    if diss: 
        string += rea_X1_diss()

    insert_string('input.dat', new_file_name, string, 70)

    if B_X:
        insert_string(new_file_name, new_file_name, H2B, 40+2*n, overwrite=True)
        n+=1
    if C_X: 
        insert_string(new_file_name, new_file_name, H2C, 40+2*n, overwrite=True)
        n+=1

gen_input('input_new.dat')

        


