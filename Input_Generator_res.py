import numpy as np

def get_energies(state):
    E = np.zeros(15)
    with open('Fantz/Table 1 Vib Eigenvalues/'+state+'_EV.txt', 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E

def species_state(state):
    energies = get_energies(state)
    t=''
    for i in range(15):
        t+='* H2(n='+state+f',v={i})\n'+\
            f'   V {energies[i]}\n'
    return t

def get_coeffs(file_name, nu_eff):

    with open(file_name, 'r') as f:
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
    
    table  = table*nu_eff

    return table

def rea_X1_B1Su(arr):
    t = ''
    for i in range(15):
        for j in range(15):
            t +=f'* MCCCDB B {i}to{j} \n'+\
                f'e + H2(n=X1,v={i}) > e + H2(n=B1,v={i})\n'+\
                f'rates/MCCC/X1Sg-excitation/vi={i}/MCCC-el-D2-B1Su_vf={j}.X1Sg_vi={i}.txt\n\n'+\
                f'* USER COEFFICIENT B1_X1_{j}to{i}\n'+\
                f'H2(n=B1,v={j}) > H2(n=X1,v={i})\n'+\
                f'{arr[i,j]}\n\n'
    return t

def rea_vibr_trans():
    t = ''
    for i in range(15):
        for j in range(15):
            t +=f'* MCCCDB trans {i}to{j}\n'+\
                f'e + H2(v={i}) > e + H2(v={j})\n'+\
                f'rates/Laporta/vibr_trans/vi={i}_vf={j}.txt\n\n'
    return t


def gen_input(new_file_name, B_X = False, vibr = False, C_X=False):

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

    #Input species ground state molecular hydrogen
    string+= species_state('X1')
        
    # Input excited state hydrogen molecules
    if B_X:
        string+=species_state('B1')
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

    # Vibrational transitions (change later to Laporta)
    if vibr:
        string+='* H2VIBR H.2 2.$v&\n' +\
                'e + H2(v=$) > e + H2(v=&)\n\n' 
    
    if B_X:
        B_X_rate = get_coeffs('Fantz/Table 2 Franck-Condon Factors/D2_B1-X1_FCF.dat', 7.7771e+08)
        string += rea_X1_B1Su(B_X_rate)

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
            'H2(n=X1,v=0)     1e10\n'+\
            '* verbose   0   # Show verbose output\n'
    
    with open(new_file_name, 'w') as f: 
        f.write(string)



gen_input('input_res.dat',B_X = True, vibr=True)


