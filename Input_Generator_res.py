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

def rea_el_exc(arr, initial_state, final_state):
    t = ''
    for i in range(15):
        for j in range(15):
            t +='* MCCCDB '+initial_state+'_'+final_state+f' {i}to{j} \n'+\
                f'e + H2(n='+initial_state[:2]+f',v={i}) > e + H2(n='+final_state[:2]+f',v={i})\n'+\
                'rates/MCCC/'+initial_state+f'-excitation/vi={i}/MCCC-el-D2-'+final_state+f'_vf={j}.'+initial_state+f'_vi={i}.txt\n\n'+\
                f'* USER COEFFICIENT '+initial_state+'_'+final_state+f'{j}to{i}\n'+\
                'H2(n='+final_state[:2]+f',v={j}) > H2(n='+initial_state[:2]+f',v={i})\n'+\
                f'{arr[i,j]}\n\n'
    return t

def rea_vibr_trans():
    t = ''
    for i in range(15):
        for j in range(15):
            t +=f'* MCCCDB trans {i}to{j}\n'+\
                f'e + H2(n=X1,v={i}) > e + H2(n=X1,v={j})\n'+\
                f'rates/Laporta/vibr_trans/vi={i}_vf={j}.txt\n\n'
    return t

def rea_ion_state(state):
    t = ''
    for i in range(15):
        t +='* MCCCDB ion '+state+f'{i} \n'+\
            'e + H2(n='+state[:2]+f',v={i}) > 2*e + H2+\n'+\
            'rates/MCCC/'+state+f'-ionization/vi={i}/MCCC-el-D2-TICS.'+state+f'_vi={i}.txt\n\n'
    return t

def rea_diss_att(state):
    t = '' 
    for i in range(15):
        t += '* MCCCDB att '+state+f'{i}\n'+\
                'e + H2(n='+state[:2]+f',v={i}) > H(n=1) + H-\n'+\
                'rates/Laporta/diss_attachment_'+state+f'/vi={i}.txt\n\n'
    return t


def gen_input(new_file_name, B_X = False, vibr_hyd = False, vibr_lap=False, C_X=False, ion_X1=False, ion_B1=False, diss_att_X1=False, diss_att_B1 = False):

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
        string+=species_state('C1')
        
    # Input negative ions
    if diss_att_X1 or diss_att_B1:
        string+='* H-\n'+\
                '   V -0.75\n' 
        
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
    if vibr_hyd:
        string+='* H2VIBR H.2 2.$v&\n' +\
                'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n'
    
    if vibr_lap:
        string+=rea_vibr_trans()
    
    # Excitation and decay from electronizally excited states
    if B_X:
        B_X_rate = get_coeffs('Fantz/Table 2 Franck-Condon Factors/D2_B1-X1_FCF.dat', 7.7771e+08)
        string += rea_el_exc(B_X_rate, 'X1Sg', 'B1Su')
    if C_X:
        C_X_rate = get_coeffs('Fantz/Table 2 Franck-Condon Factors/D2_C1-X1_FCF.dat', 1.0532e+09)
        string += rea_el_exc(C_X_rate, 'X1Sg', 'C1Pu')
    
    # Ionization 
    if ion_X1:
        string+=rea_ion_state('X1Sg')
    if ion_B1: 
        string+=rea_ion_state('B1Su')
    
    # Dissociative attachment
    if diss_att_X1:
        string+=rea_diss_att('X1Sg')
    if diss_att_B1:
        string+=rea_diss_att('B1Su')

    

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



gen_input('input_res.dat', vibr_lap=True, B_X=True)


