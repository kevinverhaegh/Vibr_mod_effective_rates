import numpy as np

def get_energies(state):
    E = np.zeros(15)
    with open('Fantz/Table 1 Vib Eigenvalues/'+state+'_EV.txt', 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1].replace('*','0')
            E[i] = float(line.split()[2])
    return E

def species_state(state):
    energies = get_energies(state)
    t=''
    for i in range(15):
        t+='* H2(n='+state+f',v={i})\n'+\
            f'   V {energies[i]}\n'
    return t

def get_coeffs(l_state, h_state):

    with open('Fantz/Aik/D2_'+h_state[:-2]+'-'+l_state[:-2]+'_Aik.dat', 'r') as f:
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

def get_coeffs_diss(l_state, h_state):

    with open('Fantz/Aik/D2_'+h_state[:-2]+'-'+l_state[:-2]+'_Aik.dat', 'r') as f:
        lines  = f.read().splitlines()
        data = []
        for line in lines:
            columns = line.split()
            data.append(columns)

    arr = np.array(data)[:,1]

    table = np.zeros(np.shape(arr))
    for i in range(len(arr)):
            table[i] = float(arr[i])

    return table

def rea_el_exc(arr, initial_state, final_state):
    t = ''
    for i in range(15):
        for j in range(15):
            # t +='* MCCCDB '+initial_state+'_'+final_state+f' {i}to{j} \n'+\
            #     f'e + H2(n='+initial_state[:2]+f',v={i}) > e + H2(n='+final_state[:2]+f',v={i})\n'+\
            #     'rates/MCCC/'+initial_state+f'-excitation/vi={i}/MCCC-el-D2-'+final_state+f'_vf={j}.'+initial_state+f'_vi={i}.txt\n\n'+\
            #     f'* USER COEFFICIENT '+initial_state+'_'+final_state+f'{j}to{i}\n'+\
            #     'H2(n='+final_state[:2]+f',v={j}) > H2(n='+initial_state[:2]+f',v={i})\n'+\
            #     f'{arr[i,j]}\n\n'
            t +=f'* USER COEFFICIENT '+final_state+'_'+initial_state+f'{i}to{j}\n'+\
                'H2(n='+final_state[:-2]+f',v={i}) > H2(n='+initial_state[:-2]+f',v={j})\n'+\
                f'{arr[i,j]}\n\n'
    return t

def rea_el_exc_diss(arr, h_state):
    t = ''
    for i in range(15):
        t +=f'* USER COEFFICIENT '+h_state+f'_v={i}\n'+\
            'e + H2(n='+h_state[:-2]+f',v={i}) > e + 2*H(n=1)\n'+\
            f'{arr[i]}\n\n'
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


def gen_input(new_file_name, B_X = False, d3_X=False, c3_X = False, vibr_hyd = False, a3Sg_X=False, EF_X=False, vibr_lap=False, C_X=False, ion_X1=False, ion_B1=False, diss_att_X1=False, diss_att_B1 = False):

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
    if EF_X:
        string+=species_state('EF1')
    if a3Sg_X:
        string+=species_state('a3')
    if c3_X:
        string+=species_state('c3')
    if d3_X:
        string+=species_state('d3')

        
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
        B_X_rate = get_coeffs('X1Sg','B1Su')
        string += rea_el_exc(B_X_rate, 'X1Sg', 'B1Su')
    if C_X:
        C_X_rate = get_coeffs('X1Sg','C1Pu')
        string += rea_el_exc(C_X_rate, 'X1Sg', 'C1Pu')
    if EF_X:
        EF_B_rate = get_coeffs('B1Su','EF1Sg')
        EF_C_rate = get_coeffs('C1Pu','EF1Sg')
        B_EF_rate = get_coeffs('EF1Sg','B1Su')
        C_EF_rate = get_coeffs('EF1Sg','C1Pu')

        string += rea_el_exc(EF_B_rate, 'B1Su', 'EF1Sg')
        string += rea_el_exc(EF_C_rate, 'C1Pu', 'EF1Sg')
        string += rea_el_exc(B_EF_rate, 'EF1Sg', 'B1Su')
        string += rea_el_exc(C_EF_rate, 'EF1Sg', 'C1Pu')        
    if a3Sg_X:
        a3_diss = get_coeffs_diss('b3Sg','a3Sg')
        string+=rea_el_exc_diss(a3_diss,'a3Sg')
    if c3_X:
        c3_a3_rate = get_coeffs('a3Sg','c3Pu')
        string+=rea_el_exc(c3_a3_rate,'a3Sg','c3Pu')
        a3_c3_rate = get_coeffs('c3Pu','a3Sg')
        string+=rea_el_exc(a3_c3_rate,'c3Pu','a3Sg')  
    if d3_X:
        d3_a3_rate = get_coeffs('a3Sg','d3Pu')
        string+=rea_el_exc(d3_a3_rate,'a3Sg','d3Pu')
     


        

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
                'H2VIBR  rates/h2vibr_custom.tex\n'+\
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



gen_input('input_false.dat',B_X=True,C_X=True,c3_X=True,d3_X=True,EF_X=True,a3Sg_X=True)


