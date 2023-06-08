#%%
from xml.etree.ElementTree import C14NWriterTarget
import numpy as np

def get_energies(state):
    E = np.zeros(15)
    with open('Fantz/Table 1 Vib Eigenvalues/'+state[:-2]+'_EV.txt', 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1].replace('*','0')
            E[i] = float(line.split()[2])
    return E

def species(state):
    energies = get_energies(state)
    t=''
    if state=='X1Sg':
        for i in range(1,15):
            t+='* H2(n='+state[:-2]+f',v={i})\n'+\
                f'   V {energies[i]}\n'
    else:
        for i in range(15):
            t+='* H2(n='+state[:-2]+f',v={i})\n'+\
                f'   V {energies[i]}\n'
    return t

def get_coeffs(i_state, f_state):

    with open('Fantz/Aik/D2_'+i_state[:-2]+'-'+f_state[:-2]+'_Aik.dat', 'r') as f:
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

def get_coeffs_diss(i_state, f_state):

    with open('Fantz/Aik/D2_'+i_state[:-2]+'-'+f_state[:-2]+'_Aik.dat', 'r') as f:
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

def rea_rad_decay(arr, i_state, f_state):
    t = ''
    if f_state=='b3Sg':
        for i in range(15):
                t +=f'* USER COEFFICIENT '+i_state+'_'+f_state+f'{i}\n'+\
                    'H2(n='+i_state[:-2]+f',v={i}) > 2*H(n=1)\n'+\
                    f'{arr[i]}\n\n'   
    else:   
        for i in range(15):
            for j in range(15):
                t +=f'* USER COEFFICIENT '+i_state+'_'+f_state+f'{i}to{j}\n'+\
                    'H2(n='+i_state[:-2]+f',v={i}) > H2(n='+f_state[:-2]+f',v={j})\n'+\
                    f'{arr[i,j]}\n\n'
    return t

def rea_el_exc_diss(arr, h_state):
    t = ''
    for i in range(15):
        t +=f'* USER COEFFICIENT '+h_state+f'_v={i}\n'+\
            'e + H2(n='+h_state[:-2]+f',v={i}) > e + 2*H(n=1)\n'+\
            f'{arr[i]}\n\n'
    return t

def rea_rad_decay_unr(arr, i_state, f_state):
    t = ''
    for i in range(15):
        t +='* USER COEFFICIENT '+i_state+'_'+f_state+f'X{i}\n'+\
            'H2(n='+i_state[:-2]+') > H2(n='+f_state[:-2]+f',v={i})\n'+\
            f'{arr[i]}\n\n'
    return t

def get_coeffs_unr(file_name, nu_eff):

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

    coeffs = nu_eff*np.sum(table, axis = 1)/15 #??? Sum or average??? 

    return coeffs


def gen_input(new_file_name, vibr_resolved=True, rad_decay=True, inter_states = True, coll_deex = True, incl_ground=False, ion = False, B1 = False, 
              d3=False, c3 = False, vibr_hyd = False, diss = False, a3=False, EF1=False, mol_cx=False, diss_ion=False,
              vibr_lap=False, C1=False, diss_att_X1=False, diss_att_B1 = False, diss_att_old=False, MA=False, MA_min=False):
    '''
    Generates an input.dat file for CRUMPET. It uses a custom file in the same format as the H2VIBR database. 

        Parameters
        ----------
        new_file_name : string
            desired path for the input file.

        vibr_resolved : bool (default: True)
            switch to determine whether the electronically excited states are 
            vibrationally resolved. If vibr_resolved = False, The model will not consider 
            the vibrational distribution in the excited states. 

        rad_decay : bool (default: True)
            switch to turn radiative transitions on and off.
        
        inter_states: bool (default: True)
            switch to determine whether collisional transitions between electronically excited states are 
            taken into account. The rates for these transitions are calculated using data for H2, in stead of D2, 
            because no such data exists for D2. 

        coll_deex: bool (default: True)
            switch to determine whether collisional deexcitation of electronically excited states is taken into account.

        States
        ------------------
        B1: bool (default: False)
            switch to include B1Su
        C1: bool (default: False)
            switch to include C1Pu
        EF1: bool (default: False)
            switch to include EF1Sg
        a3: bool (default: False)
            switch to include a3Sg
        c3: bool (default: False)
            switch to include c3Pu
        d3: bool (default: False)
            switch to include d3Pu

        
        Reactions
        ---------
        vibr_lap: (default: False)
            switch to include vibrational transitions in X1Sg, using Laporta cross sections.
        vibr_hyd: (default: False)
            switch to include vibrational transitions in X1Sg, using HYDHEL rate data. 
        ion : bool (default: False)
            switch to include ionization from every included electronic state. 
        diss: bool (default: False)
            swith to include dissociation (either collisional or radiative).
        mol_cx: bool (default: False)
            switch to include molecular charge exchange.
        diss_ion: bool (default: False)
            switch to include dissociative ionization from X1Sg.
        diss_at: bool (default: False)
            switch to include dissociative attachment.
        diss_att_old: bool (default: False)
            switch to include dissociative attachment, using H2VIBR rate data. 
    '''
    
    
    ## TITLE 
    string = '# This is an input file for the UEDGE Python CRM\n'+\
                '# Created from scratch by Stijn Kobussen\n'+\
                '# May 2023\n\n'

    ## SPECIES
    string+='** SPECIES\n\n'

    #Input species ground state molecular hydrogen (for v>0)
    if incl_ground:
        string+='* H2(n=X1,v=0)\n'+\
                '    V 0.19652\n'

    string+= species('X1Sg')

    if vibr_resolved:         
        # Input excited state hydrogen molecules
        if B1:
            string+=species('B1Su')
        if C1:
            string+=species('C1Pu')
        if EF1:
            string+=species('EF1Sg')
        if a3:
            string+=species('a3Sg')
        if c3:
            string+=species('c3Pu')
        if d3:
            string+=species('d3Pu')

            
        # Input molecular ions and H-
            # Input species molecular ions
        if ion:
            string+='* H2+\n'+\
                '    V 15.56\n'
        if diss_att_X1 or diss_att_B1:
            string+='* H-\n'+\
                    '   V -0.75\n' 
            
        string +='\n\n'

        ## BACKGROUND SPECIES
        string +='** BACKGROUND\n'+\
                    '* e\n'+\
                    '    V 0\n'+\
                    '* p\n'+\
                    '    V 15.975\n'+\
                    '* H(n=1)\n'+\
                    '    V 2.375\n'
        if not incl_ground:
              string+='* H2(n=X1,v=0)\n'+\
                        '    V 0.19652\n'

        string+='\n\n'

        ## REACTIONS
        string+='** REACTIONS\n\n'

        # Vibrational transitions
        if vibr_hyd:
            string+='* H2VIBR H.2 2.$v&\n' +\
                    'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n'
        
        if vibr_lap:
            string+='* H2VIBR H.2 2.$t&\n' +\
                    'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n'
        
        if mol_cx:
            string+='* H2VIBR H.2 2.$q6\n'+\
                    'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n' 
        
        if diss:
            string+='* H2VIBR H.2 2.$l1\n'+\
                    'e + H2(n=X1,v=$) > e + 2*H(n=1)\n\n'
        
        if ion: 
            string+='* H2VIBR H.2 2.$l2\n'+\
                    'e + H2(n=X1,v=$) > e + H2+\n\n'
            
        if diss_ion:
            string += '* AMJUEL H.4 2.2.10\n'+\
                        'e + H2(n=X1,v=0) > 2*e + p + H(n=1)\n\n'
        
        if MA:
            # Dissociative ionization
            string+='* AMJUEL H.4 2.2.11\n'+\
                    'e + H2+ > 2*e + 2*p\n\n'
            # Dissociation
            string+='* AMJUEL H.4 2.2.12\n'+\
                    'e + H2+ > e + p + H(n=1)\n\n'
            # Dissociative recombination
            string+='* AMJUEL H.4 2.2.14\n'+\
                    'e + H2+ > 2*H(n=1)\n\n'
        
        if MA_min:
            # string+='* HYDHEL H.2 7.1.1\n'+\
            #         'e + H- > 2*e + H(n=1)\n\n'
            # string+='* HYDHEL H.2 7.1.2\n'+\
            #         'e + H- > 3*e + p\n\n'
            string+='* AMJUEL H.4 7.2.3a\n'+\
                    'p + H- > 2*H(n=1)\n\n'
            string+='* AMJUEL H.4 7.2.3b\n'+\
                    'p + H- > 2*e + p + H(n=1)\n\n'

            
        # Add excited state reactions
        if B1: 
            # Vibrational transitions
            string+='* H2VIBR H.2 3.$t&\n'+\
                    'e + H2(n=B1,v=$) > e + H2(n=B1,v=&)\n\n'
            
            # Excitation between states X1 and B1
            string+='* H2VIBR H.2 1.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=B1,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 1.$d&\n'+\
                        'e + H2(n=B1,v=$) > e + H2(n=X1,v=&)\n\n'
            if ion:
                string+='* H2VIBR H.2 2.$b1\n'+\
                        'e + H2(n=B1,v=$) > e + H2+\n\n'
        
        if C1:
            # Excitation between X1 and C1
            string+='* H2VIBR H.2 3.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=C1,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 3.$d&\n'+\
                        'e + H2(n=C1,v=$) > e + H2(n=X1,v=&)\n\n'
            
            if inter_states:
                # Excitation between B1 and C1
                string+='* H2VIBR H.2 8.$u&\n'+\
                        'e + H2(n=B1,v=$) > e + H2(n=C1,v=&)\n\n'
                if coll_deex:
                    string+='* H2VIBR H.2 8.$d&\n'+\
                            'e + H2(n=C1,v=$) > e + H2(n=B1,v=&)\n\n'

            if ion:
                # Ionization from C1
                string+='* H2VIBR H.2 2.$c2\n'+\
                        'e + H2(n=C1,v=$) > e + H2+\n\n'
        
        if EF1:
            # Excitation between X1 and EF1
            string+='* H2VIBR H.2 4.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=EF1,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 4.$d&\n'+\
                        'e + H2(n=EF1,v=$) > e + H2(n=X1,v=&)\n\n'
            
            if inter_states:
                # Excitation between B1 and EF1
                string+='* H2VIBR H.2 9.$u&\n'+\
                        'e + H2(n=B1,v=$) > e + H2(n=EF1,v=&)\n\n'
                if coll_deex:
                    string+='* H2VIBR H.2 9.$d&\n'+\
                            'e + H2(n=EF1,v=$) > e + H2(n=B1,v=&)\n\n'
                # Excitation between C1 and EF1
                string+='* H2VIBR H.2 10.$u&\n'+\
                        'e + H2(n=C1,v=$) > e + H2(n=EF1,v=&)\n\n'
                if coll_deex:
                    string+='* H2VIBR H.2 10.$d&\n'+\
                            'e + H2(n=EF1,v=$) > e + H2(n=C1,v=&)\n\n'
            if ion:
                # Ionization from EF1
                string+='* H2VIBR H.2 4.$l2\n'+\
                        'e + H2(n=EF1,v=$) > e + H2+\n\n'
        if a3:
            # Excitation between states X1 and a3
            string+='* H2VIBR H.2 5.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=a3,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 5.$d&\n'+\
                        'e + H2(n=a3,v=$) > e + H2(n=X1,v=&)\n\n'
            
            if ion:
                string+='* H2VIBR H.2 5.$l2\n'+\
                        'e + H2(n=a3,v=$) > e + H2+\n\n'

        
        if c3:
            # Excitation between states X1 and c3
            string+='* H2VIBR H.2 6.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=c3,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 6.$d&\n'+\
                        'e + H2(n=c3,v=$) > e + H2(n=X1,v=&)\n\n'
            
            if ion:
                string+='* H2VIBR H.2 5.$l2\n'+\
                        'e + H2(n=c3,v=$) > e + H2+\n\n'
                
        if d3:
            # Excitation between states X1 and d3
            string+='* H2VIBR H.2 7.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=d3,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 7.$d&\n'+\
                        'e + H2(n=d3,v=$) > e + H2(n=X1,v=&)\n\n'
            if ion:
                string+='* H2VIBR H.2 5.$l2\n'+\
                        'e + H2(n=d3,v=$) > e + H2+\n\n'

            if inter_states:
                # Excitation between a3 and d3
                string+='* H2VIBR H.2 11.$u&\n'+\
                        'e + H2(n=a3,v=$) > e + H2(n=d3,v=&)\n\n'
                if coll_deex:
                    string+='* H2VIBR H.2 11.$d&\n'+\
                            'e + H2(n=d3,v=$) > e + H2(n=a3,v=&)\n\n'
                    string+='* H2VIBR H.2 11.$d&\n'+\
                            'e + H2(n=d3,v=$) > e + H2(n=a3,v=&)\n\n'
                # Excitation between c3 and d3
                string+='* H2VIBR H.2 12.$u&\n'+\
                        'e + H2(n=c3,v=$) > e + H2(n=d3,v=&)\n\n'
                if coll_deex:
                    string+='* H2VIBR H.2 12.$d&\n'+\
                            'e + H2(n=d3,v=$) > e + H2(n=c3,v=&)\n\n'                


        

    
        # Radiative decay from electronizally excited states
        if rad_decay: 
            if B1:
                B_X_rate = get_coeffs('B1Su','X1Sg')
                string += rea_rad_decay(B_X_rate,'B1Su','X1Sg')
            if C1:
                C_X_rate = get_coeffs('C1Pu','X1Sg')
                string += rea_rad_decay(C_X_rate,'C1Pu','X1Sg')
            if EF1:
                EF_B_rate = get_coeffs('EF1Sg','B1Su')
                EF_C_rate = get_coeffs('EF1Sg','C1Pu')
                B_EF_rate = get_coeffs('B1Su','EF1Sg')
                C_EF_rate = get_coeffs('C1Pu','EF1Sg')

                string += rea_rad_decay(EF_B_rate, 'EF1Sg', 'B1Su')
                string += rea_rad_decay(EF_C_rate, 'EF1Sg', 'C1Pu')
                string += rea_rad_decay(B_EF_rate, 'B1Su', 'EF1Sg')
                string += rea_rad_decay(C_EF_rate, 'C1Pu', 'EF1Sg')        
            if a3:
                if diss: 
                    a3_diss = get_coeffs_diss('a3Sg','b3Sg')
                    string+=rea_rad_decay(a3_diss,'a3Sg','b3Sg')
            if c3:
                c3_a3_rate = get_coeffs('c3Pu','a3Sg')
                string+=rea_rad_decay(c3_a3_rate,'c3Pu','a3Sg')
                a3_c3_rate = get_coeffs('a3Sg','c3Pu')
                string+=rea_rad_decay(a3_c3_rate,'a3Sg','c3Pu')  
            if d3:
                d3_a3_rate = get_coeffs('d3Pu','a3Sg')
                string+=rea_rad_decay(d3_a3_rate,'d3Pu','a3Sg')




    ## IN THE CASE THAT YOU DO NOT WANT TO VIBRATIONALLY RESOLVE THE ELECTRONICALLY EXCITED STATES

    if not vibr_resolved:
        # Input excited state hydrogen molecules
        if B1:
            string+='* H2(n=B1)\n'+\
                    '   V 11.36832\n'
        if C1:
            string+='* H2(n=C1)\n'+\
                    '   V 12.41104\n'
            
        # Input negative ions
        if diss_att_X1 or diss_att_old:
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

        #  Vibrational transitions
        if vibr_lap:
            string+='* H2VIBR H.2 2.$t&\n' +\
                        'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n' 
        if vibr_hyd:
            string+='* H2VIBR H.2 2.$v&\n' +\
                        'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n' 
        
        # Molecular charge exchange
        if mol_cx:
            string +='* H2VIBR H.2 2.$q6\n'+\
                        'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n'
        
        # Ionization from ground
        if ion: 
            string+='* H2VIBR H.2 2.$l2\n'+\
                    'e + H2(n=X1,v=$) > e + H2+\n\n'

        # Dissociation from ground
        if diss:
            string+='* H2VIBR H.2 2.$l1\n'+\
                    'e + H2(n=X1,v=$) > e + 2*H(n=1)\n\n'
        
        # Dissociative ionization
        if diss_ion: 
            string += '* HYDHEL H.2 2.2.10\n'+\
                        'e + H2(n=X1,v=0) > 2*e + p + H(n=1)\n\n'
        


        # Electronic (de-)excitation
        if B1: 
            string+='* H2VIBR H.2 13.$u1\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=B1)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 13.$d1\n'+\
                    'e + H2(n=B1) > e + H2(n=X1,v=$)\n\n'
        if C1: 
            string+='* H2VIBR H.2 14.$u1\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=C1)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 14.$d1\n'+\
                    'e + H2(n=C1) > e + H2(n=X1,v=$)\n\n'

            
        if rad_decay:
            # Excitation and emission from B1Su to ground
            if B1:
                B_X_rate = get_coeffs_unr('Fantz/Table 2 Franck-Condon Factors/D2_B1-X1_FCF.dat', 7.7771e+08)
                string += rea_rad_decay_unr(B_X_rate,'B1Su', 'X1Sg')

            # Excitation and emission from C1Pu to ground
            if C1: 
                C_X_rate = get_coeffs_unr('Fantz/Table 2 Franck-Condon Factors/D2_C1-X1_FCF.dat', 1.0532e+09)
                string += rea_rad_decay_unr(C_X_rate, 'C1Pu', 'X1Sg')



    

    # Dissociative attachment
    if diss_att_X1:
        string+='* H2VIBR H.2 2.$l4\n'+\
                'e + H2(n=X1,v=$) > H(n=1) + H-\n\n'
    if diss_att_B1:
        string+='* H2VIBR H.2 2.$z4\n'+\
                'e + H2(n=B1,v=$) > H(n=1) + H-\n\n'

    

    string +='\n\n'

    ## RATES
    string += '** RATES\n'+\
                '# Define the files for the standard inputs\n'+\
                'H2VIBR  rates/h2vibr_custom.tex\n'+\
                'HYDHEL rates/HYDHEL.tex\n'+\
                'AMJUEL rates/amjuel.tex\n\n'

    ## SETTINGS
    string +='** SETTINGS\n'+\
            '* vmax      14\n'+\
            '* n0\n'
    if incl_ground:
        string+= 'H2(n=X1,v=0)  1e10\n'
    string+='* verbose   0   # Show verbose output\n'

    with open(new_file_name, 'w') as f: 
        f.write(string)



def eval_1D(coeff, T):
    
    o = np.zeros(np.shape(T))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(T) ** i)
    return 1e-6 * np.exp(o)


def eval_2D(coeff, T, n):
    
    o = np.zeros(np.shape(T))
    for i in range(0, np.shape(coeff)[0]):
        for j in range(0, np.shape(coeff)[1]):
            o = o + coeff[i,j] * (np.log(T) ** i) * (np.log(n*1e-8) ** j) 
    return 1e-6 * np.exp(o)



def eval_1D_cs(coeff, E):
        
    o = np.zeros(np.shape(E))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(E) ** i)
    return 1e-4 * np.exp(o)


def calc_rates(E, sigma, Tev):
    import scipy.integrate as integrate
    import scipy.interpolate as interpolate

    e = 1.6e-19  #Elementary charge
    me = 9.11e-31 #electron mass
    kb = 1.38e-23 #Bolzmann constant

    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def vibr_dist(input_crm,iso_mass=1, Te_max=100,Te_reso=int(1e3),Te_min=0.1):
    import CRUMPET
    import numpy as np

    print('Running DEMO with input:')
    print(input_crm)

    indx_H2v = np.append(2,np.arange(3,17))

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    Tev = 10**np.linspace(np.log10(Te_min),np.log10(Te_max),Te_reso)
    Tiv = Tev/iso_mass # Assume Ti=Te
    ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)
    crm.source[2] = 1e-100 #add small source for numerical stability (only needed if reactions that dissociate are included)

    #compute vibrational distribution H2
    fv_H2 = np.zeros([15,len(Tev)])

    #calculate vibrational distribution using Tiv = Tev/iso_mass
    for i in range(0,len(Tev)):
        fv_H2[:,i]=crm.steady_state(Tev[i],ne[i],Ti=Tiv[i],plot=False,dt=True)[indx_H2v]
        print(i)

    #normalise vibrational distribution by dividing the distribution values to the sum of the distribution
    fv_H2 = fv_H2/(np.sum(fv_H2,axis=0)[None,:])

    return fv_H2, Tev

def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E

def fit_eval(coeffs, x):
    o = np.zeros(np.shape(x))
    for i in range(0,len(coeffs)):
        o = o + coeffs[i]*x**i
    return np.exp(o)


def err_bolzmann(fv_H2, T):
    x = 0.5+np.linspace(0,14,15)

    err = np.zeros(len(T))
    coeff = np.zeros(len(T))

    for i in range(len(T)):
        fit = np.flip(np.polyfit(x,np.log(fv_H2[:,i]), 1))
        err[i] = np.sqrt(np.sum(((fit_eval(fit,x)-fv_H2[:,i]))**2/len(T)))
        coeff[i] = fit[1]
    return err, coeff

def eff_rates(fv_H2, Tev, ne, iso_mass=2):
    indx_X1 = np.append(0,np.arange(1,14))
    indx_B1 = np.append(14,np.arange(15,29))
    indx_C1 = np.append(29,np.arange(30,44))
    indx_EF1 = np.append(44,np.arange(45,59))
    indx_a3 = np.append(59,np.arange(60,74))
    indx_c3 = np.append(74,np.arange(75,89))
    indx_d3 = np.append(89,np.arange(90,104))

    Tiv = Tev/iso_mass
    

    #Get vibrationally resolved molecular CX rates from H2VIBR
    import CRUMPET
    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr_custom.tex'})

    #Get rates as function of Te
    vibr_resolved_CX = np.zeros([15,len(Tev)])
    vibr_resolved_Diss = np.zeros([15,len(Tev)])
    vibr_resolved_DA = np.zeros([15,len(Tev)])
    vibr_resolved_DA_B1 = np.zeros([15,len(Tev)])

    states = [indx_X1, indx_B1, indx_C1, indx_EF1, indx_a3, indx_c3]
    vibr_resolved_Ion = np.zeros([15,len(Tev),len(states)])
    diss_decay_b3 = np.zeros([15, len(Tev)])

    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'Q6'],Tiv)

        vibr_resolved_Diss[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L1'],Tev)
        # diss_decay_b3[i,:] = 1/ne*crm.reactions['USER']['COEFFICIENT'][f'a3Sg_v={i}'].coeffs # Look up the units !!!

        vibr_resolved_DA[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L4'],Tev)
        vibr_resolved_DA_B1[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'Z4'],Tev)

        for j in range(len(states)):
            labels = ['2.'+str(i)+'L2', '2.'+str(i)+'B1','2.'+str(i)+'C2','4.'+str(i)+'L2', '5.'+str(i)+'L2', '7.'+str(i)+'L2']
            vibr_resolved_Ion[i,:,j] = eval_1D(X.reactions['H2VIBR']['H.2'][labels[j]],Tev)
    
    #Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
    eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2[indx_X1],axis=0)
    eff_mol_diss = np.sum(vibr_resolved_Diss*fv_H2[indx_X1], axis=0)
    eff_mol_DA = np.sum(vibr_resolved_DA*fv_H2[indx_X1]+vibr_resolved_DA_B1*fv_H2[indx_B1], axis=0)

    eff_mol_ion = np.zeros(len(Tev))
    for i, state in enumerate(states):
        eff_mol_ion += np.sum(fv_H2[state,:]*vibr_resolved_Ion[:,:,i],axis=0)


    return eff_mol_cx, eff_mol_diss, eff_mol_ion, eff_mol_DA

def get_Apq(crm, i_state, f_state):
    '''
    Get Einstein coefficients for transitions between i_state and f_state
    
    '''

    Apq = np.zeros((15,15))

    for i in range(15):
        for j in range(15):
            Apq[i,j] = crm.reactions['USER']['COEFFICIENT'][i_state+'_'+f_state+f'{i}to{j}'].coeffs
    return Apq

def get_R0p(fv,state,ne):
    '''
    Get population coefficients for state at density ne
    '''
    indx_X1 = np.append(2,np.arange(3,17))
    indx_B1 = np.append(17,np.arange(18,32))
    indx_C1 = np.append(32,np.arange(33,47))
    indx_EF1 = np.append(47,np.arange(48,62))
    indx_a3 = np.append(62,np.arange(63,77))
    indx_c3 = np.append(77,np.arange(78,92))
    indx_d3 = np.append(92,np.arange(93,107))
    indx_mol = np.append(2,np.arange(3,107))

    states = ['X1Sg','B1Su','C1Pu','EF1Sg','a3Sg','c3Pu','d3Pu']
    indx = [indx_X1,indx_B1,indx_C1,indx_EF1,indx_a3,indx_c3,indx_d3]

    for i,st in enumerate(states):
        if st==state:
            p = indx[i]

    R0p = fv[p]/(np.sum(fv[indx_mol],axis=0)[None,:]*ne)

    return R0p


def get_Xeff(crm,fv,i_state,f_state,ne):

    Apq = get_Apq(crm,i_state,f_state)
    R0p = get_R0p(fv,i_state,ne)

    E_i = get_energies(i_state)
    E_f = get_energies(f_state)

    table = np.zeros((2,15**2,np.shape(fv)[1]))

    k=0
    for i in range(15):
        for j in range(15):
            table[0,k,:] = E_i[i]-E_f[j]
            table[1,k,:] = Apq[i,j]*R0p[i,:]
            k+=1

    data = []
    for i in range(np.shape(fv)[1]):
        data.append(np.array(sorted(np.transpose(table[:,:,i]), key=lambda row: row[0])))


    return data
# %%
