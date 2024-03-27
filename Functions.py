#%%
from xml.etree.ElementTree import C14NWriterTarget
import numpy as np

def get_energies(state, vmax = 14,isotope=2):
    '''
    Retrieves vibrational eigenvalues (Energy levels) for molecular deuterium from Fantz et al.

        Parameters
        ----------
        state : string
            Determines which electronically excited state you want the energies of (e.g. X1Sg or EF1Sg).
        vmax: integer (default: 14)
            Maximum vibrational quantum number you want to retrieve.
        isotope: integer (default: 2)
            Isotope mass (1 = hydrogen, 2 = deuterium, 3 = tritium)
    '''

    E = np.zeros(vmax+1)
    with open('Fantz/Table 1 Vib Eigenvalues/'+state[:-2]+'_EV.txt', 'r') as f:
        lines = f.readlines()
        for i in range(vmax+1):
            line = lines[i+1].replace('*','0')
            E[i] = float(line.split()[isotope])
    return E


def species(state):
    '''
    Function that writes the 'species' string for the input file. 

        Parameters
        ----------
        state: string
            Determines which state you want to include in the string (e.g. X1Sg or B1Su)
    '''

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
    '''
    Retrieves radiative transition probabilities (Einstein coefficients) from Fantz et al. 

        Parameters
        ----------
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string 
            Specifies the final electronically excited state. 
    '''

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

    '''
    Retrieves radiative transition probabilities (Einstein coefficients) from Fantz et al.,
    in the case that the final state is repulsive (e.g. the b3Sg state).  

        Parameters
        ----------
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string 
            Specifies the final electronically excited state. 
    '''


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
    '''
    Function that creates the string to include radiative transitions in the input file. 

        Parameters
        ----------
        arr: numpy array
            array containing the transition probabilities between the initial and final state. 
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string
            Specifies the final electronic state. 

    '''

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
                
    if i_state=='d3Pu':
        for i in range(15):
            coeff = np.sum(arr[i,15:])
            t+=f'* USER COEFFICIENT '+i_state+'_'+f_state+f'{i}to15\n'+\
                    'H2(n='+i_state[:-2]+f',v={i}) > 2*H(n=1)\n'+\
                    f'{coeff}\n\n'
    return t

def rea_rad_decay_unr(arr, i_state, f_state):
    '''
    Function that creates a string to include radiative decay in the input file, 
    in the case that the electronically excited states are vibrationally unresolved. 

        Parameters
        ----------
        arr: numpy array
            array containing the transition probabilities between the initial and final state. 
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string
            Specifies the final electronic state. 
    '''

    t = ''
    for i in range(15):
        t +='* USER COEFFICIENT '+i_state+'_'+f_state+f'X{i}\n'+\
            'H2(n='+i_state[:-2]+') > H2(n='+f_state[:-2]+f',v={i})\n'+\
            f'{arr[i]}\n\n'
    return t

def get_coeffs_unr(i_state, f_state, nu_eff):
    '''
    Calculates transition probabilities using data from Fantz et al. This function can be used when the electronically
    excited states are vibrationally unresolved. It calculates the transition probabilities using the Franck-Condon
    factors and the effective transition frequency of the upper state. 

    Parameters
    ----------
    i_state: string
        Specifies the initial electronically excited state. 
    f_state: string
        Specifies the final electronic state. 
    nu_eff: float
        Effective transition frequency of the transition between the initial and final state. 


    '''

    with open('Fantz/Table 2 Franck-Condon Factors/D2_'+i_state[:-2]+'-'+f_state[:-2]+'_FCF.dat', 'r') as f:
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


def gen_input(new_file_name, vibr_resolved=True, rad_decay=True, inter_states = True, coll_deex = True, incl_ground=False, ion = False, ion_hyd=False, B1 = False,
              d3=False, c3 = False, vibr_hyd = False, diss = False, diss_hyd = False, a3=False, EF1=False, mol_cx=False, mol_cx_hyd=False,diss_ion=False,
              vibr_lap=False, C1=False, diss_att_X1=False, diss_att_B1 = False, diss_att_old=False, MA=False, MA_min=False,MolIonR=False):
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
            switch to include vibrational transitions in X1Sg and B1, using Laporta cross sections.
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
        if not MolIonR:
            if ion or mol_cx_hyd or mol_cx:
                string+='* H2+\n'+\
                    '    V 15.56\n'
            if diss_att_X1 or diss_att_B1 or diss_att_old:
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
        if MolIonR:
            string+='* H2+\n'+\
                    '    V 15.56\n'+\
                    '* H-\n' +\
                    '    V -0.75\n'

        string+='\n\n'

        ## REACTIONS
        string+='** REACTIONS\n\n'

        # Vibrational transitions
        if vibr_hyd:
            string+='* H2VIBR_OR H.2 2.$v&\n' +\
                    'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n'
        
        if vibr_lap:
            string+='* H2VIBR H.2 2.$t&\n' +\
                    'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n'
        
        if mol_cx:
            string+='* H2VIBR H.2 2.$q6\n'+\
                    'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n'

        if mol_cx_hyd:
            string+='* H2VIBR_OR H.2 2.$l2\n'+\
                    'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n'
        
        if diss:
            string+='* H2VIBR H.2 2.$l1\n'+\
                    'e + H2(n=X1,v=$) > e + 2*H(n=1)\n\n'

        if diss_hyd:
            string+='* H2VIBR_OR H.2 2.$l1\n'+\
                'e + H2(n=X1,v=$) > e + 2*H(n=1)\n\n'
        
        if ion: 
            string+='* H2VIBR H.2 2.$l2\n'+\
                    'e + H2(n=X1,v=$) > e + H2+\n\n'
        if ion_hyd:
            string+='* H2VIBR_OR H.2 2.$l4\n'+\
                    'e + H2(n=X1,v=$) > e + H2+\n\n'
            
        if diss_ion:
            string += '* AMJUEL H.4 2.2.10\n'+\
                        'e + H2(n=X1,v=0) > 2*e + p + H(n=1)\n\n'

            
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
                string+='* H2VIBR H.2 7.$l2\n'+\
                        'e + H2(n=c3,v=$) > e + H2+\n\n'
                
        if d3:
            # Excitation between states X1 and d3
            string+='* H2VIBR H.2 7.$u&\n'+\
                    'e + H2(n=X1,v=$) > e + H2(n=d3,v=&)\n\n'
            if coll_deex:
                string+='* H2VIBR H.2 7.$d&\n'+\
                        'e + H2(n=d3,v=$) > e + H2(n=X1,v=&)\n\n'
            if ion:
                string+='* H2VIBR H.2 12.$l2\n'+\
                        'e + H2(n=d3,v=$) > e + H2+\n\n'

            if inter_states:
                # Excitation between a3 and d3
                string+='* H2VIBR H.2 11.$u&\n'+\
                        'e + H2(n=a3,v=$) > e + H2(n=d3,v=&)\n\n'
                if coll_deex:
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
            

        if ion or ion_hyd:
            string+='* H2+\n'+\
                '    V 15.56\n'
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
                    '    V 15.975\n'+\
                    '* H(n=1)\n'+\
                    '    V 2.375\n'
        if not incl_ground:
            string+='* H2(n=X1,v=0)\n'+\
                    '    V 0.19652\n'

        string+='\n\n'

        ## REACTIONS
        string+='** REACTIONS\n\n'

        #  Vibrational transitions
        if vibr_lap:
            string+='* H2VIBR H.2 2.$t&\n' +\
                        'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n' 
        if vibr_hyd:
            string+='* H2VIBR_OR H.2 2.$v&\n' +\
                        'e + H2(n=X1,v=$) > e + H2(n=X1,v=&)\n\n' 
        
        # Molecular charge exchange
        if mol_cx:
            string +='* H2VIBR H.2 2.$q6\n'+\
                        'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n'

        if mol_cx_hyd:
            string +='* H2VIBR_OR H.2 2.$l2\n'+\
                        'p + H2(n=X1,v=$) > H2+ + H(n=1)\n\n'
        
        # Ionization from ground
        if ion: 
            string+='* H2VIBR H.2 2.$l2\n'+\
                    'e + H2(n=X1,v=$) > e + H2+\n\n'

        if ion_hyd:
            string += '* H2VIBR_OR H.2 2.$l4\n' + \
                          'e + H2(n=X1,v=$) > e + H2+\n\n'

        # Dissociation from ground
        if diss:
            string+='* H2VIBR H.2 2.$l1\n'+\
                    'e + H2(n=X1,v=$) > e + 2*H(n=1)\n\n'

        if diss_hyd:
            string+='* H2VIBR_OR H.2 2.$l1\n'+\
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
                B_X_rate = get_coeffs_unr('B1Su', 'X1Sg', 7.7771e+08)
                string += rea_rad_decay_unr(B_X_rate,'B1Su', 'X1Sg')

            # Excitation and emission from C1Pu to ground
            if C1: 
                C_X_rate = get_coeffs_unr('B1Su', 'X1Sg', 1.0532e+09)
                string += rea_rad_decay_unr(C_X_rate, 'C1Pu', 'X1Sg')


                    
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
                'AMJUEL rates/amjuel.tex\n'
    if mol_cx_hyd or vibr_hyd or diss_hyd:
        string += 'H2VIBR_OR rates/h2vibr.tex\n\n'
    else:
        string += '\n'

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
    '''
    Function to evaluate 1D polynomial logarithmic fits of reaction rate coefficients of the type in Eirene data files (e.g. AMJUEL).
    In this case, the rates only have a temperature dependence.

        Parameters
        ----------
        coeff: numpy array
            Fit constants of the rate coefficients. (For Eirene data, this is a 1x8 array)
        T: float
            Evaluation temperature in eV (may also be an array of temperatures).

        Output
        ------
        Reaction rate coefficient(s) in m^3s^-1. 
    '''
    o = np.zeros(np.shape(T))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(T) ** i)
    return 1e-6*np.exp(o)


def eval_2D(coeff, T, n):
    '''
    Function to evaluate 2D polynomial logarithmic fits of reaction rate coefficients of the type in Eirene data files (e.g. AMJUEL). 
    In this case, the rate coefficients have both a temperature and density dependence. 

        Parameters
        ----------
        coeff: numpy array
            Fit constants of the rate coefficients. (For Eirene data, this is an 8x8 array)
        T: float
            Evaluation temperature in eV (may also be an array of temperatures). 
        n: float
            Evaluation density in m^-3 (may also be an array of densities). 

        Output
        ------
        Reaction rate coefficient(s) in m^3s^-1.
    '''
    
    o = np.zeros(np.shape(T))
    for i in range(0, np.shape(coeff)[0]):
        for j in range(0, np.shape(coeff)[1]):
            o = o + coeff[i,j] * (np.log(T) ** i) * (np.log(n*1e-14) ** j) 
    return 1e-6*np.exp(o)



def eval_1D_cs(coeff, E):
    '''
    Function to evaluate 1D polynomial logarithmic fits of cross sections sigma(T) of the type in Eirene data files (e.g. AMJUEL).
    The cross sections are only dependent on temperature.

        Parameters
        ----------
        coeff: numpy array
            Fit constants of the cross sections. (For Eirene data, this is a 1x8 array)
        T: float
            Evaluation temperature in eV (may also be an array of temperatures).

        Output
        ------
        Cross sections in m^2
    '''
        
    o = np.zeros(np.shape(E))
    for i in range(0, len(coeff)):
        o = o + coeff[i] * (np.log(E) ** i)
    return 1e-4*np.exp(o)


def calc_rates(E, sigma, Tev):
    import scipy.integrate as integrate
    import scipy.interpolate as interpolat
    '''
    Calculates rate coefficients for electron collisions from cross sectional data by integrating over a Bolzmann distribution.
    Assumes stationary molecules.  

        Parameters
        ----------
        E: numpy array
            Array containing the collision energies at which the cross sections are given (in eV).  
        sigma: numpy array
            Array containing the cross sections at the aforementioned energies (in m^2).
        Tev: float
            Evaluation temperature (in eV).
        
        Output
        ------
        Returns rate coefficient in m^3s^-1.
    

    '''

    e = 1.6e-19  #Elementary charge
    me = 9.11e-31 #electron mass
    kb = 1.38e-23 #Bolzmann constant

    v = np.sqrt(2*E*e/me)
    f = (me/(2*np.pi*e*Tev))**(3/2)*4*np.pi*v**2*np.exp(-me*v**2/(2*e*Tev))  # Bolzmann velocity distribution

    rate = integrate.trapz(sigma*v*f, x=v)
    return(rate)

def vibr_dist(crm, Tev, ne, iso_mass=2, incl_ground=False):
    import CRUMPET
    import numpy as np

    # Mass rescaling of ion temperature
    Tiv = Tev/iso_mass
    if not np.shape(ne):
        ne = np.array([ne])
    if not np.shape(Tev):
        Tev = np.array([Tev])

    #compute vibrational distribution H2
    # if np.shape(ne) and (not np.shape(Tev)):
    #     fv_H2 = np.zeros([len(crm.species),len(ne)])
    #     for i in range(0,len(ne)):
    #         fv_H2[:,i]=crm.steady_state(Tev,ne[i]*1e-6,Ti=Tiv,plot=False,dt=True)

    # elif (np.shape(Tev) and (not np.shape(ne)))==True:
    #     fv_H2 = np.zeros([len(crm.species),len(ne)])
    #     for i in range(0,len(Tev)):
    #         fv_H2[:,i]=crm.steady_state(Tev[i],ne*1e-6,Ti=Tiv[i],plot=False,dt=True)

    # elif np.shape(Tev) and np.shape(ne):
    #     fv_H2 = np.zeros([len(crm.species),len(Tev),len(ne)])

    fv_H2 = np.zeros([len(crm.species),len(Tev),len(ne)])
    for i in range(0,len(Tev)):
        for j in range(0,len(ne)):
            fv_H2[:,i,j]=crm.steady_state(Tev[i],ne[j]*1e-6,Ti=Tiv[i],plot=False,dt=True)
            # print(str(i),str(j))
    # else:
    #     fv_H2 = crm.steady_state(Tev,ne*1e-6,Ti=Tiv,plot=False,dt=True)


    indx_mol = np.arange(0,105)
    
    if not incl_ground:
        fv = np.zeros((len(crm.species)+1,len(Tev),len(ne)))
        for i in range(len(Tev)):
            for j in range(len(ne)):
                fv[:,i,j]=np.append(1,fv_H2[:,i,j])
                # Renormalization
                fv[:,i,j] = fv[:,i,j]/np.sum(fv[indx_mol,i,j])


    # ind = np.shape(fv_H2).index(1)
    # fv_H2 = fv_H2.reshape([])

    return fv


def fit_eval(coeffs, x):
    '''
    Evaluates logarithmic fit. 

        Parameters
        ----------
        coeffs: numpy array
            Fit constants
        x: float 
            Independent variable at which to evaluate fit. 

        Output
        ------
        Value of the function at x
    '''
    o = np.zeros(np.shape(x))
    for i in range(0,len(coeffs)):
        o = o + coeffs[i]*x**i
    return np.exp(o)

def err_bolzmann(fv_H2, E):
    '''
    Calculates fit constants and rms error from a Bolzmann fit

        Parameters
        ----------
            fv_H2: numpy array
                Vibrational distribution as a function of temperature in a given vibrational state.
            E: numpy array
                Energies corresponding to the different vibrational levels. 
    '''
    err = np.zeros(np.shape(fv_H2)[1])
    coeff = np.zeros(np.shape(fv_H2)[1])

    for i in range(np.shape(fv_H2)[1]):
        fit = np.flip(np.polyfit(E,np.log(fv_H2[:,i]), 1))
        err[i] = np.sqrt(np.sum(((fit_eval(fit,E)-fv_H2[:,i]))**2/np.shape(fv_H2)[1]))
        coeff[i] = fit[1]
    return err, coeff

def eff_rates(fv_H2, Tev, iso_mass=2, incl_ground=False, get_vibr_res = False):
    '''
    Calculates effective rates using the calculated distribution and the rate data from the CRM. 

        Parameters 
        ----------
        fv_H2: numpy array
            Particle distribution.
        Tev: float
            Electron temperature (eV). May also be an array of different temperatures. 
        iso_mass: float (default: 2)
            Isotope mass (=2 for deuterium).
        incl_ground: bool (default: False)
            Switch to signify whether the vibrational ground state of D2 was used as a regular species (incl_ground = True), 
            or as a background species (incl_ground = False). 
        
        Output
        ------
        Returns the effective rate coefficients (in m^3s^-1) of ionization, dissociation, molecular charge exchange, and dissociative
        attachment. 
    '''


    indx_X1 = np.append(0,np.arange(1,15))
    indx_B1 = np.append(15,np.arange(16,30))
    indx_C1 = np.append(30,np.arange(31,45))
    indx_EF1 = np.append(45,np.arange(46,60))
    indx_a3 = np.append(60,np.arange(61,75))
    indx_c3 = np.append(75,np.arange(76,90))
    indx_d3 = np.append(90,np.arange(91,105))

    Tiv = Tev/iso_mass
    
    # if incl_ground:
    #     fv = fv_H2
    # else:
    #     fv = np.zeros((np.shape(fv_H2)[0]+1,np.shape(fv_H2)[1]))
    #     for i in range(np.shape(fv_H2)[1]):
    #         fv[:,i] = np.append(1,fv_H2[:,i])

    #Get vibrationally resolved molecular CX rates from H2VIBR
    import CRUMPET
    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr_custom.tex'})

    #Get rates as function of Te
    vibr_resolved_CX = np.zeros([15,np.shape(fv_H2)[1]])
    vibr_resolved_Diss = np.zeros([15,np.shape(fv_H2)[1]])
    vibr_resolved_DA = np.zeros([15,np.shape(fv_H2)[1]])
    vibr_resolved_DA_B1 = np.zeros([15,np.shape(fv_H2)[1]])

    states = [indx_X1, indx_B1, indx_C1, indx_EF1, indx_a3, indx_c3, indx_d3]
    vibr_resolved_Ion = np.zeros([15,np.shape(fv_H2)[1],len(states)])
    diss_decay_b3 = np.zeros([15, np.shape(fv_H2)[1]])

    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'Q6'],Tiv)

        vibr_resolved_Diss[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L1'],Tev)
        # diss_decay_b3[i,:] = 1/ne*crm.reactions['USER']['COEFFICIENT'][f'a3Sg_v={i}'].coeffs # Look up the units !!!

        vibr_resolved_DA[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L4'],Tev)
        vibr_resolved_DA_B1[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'Z4'],Tev)

        for j in range(len(states)):
            labels = ['2.'+str(i)+'L2', '2.'+str(i)+'B1','2.'+str(i)+'C2','4.'+str(i)+'L2', '5.'+str(i)+'L2', '7.'+str(i)+'L2', '12.'+str(i)+'L2']
            vibr_resolved_Ion[i,:,j] = eval_1D(X.reactions['H2VIBR']['H.2'][labels[j]],Tev)
    
    #Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
    eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2[indx_X1],axis=0)
    eff_mol_diss = np.sum(vibr_resolved_Diss*fv_H2[indx_X1], axis=0)
    eff_mol_DA = np.sum(vibr_resolved_DA*fv_H2[indx_X1]+vibr_resolved_DA_B1*fv_H2[indx_B1], axis=0)

    eff_mol_ion = np.zeros(fv_H2.shape[1])
    for i, state in enumerate(states):
        eff_mol_ion += np.sum(fv_H2[state,:]*vibr_resolved_Ion[:,:,i],axis=0)

    if get_vibr_res:
        return eff_mol_cx, eff_mol_diss, eff_mol_ion, eff_mol_DA, vibr_resolved_CX, vibr_resolved_Diss, vibr_resolved_Ion, vibr_resolved_DA, vibr_resolved_DA_B1
    else:
        return eff_mol_cx, eff_mol_diss, eff_mol_ion, eff_mol_DA

def get_Apq(crm, i_state, f_state):
    '''
    Get Einstein coefficients for transitions between i_state and f_state from the crm 

        Parameters
        ----------
        crm: CRUMPET CRM object
            CRM object generated using CRUMPET and a certain input file. 
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string 
            Specifies the final electronically excited state. 
        
        Output
        ------
        An array containing the Einstein coefficients for the transition between the initial and final state
    '''

    Apq = np.zeros((15,15))

    for i in range(15):
        for j in range(15):
            Apq[i,j] = crm.reactions['USER']['COEFFICIENT'][i_state+'_'+f_state+f'{i}to{j}'].coeffs
    return Apq

def get_R0p(fv, state, ne, incl_ground=False):
    '''
    Get population coefficients for state at density ne. 

        Parameters 
        ----------
        fv: numpy array
            Particle distribution.
        state: string
            Specifies which state you want the population coefficients of. 
        ne: float
            Electron density (m^-3)
        incl_ground: bool (default: False)
            Switch to signify whether the vibrational ground state of D2 was used as a regular species (incl_ground = True), 
            or as a background species (incl_ground = False). 
    
    Returns population coefficient R0p = np/(n0*ne) in m^3

    '''
    indx_X1 = np.append(0,np.arange(1,15))
    indx_B1 = np.append(15,np.arange(16,30))
    indx_C1 = np.append(30,np.arange(31,45))
    indx_EF1 = np.append(45,np.arange(46,60))
    indx_a3 = np.append(60,np.arange(61,75))
    indx_c3 = np.append(75,np.arange(76,90))
    indx_d3 = np.append(90,np.arange(91,105))

    states = ['X1Sg','B1Su','C1Pu','EF1Sg','a3Sg','c3Pu','d3Pu']
    indx = [indx_X1,indx_B1,indx_C1,indx_EF1,indx_a3,indx_c3,indx_d3]

    for i,st in enumerate(states):
        if st==state:
            p = indx[i]
    
    # if incl_ground:
    #     R0p = fv[p]/(fv[0]*ne)
    # else:
        # fv_new = np.zeros((np.shape(fv)[0]+1,np.shape(fv)[1]))
        # for i in range(np.shape(fv)[1]):
        #     fv_new[:,i] = np.append(1,fv[:,i])
        # R0p = fv_new[p]/ne
    R0p = fv[p]/ne


    

    return R0p


def get_Xeff(crm,fv,i_state,f_state,ne, vmax = 14):
    '''
    Calculates effective spontaneous emission rate coefficients.

        Parameters
        ----------
        crm: CRUMPET CRM object
            CRM object generated using CRUMPET and a certain input file. 
        fv: numpy array
            Array containing the vibrational distribution
        i_state: string
            Specifies the initial electronically excited state. 
        f_state: string 
            Specifies the final electronically excited state. 
        ne: float
            Electron density (in m^-3)
        
        Output
        ------
        An array containing the energies of the different transitions (in eV) and the corresponding
        effective spontaneous emission rate coefficients (in m^3s^-1).

    '''

    Apq = get_Apq(crm,i_state,f_state)
    R0p = get_R0p(fv,i_state,ne)

    E_i = get_energies(i_state)
    E_f = get_energies(f_state)

    table = np.zeros((2,(vmax+1)**2,np.shape(fv)[1]))

    k=0
    for i in range(vmax+1):
        for j in range(vmax+1):
            table[0,k,:] = E_i[i]-E_f[j]
            table[1,k,:] = Apq[i,j]*R0p[i,:]
            k+=1

    data = []
    for i in range(np.shape(fv)[1]):
        data.append(np.array(sorted(np.transpose(table[:,:,i]), key=lambda row: row[0])))


    return data

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
    with open(file_name_new, 'r') as f:
        lines = f.readlines()

    with open(file_name_new, 'r+') as f:
        for i, line in enumerate(lines):
            if i == line_num:
                f.write(string + '\n')
            f.write(line)

def get_teq(crm, T_final, ne, tau, delta=1e-3, T_init=None, iso_mass=2):

    if T_init==None:
        f_in=np.zeros(len(crm.species))
    else:
        f_in = crm.steady_state(T_init, ne/1e6, Ti=T_init/2, plot=False, dt=True)

    f_s = crm.steady_state(T_final, ne/1e6, Ti=T_final/2, plot=False, dt=True)

    crm.source[:]=0
    sol = crm.solve_crm(tau,T_final,ne/1e6,Ti=T_final/2, gl=False, n=f_in, densonly = True, Qres=False)
    t = sol.t
    f_fin = sol.y


    i=1
    tol=0
    try: 
        while tol<delta:
            tol = abs((np.sum(f_s)-np.sum(f_fin[:,-i]))/np.sum(f_s))
            teq=t[-i]
            i+=1
    except IndexError:
        teq = 0
    return teq

# import CRUMPET
# crm=CRUMPET.Crumpet('input_fin.dat')
# get_teq(crm,1,1e19,tau=1e-3)


# %%
