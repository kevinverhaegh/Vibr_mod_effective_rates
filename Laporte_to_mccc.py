import numpy as np


def extract_table(lines, initial_state, final_state):

    for i,line in enumerate(lines): 
        if f'PROCESS: E + D2(X,v={initial_state}) -> E + D2(X,v={final_state}), Excitation' in line:
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

def extract_table_diss_att(lines, initial_state):

    for i,line in enumerate(lines): 
        if f'PROCESS: E + D2(B,v={initial_state}) -> D(n=2) + D-(1s2), Attachment' in line:
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


def numpy_array_to_string(arr):
    # Convert the numpy array to a string in the required format
    s = ''
    for row in arr:
        s += f"{row[0]:.5E}   {row[1]:.5E}\n"
    return s

# Exctract vibrational energies of deuterium
def D2_energies(file):
    E = np.zeros(15)

    with open(file, 'r') as f:
        lines = f.readlines()
        for i in range(15):
            line = lines[i+1]
            E[i] = float(line.split()[2])
    return E



def laporta_to_mccc_vibr_trans(file, E):

    a0 = 5.29177210903e-11

    with open(file, 'r') as f:
        lines = f.readlines()

        for i in range(15):
            for j in range(i,15): 
                table = extract_table(lines,i,j)

                # Unit conversion
                table[:,1] = a0**-2*table[:,1]

                string_ex = numpy_array_to_string(table)
                filename_ex = f'rates/Laporta/vibr_trans/vi={i}_vf={j}.txt'


                table_deex = table
                table_deex[:,1] = table[:,1]*np.exp((E[j]-E[i])/table[:,0])

                string_deex = numpy_array_to_string(table_deex)
                filename_deex = f'rates/Laporta/vibr_trans/vi={j}_vf={i}.txt'

                with open(filename_ex, 'w') as dest: 
                    dest.write(string_ex)

                with open(filename_deex, 'w') as dest: 
                    dest.write(string_deex)


def laporta_to_mccc_diss_attachment(file):
    
    a0 = 5.29177210903e-11

    with open(file, 'r') as f:
        lines = f.readlines()

        for i in range(15):
                table = extract_table_diss_att(lines,i)

                # Unit conversion
                table[:,1] = a0**-2*table[:,1]

                string = numpy_array_to_string(table)
                filename = f'rates/Laporta/diss_attachment_B1Su/vi={i}.txt'

                with open(filename, 'w') as dest: 
                    dest.write(string)





E = D2_energies('Fantz/Table 1 Vib Eigenvalues/X1_EV.txt')
# laporta_to_mccc_vibr_trans('rates/Laporta/vibr_trans/Cross section.txt', E) ##Add energies from franck condon paper
laporta_to_mccc_diss_attachment('rates/Laporta/diss_attachment_B1Su/Cross section.txt')

print('Done')
