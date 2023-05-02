import numpy as np

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

def reactions(arr):
    t = ''
    for i in range(15):
        t +=f'* MCCCDB k {i} \n'+\
            f'e + H2(v={i}) > e + H2(n=B)\n'+\
            f'rates\MCCC\X1Sg-excitation\\vi={i}\MCCC-el-D2-B1Su_bound.X1Sg_vi={i}.txt\n\n'+\
            f'* USER COEFFICIENT {i}\n'+\
            f'H2(n=B) > H2(v={i})\n'+\
            f'{arr[i]}\n\n'
    return t



def FC_factors(file_name):

    nu_eff = 7.7771e08

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
    
arr = FC_factors('Table 2 Franck-Condon Factors\D2_B1-X1_FCF.dat')

insert_string('input.dat', 'input_mccc.dat', reactions(arr), 71)


        


