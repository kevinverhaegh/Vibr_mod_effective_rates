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



insert_string('input.dat', 'input_mccc.dat', reactions(arr), )


        


