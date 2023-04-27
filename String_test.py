import numpy as np

def output_string(A):
    exponent = str(f"{A:E}").split('E')[1]
    output = f"{A / 10**int(exponent):.12f}E{exponent}"
    if A>0:
        output = ' ' + output
    return output

def block_string(X):
    block = '  b0 ' + output_string(X[0]) + '  b1 ' + output_string(X[1]) + '  b2 ' + output_string(X[2])+'\n'+\
            '  b3 ' + output_string(X[3]) + '  b4 ' + output_string(X[4]) + '  b5 ' + output_string(X[5])+'\n'+\
            '  b6 ' + output_string(X[6]) + '  b7 ' + output_string(X[7]) + '  b8 ' + output_string(X[8])+'\n'
    return block


def numpy_to_string(i): 
    part_1 = "\subsection{\n" +\
                f"Reaction 2.{i}l2\n"+\
                f"$ p + H_2(v={i}) \\rightarrow H + H_2^+$ (ion conversion)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived from HYDHEL rate coeff. data.\n"+\
                "Taken at $E(H_2) = 0.1 \\approx 0.0$ eV,  and fit is for temperature $T_p=T$ with $H_2$ at rest.\n"+\
                "\n"+\
                "\\begin{small}\\begin{verbatim}\n"+\
                "\n"
                
    part_2 =   "\n"+\
                "\\end{verbatim}\end{small}\n"+\
                "\n"+\
                "\\newpage\n"
    return part_1, part_2


full_string = ''
coeffs = np.random.rand(15,9)

for i in range(15):
    str1, str2 = numpy_to_string(i)
    full_string += str1+block_string(coeffs[i,:])+str2


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

insert_string('Test.tex', 'Test_insert.tex', full_string, 4)

