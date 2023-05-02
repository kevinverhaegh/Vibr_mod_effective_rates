import numpy as np

def FC_factors(file_name):
    with open(file_name, 'r') as f:
        lines  = f.read().splitlines()
        data = []
        for line in lines:
            columns = line.split()
            data.append(columns)

    arr = np.array(data)[1:16, 1:16]

    table = np.zeros(np.shape(arr))
    for i in range(np.shape(arr)[0]):
        for j in range(np.shape(arr)[1]):
            table[i,j] = float(arr[i,j])

    coeffs = np.sum(table, axis = 1)/15

    return coeffs


coeff = FC_factors('Table 2 Franck-Condon Factors\D2_B1-X1_FCF.dat')
np.shape(coeff)
print(coeff)