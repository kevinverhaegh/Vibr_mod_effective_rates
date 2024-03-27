#This script can load in the Ichihara rates, fit it and write it into an Eirene file

#load in ichihara data
import os, mat73

ichi_tables = os.path.join(os.path.dirname(__file__), 'MolCX_Ichi.mat')
ichi_tab = mat73.loadmat(ichi_tables)

#create interpolation object in terms of T, E and vib_energy
import Functions
import numpy as np

E_H2 = np.insert(Functions.get_energies('X1sg',isotope=1),0,0)
E_D2 = np.insert(Functions.get_energies('X1sg',isotope=2),0,0)
E_T2 = np.insert(Functions.get_energies('X1sg',isotope=3),0,0)

import scipy.interpolate as interpolate

Ichi_H2 = np.dstack((ichi_tab['MolCX_Ichi_01'],ichi_tab['MolCX_Ichi_1'],ichi_tab['MolCX_Ichi_5']))
f = interpolate.RegularGridInterpolator((E_H2[0:15], np.log(ichi_tab['T']),np.log([0.1,1,5])), np.log(Ichi_H2),method='linear',bounds_error=False,fill_value=None)

#evaluate & fit rate
def ichi_fit(Tv,input_rates):
    fit = np.zeros((15,9))
    for i in range(15):
        fit[i,:] = np.flip(np.polyfit(np.log(Tv),np.log(input_rates[i,:]/1e-6),8))
    return fit

E_mol_H2 = 0.1
Tv = np.linspace(0.1,50,1000)
E_H2_m, T_m, = np.meshgrid(E_H2,Tv,indexing='ij')
Ichi_H2_v = 1e-6 * np.exp(f((np.vstack([E_H2_m.ravel(),np.log(T_m.ravel()),np.log(np.full(T_m.size,0.1))]).T))).reshape(len(E_H2),len(Tv))
p_Ichi_H2_v = ichi_fit(Tv, Ichi_H2_v)

E_D2_m, T_m, = np.meshgrid(E_D2,Tv,indexing='ij')
Ichi_D2_v = 1e-6 * np.exp(f((np.vstack([E_D2_m.ravel(),np.log(T_m.ravel()),np.log(np.full(T_m.size,0.1))]).T))).reshape(len(E_H2),len(Tv))
p_Ichi_D2_v = ichi_fit(Tv, Ichi_D2_v)

E_T2_m, T_m, = np.meshgrid(E_T2,Tv,indexing='ij')
Ichi_T2_v = 1e-6 * np.exp(f((np.vstack([E_D2_m.ravel(),np.log(T_m.ravel()),np.log(np.full(T_m.size,0.1))]).T))).reshape(len(E_H2),len(Tv))
p_Ichi_T2_v = ichi_fit(Tv, Ichi_T2_v)

#store rates

#functions for handling

def numpy_to_string_H(i):
    part_1 = "\subsection{\n" + \
             f"Reaction 2.{i}q5\n" + \
             f"$ p + H_2(v={i}) \\rightarrow H + H_2^+$ (ion conversion)\n" + \
             "}\n" + \
             "Rate coeff. for E(H_2) = 0.1 eV, obtained from Ichihara, et al. 2002 as function of the ion temperature.\n" + \
             "\n" + \
             "\\begin{small}\\begin{verbatim}\n" + \
             "\n"

    part_2 = "\n" + \
             "\\end{verbatim}\end{small}\n" + \
             "\n" + \
             "\\newpage\n"
    return part_1, part_2

def numpy_to_string_D(i):
    part_1 = "\subsection{\n" + \
             f"Reaction 2.{i}q6\n" + \
             f"$ D+ + D_2(v={i}) \\rightarrow D + D_2^+$ (ion conversion)\n" + \
             "}\n" + \
             "Rate coeff. for E(D_2) = 0.1 eV, obtained from Ichihara, et al. 2002 as function of the ion temperature.\n" + \
             "In this, the ion temperature has been mass rescaled to different isotope mass, as well as the vibrational energy.\n" + \
             "\n" + \
             "\\begin{small}\\begin{verbatim}\n" + \
             "\n"

    part_2 = "\n" + \
             "\\end{verbatim}\end{small}\n" + \
             "\n" + \
             "\\newpage\n"
    return part_1, part_2

def numpy_to_string_T(i):
    part_1 = "\subsection{\n" + \
             f"Reaction 2.{i}q7\n" + \
             f"$ T+ + T_2(v={i}) \\rightarrow T + T_2^+$ (ion conversion)\n" + \
             "}\n" + \
             "Rate coeff. for E(T_2) = 0.1 eV, obtained from Ichihara, et al. 2002 as function of the ion temperature.\n" + \
             "In this, the ion temperature has been mass rescaled to different isotope mass, as well as the vibrational energy.\n" + \
             "\n" + \
             "\\begin{small}\\begin{verbatim}\n" + \
             "\n"

    part_2 = "\n" + \
             "\\end{verbatim}\end{small}\n" + \
             "\n" + \
             "\\newpage\n"
    return part_1, part_2

def gen_full_string(coeffs_H,coeffs_D,coeffs_T):
    full_string = ''
    for i in range(15):
        str1, str2 = numpy_to_string_H(i)
        full_string += str1+Functions.block_string(coeffs_H[i,:])+str2
    for i in range(15):
        str1, str2 = numpy_to_string_D(i)
        full_string += str1+Functions.block_string(coeffs_D[i,:])+str2
    for i in range(15):
        str1, str2 = numpy_to_string_T(i)
        full_string += str1+Functions.block_string(coeffs_T[i,:])+str2
    return full_string

full_string = gen_full_string(p_Ichi_H2_v,p_Ichi_D2_v,p_Ichi_T2_v)
Functions.insert_string('rates/h2vibr.tex', 'rates/h2vibr_ichi_2.tex', full_string, 2565)
