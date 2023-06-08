import Functions as fn
import numpy as np
import CRUMPET
import matplotlib.pyplot as plt

# Get cross section coefficients
crm = CRUMPET.Crumpet('input_false.dat')
coeffs = crm.reactions['HYDHEL']['H.1']['7.2.2'].coeffs

# Calculate rates as function of T
Tev = np.linspace(0.1,100,100)
E = 10**np.linspace(np.log10(2.25e2), np.log10(1e5),10000)
sv = np.zeros(len(Tev))
for i in range(len(Tev)):
    sv[i] = fn.calc_rates(E,fn.eval_1D_cs(coeffs,E),Tev[i]/2)

sv[sv==0]=1e-300


plt.plot(Tev,sv)

# Fit rates to get coefficients
fit = np.flip(np.polyfit(np.log(Tev[20:]),np.log(sv[20:]/1e-6),8))
print(fit)

plt.plot(Tev,fn.eval_1D(fit,Tev),'--')
plt.show()

# Add to custom H2VIBR file
def numpy_to_string(): 
    part_1 = "\subsection{\n" +\
                f"Reaction 1_min\n"+\
                f"$ e + H^- \\rightarrow e + p + H $ (ion desctruction)\n" +\
                "}\n"+\
                "Rate coeff. for H2 at rest, derived by integrating HYDHEL cross sections using iso_mass = 2\n"+\
                "\n"+\
                "\\begin{small}\\begin{verbatim}\n"+\
                "\n"
                
    part_2 =   "\n"+\
                "\\end{verbatim}\end{small}\n"+\
                "\n"+\
                "\\newpage\n"
    return part_1, part_2



p1,p2=numpy_to_string()
string = p1+fn.block_string(fit)+p2

print(string)

# with open('rates/h2vibr_custom.tex','a') as file:
#     file.write(string)




