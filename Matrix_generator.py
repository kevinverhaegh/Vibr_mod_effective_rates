#%%
def matrix_generator(input_file, output_file, Te, ne):
    import CRUMPET
    crm = CRUMPET.Crumpet(input_file)
    

    import pickle
    pickle.dump(crm.getM(10,1e19)[0],open(output_file,'wb'))

matrix_generator('input_false.dat','m.pkl',10,1e19)

# %%
