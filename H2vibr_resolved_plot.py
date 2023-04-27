def run_demo(input_crm,iso_mass=1,T_request=1):
    print('Running DEMO with input:')
    print(input_crm)
    import CRUMPET
    import numpy as np 
    from scipy.optimize import curve_fit

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    Tev = np.linspace(0.2,10,100)
    Tiv = Tev # Assume Ti=Te
    ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)
    crm.source[2] = 1e-100 #add small source for numerical stability (only needed if reactions that dissociate are included)

    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr.tex', 'AMJUEL' : '/rates/amjuel.tex'})

    #Get rates as function of Te

    def eval_1D(coeff,T):
        o = np.zeros(np.shape(T))
        for i in range(0,len(coeff)):
            o = o + coeff[i]*(np.log(T)**i)
        return 1e-6*np.exp(o)

    vibr_resolved_CX = np.zeros([15,len(Tev)])
    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L2'],Tiv/iso_mass)


    # show a plot comparing different effective rates
    import matplotlib.pyplot as plt
    plt.figure()
    colors = plt.cm.rainbow(np.linspace(0, 1, 15))

    for i in range(15):
        plt.loglog(Tev,vibr_resolved_CX[i,:], c=colors[i])

    plt.xlabel('Temperature (eV)')
    plt.show()

run_demo('input.dat')