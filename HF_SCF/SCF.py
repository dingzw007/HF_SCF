def get_energy(basis):
    import numpy as np
    import psi4


    np.set_printoptions(suppress=True,precision =3)
    #Build a molecula
    mol=psi4.geometry('''
    O
    H 1 1.1
    H 1 1.1 2 104''')
    mol.update_geometry()
    mol.print_out()
    nel = 5

    #Build a basis

    bas=psi4.core.BasisSet.build(mol,target=basis)
    bas.print_out()

    #Build a Mintshelper
    mints=psi4.core.MintsHelper(bas)

    V=np.array(mints.ao_potential())
    T=np.array(mints.ao_kinetic())


    #Core Hamiltonian
    H = T + V

    S = np.array(mints.ao_overlap())
    g = np.array(mints.ao_eri())

    #print(S.shape)
    #print(I.shape)


    A = mints.ao_overlap()
    print (A)
    A.power(-0.5,1.e-14)
    A = np.array(A)

    #Fp =A.T@H@A
    #eps,cp = np.linalg.eigh(Fp)
    #C=A@cp

    #Cocc = C[:,:nel]
    #D = Cocc @ Cocc.T
    #Dagonalize Core H

    def diag(F,A):
        Fp =A.T@F@A
        eps,cp = np.linalg.eigh(Fp)
        C=A@cp
        return eps,C
    eps, C = diag(H,A)
    Cocc = C[:,:nel]
    D = Cocc @ Cocc.T
    E_old=0.0
    for iteration in range(5):

    #raise Exception("Breakpoint")
    # F = H + 2 * g_pqrs Drs-g_prqs Drs

    #
    #Jsum = np.sum(g*D,axis=(2,3))
        J = np.einsum("pqrs,rs->pq",g,D)
        K = np.einsum("prqs,rs->pq",g,D)

        F = H + 2.0*J-K
            
        grad = F@D@S-S@D@F
        grad_rms = np.mean(grad **2) **0.5
        
        E_electric =  np.sum((F+H)*D)
        E_total= E_electric+mol.nuclear_repulsion_energy()
        E_diff=E_total-E_old
        E_old =E_total
        print("% 16.2f % 8.4e % 8.4e"%( E_total,E_diff, grad_rms))
        eps, C = diag(F,A)
        Cocc = C[:,:nel]
        D = Cocc @ Cocc.T

    #print(F)
    #print(Jein) 
    print("SCF has finished!\n")
     
    return E_total
    #psi4.set_options({"scf_type":"pk"})
    #psi4_energy = psi4.energy("SCF/sto-3g",molecule=mol)

    #print("Energy matches")  
    #print("%s" % ( np.allclose(E_total,psi4_energy)))
