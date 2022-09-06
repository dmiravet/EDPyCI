import sys
sys.path.append('../') 
from QCHH_general import QuantumChemHamiltonian
import numpy as np
import matplotlib.pyplot as plt
       
def load_from_file( en_file, vee_file):
    """
    Load the coulomb matrix elements and single particle energy from files
    """
    Ei = np.loadtxt(en_file, unpack = True)
    n1,n2,n3,n4 = np.loadtxt(vee_file, unpack = True, skiprows = 1, max_rows = 1)
    print(int(n1),int(n2),int(n3),int(n4))
    Vhh = np.zeros((int(n1), int(n2), int(n3), int(n4)), dtype=np.float)
    #print( Vhh.shape())
    iv,jv,kv,lv,vv,_ = np.loadtxt(vee_file, unpack = True, skiprows = 2)
    for i in range(len(vv)):
        Vhh[int(iv[i])-1,int(jv[i])-1,int(kv[i])-1,int(lv[i])-1] =  vv[i] #using only real part
    return Ei,Vhh

L = 12 # number of orbitals
npart = 3 # number of particles

Ei, Vhh = load_from_file("data/En.dat", "data/Vhh.dat") # loading coulomb parameters from file
HHam = QuantumChemHamiltonian(L, npart, Ei , Vhh)
w, v = HHam.Diagonalize()
output_file = "data/spectrum_" + str(npart) + "holes_" + str(L) + "orbitals.dat"
np.savetxt(output_file, w)
plt.plot(w,'o')