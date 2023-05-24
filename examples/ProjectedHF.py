#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:39:39 2023

@author: dmiravet
"""

import sys
sys.path.append('../') 
from QCHH_withHopping import QuantumChemHamiltonianWithHopping
import numpy as np
import matplotlib.pyplot as plt
       
def load_from_file( t_file, vee_file,nF,nD):
    """
    Load the coulomb matrix elements and single particle energy from files
    """
    iv,jv,tv = np.loadtxt(t_file, unpack = True, skiprows = 0)
    nD1 = int(np.sqrt(len(iv)))
    
    
    t_matrix = np.zeros((nD,nD), dtype=np.float)
    for m in range(len(tv)):
        i = int(iv[m])-1 -nF
        j = int(jv[m])-1 -nF
        t_matrix[i,j] =  tv[m]
        
        
       
    iv,jv,kv,lv,  vv = np.loadtxt(vee_file, unpack = True, skiprows = 0)
    Vhh = np.zeros((nD,nD,nD,nD), dtype=np.float)
    
    
    for m in range(len(vv)):
        i = int(iv[m])-1 - nF
        j = int(jv[m])-1 - nF
        k = int(kv[m])-1 - nF
        l = int(lv[m])-1 - nF
        Vhh[i,j,k,l] =  vv[m]
        
    return t_matrix,Vhh
nF = 99
nD = 4
npart = 6
L = 2*nD
t, Vee = load_from_file("data/tCarbonTriangle.dat", "data/VCarbonTriangle.dat", nF, nD) # loading coulomb parameters from file

HHam = QuantumChemHamiltonianWithHopping(L, npart, t, Vee)
w, v = HHam.Diagonalize(50)
output_file = "data/spectrum_" + str(npart) + "holes_" + str(nD) + "orbitals.dat"
np.savetxt(output_file, w)
plt.plot(w,'o')