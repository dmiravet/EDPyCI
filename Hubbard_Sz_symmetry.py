#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:23:47 2022

@author: dmiravet
"""
import numpy as np
import copy


class  HubbardHamiltonian:
      """
      Define and solve the Hubbard Hamiltonian for L sites 
      H = \sum_i c_{i,\sigma}^dagger c_{i+1,sigma} + \sum_i n_{i,uparray} n_{i,downarray} 
       """
     
      
      def __init__(self, L, nup, ndn, t, U):
          #self.state = [0] * L
          
          #self.state = np.zeros(L,'int')
          self.L = L
          self.b = []
          self.dic = {}
          npart = nup + ndn
          self.state = list(0 for _ in range(2*L))
          self.generateFockSpace(self.state, 0, npart, nup, ndn)
          self.M = np.zeros((len(self.b),len(self.b)))
          self.Hopping(t)
          self.HubbardU(U)
          
          
      def generateFockSpace(self,state, i, npart, nup, ndn):
    
          if i + npart > self.L*2:
             return
 
          if npart == 0:
               if sum(state[:L]) == nup and sum(state[L:]) == ndn:
                   self.b.append(list(state))
                   self.dic[tuple(state)] =  len(self.b)-1
                   print(tuple(state))
          else:
              self.generateFockSpace(state, i + 1, npart,nup, ndn)
              state[i] = 1
              self.generateFockSpace(state, i + 1, npart-1,nup, ndn)
              state[i] = 0
              
      def Hopping(self, t):
       for s in range(2):
        shift = s*self.L
        for i in range(len(self.b)):          
          for j in range(shift,shift+L-1):
            if self.b[i][j] + self.b[i][j+1] == 1:
              aux = copy.copy(self.b[i])
              if self.b[i][j] == 0:
                aux[j] = 1
                aux[j+1] = 0                      
              else:
                aux[j+1] = 1
                aux[j] = 0             
              self.M[i][self.dic[tuple(aux)]] += t
        
      def HubbardU(self,U):
          for i in range(len(self.b)):
              for j in range(0,L):
                  self.M[i][i] += U*self.b[i][j]*self.b[i][j+self.L]
      def Diagonalize(self):
          return np.linalg.eigh(self.M)

 
L = 2   
npart = 2
nup = 2;
ndn = npart - nup
t = 1
U = 1;
HHam = HubbardHamiltonian(L,nup, ndn,t,U)
w, v = HHam.Diagonalize()
print("Energies: ", w)







#def generateFullFockSpace(n, state, i):
# 
#    if i == n:
#        printTheArray(state, n)
#        return
#    
#    generateFullFockSpace(n, state, i + 1)
#    state[i] = 1
#    generateFullFockSpace(n, state, i + 1)
#    state[i] = 0
#    
