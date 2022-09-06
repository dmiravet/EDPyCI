#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 10:23:47 2022

@author: dmiravet
"""
import numpy as np
import copy
import matplotlib.pyplot as plt
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from sys import getsizeof


class  HubbardHamiltonian:
      """
      Define and solve the Hubbard Hamiltonian for L sites 
      H = \sum_i c_{i,\sigma}^dagger c_{i+1,sigma} + \sum_i n_{i,uparray} n_{i,downarray} 
       """
     
      
      def __init__(self, L, npart, t, U,periodic = 0):
          #self.state = [0] * L
          
          #self.state = np.zeros(L,'int')
          self.L = L
          self.b = []
          self.dic = {}
          self.row = np.array((0), dtype=np.int)
          self.col = np.array((0), dtype=np.int)
          self.data = np.array((0), dtype=np.float)
          self.state = list(0 for _ in range(2*L))
          self.generateFockSpace(self.state, 0, npart)
          self.M = np.zeros((len(self.b),len(self.b)))
          self.Hopping(t, periodic)
          self.HubbardU(U)
          self.M_sparse = coo_matrix((self.data, (self.row, self.col)), shape=(len(self.b),len(self.b)))
          
          
      def generateFockSpace(self,state, i, npart):
    
          if i + npart > self.L*2:
             return
 
          if npart == 0:
               self.b.append(list(state))
               self.dic[tuple(state)] =  len(self.b)-1
               #print(tuple(state))
          else:
              self.generateFockSpace(state, i + 1, npart)
              state[i] = 1
              self.generateFockSpace(state, i + 1, npart-1)
              state[i] = 0
              
   
              
      def Hopping(self, t, periodic = 0):
         L = self.L;
         for s in range(2):
             shift = s*self.L
             for i in range(len(self.b)):          
                 for j in range(shift, shift + L-1):
                     if self.b[i][j] + self.b[i][j+1] == 1:
                         aux = self.b[i].copy()
                         if self.b[i][j] == 0:
                             aux[j] = 1
                             aux[j+1] = 0                      
                         else:
                                 aux[j+1] = 1
                                 aux[j] = 0             
                         self.M[i][self.dic[tuple(aux)]] += t
                         self.row = np.append( self.row,i ) 
                         self.col = np.append( self.col,self.dic[tuple(aux)])
                         self.data = np.append( self.data, t)
         if periodic == 1 and  L > 2:
            for s in range(2):
                shift = s*self.L
                for i in range(len(self.b)):
                    if self.b[i][shift + L - 1] + self.b[i][shift] == 1:
                        aux = self.b[i].copy()
                        phase = (-1)**(sum(self.b[i][shift:shift + L - 1])%2)
                        if self.b[i][shift] == 0:
                          aux[shift] = 1
                          aux[shift + L - 1] = 0
                        else:
                          aux[shift + L - 1] = 1
                          aux[shift] = 0             
                          phase = phase*(-1)
                        self.M[i][self.dic[tuple(aux)]] += t*phase
                        self.row = np.append( self.row,i ) 
                        self.col = np.append( self.col,self.dic[tuple(aux)])
                        self.data = np.append( self.data, t*phase)
                       
            
        
      def HubbardU(self,U):
          for i in range(len(self.b)):
              for j in range(0,L):
                  self.M[i][i] += U*self.b[i][j]*self.b[i][j+self.L]
                  self.row = np.append( self.row,i ) 
                  self.col = np.append( self.col,i)
                  self.data = np.append( self.data, U*self.b[i][j]*self.b[i][j+self.L])
                  
      def Diagonalize(self):
          print("Sparse size :",getsizeof(self.M_sparse))
          print("Dense size:",getsizeof(self.M))
          #return np.linalg.eigh(self.M)
          return eigsh(self.M_sparse,len(self.b)-2, which='SA')

L = 2 
npart = 2
t = 1
U = 1
HHam = HubbardHamiltonian(L, npart, t, U, 1)
w, v = HHam.Diagonalize()
print("Energies: ", w) 




#comparing with the analitical solution for U = 0, single partile case
L = 20 
npart = 1
t = 1
U = 0
HHam = HubbardHamiltonian(L, npart, t, U, 1)
w, v = HHam.Diagonalize()
#print("Energies: ", w)
k =  np.arange(0,L)
Ek = 2*t*np.cos(2*np.pi/L*k)
Ek2 = np.sort(np.repeat(Ek,2))
index =  np.arange(0,2*L)
plt.plot(index[:-2], w, 'ro')
plt.plot(index, Ek2, '>')
plt.xlabel("Index")
plt.ylabel("Energy")
plt.show()




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
