#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 11:42:45 2022

@author: dmiravet
"""

import numpy as np
import copy

from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from sys import getsizeof


class  QuantumChemHamiltonian:
      """
     Quantum Chemical Hamiltonian definition in terms of creation and anihilation operators
     H = sum_i En[i] C^\dagger_i C_i + \sum_{i,j,k,l} V_{i,j,k,l} C^\dagger_i C^\dagger_j C_k C_l
       """
     
      
      def __init__(self, L, npart, En, Vee):
          #self.state = [0] * L
          
          #self.state = np.zeros(L,'int')
          self.L = L
          self.b = []
          self.dic = {}
          self.row = np.array((0), dtype=np.int)
          self.col = np.array((0), dtype=np.int)
          self.data = np.array((0), dtype=np.float)
          self.state = list(0 for _ in range(L))
          self.generateFockSpace(self.state, 0, npart)
                   
          dim = len(self.b)
          print("Hilbert space size:",dim)
          print("Generating matrix:")
          self.kinetic(En)
          self.interacting(Vee)
          print("Matrix generated")
          self.M_sparse = coo_matrix((self.data, (self.row, self.col)), shape=(dim,dim))          
          
      def generateFockSpace(self,state, i, npart):
    
          if i + npart > self.L:
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
              
   
        

      def create(self, state, i):
          if state[i] == 1:
              return 0,0
          else :
              aux = state.copy()
              aux[i] = 1
              if sum(aux[:i])%2 == 0:
                  sg = 1
              else:
                  sg = -1;
              return sg,aux
          
      def destroy(self, state, i):
          if state[i] == 0:
              return 0,0
          else :
              aux = state.copy()
              aux[i] = 0
              if sum(aux[:i])%2 == 0:
                  sg = 1
              else:
                  sg = -1;
              return sg,aux
          
              
      def interacting(self,Vee):
         L = self.L;
         tol = 1e-8*Vee[0,0,0,0]
         for s in range(len(self.b)):  
             for l in range(L):
                 sgl,auxl = self.destroy(self.b[s],l)
                 if sgl != 0:
                     for k in range(l):
                         sgk,auxk = self.destroy(auxl,k)
                         if sgk != 0:
                             for i in range(L):
                                 sgi,auxi = self.create(auxk,i)
                                 if sgi != 0:
                                     for j in range(i):
                                         V = Vee[i,j,k,l]- Vee[i,j,l,k]
                                         if abs(V) > tol:
                                             sgj,auxj = self.create(auxi,j)
                                             if sgj != 0:
                                                 sgT = -sgi*sgj*sgk*sgl
                                                 #print(self.dic[tuple(auxj)],self.dic[tuple(self.b[s])],V*sgT)
                                                 self.row = np.append( self.row,s) 
                                                 self.col = np.append( self.col,self.dic[tuple(auxj)])
                                                 self.data = np.append( self.data, V*sgT)         
       
                       
            
        
      def kinetic(self,En):
          for s in range(len(self.b)):#
              for j in range(self.L):
                  sg1,aux1 = self.destroy(self.b[s],j)
                  if sg1 != 0:
                    sg2,aux2 = self.create(aux1,j)
                    self.row = np.append( self.row,s ) 
                    self.col = np.append( self.col,s)
                    self.data = np.append( self.data, sg1*sg2*En[j])

                  
      def Diagonalize(self):
          print("Diagonalization...")
          print("Sparse size :", getsizeof(self.M_sparse))
          return eigsh(self.M_sparse,len(self.b)-2, which='SA')
          
      

