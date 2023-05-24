#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 13:29:58 2023

@author: dmiravet
"""

import numpy as np
import copy

from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh
from sys import getsizeof


class  QuantumChemHamiltonianWithHopping:
      """
     Quantum Chemical Hamiltonian definition in terms of creation and anihilation operators
     H = sum_{i,j,s} t[i,j] C^\dagger_{i,s} C_{j,s} + \sum_{i,j,k,l,s1,s2} V_{is1, js2, ks2, ls1} C^\dagger_{is1} C^\dagger_{js2} C_{ks2} C_{ls1}
       """
     
     
      def __init__(self, L, npart, t, Vee):
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
          self.Vee = Vee
          self.t = t
          dim = len(self.b)
          print("Hilbert space size:",dim)
          print("Generating matrix:")
          self.kinetic()
          self.interacting()
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
          
      def destroy(self, state,i):
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
          
              
      def interacting(self):
         L = self.L;
         tol = 1e-8*self.V_Coulomb(0,0,0,0)
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
                                         V = self.V_Coulomb(i,j,k,l)- self.V_Coulomb(i,j,l,k)
                                         if abs(V) > tol:
                                             sgj,auxj = self.create(auxi,j)
                                             if sgj != 0:
                                                 sgT = -sgi*sgj*sgk*sgl
                                                 #print(self.dic[tuple(auxj)],self.dic[tuple(self.b[s])],V*sgT)
                                                 self.row = np.append( self.row,s) 
                                                 self.col = np.append( self.col,self.dic[tuple(auxj)])
                                                 self.data = np.append( self.data, V*sgT)         
       
      
          
                
            
        
      def kinetic(self):
          for s in range(len(self.b)):#
              for i in range(self.L):
                  for j in range(self.L):
                      sg1,aux1 = self.destroy(self.b[s],j)
                      if sg1 != 0:
                          sg2,aux2 = self.create(aux1,i)
                          if sg2 != 0:
                              self.row = np.append( self.row,s ) 
                              self.col = np.append( self.col,self.dic[tuple(aux2)])
                              self.data = np.append( self.data, sg1*sg2*self.t_Hopping(i,j))
                    
      def get_quantum_numbers(self,i):
          if i%2==0:
              s = 0
          else:
              s = 1
          pos = int(i/2)
          
          return pos,s
      
      def V_Coulomb(self,i,j,k,l):
          qn1 = self.get_quantum_numbers(i)
          qn2 = self.get_quantum_numbers(j)
          qn3 = self.get_quantum_numbers(k)
          qn4 = self.get_quantum_numbers(l)
          n1 = qn1[0]
          n2 = qn2[0]
          n3 = qn3[0]
          n4 = qn4[0]
          
          s1 = qn1[1]
          s2 = qn2[1]
          s3 = qn3[1]
          s4 = qn4[1]
          
          
          if s1 != s4:
              return 0;
          if s2 != s3:
              return 0;

          return self.Vee[n1,n2,n3,n4];
      
      def t_Hopping(self,i,j):
          qn1 = self.get_quantum_numbers(i)
          qn2 = self.get_quantum_numbers(j)
          
          n1 = qn1[0]
          n2 = qn2[0]
          
          
          s1 = qn1[1]
          s2 = qn2[1]
          
          
          
          if s1 != s2:
              return 0;
         

          return self.t[n1,n2];
          

                  
      def Diagonalize(self,n_states):
          print("Diagonalization...")
          print("Sparse size :", getsizeof(self.M_sparse))
          return eigsh(self.M_sparse,min(len(self.b)-1,n_states), which='SA')
          