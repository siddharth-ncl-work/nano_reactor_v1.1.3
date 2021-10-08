import numpy as np
import math
import networkx as nx
import shared_methods as shm

class molecule_class:
  atoms=-1
  cords=-1
  mol_adj_mat=-1
  mol_labels=-1
  mol_graph=-1
  mol_name=-1

  def __init__(self,file_name):
    self.mol_name=file_name.split('/')[-1].split('.')[0]
    file=open(file_name,'r')
    self.cords=self.readFile(file)
    mol_label_file=open('output/mol_labels.csv','w')
    self.mol_labels=shm.getGraphLabels(self.cords,mol_label_file)
    self.mol_adj_mat=shm.getAdjMat(self.cords)
    self.mol_graph=shm.getGraph(self.mol_adj_mat,self.mol_labels)
    file=open('output/mol_adjancency_matrix.csv','w')
    shm.writeCsvFile(file,self.cords,self.mol_adj_mat,0)
    print "molecule class '%s' initialized"%file_name 
  
  def readFile(self,file):
    cords=-1
    self.atoms=int(file.readline())
    file.readline()
    cords=[file.readline().split() for i in range(self.atoms)]
    for i in range(self.atoms):
      cords[i][1:]=map(float,cords[i][1:])
    return cords

