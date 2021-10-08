import networkx as nx
import numpy as np
import math 
import shared_methods as shm
from param.gen_params import coll_steps,non_coll_steps,include_coll_step,include_hmm
import hmm

#frame_labels={}

class frame_class:
  frame_no=-1
  atoms=-1
  cords=-1
  _orig_frame_adj_mat=-1
  frame_adj_mat=-1
  _frame_labels=-1
  frame_graph=-1

  def __init__(self,file,frame_no=-1):
    global frame_labels
    is_found=-1
    file_type=shm.getFileExtension(file.name)
    if file_type=='mol':
      is_found=self.readFile_mol(file,frame_no)
    elif file_type=='xyz':
      is_found=self.readFile_xyz(file,frame_no)

    if is_found==1:
      pass
      #print 'Successfully Found Frame '+str(frame_no)
    else:
      print 'Error:Frame %d is not present '%frame_no
      return
  
    self._orig_frame_adj_mat=self.frame_adj_mat 
    self.frame_no=frame_no
    if shm.isCollStep(frame_no):
      if include_coll_step:
        if include_hmm:
          self.frame_adj_mat=getCollAdjMat(frame_no)
      else:
        self.frame_no=-1
        return
    frame_label_file=open('output/frame_labels.csv','w') 
    self._frame_labels=shm.getGraphLabels(self.cords,frame_label_file)
    self.frame_graph=shm.getGraph(self.frame_adj_mat,self._frame_labels)
    self._orig_frame_graph=shm.getGraph(self._orig_frame_adj_mat,self._frame_labels)
    wfile=open('output/frame_adjancency_matrix.csv','w')
    shm.writeCsvFile(wfile,self.cords,self.frame_adj_mat,self.frame_no)
    #print 'Frame %d class initialized\r'%self.frame_no,

  def readFile_mol(self,file,frame_no):
    is_found=-1
    f=file.readline()
    while f!='':
      f=int(f.split()[1])
      file.readline()
      file.readline()
      self.atoms=int(file.readline().split()[0])
      if f==frame_no:
        self.frame_adj_mat=np.zeros((self.atoms,self.atoms))
        self.cords=[file.readline().split()[:4] for i in range(self.atoms)]
        self.cords=map(swap,self.cords)
        for i in range(self.atoms):
          self.cords[i][1:]=map(float,self.cords[i][1:])
  
        edge=file.readline()
        while edge.split()[0]!='M':
          n1=int(edge[:3])-1
          n2=int(edge[3:6])-1
          #print str(n1)+" "+str(n2)
          self.frame_adj_mat[n1][n2]=self.frame_adj_mat[n2][n1]=1
          edge=file.readline()
          
        end='ssss'
        while end!='' and end[0]!='$':
          end=file.readline().strip()              
        is_found=1
        break
      else:
        end='ssss'
        while end!='' and end[0]!='$':
          end=file.readline().strip()
      f=file.readline()
    return is_found


  def readFile_xyz(self,file,frame_no):
    is_found=-1
    _atoms=self.atoms=int(file.readline())
    while _atoms!='':
      self.atoms=int(_atoms)
      f=int(file.readline().split()[1])
      if f==frame_no:
        self.cords=[file.readline().split() for i in range(self.atoms)]
        for i in range(self.atoms):
          self.cords[i][1:]=map(float,self.cords[i][1:])
        is_found=1
        #Adj_mat
        self.frame_adj_mat=shm.getAdjMat(self.cords)       

        break
      else:
        for i in range(self.atoms):file.readline()
      _atoms=file.readline()
    return is_found

def swap(list):
  temp=list[-1]
  list[-1]=list[0]
  list[0]=temp
  return list

def getCollAdjMat(frame_no):
  return hmm.getAdjMat(frame_no)

"""
def initFrameLabels():
  global frame_labels
  cords
  frame_label_file=open('output/frame_labels.csv','w')
  frame_labels=shm.getGraphLabels(self.cords,frame_label_file)
"""
# NEED TO SEPARATE READING FROM FRAME CLASS
