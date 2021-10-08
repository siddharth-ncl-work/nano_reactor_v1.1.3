import numpy as np
import math
from param.bond_lengths import bl
import networkx as nx
import networkx.algorithms.isomorphism as iso
from param.gen_params import coll_steps,non_coll_steps

seek_dict={}

def getAdjMat(cords):
  atoms=len(cords)
  adj_mat=np.zeros((atoms,atoms))
  for i in range(atoms):
    for j in range(i+1,atoms):
      d=getDist(cords[i][1:],cords[j][1:])
      if d<=bl[cords[i][0]][cords[j][0]]:
        adj_mat[i][j]=adj_mat[j][i]=1
        #print '%s %s %d %d %f'%(cords[i][0],cords[j][0],i,j,d)
  return adj_mat

def getGraph(adj_mat,labels):
  G=nx.from_numpy_matrix(adj_mat)
  nx.set_node_attributes(G,name='element',values=labels)
  return G

def getGraphLabels(cords,label_file):
  atoms=len(cords)
  for i in range(atoms):
    label_file.write(str(i)+' '+cords[i][0]+'\n')
  return {i:cords[i][0] for i in range(atoms)}

def getDist(v1,v2):
  d=0.0
  for i in range(len(v1)):
    d+=(v2[i]-v1[i])**2
  return math.sqrt(d)

def writeCsvFile(wfile,cords,adj_mat,frame_no):
  atoms=len(cords)
  wfile.write('elements,')
  for i in range(atoms):
    if i!=atoms-1:
      wfile.write(str(i)+',')
    else:
      wfile.write(str(i)+" "+str(frame_no)+'\n')
  for row in range(atoms):
    wfile.write(str(row)+','+cords[row][0]+',')
    for col in range(atoms):
      if col!=atoms-1:
        wfile.write(str(adj_mat[row][col])+',')
      else:
        wfile.write(str(adj_mat[row][col])+'\n')

def addUnique(final_list,found):
  if len(final_list)==0:
    final_list.append(found)
    return True

  if not final_list[-1].isEqual(found):
    final_list.append(found)
    return True
  else:
   return False

def getFileExtension(file_name):
  return file_name.strip().split('/')[-1].split('.')[-1]

def includesFrame(file_name):
  file=open(file_name,'r')
  line1=file.readLine().strip()[0].lower()
  line2=file.readLine().strip()[0].lower()
  return 'frame'==line1 or 'frame'==line2

def isCollStep(frame_no):
  return frame_no%(coll_steps+non_coll_steps)>=non_coll_steps 

class found:
  _frame_no=-1
  _no_of_insts=-1
  _insts=-1
  _found_products=-1

  def __init__(self,frame_no,insts):
    self._frame_no=frame_no
    self._no_of_insts=len(insts)
    self._insts=insts

  def __str__(self):
    s=str(self._frame_no)+" "+str(self._no_of_insts)+" "
    for G in self._insts:
      s+=str(G.nodes())+" "    
    return s

  def isEqual(self,f):
    if self._no_of_insts!=f._no_of_insts: 
      return False
   
    if self._no_of_insts==0 and f._no_of_insts==0:
      return True
  
    nm=iso.categorical_node_match('element','C')
    present=False
    for G in self._insts:
      present=False
      for H in f._insts:
        if set(G.nodes())==set(H.nodes()):
          GM=iso.GraphMatcher(H,G,node_match=nm)
          present=GM.is_isomorphic()
          break
      if not present:
        return present

    return present


