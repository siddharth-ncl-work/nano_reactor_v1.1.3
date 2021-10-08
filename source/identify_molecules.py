from source.frame_class import frame_class
import networkx as nx
import networkx.algorithms.isomorphism as iso
from param.gen_params import coll_steps,non_coll_steps,inc
from source.shared_methods import found,addUnique

def findInFrame(molecule,frame):
  insts=[]
  G=frame.frame_graph
  H=molecule.mol_graph
  nm=iso.categorical_node_match('element','C')
  for sg in list(nx.connected_component_subgraphs(G)):
    GM=iso.GraphMatcher(H,sg,node_match=nm)
    if GM.is_isomorphic():
      insts.append(sg)
  
  return found(frame.frame_no,insts)

def identifyMolecules(molecule,file_name,start_frame_no=0,end_frame_no=100000,inc=1):
  found_molecules_list=[]
  curr_list=[]
  file=open(file_name,'r')
  for frame_no in range(start_frame_no,end_frame_no+1,inc):
    frame=frame_class(file,frame_no)
    if frame.frame_no==-1:
      continue
    #print "At frame "+str(frame_no)
    found_molecules=findInFrame(molecule,frame)
    is_added=addUnique(found_molecules_list,found_molecules)
    if is_added:
      print found_molecules    
  return found_molecules_list

def getUniqueFoundObj(found_objs):

  """
  first occurence of non-zero unique found objects
  """
  filtered_found_objs=[]
    
  for f in found_objs:
    if f._no_of_insts!=0:   
     filtered_found_objs.append(f)
     break

  for f in found_objs:
    if f._no_of_insts!=0:
     present=False
     for ff in filtered_found_objs:
       if ff.isEqual(f):
         present=True
         break
     if not present:
       filtered_found_objs.append(f)
 
  return filtered_found_objs

def getUniqueMolecules(found_objs):
  unique_mols=[]

  for f in found_objs:
    if f._no_of_insts!=0:
     for G in f._insts:
       unique_mols.append((f._frame_no,G))
     break

  for f in found_objs:
    for G in f._insts:
      present=False
      for H in unique_mols:
        if set(G.nodes())==set(H[1].nodes()):
          present=True
          break
      if not present:
        unique_mols.append((f._frame_no,G))

  return unique_mols


"""
(13/9/18)
#count number of unique molecules
"""

def countMolecules(file,start_frame_no,end_frame_no):
  no_of_mols={}
  #file=open(file_name,'r')
  for frame_no in range(start_frame_no,end_frame_no+1,inc):
    frame=frame_class(file,frame_no)
    if frame.frame_no==-1:
      continue
    #print "At frame "+str(frame_no)
    unique_mols=countInFrame(frame)
    no_of_mols[frame.frame_no]=len(unique_mols)
    print frame_no
    for i in unique_mols:
      print getNodeName(i)+" "+str(unique_mols[i])
    print ''
  for i in no_of_mols:
    print str(i)+' '+str(no_of_mols[i])


def countInFrame(frame):
  unique_mols={}#myDict()
  G=frame.frame_graph
  nm=iso.categorical_node_match('element','C')
  sg_list=list(nx.connected_component_subgraphs(G))
  for sg in sg_list:
    is_present=False
    for H in unique_mols:
      GM=iso.GraphMatcher(H,sg,node_match=nm)
      if GM.is_isomorphic():
        unique_mols[H]+=1
        is_present=True
        break
    if not is_present:
       unique_mols[sg]=1

  return unique_mols
 

def getNodeName(H):
  #s=frame='+str(frame_no)
  #nodes=''
  #for i in list(H.nodes()):
    #nodes+=str(i)+'_'
  elements=''
  for i in list(H.nodes()):
    elements+=str(H.node[i]['element'])
  return elements


class myDict(dict):
  nm=iso.categorical_node_match('element','C')
  def __contains__(self, key):
    for sg in self:
      GM=iso.GraphMatcher(key,sg,node_match=nm)
      return GM.is_isomorphic()
      #insts.append(sg)  
    

  
