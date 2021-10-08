from param.gen_params import file_name,coll_steps,non_coll_steps,start_frame,end_frame,molecule_name,v_no,inc,filter,tracking_start_frame,include_coll_step,include_hmm,output_file_prefix
from frame_class import frame_class
import networkx as nx
from operator import itemgetter
from datetime import datetime
from shared_methods import found
from backward_reading import seekIni,getSeekPos

import networkx as nx
import graphviz
from networkx.drawing.nx_agraph import to_agraph,graphviz_layout

##########################
#___PATH___
PATH=[]
path_graph=nx.DiGraph()
found_paths=dict()
########################

timestamp=-1
outfile_name=''

def findPaths(unique_molecule_list):
    global PATH,timestamp      
    timestamp=datetime.now().strftime("%d_%m_%Y")
    print "Tracking Started!"  
    file=open(file_name,'r')
    seekIni(file,unique_molecule_list[0][0])
  #for found in found_frame:
    frame_no,molecule=unique_molecule_list[0]
    trackAtoms(file,tracking_start_frame,frame_no,molecule)
    #PATH=sorted(PATH,key=itemgetter(0))
    post_track()
    file.close()

def post_track():
  global timestamp,outfile_name
  #timestamp=datetime.now().strftime("%d_%m_%Y")
  outfile_name=output_file_prefix+'path_'+molecule_name+'_'+str(timestamp)+'_with'
  if not include_coll_step:
    outfile_name+='out'
  outfile_name+='_collsion_'+'with'
  if include_coll_step:
    if not include_hmm:
      outfile_name+='out'
    outfile_name+='_hmm_'+'with'
  if not filter:
    outfile_name+='out'
  outfile_name+='_filter_'+v_no
  nx.write_gexf(path_graph,'output/'+outfile_name+'.gexf')
  out_file=open('output/'+outfile_name+'.dat','w')
  for i in PATH:
    out_file.write(str(i._frame_no)+',')
    out_file.write(str(i._no_of_insts)+',')
    for j in i._insts:
      out_file.write(str(j.nodes(data=False)))
    out_file.write('\n')
    out_file.write('---->')
    out_file.write(str(i._found_products._frame_no)+',')
    out_file.write(str(i._found_products._insts[0].nodes(data=False)))
    out_file.write('\n')
  #nx.write_gexf(path_graph,outfile_name+'.gexf')
  #nx.write_gml(path_graph,outfile_name+'.gml')
  #nx.write_gpickle(path_graph,outfile_name'+.gpickle')
  #nx.write_graphml(path_graph,outfile_name+'.graphml')
  #nx.write_yaml(path_graph,outfile_name+'.yaml') #library problem
  #nx.write_pajek(path_graph,outfile_name+'.net')
  createPathImage()

def createPathImage():
  global outfile_name
  G=nx.read_gexf('output/'+outfile_name+'.gexf')
  aG=to_agraph(G)
  U=graphviz.Source(aG,filename='output/'+outfile_name,format='png')
  U.render(cleanup=True)

def trackAtoms(file,start_frame_no,end_frame_no,H):
  global PATH
  global found_paths

  print "\nTracking:("+str(end_frame_no)+','+ str(list(H.nodes()))+")"

  if end_frame_no==start_frame_no:
    return 
  is_found=isPathFound(end_frame_no,H)
  print 'is_found:'+str(is_found)  
  if is_found:
    return

  last_found_fragments=getLastFragments(file,start_frame_no,end_frame_no,H)
  print 'Last Fragments:'+str(last_found_fragments)

  PATH.append(last_found_fragments)
  generatePathGraph(end_frame_no,last_found_fragments,H)
  #post_track()
  if end_frame_no in found_paths.keys():
    found_paths[end_frame_no].append(set(H.nodes()))
  else:
    found_paths[end_frame_no]=[set(H.nodes())]

  for inst in last_found_fragments._insts:
    trackAtoms(file,start_frame_no,last_found_fragments._frame_no,inst)

def getLastFragments(file,start_frame_no,end_frame_no,H):
  last_found_fragments=found(-1,[])
  curr_frame_no=end_frame_no-1
  while curr_frame_no>=start_frame_no:
    file.seek(getSeekPos(curr_frame_no))
    frame=frame_class(file,curr_frame_no)
    if frame.frame_no==-1:
      curr_frame_no-=1
      continue
    found_fragments=findInFrameTracking(frame,H)
    found_fragments._found_products=found(end_frame_no,[H])     ## for path.dat ##
    #filter here
    if not filter or not filter_mols(found_fragments):
       if last_found_fragments._frame_no==-1:
         last_found_fragments=found_fragments
       else:
         if last_found_fragments.isEqual(found_fragments):
           last_found_fragments=found_fragments
         else:
           break
    curr_frame_no-=1 
  print last_found_fragments
  return last_found_fragments

def generatePathGraph(frame_no,found,H):
  global path_graph
  parent_name=getNodeName(H,frame_no)
  path_graph.add_node(parent_name)
  child_nodes=[getNodeName(i,found._frame_no) for i in found._insts]
  path_graph.add_nodes_from(child_nodes)
  for i in child_nodes:
    path_graph.add_edge(parent_name,i)

def getNodeName(H,frame_no):
  s='frame='+str(frame_no)
  nodes=''
  for i in list(H.nodes()):
    nodes+=str(i)+'_'
  elements=''
  for i in list(H.nodes()):
    elements+=str(H.node[i]['element'])+'_'
  return s+'\n'+nodes+'\n'+elements+'\n'+str(list(H.edges()))  

def isPathFound(frame_no,H):
  global found_paths

  if frame_no not in found_paths.keys():
    return False
   
  found_mols=found_paths[frame_no]
  present=False
  for i in found_mols:
    if i==set(H.nodes()):
      present=True  
      break

  return present

def findInFrameTracking(frame,molecule):
  insts=[]
  G=frame.frame_graph
  H=molecule
  H_nset=set(H.nodes())
  for sg in list(nx.connected_component_subgraphs(G)):
    sg_nset=set(sg.nodes())
    if not H_nset.isdisjoint(sg_nset):
      insts.append(sg)
  
  return found(frame.frame_no,insts) 

#filter conditions
def filter_mols(found):
  for inst in found._insts:
    if nx.number_of_nodes(inst)<3:
      print '1 or 2 atom mol found....not adding'
      return True
  return False




   
