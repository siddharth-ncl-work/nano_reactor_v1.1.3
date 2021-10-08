##_____main_file_____##

from source.molecule_class import molecule_class #(was import loop)
from source.tracking import findPaths
from param.gen_params import file_name,coll_steps,non_coll_steps,start_frame,end_frame,molecule_name,inc,_tracking
from source.identify_molecules import identifyMolecules,getUniqueFoundObj,getUniqueMolecules
from source.frame_class import frame_class

molecule=molecule_class("input/molecules/"+molecule_name+".xyz") 
found_molecules=identifyMolecules(molecule,file_name,start_frame_no=start_frame,end_frame_no=end_frame,inc=inc)

for i in found_molecules:
  print i

ufo=getUniqueFoundObj(found_molecules)
print ""
for i in ufo:
  print i

unique_molecules_list=getUniqueMolecules(found_molecules)
print ""
for i,G in unique_molecules_list:
  print str(i)+" "+str(G.nodes())
print ""

if _tracking:
  paths=findPaths(unique_molecules_list)





