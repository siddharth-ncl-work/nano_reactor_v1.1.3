from param.gen_params import file_name,start_frame,end_frame
from source.identify_molecules import countMolecules
from source.backward_reading import seekIni,getSeekPos

file=open(file_name,'r')
seekIni(file,3000)
print 'dfd'
file.close()
file=open(file_name,'r')
countMolecules(file,0,2000)

