from source import hmm
from source.frame_class import frame_class
from param.gen_params import file_name,coll_steps,non_coll_steps
from source.backward_reading import seekIni,getSeekPos
import numpy as np

X=[0,1,0,1,0,1,0,1,0]
X = np.atleast_2d(X).T
hmm.test_hmm1(X)

X=[0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
X = np.atleast_2d(X).T
hmm.test_hmm1(X)

file=open(file_name,'r')
seekIni(file,78860)

def _count_nonzero(a):
  c=0
  for i in a:
    for j in i:
      if j!=0:
        c+=1
  return c

def acc(a,b):
  acc=0
  for i in range(len(a)):
    for j in range(i+1,len(a[0])):
      if a[i][j]==b[i][j]:
        acc+=1
  total_ubonds=(len(a)*(len(a)-1)/2)
  #print str(_count_nonzero(b))+" "+str(_count_nonzero(a))+" "+str(total_ubonds-acc)
  return (total_ubonds-acc,acc*100/total_ubonds)

def getStats(l):
  if len(l)==0:return
  return [min(l),max(l),round(np.mean(l),2),round(np.std(l),2)]

def getExtraStats(l):
  if len(l)==0:return
  a=len(filter(lambda x:x>=0,l))
  b=len(filter(lambda x:x>0,l))
  return (a*100/len(l),b*100/len(l))

def analyse(start_frame_no,end_frame_no):
  ubonds_diff=[]
  similarity=[]
  for frame_no in range(start_frame_no,end_frame_no+1):
    print frame_no
    a=hmm.getAdjMat(frame_no)
    file.seek(getSeekPos(frame_no))
    frame=frame_class(file,frame_no)
    z,t=acc(a,frame._orig_frame_adj_mat)
    ubonds_diff.append(z)
    similarity.append(t)
  print [getStats(ubonds_diff),getStats(similarity)]


  
def analyseCollStep(start_frame_no,end_frame_no):
  bonds_diff=[]
  hmm_bonds_diff=[]
  
  
                                                                   
  file.seek(getSeekPos(start_frame_no-1))
  ref_frame=frame_class(file,start_frame_no-1)
  #ref_no_of_bonds=np.count_nonzero(ref_frame._orig_frame_adj_mat)
  ref_no_of_edges=ref_frame.frame_graph.number_of_edges()
  for frame_no in range(start_frame_no,end_frame_no+1):
    file.seek(getSeekPos(frame_no))
    frame=frame_class(file,frame_no)
    #no_of_bonds=np.count_nonzero(frame._orig_frame_adj_mat)
    no_of_edges=frame._orig_frame_graph.number_of_edges()
    #bonds_diff.append(no_of_bonds-ref_no_of_bonds)
    bonds_diff.append(no_of_edges-ref_no_of_edges)
    #print str(frame_no)+" "+str(ref_no_of_bonds)+" "+str(no_of_bonds)+" "+str(no_of_bonds-ref_no_of_bonds)
    print ref_no_of_edges,no_of_edges
    #a=hmm.getAdjMat(frame_no)
    #hmm_no_of_bonds=np.count_nonzero(frame.frame_adj_mat)
    hmm_no_of_edges=frame.frame_graph.number_of_edges()
    #hmm_bonds_diff.append(no_of_bonds-hmm_no_of_bonds)
    hmm_bonds_diff.append(no_of_edges-hmm_no_of_edges)
    #print str(frame_no)+" "+str(no_of_bonds)+" "+str(hmm_no_of_bonds)+" "+str(no_of_bonds-hmm_no_of_bonds)
    print no_of_edges,hmm_no_of_edges
  bond_diff_stats=getStats(bonds_diff)
  bond_diff_extra_stats=getExtraStats(bonds_diff)
  print [bond_diff_stats,bond_diff_extra_stats,getStats(hmm_bonds_diff),getExtraStats(hmm_bonds_diff)]
  wl=filter(lambda x:x[0]>0 and x[1]>0,zip(bonds_diff,hmm_bonds_diff))
  e=len(wl)
  f=len(filter(lambda x:x>0,bonds_diff))
  g=e*100/f
  print e,f,g

  m=max(bonds_diff)
  power_stats=getStats(map(lambda x:x[1]*100/m,wl))
  print power_stats
  if power_stats==None:
    power_stats=[0,0,0,0]
  elif power_stats[2]>100:
    power_stats[2]=100
  stats=(bond_diff_extra_stats[1],bond_diff_stats[1],g,power_stats[2])
  return stats


def cycleAnalysis(cycle):
  global coll_steps,non_coll_steps
  performance=[]
  cycle_frames=coll_steps+non_coll_steps
  start_frame_no=-1
  end_frame_no=-1
  
  """
  #non-collsion
  start_frame_no=cycle_frames*cycle
  end_frame_no=start_frame_no+non_coll_steps-1
  #print (start_frame_no,end_frame_no) 
  #analyse(start_frame_no,end_frame_no)
  """
  
  #collision  
  start_frame_no=cycle_frames*cycle+non_coll_steps
  end_frame_no=cycle_frames*(cycle+1)-1#start_frame_no+coll_steps-1
  print (start_frame_no,end_frame_no)
  return analyseCollStep(start_frame_no,end_frame_no)


def writeCsvFile(wfile,p,cycle):
  wfile.write(str(cycle)+",")
  for i in range(len(p)-1):
    wfile.write(str(p[i])+",")
  wfile.write(str(p[-1])+"\n")



max_cycle=int(78860/coll_steps+non_coll_steps)
performance=[[],[],[],[]]
wfile=open("analysis.dat","w")

for cycle in range(0,max_cycle-1):
  #cycle_frames=coll_steps+non_coll_steps
  #cycle=int(i/cycle_frames)
  print cycle
  p=cycleAnalysis(cycle)
  print p

  for i in range(len(p)):
    performance[i].append(p[i])

  print "\n"
  for i in range(len(p)):
    print getStats(performance[i])
  print "==================================================="
  writeCsvFile(wfile,p,cycle)





"""
file=open(file_name,'r')
seekIni(file,15000)
a=hmm.getAdjMat(986)
file.seek(getSeekPos(986))
frame=frame_class(file,986)
print np.array_equal(a,frame.frame_adj_mat)
"""

#=========== Test Examples ===============#

# 1) toy example
"""
Observed Signal =  GGCACTGAA = 221013200
ACGT=0123
Expected Output = HHHLLLLLL = 000111111
HL = 01

Hidden States = 2,4
Start Prob = [0.5,0.5]
Transition Prob = [[0.5, 0.5],
                   [0.4, 0.6]]
Emission Prob = [[0.2,0.3,0.3,0.2],
                 [0.3,0.2,0.2,0.3]]

"""
#hmm.toy_hmm()
#------------------------------------

# 2) fever example
"""
Observed Signal = Normal,Cold,Dizzy = 012
Normal,Cold,Dizzy = 012
Expected Output = HHF = 001
HF=01

Hidden States = 2,3
Start Prob = [0.6,0.4]
Transition Prob = [[0.7,0.3],
                   [0.4,0.6]]
Emission Prob = [[0.5,0.4,0.1],
                 [0.1.0.3.0.6]]
"""
#hmm.fever_hmm()
#----------------------------------------


# 3) mini_example
# NOT WORKING
"""
Observed Signal = SxyyE = 01223
SxyE=0123
Expected Output = SABBE = 01223 
S AB E = 0 12 3

Hidden and Observed States = 2,2 
Start Prob = [1,0,0,0]
Transition Prob = [[0,0.7,0.3,0],
                   [0,0.2,0.7,0.1],
                   [0,0.7,0.2,0.1],
                   [0,0,0,0]]
Emmision Prob = [[1,0,0,0],
                 [0,0.4,0.6,0],
                 [0,0.3,0.7,0],
                 [0,0,0,0]]                         
"""
# NOT WORKING
#hmm.mini_ex_hmm()
#------------------------------------------


# 4) lazy programmer example
#for i in range(50):
 # hmm.lazy_prog_hmmd(i)

#--------------------------------------------





