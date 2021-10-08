import numpy as np
from hmmlearn import hmm
from param.gen_params import file_name,coll_steps,non_coll_steps
import shared_methods as shm
from backward_reading import seekIni,getSeekPos

adj_mat_dict={}
atoms=-1
symm=False



def getAdjMat(frame_no=1000,type="determine",start_frame_no=0,end_frame_no=1249):
  global adj_mat_dict

  if type=="custom" and frame_no not in range(start_frame_no,end_frame_no+1):
    print "wrong frame range"
    return
  if frame_no in adj_mat_dict.keys():
    #print "Already found"
    return adj_mat_dict[frame_no]
  
  #adj_mat_dict.clear() #keep all collision adj mat
  if type!="custom":
    start_frame_no,end_frame_no=getTopBottomFrameNos(frame_no)   
    if frame_no not in range(start_frame_no,end_frame_no+1):  
      print "wrong frame range from getTopBottomFrameNos"
      return
   
  print (start_frame_no,end_frame_no)
  createAllAdjMat(start_frame_no,end_frame_no)
  return adj_mat_dict[frame_no]




def getTopBottomFrameNos(frame_no):
  global coll_steps,non_coll_steps
  cycle_frames=coll_steps+non_coll_steps
  cycle=int(frame_no/cycle_frames)
  start_frame_no=-1
  end_frame_no=-1

  start_frame_no=cycle_frames*cycle
  #end_frame_no=cycle_frames*(cycle+1)-1#start_frame_no+coll_steps-1
  end_frame_no=cycle_frames*(cycle+1)+non_coll_steps-1
  """
  if shm.isCollStep(frame_no):
    start_frame_no=cycle_frames*cycle+non_coll_steps
    end_frame_no=cycle_frames*(cycle+1)-1#start_frame_no+coll_steps-1
  else:
    start_frame_no=cycle_frames*cycle
    end_frame_no=start_frame_no+non_coll_steps-1
  """
  return(start_frame_no,end_frame_no)
  



def createAllAdjMat(start_frame_no,end_frame_no):
  global adj_mat_dict
  ids=readData(file_name,start_frame_no,end_frame_no)
  pds=predictSeries(ids,start_frame_no,end_frame_no)
  for frame_no in range(start_frame_no,end_frame_no+1):
    adj_mat_dict[frame_no]=createAdjMat(pds,frame_no)


def createAdjMat(pds,frame_no):
  frame_adj_mat=np.zeros((atoms,atoms),dtype=np.int8)
  for i in range(atoms):
     for j in range(i,atoms):
       frame_adj_mat[i][j]=frame_adj_mat[j][i]=pds.get(i,j,frame_no)
  return frame_adj_mat



def readData(file_name,start_frame_no,end_frame_no):
  global atoms,symm
  file=open(file_name,'r')
  #seekIni(file,end_frame_no) # remove afterwards
  file.readline()  
  file.readline()
  file.readline()
  atoms=int(file.readline().split()[0])
  ds=dataStorage(atoms,start_frame_no,end_frame_no,symm)
  file_type=shm.getFileExtension(file.name)
  if file_type=="mol":
    for frame_no in range(start_frame_no,end_frame_no+1):
       frame_adj_mat=np.zeros((atoms,atoms))
       file.seek(getSeekPos(frame_no))
       if readFile_mol(file,frame_no,frame_adj_mat):
         for i in range(atoms):
           for j in range(i,atoms):
             ds.put(i,j,frame_no,frame_adj_mat[i][j])    
       else:
         print "cannot read file"
  else:
    pass #for .xyz file
  file.close()
  return ds



def readFile_mol(file,frame_no,frame_adj_mat):
  global atoms
  is_found=False
  f=file.readline()
  while f!='':
    f=int(f.split()[1])
    file.readline()
    file.readline()
    atoms=int(file.readline().split()[0])
    if f==frame_no:
      [file.readline().split()[:4] for i in range(atoms)]

      edge=file.readline()
      while edge.split()[0]!='M':
        n1=int(edge[:3])-1
        n2=int(edge[3:6])-1
        frame_adj_mat[n1][n2]=frame_adj_mat[n2][n1]=1
        edge=file.readline()

      end='ssss'
      while end!='' and end[0]!='$':
        end=file.readline().strip()
      is_found=True
      break
    else:
      end='ssss'
      while end[0]!='$':
        end=file.readline().strip()
    f=file.readline()
  return is_found





def predictSeries(ids,start_frame_no,end_frame_no):
  pds=dataStorage(atoms,start_frame_no,end_frame_no,symm)
  for i in range(atoms):
    for j in range(i+1,atoms):
      _,series=_hmm1(ids.getSeries(i,j))
      series=np.reshape(series,(end_frame_no-start_frame_no+1,1))#atleast_2d(series).T----doing it
      pds.putSeries(i,j,series)
  return pds





def _hmm1(X):
  model = hmm.MultinomialHMM(n_components=2)
  model.startprob_ = np.array([1-X[0][0], X[0][0]])
  model.transmat_ = np.array([[0.999, 0.001],
                              [0.001, 0.999]])
  model.emissionprob_ = np.array([[0.6, 0.4],
                                  [0.4, 0.6]])

  # Predict the optimal sequence of internal hidden state
  #X=map(lambda x:[x[0]],X)
  #print len(X)
  predict = model.decode(X,algorithm="viterbi")
  return predict
  

 

class dataStorage:
  _data=[]
  _atoms=-1
  _no_of_frames=-1
  _symm=False
  _start_frame_no=-1

  def  __init__(self,atoms,start_frame_no,end_frame_no,symm=False):
    self._atoms=atoms
    self._no_of_frames=end_frame_no-start_frame_no+1
    self.createStorage()
    self._symm=symm
    self._start_frame_no=start_frame_no

  def createStorage(self):
    if self._symm:
      self._data=np.zeros((self._atoms*(self._atoms+1)/2,self._no_of_frames,1),dtype=np.int8)
    else:
      self._data=np.zeros((self._atoms,self._atoms,self._no_of_frames,1),dtype=np.int8)

  def put(self,i,j,k,data):
    k=k-self._start_frame_no
    if self._symm:
      self._data[i*(i+1)/2+j][k][0]=data
    else:
      self._data[j][i][k][0]=self._data[i][j][k][0]=data

  def get(self,i,j,k):
    k=k-self._start_frame_no
    if self._symm:
      return self._data[i*(i+1)/2+j][k][0]
    else:
      return self._data[i][j][k][0]


  def getSeries(self,i,j):
    if self._symm:
      return self._data[i*(i+1)/2+j]
    else:
      return self._data[i][j]

  def putSeries(self,i,j,series):
    if self._symm:
      self._data[i*(i+1)/2+j]=series
    else:
      self._data[i][j]=self._data[j][i]=series





"""
def _predictSeries(ids,start_frame_no,end_frame_no):
  pds=dataStorage(atoms,start_frame_no,end_frame_no,symm)
  p,series=_hmm1(ids.getSeries(10,5))
  print (p,series)
  p,series=_hmm1(ids.getSeries(12,40))
  print (p,series)
  p,series=_hmm1(ids.getSeries(40,52))
  print (p,series)
  p,series=_hmm1(ids.getSeries(20,25))
  print (p,series)
  p,series=_hmm1(ids.getSeries(4,8))
  print (p,series)
"""





#================ Test Examples ====================#

#toy example
#observed states should start from 0
def toy_hmm():
  model = hmm.MultinomialHMM(n_components=2)
  model.startprob_ = np.array([0.5, 0.5])
  model.transmat_ = np.array([[0.5, 0.5],
                              [0.4, 0.6]])
  model.emissionprob_ = np.array([[0.2,0.3,0.3,0.2],
                                  [0.3,0.2,0.2,0.3]])

  # Predict the optimal sequence of internal hidden state
  X = np.atleast_2d([2,2,1,0,1,3,2,0,0]).T
  #X=np.random.rand(10,4)
  #print X
  print len(X)
  print X.shape
  print model.decode(X,algorithm="viterbi")
  #print model.sample(10)


#fever example
#observed states should start from 0
def fever_hmm():
  model = hmm.MultinomialHMM(n_components=2)
  model.startprob_ = np.array([0.6, 0.4])
  model.transmat_ = np.array([[0.7, 0.3],
                              [0.4, 0.6]])
  model.emissionprob_ = np.array([[0.5,0.4,0.1],
                                  [0.1,0.3,0.6]])

  # Predict the optimal sequence of internal hidden state
  X = np.atleast_2d([0,1,2]).T
  #X=np.random.rand(10,4)
  #print X
  print len(X)
  print X.shape
  print model.decode(X,algorithm="viterbi")
  #print model.sample(10)


#mini example (NOT WORKING)
#observed states should start from 0
def mini_ex_hmm():
  model = hmm.MultinomialHMM(n_components=4)
  model.startprob_ = np.array([1,0,0,0])
  model.transmat_ = np.array([[0,0.7,0.3,0],
                   	      [0,0.2,0.7,0.1],
                   	      [0,0.7,0.2,0.1],
                   	      [0,0,0,1]])
  model.emissionprob_ = np.array([[1,0,0,0],
                 		  [0,0.4,0.6,0],
                 		  [0,0.3,0.7,0],
                 		  [0,0,0,1]])

  # Predict the optimal sequence of internal hidden state
  X = np.atleast_2d([0,1,2,2,3]).T
  #X=np.random.rand(10,4)
  #print X
  print len(X)
  print X.shape
  print model.decode(X,algorithm="viterbi")
  #print model.sample(10)


#lazy_programmer
def lazy_prog_hmmd(i):
  model = hmm.MultinomialHMM(n_components=2)
  model.startprob_ = np.array([0.5, 0.5])
  model.transmat_ = np.array([[0.1, 0.9], [0.8, 0.2]])
  model.emissionprob_ = np.array([[0.6, 0.4], [0.3, 0.7]])

  X = []
  for line in open('/home/vanka/siddharth/python/HMM/coin_data.txt'):
     # 1 for H, 0 for T
     x = [1 if e == 'H' else 0 for e in line.rstrip()]
     X.append(x)

  
  # Predict the optimal sequence of internal hidden state
  Xt = np.atleast_2d(X[i]).T
  #X=np.random.rand(10,4)
  #print X
  print len(Xt)
  print Xt.shape
  print model.decode(Xt,algorithm="viterbi")
  #print model.sample(10)



def test_hmm1(X):
  model = hmm.MultinomialHMM(n_components=2)
  model.startprob_ = np.array([1-X[0][0], X[0][0]])
  model.transmat_ = np.array([[0.999, 0.001],
                              [0.001, 0.999]])
  model.emissionprob_ = np.array([[0.6, 0.4],
                                  [0.4, 0.6]])

  # Predict the optimal sequence of internal hidden state
  #X = np.atleast_2d(X[i]).T
  #print len(X)
  predict = model.decode(X,algorithm="viterbi")
  print predict
  return predict



