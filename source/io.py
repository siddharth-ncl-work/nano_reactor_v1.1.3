import shared_meths as shm
from param.gen_params import coll_steps,non_coll_steps,include_coll_step,backward_reading



class readFile:

  def __init__(self,file=None,frame_no=-1,file_name=None):

    if file!=None:
      self.file=file
      self.file_name=file.name
    elif file_name!=None:
      self.file_name=file_name
      self.file=open(file_name,'r')
    else:
      print "file either file object or file name is not given"
      return

    self.file_type=shm.getFileExtension(self.file_name)
    self.includes_frames=shm.includesFrame(self.file_name)
    if self.file_type=='mol':
      if self.includes_frames:
        self.readFrame_mol(file,frame_no)
      else:
        self.readNoFrame_mol(file)
    elif self.file_type=='xyz':
      if self.includes_frames:
        self.readFrame_xyz(file,frame_no)
      else:
        self.readNoFrame_xyz(file)

  def readFrame_mol(self):
    pass

  def readFrame_xyz(self):
    pass

  def readNoFrame_mol(self):
    pass
  
  def readNoFrame_xyz(self):
    pass
  

