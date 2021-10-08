seek_dict={}

def seekReadFile_mol(file,frame_no):
    seek_pos=file.tell()
    f=file.readline()
    while f!='':
      f=int(f.split()[1])
      file.readline()
      file.readline()
      atoms=int(file.readline().split()[0])
      if f==frame_no:
        end='ssss'
        while end!='' and end[0]!='$':
          end=file.readline().strip()
        return seek_pos
      else:
        end='ssss'
        while end!='' and  end[0]!='$':
          end=file.readline().strip()
      seek_pos=file.tell()
      f=file.readline()
    return -1

def seekIni(file,end_frame):
  global seek_dict
  for frame_no in range(end_frame+1):
    seek_dict[frame_no]=seekReadFile_mol(file,frame_no)

def getSeekPos(frame_no):
  return seek_dict[frame_no]

