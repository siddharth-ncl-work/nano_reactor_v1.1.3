import networkx as nx
import networkx.algorithms.isomorphism as iso

def fbEquality(fpath,bpath):
  fgraph=nx.read_gexf(fpath)
  bgraph=nx.read_gexf(bpath)
 
  fnodes=set(fgraph.nodes())
  bnodes=set(bgraph.nodes())

  fedges=set(fgraph.edges())
  bedges=set(bgraph.edges())
  
  GM=iso.GraphMatcher(fgraph,bgraph)

  print "nodes:"
  print type(next(iter(fnodes)))
  print str(len(fnodes))+","+str((len(bnodes)))
  #print fnodes
  #print bnodes
  print "edges:"
  print str(len(fedges))+","+str((len(bedges)))
  #print fedges
  #print bedges

  print fnodes==bnodes
  print fedges==bedges
  print GM.is_isomorphic()

print "nhchoh"
fbEquality("/home/vanka/siddharth/tamal/nano_reactor_v1.0.6f/output/TEST_after_bugfix_path_nhchoh_02_06_2018_without_filter_v1.0.6f.gexf","output/TEST_path_nhchoh_02_06_2018_without_filter_v1.0.6b.gexf")

print "\nnh2cho"
fbEquality("/home/vanka/siddharth/tamal/nano_reactor_v1.0.6f/output/TEST_after_bugfix_path_nh2cho_02_06_2018_without_filter_v1.0.6f.gexf","output/TEST_path_nh2cho_02_06_2018_without_filter_v1.0.6b.gexf")

print "\nnhchoh,nh2cho"
fbEquality("/home/vanka/siddharth/tamal/nano_reactor_v1.0.6f/output/TEST_after_bugfix_path_nhchoh_02_06_2018_without_filter_v1.0.6f.gexf","output/TEST_path_nh2cho_02_06_2018_without_filter_v1.0.6b.gexf")

