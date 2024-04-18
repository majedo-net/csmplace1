import numpy as np

def argmax(listvals):
    return max(enumerate(listvals), key=lambda x:x[1])[0]

class Vertex:
    '''
    Class to represent a vertex in clustering level graphs
    '''
    neighbors=set()
    affinities=list()
    cell_idxes = list()

    def __init__(self,cell_idxes_,neighbors_=set(),affinities_=[]):
        self.neighbors=neighbors_
        self.affinities_=affinities_
        self.cell_idxes = cell_idxes_

class LevelGraph:
    '''
    Class to represent hypergraph of clustering levels in coarsening/relaxation process
    '''
    edges = list() # list of lists [i][j][k] where i=level index, j=edge index, k=vertex index
    master_cell_array = list()
    level_index_map = list()
    Ncells = 0 # store the number of master cells for convenience
    Nverts = 0 # number of vertices changes each level
    verts = list() # only need the vertice objects for the current level

    def __init__(self,edges_,cell_array_):
        self.edges = [edges_]
        self.master_cell_array = cell_array_
        self.Ncells = len(cell_array_)
        self.Nverts = len(cell_array_)
        self.level_index_map = np.linspace(0,self.Ncells,self.Ncells,dtype=np.int16)
        self.current_level = 0

    def nextLevel(self):
        self.current_level += 1
        if self.current_level>1: 
            self.level_index_map=np.vstack((self.level_index_map,np.zeros_like(self.level_index_map[0,:])))
        else:
            self.level_index_map=np.vstack((self.level_index_map,np.zeros_like(self.level_index_map)))
        for idx in range(self.Nverts):
            self.verts.append(Vertex([idx]))
            for edge in self.edges[-1]:
                if idx in edge:
                    self.verts[-1].neighbors=self.verts[-1].neighbors.union(edge)
            # after using the set type to eliminate duplicates, convert neighbors to a list so we can index it later
            self.verts[-1].neighbors = list(self.verts[-1].neighbors)



    def sortDegree(self):
        self.degree = dict(zip(range(self.Nverts), (len(vert.neighbors) for vert in self.verts)))
        self.degree = dict(sorted(self.degree.items(),key=lambda item: item[1],reverse=True))

    def calcAffinity(self):
        for cidx in range(self.Nverts):
            self.verts[cidx].affinities = [0]*len(self.verts[cidx].neighbors)
            for nidx,neighb in enumerate(self.verts[cidx].neighbors):
                for edge in self.edges[-1]:
                    if (cidx in edge) and (neighb in edge) and (cidx!=neighb):
                        weight = 1 # TODO: What do the edge weights mean?
                        mage = len(edge)
                        areaE = 0.0
                        for vert in edge:
                            for cell in self.verts[vert].cell_idxes:
                                areaE += self.master_cell_array[cell].area 
                        self.verts[cidx].affinities[nidx] += weight/((mage-1)*areaE)

    def cluster(self):
        mergedverts = [] # vertices that have already been merged
        nextcluster=0 # running index of clusters
        for cidx in self.degree.keys():
            if cidx in mergedverts: continue
            # find cell with max affinity
            nidx = argmax(self.verts[cidx].affinities)
            mergevertidx = self.verts[cidx].neighbors[nidx]
            if mergevertidx in mergedverts: 
                self.level_index_map[self.current_level,cidx] = self.level_index_map[self.current_level,mergevertidx]
                mergedverts.append(cidx)

            else:
                self.level_index_map[self.current_level,cidx] = nextcluster
                self.level_index_map[self.current_level,mergevertidx] = nextcluster
                nextcluster += 1
                mergedverts.append(cidx)
                mergedverts.append(mergevertidx)
        self.Nverts = np.max(self.level_index_map[-1,:])+1






