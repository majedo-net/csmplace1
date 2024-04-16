import numpy as np

def argmax(listvals):
    return max(enumerate(listvals), key=lambda x:x[1])[0]

class LevelGraph:
    '''
    Class to represent hypergraph of each level in coarsening/relaxation process
    '''
    cells = list()
    edges = list()
    Ncells = 0

    def __init__(self, cells_,edges_):
        self.cells = cells_
        self.edges = edges_
        self.Ncells = len(self.cells)


    def findNeighbors(self):
        for cidx in range(self.Ncells):
            for edge in self.edges: 
                if cidx in edge:
                    self.cells[cidx].neighbors=self.cells[cidx].neighbors.union(edge)

    def sortDegree(self):
        self.degree = dict(zip(range(self.Ncells), (len(cell.neighbors) for cell in self.cells)))
        self.degree = dict(sorted(self.degree.items(),key=lambda item: item[1],reverse=True))

    def calcAffinity(self):
        for cidx in range(self.Ncells):
            self.cells[cidx].affinities = [0]*len(self.cells[cidx].neighbors)
            for nidx,neighb in enumerate(self.cells[cidx].neighbors):
                for edge in self.edges:
                    if (cidx+1 in edge) and (neighb in edge) and (cidx!=neighb):
                        weight = 1 # TODO: What do the edge weights mean?
                        mage = len(edge)
                        areaE = 0.0
                        for vert in edge:
                            areaE += self.cells[vert].area 
                        self.cells[cidx].affinities[nidx] = weight/((mage-1)*areaE)

    def cluster(self):
        lplus1cells  = [] # cell array for the next level
        mergedverts = [] # cell vertices that have already been merged
        for cidx in self.degree.keys():
            if cidx in mergedverts: continue
            # find cell with max affinity
            nidx = argmax(self.cells[cidx].affinities)
            mergecell = self.cells[cidx].neighbors[nidx]





