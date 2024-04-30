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
        self.x = None
        self.y = None

    def getVertexPosition(self):
        '''
        The position of each vertex is the average of the constituent cell positions
        '''
        x = 0.0
        y = 0.0
        for cidx in self.cell_idxes: # TODO need a global cell array
            cell = master_cell_array[cidx]
            x += cell.cx
            y += cell.cy

        x = x / len(self.cell_idxes)
        y = y / len(self.cell_idxes)
        return x , y

    def getVertexArea(self):
        a = 0.0
        for cidx in self.cell_idxes: # TODO need a global cell array
            cell = master_cell_array[cidx]
            a += cell.area
        return a

    def setVertexPosition(self,x,y):
        '''
        The vertex position is applied to the cells it contains
        '''
        self.x = x
        self.y = y
        for cidx in self.cell_idxes:# TODO need a global cell array
            cell = master_cell_array[cidx]
            cell.cx = x
            cell.cy = y
        
class LevelGraph:
    '''
    Class to represent hypergraph of clustering levels in coarsening/relaxation process
    '''
    edges = list() # list of lists [i][j][k] where i=level index, j=edge index, k=vertex index
    master_cell_array = list()
    level_index_map = list()
    Ncells = 0 # store the number of master cells for convenience
    Nverts = 0 # number of vertices changes each level
    verts = list() # only need the vertex objects for the current level

    def __init__(self,edges_,cell_array_):
        self.edges = [edges_]
        self.master_cell_array = cell_array_
        self.Ncells = len(cell_array_)
        self.Nverts = len(cell_array_)
        self.level_index_map = np.linspace(0,self.Ncells,self.Ncells,dtype=np.int16)
        self.current_level = 0
        self.average_cluster_area = 0.0
        for idx in range(self.Ncells):
            self.average_cluster_area += self.master_cell_array[idx].area
        self.average_cluster_area = self.average_cluster_area / self.Ncells

    def xvec(self):
        '''
        helper function to generate a numpy vector of x values for vertices in the 
        current level
        '''
        vx = []
        for vert in self.verts:
            vx.append(vert.x)
        return np.asarray(vx)

    def yvec(self):
        '''
        helper function to generate a vector of y values for vertices in the 
        current level
        '''
        vy = []
        for vert in self.verts:
            vy.append(vert.y)
        return np.asarray(vy)

    def vvec(self):
        '''
        helper function to generate a vector of concatenated x and y values 
        for vertices in the current level
        '''
        return np.concatenate((self.xvec(),self.yvec()))

    def avec(self):
        '''
        helper function to generate of vector of cell/vertex areas
        '''
        for vert in self.verts:


    def nextLevel(self):
        '''
        Initialize the next level of coarsening with vertex objects and their neighbors
        '''
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

    def getClusterArea(self,idx):
        '''
        calculate the area of the cells in a cluster with index=idx
        '''
        area = 0.0
        for cell in range(self.Ncells):
            if self.level_index_map[self.current_level,cell] == idx:
                area += self.master_cell_array[cell].area

        return area

    def getClusterCenter(self,idx):
        '''
        compute the center of the cluster with index=idx
        '''
        xc = 0.0
        yc = 0.0
        nc = 0 
        for cell in range(self.Ncells):
            if self.level_index_map[self.current_level,cell] == idx:
                nc += 1
                xc += self.master_cell_array[cell].tlx
                yc += self.master_cell_array[cell].tly
        xc = xc / nc
        yc = yc / nc
        return xc, yc

    def doInitialPlace(self):
        '''
        topological sort the current coarsening level and arrange in grid
        '''
        grid_nw = int(np.sqrt(self.Nverts))#Make grid of bins square by default
        grid_nh = grid_nw
        self.sortDegree()
        row=0
        column = 0
        x = 0.0
        y = 0.0
        for cidx in self.degree.keys():
            self.verts[cidx].setVertexPosition(x,y)
            t_area = self.getClusterArea(cidx)
            t_wh = np.sqrt(t_area)
            if column >= grid_nw:
                row += 1
                column = 0
                y += t_wh
                x = 0.0
            else:
                column += 1
                x += t_wh
            
    def sortDegree(self):
        '''
        sort the vertices in the current coarsening level by degree of connectivity. 
        Clustering starts with the highest degree vertex
        '''
        self.degree = dict(zip(range(self.Nverts), (len(vert.neighbors) for vert in self.verts)))
        self.degree = dict(sorted(self.degree.items(),key=lambda item: item[1],reverse=True))

    def calcAffinity(self):
        '''
        Calculates the affinity matrix between vertices in the graph of the current coarsening level
        '''
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
        '''
        Cluster the vertices by affinity
        '''
        mergedverts = [] # vertices that have already been merged
        nextcluster=1 # running index of clusters
        self.edges.append([]) # initialize next level edges
        for cidx in self.degree.keys():
            if cidx in mergedverts: continue
            while cidx not in mergedverts:
                # check if we have run out of neighbors
                if len(self.verts[cidx].affinities) == 0:
                    self.level_index_map[self.current_level,cidx] = nextcluster
                    nextcluster += 1
                    mergedverts.append(cidx)
                    continue
                # find cell with max affinity
                nidx = argmax(self.verts[cidx].affinities)
                mergevertidx = self.verts[cidx].neighbors[nidx]
                # check if the maximum affinity neighbor has already been assigned a new cluster
                if mergevertidx in mergedverts: 
                    self.level_index_map[self.current_level,cidx] = self.level_index_map[self.current_level,mergevertidx]
            
                    # check if the cluster is too big
                    cluster_idx = self.level_index_map[self.current_level,cidx]
                    if self.getClusterArea(cluster_idx) > 1.5*self.average_cluster_area:
                        self.verts[cidx].affinities.pop(nidx)
                        self.verts[cidx].neighbors.pop(nidx)
                        continue


                    mergedverts.append(cidx)
                # if not, create a new cluster and add both the current vertex and neighbor
                else:
                    self.level_index_map[self.current_level,cidx] = nextcluster
                    self.level_index_map[self.current_level,mergevertidx] = nextcluster
                    

                    # check if the cluster is too big
                    cluster_idx = self.level_index_map[self.current_level,cidx]
                    if self.getClusterArea(cluster_idx) > 1.5*self.average_cluster_area:
                        self.verts[cidx].affinities.pop(nidx)
                        self.verts[cidx].neighbors.pop(nidx)
                        continue

                    nextcluster += 1
                    mergedverts.append(cidx)
                    mergedverts.append(mergevertidx)

        self.Nverts = np.max(self.level_index_map[-1,:])
        # update average cluster area for this level
        self.average_cluster_area = 0.0
        for idx in range(self.Nverts):
            self.average_cluster_area += self.getClusterArea(idx)
        self.average_cluster_area = self.average_cluster_area / self.Nverts
        # generate edges for new graph
        for edge in self.edges[-2]:
            newverts = []
            for vert in edge: newverts.append(self.level_index_map[self.current_level,vert]) 
            newverts = list(set(newverts)) # extract unique elements by converting to set and back again
            if len(newverts) < 2: continue # drop unitary edges
            else:
                self.edges[-1].append(newverts)

    def doFCCluster(self):
        self.nextLevel()
        self.calcAffinity()
        self.sortDegree()
        self.cluster()







