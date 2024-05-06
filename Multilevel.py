import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

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
        self.w = None
        self.h = None
        self.area = None

    def setVertexPosition(self,x,y):
        '''
        The vertex position is applied to the cells it contains
        '''
        self.x = x
        self.y = y

    def updateVertex(self,mca):
        '''
        Updates the size of the vertex based on constinuent cells.
        Assumes that cell_idxes parameter has already been updated.
        '''
        w = 0.0
        h = 0.0
        a = 0.0
        for cidx in self.cell_idxes:
            w += mca[cidx].w
            h += mca[cidx].h
            a += mca[cidx].w * mca[cidx].h
        self.w = w
        self.h = h
        self.area = a
        
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
    OVR_H=0.0
    OVR_W=0.0


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
            self.verts.append(Vertex([idx]))
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

    def wvec(self):
        '''
        helper function to generate vector of cell/vertex widths 
        '''
        vw = []
        for vert in self.verts:
           vw.append(vert.w)
        return np.asarray(vw)

    def hvec(self):
        '''
        helper function to generate a vector of cell/vertex heights
        '''
        vh = []
        for vert in self.verts:
            vh.append(vert.h)
        return np.asarray(vh)

    
    def getMaxArea(self):
        ma = 0.00
        for vert in self.verts:
            if vert.area > ma: ma = vert.area
        return ma
    
    def deCluster(self,xnew):
        '''
        Function to update vertex positions and decluster to the next finest 
        level

        Parameters:
            xnew: new positions of 2*Nverts with x values followed by y values
        '''
        self.updatePositions(xnew)
        self.updateCellPositions()
        self.current_level -= 1
        self.Nverts = np.max(self.level_index_map[self.current_level,:])+1
        # generate vertex objects for new graph
        self.verts.clear()
        for idx in range(self.Nverts):
            cell_list = []
            for cidx in range(self.Ncells):
                if self.level_index_map[self.current_level,cidx]==idx:
                    cell_list.append(cidx)
            self.verts.append(Vertex(cell_list))
            (self.verts[idx].x, self.verts[idx].y) = self.getClusterCenter(idx)
            #print("declustered x: " + str(self.verts[idx].x) + "  declustered y: " + str(self.verts[idx].y))
        self.updateVerts()

    def updatePositions(self,xnew):
        '''
        update the positions of the vertices with new iteration of optimized 
        positions
        Parameters:
            xnew: new positions of 2*Nverts with x values followed by y values
        '''
        for vidx in range(self.Nverts):
            self.verts[vidx].x = xnew[vidx]
            self.verts[vidx].y = xnew[self.Nverts+vidx]
    
    def updateCellPositions(self):
        '''
        Update cell objects in master array with vertex positions
        '''
        for idx in range(self.Nverts):
            for cidx in range(self.Ncells):
                if idx == 0:
                    print("cell x: " + str(self.master_cell_array[cidx].x) + "  cell y: " + str(self.master_cell_array[cidx].y))
                if self.level_index_map[self.current_level,cidx] == idx:
                    self.master_cell_array[cidx].x = self.verts[idx].x
                    self.master_cell_array[cidx].y = self.verts[idx].y
        
        #For cells that escape clustering, place them like in initial placement
        row =0
        column = 0
        x0 = self.OVR_W /2
        y0 = self.OVR_H/2
        x_ = x0
        y_ = y0
        grid_nw = int(np.sqrt(self.Nverts))#Make grid of bins square by default
        for cidx in range(self.Ncells):
            if (self.master_cell_array[cidx].x >= self.OVR_W) or (self.master_cell_array[cidx].x <= 0.0) or (self.master_cell_array[cidx].y >= self.OVR_H) or (self.master_cell_array[cidx].y <= 0.0): 
                self.master_cell_array[cidx].x = x_
                self.master_cell_array[cidx].y = y_
                t_area = self.master_cell_array[cidx].area
                t_wh = np.sqrt(t_area)
                if column >= grid_nw:
                    row += 1
                    column = 0
                    y_ += t_wh
                    x_ = x0
                else:
                    column += 1
                    x_ += t_wh


    def plotVerts(self,filename=None):
        '''
        Plot the current positions and sizes of vertices 
        Parameters:
            filename: string to name the plot to save to file
        '''
        rect_array = []
        for vert in self.verts:
            rect_array.append(Rectangle((vert.x-0.5*vert.w , vert.y - 0.5*vert.h),vert.w,vert.h))
        pc = PatchCollection(rect_array, facecolor="blue", edgecolor="black")
        fig, ax = plt.subplots(1)
        ax.set_xlim([0.0, self.OVR_W])
        ax.set_ylim([0.0, self.OVR_H])
        ax.add_collection(pc)
        plt.show()
        # saving png instead of showing because remote machine
        if filename is not None:
            plt.savefig(f'{filename}.png',bbox_inches='tight')

    def nextLevel(self):
        '''
        Initialize the next level of coarsening with vertex objects and their neighbors
        '''
        self.current_level += 1
        if self.current_level>1: 
            self.level_index_map=np.vstack((self.level_index_map,-1*np.ones_like(self.level_index_map[0,:])))
        else:
            self.level_index_map=np.vstack((self.level_index_map,-1*np.ones_like(self.level_index_map)))

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
                xc += self.master_cell_array[cell].x
                yc += self.master_cell_array[cell].y

        if nc == 0:
            xc = 0.5*self.OVR_W
            yc = 0.5*self.OVR_H
            nc += 1


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
        x0= self.OVR_W /2
        y0= self.OVR_H/2
        x = x0
        y = y0
        for cidx in self.degree.keys():
            self.verts[cidx].setVertexPosition(x,y)
            t_area = self.getClusterArea(cidx)
            t_wh = np.sqrt(t_area)
            if column >= grid_nw:
                row += 1
                column = 0
                y += t_wh
                x = x0
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
        nextcluster=0 # running index of clusters
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

        self.Nverts = np.max(self.level_index_map[-1,:])+1
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
        # generate vertex objects for new graph
        self.verts.clear()
        for idx in range(self.Nverts):
            cell_list = []
            for cidx in range(self.Ncells):
                if self.level_index_map[self.current_level,cidx]==idx:
                    cell_list.append(cidx)
            self.verts.append(Vertex(cell_list))
            for edge in self.edges[-1]:
                if idx in edge:
                    self.verts[-1].neighbors=self.verts[-1].neighbors.union(edge)
            # after using the set type to eliminate duplicates, convert neighbors to a list so we can index it later
            self.verts[-1].neighbors = list(self.verts[-1].neighbors)

    def updateVerts(self):
        for vert in self.verts:
            vert.updateVertex(self.master_cell_array)


    def doFCCluster(self):
        self.nextLevel()
        self.calcAffinity()
        self.sortDegree()
        self.cluster()
        self.updateVerts()







