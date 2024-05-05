import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import heapq
from Node import Cell

MIN_SPACING = 0.500 #Microns, minimum spacing between cells for them to be considered a 
                    #legal placement

def lgOrdering(cell_arr):
    """
    Function to determine the order in which cells will be placed in LG

    Formula is based on priority, implementation is with Python's heapq
    priority queue

    Parameters:

        cell_arr: master_cell_array from LevelGraph object

    Return:

        ord: Array of cell indices in order of decreasing priority 
        max_height: Maximum height of all cells
    """
    k1 = 1000.0
    k2 = 1.0
    k3 = 1.0
    q = []#Min. priority queue of priorities
    L = len(cell_arr)
    max_height = -1.0
    #sum_w = 0.0#Sum of all cell widths
    for i in np.arange(L):
        priority = k1*(cell_arr[i].x - 0.5*cell_arr[i].w) + (k2*cell_arr[i].w) + (k3*cell_arr[i].h)
        #sum_w += cell_arr[i].w

        #Find max. height of all cells
        if cell_arr[i].h > max_height:
            max_height = cell_arr[i].h
    
        heapq.heappush(q, (priority, i))

    ord = np.zeros(L, dtype=np.intc)
    for i in np.arange(L-1,-1,-1):#heapq is a min-priority queue, but we want a max. priority queue
        ord[i] = heapq.heappop(q)[1]

    #avg_w_ = sum_w/L#Average cell width
    return ord, max_height

def lgMain(Hfgp, OVR_W, OVR_H):
    """
    Function that is the main for the legalization (LG) phase
    of NTUPlace3 algorithm

    Parameters:

        Hfgp_arr: Hypergraph of final resulting hypergraph after global placement(GP)
              phase of NTUPlace3
        OVR_W: Overall width of layout region
        OVR_H: Overall height of layout region

    Return:
        Hfgp with cell locations modified such that the placement is legal
        (non-overlapping)
        lg_layout: Numpy array of lists of indices of legal cell locations for 
                   visualization
    """
    sort_ind, max_h = lgOrdering(Hfgp.master_cell_array)
    n_r = int(OVR_H / (max_h + MIN_SPACING))#Num rows
    rowh = OVR_H / n_r
    row_c = np.linspace(0.0, OVR_H, num=n_r) + 0.5*rowh
    row_c = row_c[:(n_r-1)]
    print(row_c)
    lg_layout = [None]*len(row_c)#numpy array of python lists representing 
                                                       #legal layout rows and the indices of cells
                                                       #in each row, going from left to right 
                                                       #(x=0.0 to x=OVR_W)
    for idx in sort_ind:
        #Start Subroutine 330 from Patent 
        valid_locations = []#List of lists of length 3 [lg_layout_row_index, x, y] of valid placement locations
        for i in np.arange(len(lg_layout)):
            #print(row_c[i])
            if lg_layout[i] is None:
                valid_locations.append([i, MIN_SPACING + 0.5*Hfgp.master_cell_array[idx].w, row_c[i]])
            else:
                #L = len(lg_layout[i])
                left_bound = Hfgp.master_cell_array[lg_layout[i][-1]].x + 0.5*Hfgp.master_cell_array[lg_layout[i][-1]].w + MIN_SPACING
                if (left_bound + Hfgp.master_cell_array[idx].w) < OVR_W:
                    #print("Left bound: " + str(left_bound) + " Valid Location: " + str(left_bound + 0.5*Hfgp[idx].w) + " Current Cell Width: " + str(Hfgp[idx].w))
                    valid_locations.append([i, left_bound + 0.5*Hfgp.master_cell_array[idx].w, row_c[i]])
        
        #Determine best valid location by finding location closest in l1 distance to cell's (x,y) location 
        #at end of GP
        c_x = Hfgp.master_cell_array[idx].x
        c_y = Hfgp.master_cell_array[idx].y
        lg_x = 2.0*OVR_H
        lg_y = 2.0*OVR_W
        lg_i = 0
        for loc in valid_locations:
            if (np.abs(loc[1] - c_x) + np.abs(loc[2] - c_y)) < (np.abs(lg_x - c_x) + np.abs(lg_y - c_y)):
                lg_x = loc[1]
                lg_y = loc[2]
                lg_i = loc[0]
        #Now, we have the best (legal_x, legal_y) placement location of the cell
        if lg_layout[lg_i] is None:
            lg_layout[lg_i] = []
        
        lg_layout[lg_i].append(idx)
        Hfgp.master_cell_array[idx].x = lg_x
        Hfgp.master_cell_array[idx].y = lg_y

    return Hfgp, lg_layout

def visualizeLG(Hfgp_arr, lg_layout, OVR_W, OVR_H):
    """
    Function to visualize the legalized layout using matplotlib

    Parameters:
        Hfgp_arr: Master cell array of LevelGrid Object after LG
        lg_layout: List of lists of cell indices in the legalized layout
                   grid
        OVR_W: Overall width of layout region
        OVR_H: Overall heigfht of layout region
    """
    rect_arr = []
    for i in np.arange(len(lg_layout)):
        if lg_layout[i] is not None:
            for j in np.arange(len(lg_layout[i])):
                rect_arr.append(Rectangle((Hfgp_arr[lg_layout[i][j]].x - 0.5*Hfgp_arr[lg_layout[i][j]].w, Hfgp_arr[lg_layout[i][j]].y - 0.5*Hfgp_arr[lg_layout[i][j]].h),
                                           Hfgp_arr[lg_layout[i][j]].w, Hfgp_arr[lg_layout[i][j]].h))
    pc = PatchCollection(rect_arr, facecolor="blue", edgecolor="black")
    fig, ax = plt.subplots(1)
    ax.set_xlim([0.0, OVR_W])
    ax.set_ylim([0.0, OVR_H])
    ax.add_collection(pc)
    plt.show()

def main():
    OVR_W = 500.0
    OVR_H = 500.0
    N = 100
    x = np.random.uniform(50.0, 450.0, N)
    y = np.random.uniform(50.0, 450.0, N)
    w = np.random.uniform(10.0, 50.0, N)
    h = np.random.uniform(10.0, 50.0, N)
    cell_x = np.concatenate((x,y))
    Hfgp = [None]*N
    for idx in np.arange(N):
        Hfgp[idx] = Cell(cell_x[idx], cell_x[idx+N], w[idx], h[idx], w[idx]*h[idx])
    Hfgp, lg_lay = lgMain(Hfgp, OVR_W, OVR_H)
    visualizeLG(Hfgp, lg_lay, OVR_W, OVR_H)

if __name__ == "__main__":
    main()