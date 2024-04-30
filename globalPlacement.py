import numpy as np
from Node import Node, Cell, Partition
from Multilevel import LevelGraph
from verilog2hgr import parse_func
from spef_extractor.lef_def_parser import LefParser, DefParser

def log_sum_exp(x_k):
    return np.log(np.sum(np.exp(x_k/gamma)))

def W(LG):
    '''
    Compute the smooth wirelength function
    Parameters:
        LG: level graph object
    Returns:
        W: scalar wirelength
    '''
    xvec = LG.xvec() # TODO make sure this comes from the current level
    yvec = LG.yvec()
    lidx = LG.current_level # TODO update the current level during coarsening/relaxing
    
    w = np.zeros(len(xvec))
    
    for edge in LG.edges[current_level]:
            txvec = xvec[edge]
            tyvec = yvec[edge]
            t = log_sum_exp(txvec)
            t2 = log_sum_exp(-txvec)
            t3 = log_sum_exp(tyvec)
            t4 = log_sum_exp(-tyvec)
            w += t + t2 + t3 + t4     
    return gamma*(w)

def grad_W(LG):
    '''
    Compute the gradient of the smooth wirelength wrt vertex positions
        Parameters:
            LG: level graph object
        Returns:
            gradient vector in 2*Nverts x 1
    '''
    xvec = LG.xvec()
    yvec = LG.yvec()
    lidx = LG.current_level

    grad_w = np.zeros((len(xvec),2))

    for edge in LG.edges[current_level]:
        txvec = xvec[edge]
        tyvec = yvec[edge]
        pos_x_den = np.sum(np.exp(txvec))
        neg_x_den = np.sum(np.exp(-txvec))
        pos_y_den = np.sum(np.exp(tyvec))
        neg_y_den = np.sum(np.exp(-tyvec))
        for vidx in edge:
            grad_w[vidx,0] += (np.exp(xvec[vidx])/pos_x_den) - (np.exp(-xvec)/neg_x_den)
            grad_w[vidx,1] += (np.exp(yvec[vidx])/pos_y_den) - (np.exp(-yvec)/neg_y_den)

    return np.flatten(grad_w)
    

def calcCellPotential(vx, vy, bx, by, wv, hv, wb, hb):
    """
    Function to calculate the total potential movable area for a specific 
    cell/vertex inside a specific bin b

    Parameters:
        v_x, v_y: Center x- and y- coordinate, respectively, of this vertex
        b_x, b_y: Center x- and y- coordinate, respectively, of this bin
        wv, hv: Width and height, respectively, of this vertex
        wb, hb: Width and height, respectively, of this bin

    Return:
        pot: potential of this cell in this bin
    """
    pot = 0.0
    #Calculate p_x(b,v)
    p_x = 0.0
    p_y = 0.0
    abs_dx = np.abs(vx - bx)
    a = 4.0/((wv + 2.0*wb)*(wv + 4.0*wb))
    b = 2.0/(wb*(wv + 4.0*wb))
    if abs_dx <= (0.5*wv + wb):
        p_x = 1.0 - a*(abs_dx**2)
    elif (abs_dx >= (0.5*wv + wb)) and (abs_dx <= (0.5*wv + 2.0*wb)):
        p_x = b*((abs_dx - 0.5*wv - 2.0*wb)**2)
    else:
        p_x = 0.0
    #Calculate p_y(b,v)
    abs_dy = np.abs(vy - by)
    a = 4.0/((hv + 2.0*hb)*(hv + 4.0*hb))
    b = 2.0/(hb*(hv + 4.0*hb))
    if abs_dy <= (0.5*hv + hb):
        p_y = 1.0 - a*(abs_dy**2)
    elif (abs_dy >= (0.5*hv + hb)) and (abs_dy <= (0.5*hv + 2.0*hb)):
        p_y = b*((abs_dy - 0.5*hv - 2.0*hb)**2)
    else:
        p_y = 0.0
            
    pot = p_x*p_y
        
    return pot

def calcOvrPotential(x, w, h, bin_x, bin_y):
    """
    Function to calculate the total potential movable area for each cell
    over all bins

    Parameters:
        x: Vector of x,y coordinates of vertex/cell/cluster centers, x-coordinates first
        w: Array of vertex/cell widths. w[i] is the width of vertex/cell with center (x[i],x[i+N])
        h: Array of vertex/cell heights. h[i] is the height of vertex/cell with center (x[i],x[i+N])
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'

    Return:
        ovr_pots: Array of overall potential movable area for each cell
    """
    N = x.shape[0] // 2
    Nx = bin_x.shape[0]#Grid should be square
    Ny = bin_y.shape[1]#Grid should be square
    wb = bin_x[1][0] - bin_x[0][0]
    hb = bin_y[0][1] - bin_y[0][0]
    ovr_pots = np.zeros(N, dtype=np.double)
    for k in np.arange(N):
        pot = 0.0#Potential of this cell
        for j in np.arange(Ny - 1):
            for i in np.arange(Nx - 1):
                cbin_x = bin_x[i][j] + 0.5*wb
                cbin_y = bin_y[i][j] + 0.5*hb
                if (np.abs(x[k] - cbin_x) < (0.5*w[k] + 2.0*wb)) or (np.abs(x[k+N] - cbin_y) < (0.5*h[k] + 2.0*hb)):
                    pot += calcCellPotential(x[k], x[k+N], cbin_x, cbin_y, w[k], h[k], wb, hb)
        ovr_pots[k] = pot

    return ovr_pots

def Db(vx, vy, bx, by, wv, hv, wb, hb, ov_pots):
    """
    Function to calculate the bin density Db contribution of
    a particular vertex v in a particular bin b

    Parameters:
        v_x, v_y: Center x- and y- coordinate, respectively, of this vertex v
        b_x, b_y: Center x- and y- coordinates, respectively, of this bin b
        wv, hv: Widths and height, respectively, of this vertex v
        wb, hb: Width and height, respectively, of this bin b
        ov_pots: Total potential movable area for this vertex v

    Return:
        Term in sum for Db for this vertex v inside this bin b
    """
    val = 0.0
    #Calculate p_x(b,v)
    p_x = 0.0
    p_y = 0.0
    abs_dx = np.abs(vx - bx)
    a = 4.0/((wv + 2.0*wb)*(wv + 4.0*wb))
    b = 2.0/(wb*(wv + 4.0*wb))
    if abs_dx <= (0.5*wv + wb):
        p_x = 1.0 - a*(abs_dx**2)
    elif (abs_dx >= (0.5*wv + wb)) and (abs_dx <= (0.5*wv + 2.0*wb)):
        p_x = b*((abs_dx - 0.5*wv - 2.0*wb)**2)
    else:
        p_x = 0.0
    #Calculate p_y(b,v)
    abs_dy = np.abs(vy - by)
    a = 4.0/((hv + 2.0*hb)*(hv + 4.0*hb))
    b = 2.0/(hb*(hv + 4.0*hb))
    if abs_dy <= (0.5*hv + hb):
        p_y = 1.0 - a*(abs_dy**2)
    elif (abs_dy >= (0.5*hv + hb)) and (abs_dy <= (0.5*hv + 2.0*hb)):
        p_y = b*((abs_dy - 0.5*hv - 2.0*hb)**2)
    else:
        p_y = 0.0
    cv = (wv*hv) / ov_pots#Normalization coefficient
    val += cv*p_x*p_y
        
    return val

def fDb(x, w, h, bin_x, bin_y, ovr_pots, td=0.6):
    """
    Function to calculate the sum of the Db's over all the bins

    Parameters:
        x: Array of x,y coordinates of vertex/cell/cluster centers, x-coordinates first
        w: Array of vertex/cell widths. w[i] is the width of vertex/cell with center (x[i],x[i+N])
        h: Array of vertex/cell heights. h[i] is the height of vertex/cell with center (x[i],x[i+N])
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pots: Array of total potential movable area for each cell over all bins
        td: Target density (0.6 by default)

    Return:
        Value of sum of the Db's over all bins
    """
    N = x.shape[0] // 2
    Nx = bin_x.shape[0]#Grid should be square
    Ny = bin_y.shape[1]#Grid should be square
    wb = bin_x[1][0] - bin_x[0][0]
    hb = bin_y[0][1] - bin_y[0][0]
    Mb = td*wb*hb#All area inside the bin can be moved - there are no pre-placed cells
    fval = 0.0
    for j in np.arange(Ny - 1):
        for i in np.arange(Nx - 1):
            #Potential movable area calculation requires all verticies
            sumDb = 0.0
            for k in np.arange(N):
                cbin_x = bin_x[i][j] + 0.5*wb
                cbin_y = bin_y[i][j] + 0.5*hb
                if (np.abs(x[k] - cbin_x) < (0.5*w[k] + 2.0*wb)) or (np.abs(x[k+N] - cbin_y) < (0.5*h[k] + 2.0*hb)):
                    sumDb += Db(x[k], x[k+N], cbin_x, cbin_y, w[k], h[k], wb, hb, ovr_pots[k])

            fval += (sumDb - Mb)**2#After computing sum over this bin, calculate penalty term for this bin
        
    return fval

def f(x):
    return x

def grad_f(x):
    return x

def armijo(a_0, pk, xk):
    """
    Auxiliary function to perform Armijo backtracking to obtain a step length alpha_k that leads to
    convergence, i.e. satisfies the Wolfe Conditions

    Parameters:
        a_0: Initial guess for step length (1.0)
        pk: Descent direction
        xk: Value of x at current iteration

    Return:
        alpha: Step length that satisfies Wolfe Conditions
    """
    alpha = a_0
    rho = 0.5#Needs to be between 0 and 1
    c1 = 0.25
    while (f(xk + alpha*pk) > f(xk) + alpha*c1*np.dot(grad_f(xk), pk)):
        alpha = rho*alpha
        
    return alpha
    
def bfgs(x0, H0, eps=1.0e-3):
    """

    Function 
    Parameters:
        x0: Initial guess for trajectory being minimized x
        H0: Initial guess for approximate Hessian

    Return:
        x_new: Minimized trajectory x(t)
        f_vals: Values of f(x(t)) at each iteration
        grad_norms: values of del(f(x(t))) at each iteration
    """
    #Re-use x_prev is x_k in BFGS iteration, x_new is x_k+1, Hk is approximate inverse of del^2 f(x)
    k = 0
    x_prev = x0
    Hk = H0
    Hk1 = H0
    x_new = x0
    f_vals = []
    f_vals.append(f(x_prev))
    grad_norms = []
    grad_norms.append(np.linalg.norm(grad_f(x_prev)))
    while(np.linalg.norm(grad_f(x_prev)) >= eps):
        grad_f_p = grad_f(x_prev)
        print("Iteration " + str(k) + ":        f(x_k): " + str(f(x_prev)) + "        ||grad_f(x_k)||_2: " + str(np.linalg.norm(grad_f_p)))
        pk = -Hk@grad_f_p#Descent direction
        ak = armijo(1.0, pk, x_prev)#Armijo backtracking routine to find alpha_k that will lead to convergence
        x_new = x_prev + ak*pk
        sk = x_new - x_prev
        grad_f_n = grad_f(x_new)
        yk = grad_f_n - grad_f_p
        rho_k = 1.0/np.dot(yk,sk)
        #Expand BFGS update formula 6.17 in Ch.6 of Nocedal and write to calculate H_k+1 using only matrix addition,
        #matrix-vector products and outer products of vectors so total time complexity of update is O(n^2)
        uk = Hk@yk
        vk = (Hk.T)@yk
        Hk1 = Hk - rho_k*np.outer(sk,vk) - rho_k*np.outer(uk,sk) + ((rho_k**2)*np.dot(yk,uk))*np.outer(sk,sk) + rho_k*np.outer(sk,sk)
        f_vals.append(f(x_new))
        grad_norms.append(np.linalg.norm(grad_f(x_new)))
        x_prev = x_new
        Hk = Hk1
        k += 1
    return x_new

def isCellInside(bin_tlx, bin_tly, bin_w, bin_h, cell):
    """
    Auxiliary function to determine if a cell is inside a bin
    """
    if (cell.cx >= bin_tlx) and (cell.cx < (bin_tlx + bin_w)) and (cell.cy >= bin_tly) and (cell.cy < (bin_tly + bin_h)):
        return True
    return False

def calcDb(bin_tlx, bin_tly, bin_w, bin_h, level, cell_array, i_map):
   return 

def calcOverflowRatio(bin_x, bin_y, level, cell_nums, i_map, cell_arr, Mb, OVR_AREA):
    """
    Function to calculate the total overflow ratio at each level in the
    V-Cycle iteration

    Parameters:
        bin_x: numpy meshgrid of x-coordinates of bin grid
        bin_y: numpy meshgrid of y-coordinates of bin grid
        level: Current level in V-Cycle iteration
        cell_nums: Array of cell numbers at this level
        i_map: Index map of cell indices at all different levels
        cell_arr: Array of all cell objects
        Mb: Max. allowable area in each bin
        OVR_AREA: Overall area of layout region

    Return:
        Overflow ratio of cell layout at this level in V-Cycle
    """
    nx = bin_x.shape[0]
    ny = bin_y.shape[1]
    max_overflow = 0.0
    bin_w = bin_x[1][0] - bin_x[0][0]
    bin_h = bin_y[0][1] - bin_y[0][0]
    for j in np.arange(ny-1):
        for i in np.arange(nx-1):
            Db = calcDb(bin_x[i][j], bin_y[i][j], bin_w, bin_h, level, cell_nums, i_map, cell_arr)
            max_overflow += np.max(np.array([Db - Mb,0.0], dtype=np.double))

    return (max_overflow / OVR_AREA)

def read_hgr(hgr_file):
    '''
    Reads the text file produced by verilog2hgr.py and returns the header and list of edges
    Each edge is also a list of cell indexes
    '''
    with open(hgr_file) as file:
        edges = [line.split() for line in file]
    header=edges.pop(0)

    for idx, edge in enumerate(edges):
        edges[idx] = [int(vertex)-1 for vertex in edge]
        
    return header, edges

def populateCells(verilog_netlist, hgr_filename):
    """
    Function to read a gate level verilog netlist and generate the hypergraph
    representation and get information about instantiated cells from standard
    cell library.

    Parameters:
        verilog_netlist: file name string of gate level verilog netlist
        hgr_filename: file name string where the hgr file will be found, default 
        to the verilog netlist name with extension changed from *.vg to *.hgr

    Return:
        master_cell_array: list of cell objects instantiated from netlist
        master_hg: master hypergraph as a list of edges with vertex indices 
        referring to the master_cell_array

    """
    # convert verilog netlist to hypergraph and write to file  using verilog2hgr
    parse_func(verilog_netlist)
    header,master_hg = read_hgr(hgr_file)
    num_nets = int(header[0])
    master_num_cells = int(header[1])
    master_cell_array = []
    for cidx in range(master_num_cells):
        width = 0
        height = 0
        master_cell_array.append(Cell(1.0+5*cidx,1, width,height,width*height)) 

    # parse def file for cell information
    def_file = 'csmplace1/KSA16_yosys.def'
    defparser = DefParser(def_file)
    defparser.parse()
    # build the indexed list of cell types (macros) used in the design
    ctypes = []
    for comp in defparser.components:
        if 'tap' in comp.macro: # tap cells not placed
            continue
        else:
            ctypes.append(comp.macro.replace('hd','hs'))
    if len(ctypes) != len(master_cell_array):
        print('Something is wrong with macro list')

    for cidx in range(master_num_cells):
        lef_file = f'csmplace1/cell_lefs/{ctypes[cidx]}.lef'
        lefparser = LefParser(lef_file)
        lefparser.parse()
        master_cell_array[cidx].w = lefparser.macro_dict[ctypes[cidx]].info['SIZE'][0]
        master_cell_array[cidx].h = lefparser.macro_dict[ctypes[cidx]].info['SIZE'][1]
        master_cell_array[cidx].area = master_cell_array[cidx].w * master_cell_array[cidx].h

    return master_cell_array, master_hg


def gpMain(H_0, N_MAX, OVR_W, OVR_H):
    """
    Main function for GP phase of placement algorithm. All other
    subroutines used in GP are called from this function

    Parameters:
        H0: Hypergraph that represents circuit at the closest level
        N_MAX: Max. number of cells allowed at coarsest level
        OVR_W: Overall width, in microometers, of layout region
        OVR_H: Overall height, in micrometers, of layout region

    Return:
        H0 with cell positions modified so the placement is
        optimal, but may have intersects. This will be passed to
        lgMain() for legalization (LG) phase
    """
    H_current = H_0#Current level
    level = 0#Current level, 0 is the finest/least clustered

    #Descent/Coarsening/Decimation Stage of V-Cycle Multigrid
    while (H_current.Nverts > N_MAX):
        level += 1
        H_current.doFCCluster()
        
    # Do initial placement at coarsest level. Topologically sort cells at closest level and arrange in grid
    H_current.doInitialPlace()
    #Ascent/Refining/Interpolation Stage of V-Cycle Multigrid
    for i in range(level,-1,-1):
        grid_nw = int(np.sqrt(H_current.Nverts))#Make grid of bins square by default
        grid_nh = grid_nw
        bins_x, bins_y = np.meshgrid(np.linspace(0.0, OVR_W, grid_nw, dtype=np.double), np.linspace(0.0, OVR_H, grid_nh, dtype=np.double), indexing='ij')
        bin_area = (bins_x[1][0] - bins_x[0][0])*(bins_y[0][1] - bins_y[0][0])#Mb
        #No need for base potentials as we have no pre-placed blocks
        vvec = H_current.vvec()
        lambda_m = np.linalg.norm(W(vvec), 1) / np.linalg.norm(Db(vvec), 1)#Initialize lambda to be 1-norm of gradient
        lambda_m1 = lambda_m
        prev_overflow_ratio = calcOverflowRatio(bins_x, bins_y, i, H_current.verts, H_current.level_index_map, H_current.master_cell_array, bin_area, OVR_W*OVR_H)
        new_overflow_ratio = 100.0
        m = 0

        #do-while optimization loop
        while(True):
            x_new = bfgs(x, y)
            #Update this level's cluster/cell/vertex coordinates using x_new
            m += 1
            lmabda_m1 = 2.0*lambda_m
            #No need for look-ahead LG
            new_overflow_ratio = calcOverflowRatio(bins_x, bins_y, i, H_current.verts, H_current.level_index_map, H_current.master_cell_array, bin_area, OVR_W*OVR_H)
            if (new_overflow_ratio - prev_overflow_ratio >= -0.01):#If reduction in overflow ratio is >= 1%, keep going with bfgs
                continue
            else:#If reduction in overflow ratio is < 1%, stop with bfgs
                break
        
        #Post-processing after optimization loop to spread cells
        part = Partition(bin_area, H_current.verts)
        part.construct(0.0, 0.0, OVR_W, OVR_H)
        ws = part.calcWhiteSpace(part.topNode)
        part.WSA(part.topNode)
        part.printLeaves(part.topNode)
        part.updateCells(part.topNode)

        #De-cluster and update cell positions after WSA()

    H_0 = H_current#Return original hypergraph at finest/least clustered level with final GP cell positions
    return H_0


# Testing GP
if (__name__ == '__main__'):
    verilog_netlist = 'csmplace1/KSA16_yosys.vg'
    hgr_file = 'csmplace1/KSA16_yosys.hgr'
    master_cell_array, master_hg = populateCells(verilog_netlist,hgr_file)
    H_0 = LevelGraph(master_hg,master_cell_array)

    #OVR_W and OVR_H taken from openROAD floorplan
    OVR_W = 89.395
    OVR_H = 89.395

    # require number of vertices at coarsest level to be half the number of cells
    N_MAX = H_0.Nverts/2

    H_0 = gpMain(H_0,N_MAX,OVR_W,OVR_H)



