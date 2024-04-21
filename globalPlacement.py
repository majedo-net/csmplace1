import numpy as np
from Node import Node, Cell, Partition
from Multilevel import LevelGraph

def W(x,y):
    return x

def Db(x,y):
    return x

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
    print()

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

def gpMain(H_0, N_MAX, OVR_W, OVR_H):
    """
    Main function for GP phase of placement algorithm. All other
    subroutines used in GP are called from this function

    Parameters:
        H0: Hypergraph that represents circuit at the closest level
        N_MAX: Max. number of cells allowed at coarsest level
        OVR_W: Overall width, in nanometers, of layout region
        OVR_H: Overall height, in nanometers, of layout region

    Retrurn:
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

    #Do initial placement at coarsest level - a few options for this:
    #1. Use texternal library LEF/DEF floorplanner
    #2. Topologically sort cells at closest level and arrange in grid
    #3. Run BFGS with initial random cell coordinates

    #Ascent/Refining/Interpolation Stage of V-Cycle Multigrid
    for i in range(level,-1,-1):
        grid_nw = int(np.sqrt(H_current.Nverts))#Make grid of bins square by default
        grid_nh = grid_nw
        bins_x, bins_y = np.meshgrid(np.linspace(0.0, OVR_W, grid_nw, dtype=np.double), np.linspace(0.0, OVR_H, grid_nh, dtype=np.double), indexing='ij')
        bin_area = (bins_x[1][0] - bins_x[0][0])*(bins_y[0][1] - bins_y[0][0])#Mb
        #No need for base potentials as we have no pre-placed blocks
        x = np.zeros(1, dtype=np.double)#Center x-coordinates of cells
        y = np.zeros(1, dtype=np.double)#Center y-coordinates of cells
        lambda_m = np.linalg.norm(W(x,y), 1) / np.linalg.norm(Db(x,y), 1)#Initialize lambda to be 1-norm of gradient
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





