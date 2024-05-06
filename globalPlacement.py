import numpy as np
from Node import Node, Cell, Partition
from Multilevel import LevelGraph
from verilog2hgr import parse_func
from spef_extractor.lef_def_parser import LefParser, DefParser

def log_sum_exp(x_k):
    gamma = 0.89
    x_k = x_k.astype(float)
    #print(np.exp(x_k/gamma))
    return np.log(np.sum(np.exp(x_k/gamma)))

def W(LG,pk=None):
    '''
    Compute the smooth wirelength function
    Parameters:
        LG: level graph object
        pk: indicates if a descent step should be added x
    Returns:
        W: scalar wirelength
    '''
    xvec = LG.xvec() 
    yvec = LG.yvec()
    if pk is not None:
        N = len(xvec)
        xvec += pk[0:N]
        yvec += pk[N:]
    lidx = LG.current_level # TODO update the current level during coarsening/relaxing
    
    w = 0.0
    
    for edge in LG.edges[lidx]:
            txvec = xvec[edge]
            tyvec = yvec[edge]
            t = log_sum_exp(txvec)
            t2 = log_sum_exp(-txvec)
            t3 = log_sum_exp(tyvec)
            t4 = log_sum_exp(-tyvec)
            w += t + t2 + t3 + t4     
    return gamma*(w)

def grad_W(LG,pk=None):
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
    if pk is not None:
        N = len(xvec)
        xvec += pk[0:N]
        yvec += pk[N:-1]

    grad_w = np.zeros((len(xvec),2))

    for edge in LG.edges[lidx]:
        txvec = xvec[edge]
        tyvec = yvec[edge]
        pos_x_den = np.sum(np.exp(txvec))
        neg_x_den = np.sum(np.exp(-txvec))
        pos_y_den = np.sum(np.exp(tyvec))
        neg_y_den = np.sum(np.exp(-tyvec))
        for vidx in edge:
            grad_w[vidx,0] += (np.exp(xvec[vidx])/pos_x_den) - (np.exp(-xvec[vidx])/neg_x_den)
            grad_w[vidx,1] += (np.exp(yvec[vidx])/pos_y_den) - (np.exp(-yvec[vidx])/neg_y_den)

    return grad_w.flatten()
    

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

def calcOvrPotential(LG, bin_x, bin_y):
    """
    Function to calculate the total potential movable area for each cell
    over all bins

    Parameters:
        LG: LevelGraph object containing position and size information about the cells/vertices
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'

    Return:
        ovr_pots: Array of overall potential movable area for each cell
    """
    x = LG.vvec()
    w = LG.wvec()
    h = LG.hvec()
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

def dDb(vx, vy, bx, by, wv, hv, wb, hb, ov_pots):
    """
    Function to calculate the derivative of the contribution of a specific
    vertex v to the bin density Db of a specific bin b

    Parameters:
        v_x, v_y: Center x- and y- coordinate, respectively, of
                  this vertex v
        b_x, b_y: Center x- and y- coordinates, respectively, of this bin b
        wv, hv: Widths and height, respectively, of this vertex v
        wb, hb: Width and height, respectively, of this bin b
        ov_pots: Total potential of this vertex v

    Return:
        res_x: Gradient component with respect to x of this vertex
        res_y: Gradient component with respect to y of this vertex
    """
    res_x = 0.0
    res_y = 0.0
    #Calculate p_x(b,v)
    p_x = 0.0
    p_x_dx = 0.0
    p_y = 0.0
    p_y_dy = 0.0
    abs_dx = np.abs(vx - bx)
    dx = vx - bx
    a = 4.0/((wv + 2.0*wb)*(wv + 4.0*wb))
    b = 2.0/(wb*(wv + 4.0*wb))
    if abs_dx <= (0.5*wv + wb):
        p_x = 1.0 - a*(abs_dx**2)
        p_x_dx = -2.0*a*dx
        if dx < 0:
            p_x_dx = 2.0*a*(bx - vx)
    elif (abs_dx >= (0.5*wv + wb)) and (abs_dx <= (0.5*wv + 2.0*wb)):
        p_x = b*((abs_dx - 0.5*wv - 2.0*wb)**2)
        p_x_dx = 2.0*b*(dx - 0.5*wv - 2.0*wb)
        if dx < 0:
            p_x_dx = -2.0*b*(bx - vx - 0.5*wv - 2.0*wb)
    else:
        p_x = 0.0
        p_x_dx = 0.0
    #Calculate p_y(b,v)
    abs_dy = np.abs(vy - by)
    dy = vy - by
    a = 4.0/((hv + 2.0*hb)*(hv + 4.0*hb))
    b = 2.0/(hb*(hv + 4.0*hb))
    if abs_dy <= (0.5*hv + hb):
        p_y = 1.0 - a*(abs_dy**2)
        p_y_dy = -2.0*a*dy
        if dy < 0:
            p_y_dy = 2.0*a*(by - vy)
    elif (abs_dy >= (0.5*hv + hb)) and (abs_dy <= (0.5*hv + 2.0*hb)):
        p_y = b*((abs_dy - 0.5*hv - 2.0*hb)**2)
        p_y_dy = 2.0*b*(dy - 0.5*hv - 2.0*hb)
        if dy < 0:
            p_y_dy = -2.0*b*(by - vy - 0.5*hv - 2.0*hb)
    else:
        p_y = 0.0
        p_y_dy = 0.0
    cv = (wv*hv) / ov_pots
    #Update res_x, res_y
    res_x = cv*p_x_dx*p_y
    res_y = cv*p_x*p_y_dy
        
    return res_x, res_y

def maxDb(LG, bin_x, bin_y, ovr_pots, td=0.6):
    """
    Function to calculate maximum cell Db for each bin

    Parameters:
        LG: LevelGrid object, contains all information about cell positions and sizes
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pots: Array of total potential movable area for each cell over all bins
        td: Target density (0.6 by default)

    Return:
        Value of sum of the Db's over all bins
    """
    x = LG.vvec()
    w = LG.wvec()
    h = LG.hvec()
    N = x.shape[0] // 2
    Nx = bin_x.shape[0]#Grid should be square
    Ny = bin_y.shape[1]#Grid should be square
    wb = bin_x[1][0] - bin_x[0][0]
    hb = bin_y[0][1] - bin_y[0][0]
    Mb = td*wb*hb#All area inside the bin can be moved - there are no pre-placed cells
    bin_dB = np.zeros_like(bin_x)
    for j in np.arange(Ny - 1):
        for i in np.arange(Nx - 1):
            #Potential movable area calculation requires all verticies
            maxDb = 0.0
            for k in np.arange(N):
                cbin_x = bin_x[i][j] + 0.5*wb
                cbin_y = bin_y[i][j] + 0.5*hb
                if (np.abs(x[k] - cbin_x) < (0.5*w[k] + 2.0*wb)) or (np.abs(x[k+N] - cbin_y) < (0.5*h[k] + 2.0*hb)):
                    this_Db = Db(x[k], x[k+N], cbin_x, cbin_y, w[k], h[k], wb, hb, ovr_pots[k])

                    if this_Db > maxDb: maxDb = this_Db
                        
            bin_dB[i][j] = maxDb
    return bin_dB
def fDb(LG, bin_x, bin_y, ovr_pots, pk=None, td=0.6):
    """
    Function to calculate the sum of the Db's over all the bins

    Parameters:
        LG: LevelGrid object, contains all information about cell positions and sizes
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pots: Array of total potential movable area for each cell over all bins
        pk: descent step to be added to x during linesearch
        td: Target density (0.6 by default)

    Return:
        Value of sum of the Db's over all bins
    """
    x = LG.vvec()
    if pk is not None:
        x += pk
    w = LG.wvec()
    h = LG.hvec()
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

def grad_fDb(LG, bin_x, bin_y, ovr_pots, pk=None, td=0.6):
    """
    Function to calculate the sum of the grad_Db's over all the bins

    Parameters:
        LG: LevelGraph object containing information about the cells
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pots: Array of total potential movable area for each cell over all bins
        pk: descent step to be added to x during linesearch
        td: Target density (0.6 by default)

    Return:
        Value of sum of the Db's over all bins
    """
    x = LG.vvec()
    if pk is not None:
        x += pk
    w = LG.wvec()
    h = LG.hvec()
    N = x.shape[0] // 2
    del_f = np.zeros(x.shape[0], dtype=np.double)
    Nx = bin_x.shape[0]#Grid should be square
    Ny = bin_y.shape[1]#Grid should be square
    wb = bin_x[1][0] - bin_x[0][0]
    hb = bin_y[0][1] - bin_y[0][0]
    Mb = td*wb*hb#All area inside the bin can be moved - there are no pre-placed cells
    for j in np.arange(Ny - 1):
        for i in np.arange(Nx - 1):
            #Potential movable area calculation uses all verticies
            res = np.zeros(2*N, dtype=np.double)
            sumDb = 0.0
            for k in np.arange(N):
                cbin_x = bin_x[i][j] + 0.5*wb
                cbin_y = bin_y[i][j] + 0.5*hb
                if (np.abs(x[k] - cbin_x) < (0.5*w[k] + 2.0*wb)) or (np.abs(x[k+N] - cbin_y) < (0.5*h[k] + 2.0*hb)):
                    res_x, res_y = dDb(x[k], x[k+N], cbin_x, cbin_y, w[k], h[k], wb, hb, ovr_pots[k])
                    res[k] = res_x
                    res[k+N] = res_y
                    sumDb += Db(x[k], x[k+N], cbin_x, cbin_y, w[k], h[k], wb, hb, ovr_pots[k])
            
            del_f = del_f + (2.0*(sumDb - Mb)*res)
            
    return del_f

def f(LG,bins_x, bins_y, over_pots,pk=None):
    return W(LG,pk) + fDb(LG,bins_x,bins_y,over_pots,pk)

def grad_f(LG,bins_x,bins_y,over_pots,pk=None):
    return grad_W(LG,pk) + grad_fDb(LG,bins_x,bins_y,over_pots,pk)

def armijo(a_0, pik, LG,bins_x,bins_y,ovr_pot):
    """
    Auxiliary function to perform Armijo backtracking to obtain a step length alpha_k that leads to
    convergence, i.e. satisfies the Wolfe Conditions

    Parameters:
        a_0: Initial guess for step length (1.0)
        pik: Descent direction
        LG: LevelGraph object
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pot: Array of total potential movable area for each cell over all bins
    Return:
        alpha: Step length that satisfies Wolfe Conditions
    """
    alpha = a_0
    rho = 0.8#Needs to be between 0 and 1
    c1 = 0.01
    while (f(LG,bins_x,bins_y,ovr_pot,pk=alpha*pik) > 
           f(LG, bins_x,bins_y,ovr_pot) + alpha*c1*np.dot(grad_f(LG,bins_x,bins_y,ovr_pot), pik)):
        alpha = rho*alpha
        
    return alpha

def gsLS(a0, pik, LG, bins_x, bins_y, ovr_pot):
    """
    Golden section line search
    """
    a = 0.0#Initial values for endpoints of search brackets
    prevfa = f(LG,bins_x,bins_y,ovr_pot)
    b = a0
    prevfb = f(LG,bins_x,bins_y,ovr_pot,pk=b*pik)
    rho = (3.0 - np.sqrt(5)) / 2.0
    newa = a + (rho*(b-a))
    newfa = f(LG,bins_x,bins_y,ovr_pot,pk=newa*pik)
    newb = b - (rho*(b-a))
    newfb = f(LG,bins_x,bins_y,ovr_pot,pk=newb*pik)

    while (b - a) > 1.0e-6 and (np.abs(newfb - newfa) > 1.0e-5):
        if newfa < newfb: #Next search interval is [a, newb]
            b = newb
            newb = newa
            newfb = newfa
            newa = a + (rho*(b-a))
            newfa = f(LG,bins_x,bins_y,ovr_pot,pk=newa*pik)
        else: #Next search interval is [newa, b]
            a = newa
            newa = newb
            newfa = newfb
            newb = b - (rho*(b-a))
            newfb = f(LG,bins_x,bins_y,ovr_pot,pk=newb*pik)
    
    return (0.5*(b+a))


def lineSearch(a0, pik, LG, bins_x, bins_y, ovr_pot):
    """
    Auxiliary function to perform line search using low order polynomial interpolation 
    to obtain a step length alpha_k that leads to convergence, i.e. satisfies the Wolfe 
    Conditions

    Parameters:
        a_0: Initial guess for step length (1.0)
        pik: Descent direction
        LG: LevelGraph object
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pot: Array of total potential movable area for each cell over all bins
    Return:
        alpha: Step length that satisfies Wolfe Conditions
    """
    c1 = 1e-4#Nocedal and Wright specify for this to be small
    alpha = a0
    phi_0 = f(LG, bins_x,bins_y,ovr_pot)
    phi_prime_0 = np.dot(grad_f(LG,bins_x,bins_y,ovr_pot), pik)
    if (f(LG,bins_x,bins_y,ovr_pot,pk=alpha*pik) <= phi_0 + alpha*c1*phi_prime_0):
        return alpha
    else:
        #Start with quadratic
        alpha1 = (-phi_prime_0*(alpha**2)) / (2.0*(f(LG,bins_x,bins_y,ovr_pot,pk=alpha*pik) - phi_0 - alpha*phi_prime_0))
        if (f(LG,bins_x,bins_y,ovr_pot,pk=alpha1*pik) <= phi_0 + alpha1*c1*phi_prime_0):
            alpha = alpha1
            return alpha
        else:
            #Quadratic interpolation does not give precise enough - go cubic with Hermite Interpolation as
            #described in Section 3.5
            alpha0 = alpha
            while (f(LG,bins_x,bins_y,ovr_pot,pk=alpha1*pik) > phi_0 + alpha1*c1*phi_prime_0):
                v1 = f(LG,bins_x,bins_y,ovr_pot,pk=alpha1*pik) - phi_0 - alpha1*phi_prime_0
                v2 = f(LG,bins_x,bins_y,ovr_pot,pk=alpha0*pik) - phi_0 - alpha0*phi_prime_0
                a = ((alpha0**2)*v1 - (alpha1**2)*v2) / ((alpha0**2)*(alpha1**2)*(alpha1 - alpha0))
                b = (-1.0*(alpha0**3)*v1 + (alpha1**3)*v2) / ((alpha0**2)*(alpha1**2)*(alpha1 - alpha0))
                alpha2 = (-1.0*b + np.sqrt(b**2 - (3.0*a*phi_prime_0))) / (3.0*a)
                alpha0 = alpha1
                alpha1 = alpha2
                #If successive alphas are too close together or alpha 1 << alpha0, set alpha1 = alpha0/2
                if (np.abs(alpha1 - alpha0) < 1.0e-6):
                    alpha1 = 0.5*alpha0
            alpha = alpha1
        return alpha

def bfgs(LG,bins_x, bins_y, ovr_pot, H0, lambda_m, eps=1.0e-3):
    """

    Function 
    Parameters:
        LG: LevelGraph object
        bin_x: Numpy meshgrid of bin top-left x-coordinates, with indexing='ij'
        bin_y: Numpy meshgrid of bin top-left y-coordinates, with indexing='ij'
        ovr_pots: Array of total potential movable area for each cell over all bins
        H0: Initial guess for approximate Hessian
        lambda_m: Value of scaling factor for crowdedness term involving Db
        eps: Value of epsilon for stopping condition

    Return:
        x_new: Minimized trajectory x(t)
        f_vals: Values of f(x(t)) at each iteration
        grad_norms: values of del(f(x(t))) at each iteration
    """
    #Re-use x_prev is x_k in BFGS iteration, x_new is x_k+1, Hk is approximate inverse of del^2 f(x)
    k = 0
    Hk = H0
    Hk1 = H0
    f_vals = []
    f_vals.append(f(LG,bins_x,bins_y,ovr_pot))
    grad_norms = []
    grad_norms.append(np.linalg.norm(grad_f(LG,bins_x,bins_y,ovr_pot)))
    while(np.linalg.norm(grad_f(LG,bins_x,bins_y,ovr_pot)) >= eps):
        grad_f_p = grad_f(LG,bins_x,bins_y,ovr_pot)
        print("Iteration " + str(k) + ":        f(x_k): " + str(f(LG,bins_x,bins_y,ovr_pot)) + "        ||grad_f(x_k)||_2: " + str(np.linalg.norm(grad_f_p)))
        pk = -Hk@grad_f_p#Descent direction
        ak = gsLS(1.0, pk, LG,bins_x,bins_y,ovr_pot)#Armijo backtracking routine to find alpha_k that will lead to convergence
        print("step length ak: " + str(ak))
        x_prev = LG.vvec()
        x_new = x_prev + ak*pk
        LG.updatePositions(x_new)
        sk = x_new - x_prev
        grad_f_n = grad_f(LG,bins_x,bins_y,ovr_pot)
        yk = grad_f_n - grad_f_p
        rho_k = 1.0/np.dot(yk,sk)
        #Expand BFGS update formula 6.17 in Ch.6 of Nocedal and write to calculate H_k+1 using only matrix addition,
        #matrix-vector products and outer products of vectors so total time complexity of update is O(n^2)
        uk = Hk@yk
        vk = (Hk.T)@yk
        Hk1 = Hk - rho_k*np.outer(sk,vk) - rho_k*np.outer(uk,sk) + ((rho_k**2)*np.dot(yk,uk))*np.outer(sk,sk) + rho_k*np.outer(sk,sk)
        f_vals.append(f(LG,bins_x,bins_y,ovr_pot))
        grad_norms.append(np.linalg.norm(grad_f(LG,bins_x,bins_y,ovr_pot)))
        Hk = Hk1
        k += 1
    return LG

def isCellInside(bin_tlx, bin_tly, bin_w, bin_h, cell):
    """
    Auxiliary function to determine if a cell is inside a bin
    """
    if (cell.cx >= bin_tlx) and (cell.cx < (bin_tlx + bin_w)) and (cell.cy >= bin_tly) and (cell.cy < (bin_tly + bin_h)):
        return True
    return False

def calcOverflowRatio(LG,bin_x, bin_y,ovr_pot, OVR_AREA):
    """
    Function to calculate the total overflow ratio at each level in the
    V-Cycle iteration

    Parameters:
        LG: LevelGraph object
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
    Mb = bin_w * bin_h
    Db = maxDb(LG,bin_x, bin_y,ovr_pot)
    max_overflow = np.max(Db - Mb)

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
    def_file = 'KSA16_yosys.def'
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
        lef_file = f'cell_lefs/{ctypes[cidx]}.lef'
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
        print()
        print()
        print("######################### Level " + str(i) + " #########################")
        print()
        print()
        grid_nw = int(np.sqrt(H_current.Nverts))#Make grid of bins square by default
        grid_nh = grid_nw
        bins_x, bins_y = np.meshgrid(np.linspace(0.0, OVR_W, grid_nw, dtype=np.double), np.linspace(0.0, OVR_H, grid_nh, dtype=np.double), indexing='ij')
        bin_area = (bins_x[1][0] - bins_x[0][0])*(bins_y[0][1] - bins_y[0][0])#Mb
        #No need for base potentials as we have no pre-placed blocks
        ovr_pots = calcOvrPotential(H_current, bins_x, bins_y)
        lambda_m = np.linalg.norm(grad_W(H_current), 1) / np.linalg.norm(grad_fDb(H_current, bins_x, bins_y, ovr_pots), 1)#Initialize lambda to be 1-norm of gradient

        prev_overflow_ratio = calcOverflowRatio(H_current,bins_x, bins_y,ovr_pots, OVR_W*OVR_H)
        new_overflow_ratio = 100.0
        m = 0

        #do-while optimization loop
        while(True):
            Nx = H_current.Nverts * 2
            y0 = grad_f(H_current, bins_x, bins_y,ovr_pots)
            x0 = H_current.vvec()
            H_guess = (np.dot(y0,x0)/np.dot(y0,y0)) * np.eye(Nx)
            H_current = bfgs(H_current, bins_x, bins_y, ovr_pots, H_guess, lambda_m)
            #Update this level's cluster/cell/vertex coordinates using x_new
            m += 1
            lambda_m *= 2.0
            #No need for look-ahead LG
            new_overflow_ratio = calcOverflowRatio(H_current,bins_x,bins_y,ovr_pots,OVR_W*OVR_H)
            if (new_overflow_ratio - prev_overflow_ratio >= -0.01):#If reduction in overflow ratio is >= 1%, keep going with bfgs
                continue
            else:#If reduction in overflow ratio is < 1%, stop with bfgs
                break
        
        #Post-processing after optimization loop to spread cells
        part = Partition(bin_area, H_current.verts)
        part.construct(0.0, 0.0, OVR_W, OVR_H)
        ws = part.calcWhiteSpace(part.topNode)
        part.WSA(part.topNode)
        #part.printLeaves(part.topNode)
        part.updateCells(part.topNode)

        #De-cluster and update cell positions after WSA()
        H_current.deCluster(part.vvec())

    H_0 = H_current#Return original hypergraph at finest/least clustered level with final GP cell positions
    return H_0


# Testing GP
if (__name__ == '__main__'):
    verilog_netlist = 'KSA16_yosys.vg'
    hgr_file = 'KSA16_yosys.hgr'
    global master_cell_array
    master_cell_array, master_hg = populateCells(verilog_netlist,hgr_file)
    H_0 = LevelGraph(master_hg,master_cell_array)

    #OVR_W and OVR_H taken from openROAD floorplan
    OVR_W = 89.395
    OVR_H = 89.395

    gamma = 0.01 * OVR_W

    # require number of vertices at coarsest level to be half the number of cells
    N_MAX = H_0.Nverts/2

    H_0 = gpMain(H_0,N_MAX,OVR_W,OVR_H)



