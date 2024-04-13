import numpy as np

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
    return x_new, f_vals, grad_norms