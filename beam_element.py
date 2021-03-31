from numpy import array, arctan2, zeros, ix_
from scipy.linalg import norm 

def beam_element(xy , properties):
    
    E = properties['E']
    I = properties['I']
    A = properties['A']
    qx = properties['qx']
    qy = properties['qy']

    xi = xy[0,:]
    xj = xy[1,:]
    
    L  = norm(xj -xi)
    #θ = arctan2(xj[1] - xi[1], xj[0] - xi[0])
    
    cosθ = (xj[0] -xi[0])/L
    sinθ = (xj[1] -xi[1])/L
    print (f'sinθ ={sinθ}')

    ke = zeros((6,6))
    fe = zeros((6,1))
    ke_tilde = zeros((6,6))

    ke_tilde[0,0] = A*E/ L 
    ke_tilde[3,3] = A*E/ L 
    ke_tilde[0,3] = -A*E/ L 
    ke_tilde[3,0] = -A*E/ L 
    
    bending_dofs = ix_([1,2,4,5],[1,2,4,5])
    
    ke_tilde[bending_dofs] = E*I*array(
         [[12/L**3, 6/L**2, -12/L**3, 6/L**2],
         [6/L**2, 4/L, -6/L**2, 2/L],
         [-12/L**3, -6/L**2, 12/L**3, -6/L**2],
         [6/L**2, 2/L, -6/L**2, 4/L]])
    
    T = zeros((6,6))
    
    T[0:2,0:2] = array([[cosθ, - sinθ],  [sinθ, cosθ]])
    T[3:5,3:5] = array([[cosθ, - sinθ],  [sinθ, cosθ]])
    T[2,2] = 1.0
    T[5,5] = 1.0
    
    # Aqui debo computar el qx y el qy
    
    
    fe = zeros((6, 1))
    fe_local = zeros((6, 1))
    
    fe_local[0] =  qx*L/2
    fe_local[1] =  qy*L/2
    fe_local[2] =  qy*L**2/12
    fe_local[3] =  qx*L/2
    fe_local[4] =  qy*L/2
    fe_local[5] = -qy*L**2/12

    fe =T@fe_local

    # print (T)
    
    ke = T @ ke_tilde @ T.T
    # print (ke_tilde)
    # print (ke)

    return ke ,fe


# xy = array([[0,0],
#             [1,1]])
    
# properties ={}
# properties['E'] = 100
# properties['I'] = 1
# properties['A'] = 1


# beam_element(xy, properties)

