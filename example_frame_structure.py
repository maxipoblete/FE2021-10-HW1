from numpy import array, pi, zeros, ix_
from beam_element import beam_element
import matplotlib.pylab as plt


#Vamos a definir nuestra estructura

#------------NODOS----------
xy = array([ 
    [0.,0.],    #0
    [0.,3.],    #1
    [0.,5.],    #2
    [0.,6.],    #3
    [6.,6.5],   #4
    [6.,5.],    #5
    [6.,3.],    #6
    [6.,0.],    #7
    ])

#------------CONECCIONES- ELEMENTOS]----------
conec = array([
    [0,1],      #0
    [1,2],      #1
    [2,3],      #2
    [7,6],      #3
    [6,5],      #4
    [5,4],      #5
    [1,6],      #6
    [2,5],      #7
    [3,4],      #8
    ], dtype = int
    )

#---------------------------------
#------------PROPIEDADES----------
#---------------------------------

a = 30e-2       #30 cm
b = 20e-2       #20 cm
c = 40e-2       #40 cm


properties_0 ={}  #Propiedades del elementos Azul
properties_1 ={}  #Propiedades del elementos Verde
properties_2 ={}  #Propiedades del elementos Celeste

# qy: Para arriba positivo
# qx: Para derecha positivo

#AZUL
properties_0['E']   = 200e9 
properties_0['I']   = a*a**3/12
properties_0['A']   = a*a
properties_0['qx']  = -2206.50
properties_0['qy']  =  0            #N Carga distribuida

#VERDE
properties_1['E']   = 200e9 
properties_1['I']   = b*c**3/12
properties_1['A']   = c*b
properties_1['qx']  =  0
properties_1['qy']  = -1961.30      #N Carga distribuida

#CELESTE
properties_2['E']   =  200e9 
properties_2['I']   =  b*b**3/12
properties_2['A']   =  b*b
properties_2['qx']  = -81.44
properties_2['qy']  = -977.31       #N Carga distribuida




properties = [properties_0, properties_0, properties_0, 
              properties_0, properties_0, properties_0, 
              properties_1, properties_1, properties_2]

#------------PARAMETROS----------

Nnodes =xy.shape[0]
Nelems = conec.shape[0]

NDOFs_per_node = 3
NDOFs = NDOFs_per_node * Nnodes
K = zeros((NDOFs, NDOFs))
f = zeros((NDOFs, 1))



for e in range(Nelems):
    
    #Para el elemento e 
    ni = conec[e,0]                         #Nodo i del elemento e
    nj = conec[e,1]                         #Nodo j del elemento e

    #print(f'e = {e} ni = {ni} nj ={nj}')
    xy_e = xy[[ni, nj],:]                    #Posiciones de los nodos del elemento e 
    
    #print (f'xy_e = {xy_e}')
    
    ke, fe = beam_element(xy_e, properties[e])
    
    #print (f'ke = {ke}')
    
    
    d = [3*ni, 3*ni+1,3*ni+2 , 3*nj, 3*nj+1,3*nj+2 ]
    
    #Metodo de rigidez directa
    for i in range(2*NDOFs_per_node):
        p =d[i]
        for j in range(2*NDOFs_per_node):
            q= d[j]
            K[p, q] += ke[i,j]
        f[p]+= fe[i]
    


free_DOFs =[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
constrained_DOFs = [0,1,2,21,22,23]

Kff = K[ix_ (free_DOFs, free_DOFs) ]
Kfc = K[ix_ (free_DOFs, constrained_DOFs) ]
Kcf = K[ix_ (constrained_DOFs,free_DOFs) ]
Kcc = K[ix_ (constrained_DOFs,constrained_DOFs) ]

ff = f[free_DOFs]
fc = f[constrained_DOFs]

from scipy.linalg import solve
u = zeros((NDOFs,1))

u[free_DOFs] = solve(Kff, ff)


R = Kcf @ u[free_DOFs] + Kcc@u[constrained_DOFs] - fc

print (f'u = {u}')
print (f'R = {R}')
print (f'f = {f}')



# dd = 0
# for i in f:
#     print (f[dd][0])
#     dd+= 1



