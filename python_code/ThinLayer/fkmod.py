import numpy as np
from numpy.linalg import inv
import hankel
from hankel import HankelTransform
import matplotlib.pyplot as plt
import scipy.special as scs

# progress tracker
import time, sys
from IPython.display import clear_output

def update_progress(progress,msg=''):
    bar_length = 20
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
    if progress < 0:
        progress = 0
    if progress >= 1:
        progress = 1
            
    block = int(round(bar_length * progress))
    clear_output(wait = True)
    text = msg+"Progress: [{0}] {1:.1f}%".format( "#" * block + "-" * (bar_length - block), progress * 100)
    print(text)

    ## vertical slowness and transformation matrices
def QLmatrices(layers_dict, layer, p):
    # based on Ursin and Stovas, Geophysics (2003)
    [c11, c13, c33, c44, d] = layers_dict[layer]
    a0 = np.sqrt(c33/d)
    b0 = np.sqrt(c44/d)
    s0 = 1-c44/c33
    dD = (c13 - c33 + 2*c44)/c33
    eta = (c11*c33 - (c13 + 2*c44)**2)/(2*c33**2)
    R1 = 2*(1 - p**2*b0**2)*(dD + 2*p**2*a0**2*eta)**2
    R2 = s0 + 2*p**2*b0**2*dD - 2*p**2*a0**2*(1 - 2*p**2*b0**2)*eta
    R = R1/(R2 + np.sqrt(R2**2 + 2*p**2*b0**2*R1))
    Sa = 2*dD + 2*p**2*a0**2*eta + R
    Sb = 2*(1 - p**2*b0**2)*a0**2/b0**2*eta - R
    q1 = 1/a0**2 - p**2 - p**2*Sa
    q2 = 1/b0**2 - p**2 - p**2*Sb
    qa = np.sqrt(q1.real + 1j*abs(q1.imag))
    qb = np.sqrt(q2.real + 1j*abs(q2.imag))
    
    d2 = np.sqrt((s0 + dD)/(s0 + Sa))
    d3 = 2*b0**2*(s0 + 
     1/2*(Sa + dD))/(s0 + dD)
    d4 = np.sqrt((s0 - p**2*b0**2*(s0 + Sb))/((1 - p**2*b0**2*(1 + Sb))*(s0 + dD)))
    d5 = (s0 - 2*p**2*b0**2*(s0 + 1/2*(Sb + dD)))/(s0 + dD)
    d1 = 1/np.sqrt(p**2*d3 + d5)
    L1 = d1*np.array([[d2*np.sqrt(qa/d), 1/d4*p/np.sqrt(d*qb)],
                      [d3*d2*p*np.sqrt(d*qa), -d5/d4*np.sqrt(d/qb)]
                     ], dtype = np.complex_)
    L2 = d1*np.array([[d5/d2*np.sqrt(d/qa), d3*d4*p*np.sqrt(d*qb)],
                      [1/d2*p/np.sqrt(d*qa), -d4*np.sqrt(qb/d)]
                     ], dtype = np.complex_)
    return(np.array([np.diagflat([qa,qb]),L1,L2]))


## reflection and transmission coefficients for a single interface
def LayerReflectionTransmission(layers_dict, layer1, layer2, p):
    # based on Ursin and Stovas, Geophysics (2002)
    # no need to use slownesses. Get the L matrices only
    # and use equations B-1, B-2 to calculate T, R
    [L1t, L2t] = QLmatrices(layers_dict, layer1, p)[1:5]
    [L1b, L2b] = QLmatrices(layers_dict, layer2, p)[1:5]
    
    # C and D as defined in eq. B-2
    C = np.dot(L2t.T, L1b)
    D = np.dot(L1t.T, L2b)
    
    # Reflection and transmission matrices of a single interface
    trans = 2*inv(C + D)
    refl = np.dot((C-D),trans/2)
    
    return(np.array([trans,refl]))

## stack reflectivity
def Reflectivity(layers_dict, layers, thicknesses, p, w):
    # based on Ursin and Stovas, Geophysics (2002)
    # Recursion is eq. 20
    n_of_layers = len(layers)-1
    
    # set thickness of the upper half-space to 0
    thicknesses[0] = 0
    
    # QL matrices for layer stack
    QL = [QLmatrices(layers_dict,i,p) for i in layers]
    
    # Pairwise reflection transmission for layer stack
    RT = [LayerReflectionTransmission(layers_dict,layers[i],layers[i+1], p) for i in range(n_of_layers)]
    
    # Phase propagator for layer stack
    phase = [np.diagflat(np.exp(1j*w*i[0]*np.diag(i[1]))) for i in zip(thicknesses,[item[0] for item in QL])]
    
    # Recursive calculation of reflection response
    def rec(n):
        if n == n_of_layers-1:
            return(np.zeros(shape=(2,2), dtype=np.complex_))
        else:
            inv_mat = inv(np.identity(2, dtype=np.complex_) + np.dot(RT[n+1][1],rec(n+1)))
            inv_t = np.dot(np.transpose(RT[n+1][0]) , np.dot(rec(n+1) , np.dot(inv_mat , RT[n+1][0])))
            return(np.dot(np.dot(phase[n + 1],RT[n + 1][1] + inv_t),phase[n + 1]))
    return(rec(-1))

## response due to a vertical source
def ResponseVforce(model, p, w):

    layers_dict = model["layers"]
    layers = model["model_layers"]
    thicknesses = model["model_thickness"]
    zSource = model["zSource"]
    zReceiver = model["zReceiver"]

    # QL matrices of the upper half-space
    QLtop = QLmatrices(layers_dict,layers[0],p)
    
    # phase due to the source
    phaseSource = np.diagflat(np.exp(-1j*w*zSource*np.diag(QLtop[0])))
    phaseReceiver = np.diagflat(np.exp(1j*w*zReceiver*np.diag(QLtop[0])))
    
    # vertical force
    vforce = np.array([1/w,0])
    
    # source discontinuity vector
    Sigma = 1/np.sqrt(2) * np.dot(QLtop[1].T, vforce)
    
    # receiver-adjusted reflectivity
    refReceiver = np.dot(np.dot(phaseReceiver,Reflectivity(layers_dict, layers, thicknesses, p, w)), phaseReceiver)
    
    # wavefield vector
    u = np.dot(np.dot(refReceiver, phaseSource), Sigma )
    
    # displacement vector
    b = 1/np.sqrt(2) * w * np.array([1j * np.dot(QLtop[1], u)[0], np.dot(QLtop[2], u)[1]])
    
    return b
    

## critical horizontal slowness
def getpmax(layers_dict,layers):
    phor = [np.real(np.sqrt(layers_dict[x][4]/layers_dict[x][3])) for x in layers]
    return max(phor)

## discrete Hankel transfrom
def DHTR(f, nu, offsets, rmax, jzero):
    vec = jzero * np.pi / jzero[-1]
    scorrection = np.sin(vec)/vec
    dhtlist = 2/(rmax **2 * scs.jv(nu+1, jzero) **2)
    return [np.dot(scorrection * dhtlist * f, scs.jv(nu, jzero * r/rmax)) for r in offsets]

## Ricker wavelet
def ricker(w,wp,dt=0):
    return (2* w**2)/(np.sqrt(np.pi)* wp**3)* np.exp(-w**2/wp**2) * np.exp(-1j * dt* w)

## f-t trace transfrom
def ftot(tresp,omega,w0 = 25):
    src = ricker(omega,2*np.pi*w0)
    taxis = np.arange(2*len(omega)+1)/(omega[1]/2/np.pi)/(2*len(omega))
    traceSrc = tresp * src #* (1j* omega) ** 2
    tracef = np.concatenate(([0],traceSrc,np.conj(np.flip(traceSrc,0))))
    return taxis, np.real(np.fft.fft(tracef))


def FKmodeling(model, pI=0.01, msg=''):
    
    # default values if not supplied
    try:
        omega = model['omega']
    except:
        omega = 2*np.pi*np.arange(1,101)/2
        
    try:
        r = model['r']
    except:
        r = np.arange(0.05,2.1,0.05)
        
    try:
        nslow = model['nslow']
    except:
        nslow = 1000
    
    # get critical slowness (integration limit)
    pmax = 1.05*getpmax(model['layers'],model['model_layers'])

    # Zeroes of the Bessel function (components are 0: Z, 1: R )
    jzeros = np.array([scs.jn_zeros(0,nslow), scs.jn_zeros(1,nslow)])

    # functions to transform (u[0] - Z, u[1] - R)
    u = lambda k,w: ResponseVforce(model,k/w - pI*1j,w)

    # transfrom Z component
    tresp = []
    
    start = time.time()

    # for i, w in enumerate(omega):
    #     inp = [u(k,w)[0] for k in jzeros[0]*pmax*w/jzeros[0,-1]]
    #     rmax = jzeros[0,-1] /pmax /w
    #     tranf = DHTR(inp, 0, r, rmax, jzeros[0])
    #     tresp.append(tranf)
        
    #     update_progress(i / len(omega),msg)
        
    for w in omega:
        inp = [u(k,w)[0] for k in jzeros[0]*pmax*w/jzeros[0,-1]]
        rmax = jzeros[0,-1] /pmax /w
        tranf = DHTR(inp, 0, r, rmax, jzeros[0])
        tresp.append(tranf)

    tresp = np.array(tresp).T

    end = time.time()
    elapsed = end - start

    return tresp, omega, r, elapsed

def FKmodeling_par(DV, model, pI=0.01):
    
    # default values if not supplied
    try:
        omega = model['omega']
    except:
        omega = 2*np.pi*np.arange(1,101)/2
        
    try:
        r = model['r']
    except:
        r = np.arange(0.05,2.1,0.05)
        
    try:
        nslow = model['nslow']
    except:
        nslow = 1000
    
    # get critical slowness (integration limit)
    pmax = 1.05*getpmax(model['layers'],model['model_layers'])

    # Zeroes of the Bessel function (components are 0: Z, 1: R )
    jzeros = np.array([scs.jn_zeros(0,nslow), scs.jn_zeros(1,nslow)])

    # functions to transform (u[0] - Z, u[1] - R)
    u = lambda k,w: ResponseVforce(model,k/w - pI*1j,w)[0]

    # transfrom Z component
    DV.use_cloudpickle()

    def resp_calc(w):
        inp = [u(k,w) for k in jzeros[0]*pmax*w/jzeros[0,-1]]
        rmax = jzeros[0,-1] /pmax /w
        return DHTR(inp, 0, r, rmax, jzeros[0])
    
    start = time.time()

    tresp = DV.map_sync(resp_calc,omega)
    tresp = np.array(tresp).T
    
    end = time.time()
    elapsed = end - start

    return tresp, omega, r, elapsed