import numpy as np
from numpy.linalg import inv

class Fluid:
    def __init__(self, name, modulus, viscosity, density):
        self.name = name
        self.modulus = modulus
        self.viscosity = viscosity
        self.density = density
        
class MultiFluid:
    def __init__(self, fluid1, fluid2, saturation, patch, lam = 1, residualSw = 0, residualSnw = 0):
        def bulkmix(K1, K2, s1, q):
            return((q*(1 - s1) + s1)/((q*(1 - s1))/K2 + s1/K1))
        def BrooksCorey(l, sw, snw,s):
            def seff(s):
                return(np.clip((s-sw)/(1-sw-snw),0,1))
            def kw(s):
                return(np.clip(seff(s)**((2+3*l)/l),0,1))
            def knw(s):
                return(np.clip(((1-seff(s))**2*(1-seff(s)**((2+l)/l))),0,1))
            return([kw(s),knw(s)])
        def ViscosityMix(h1, h2, q, l, sw, snw,s):
            k = BrooksCorey(l, sw, snw,s)
            return((s + q*(1 - s))/(k[0]/h1 + q*(k[1]/h2)))
        def Density(rs, rf, f):
            return((1-f)*rs + f*rf)
        self.modulus = bulkmix(fluid1.modulus, fluid2.modulus, saturation, patch)
        self.relativePermeability = BrooksCorey(lam, residualSw, residualSnw, saturation)
        self.viscosity = ViscosityMix(fluid1.viscosity, fluid2.viscosity, patch, lam, residualSw, residualSnw, saturation)
        self.tau = fluid1.viscosity/self.viscosity
        self.density = Density(fluid1.density, fluid2.density, saturation)

class Rock:
    def __init__(self, name, grain_modulus, dry_bulk_modulus, dry_shear_modulus, porosity, mineral_density, crack_density = 0):
        self.name = name
        self.grain_modulus = grain_modulus
        self.dry_bulk_modulus = dry_bulk_modulus
        self.dry_shear_modulus = dry_shear_modulus
        self.porosity = porosity
        self.mineral_density = mineral_density
        self.crack_density = crack_density
        
        
class VTIRock:
    def __init__(self, name, frame_lambda, frame_mu, porosity, mineral_density, crack_density = 0, fracture_density = 0, ratio = 100):
        self.name = name
        self.lam = frame_lambda
        self.mu = frame_mu
        self.porosity = porosity
        self.mineral_density = mineral_density
        self.crack_density = crack_density
        self.fracture_density = fracture_density
        self.ratio = ratio

def GassmannModel(Kd, mu, Km, phi, Kf):
    return Kd + (1 - Kd/Km)**2/(phi/Kf - Kd/Km**2 + (1 - phi)/Km)

def SquirtModel(Kd, mu, Km, phi, Kf, ee, aspectRatio, tau, w, eta):
    fc = (4.*aspectRatio*np.pi*ee)/3.
    fp = phi - fc
    m = -(16*Km**2*fc +3*aspectRatio*Km*np.pi*(4*Kd + Km*(-4 + 5*fp)) + 
          np.sqrt(Km**2*(256*Km**2*fc**2 + 9*aspectRatio**2*np.pi**2*
                      (4*Kd - 4*Km + 3*Km*fp)**2 + 
                      96*aspectRatio*Km*np.pi*fc*(2*Kd - 2*Km + 3*Km*fp))))/(8.*aspectRatio*np.pi*(Kd + Km*(-1 + fp)))
    nu = (3*Km - 2*m)/(2.*(3*Km + m))
    sc =(aspectRatio*m*np.pi)/(2.*(1 - (3*Km - 2*m)/(2.*(3.*Km + m))))
    Kc = sc/Kf
    Kp = (4*m)/(3*Kf)
    gam = (3*aspectRatio*fp*(1 + Kp)*np.pi)/(8.*fc*(1 + Kc)*(1 - (3*Km - 2*m)/(2.*(3*Km + m))))
    gam2 = (3*Km*gam + 4*m*gam)/(9*Km + 9*Km*Kp)
    mhf =(4*fc*m*((2*(1 - (3*Km - 2*m)/(2.*(3*Km + m))))/(aspectRatio*np.pi) + (6 - 6*nu)/(2*aspectRatio*np.pi - aspectRatio*nu*np.pi) - 
                  (6*(-1 + nu)*(1j*m + eta*w))/(-2*eta*(-1 + nu)*w + aspectRatio*(-2 + nu)*np.pi*(1j*m + eta*w)) - 
                  (2*(1 - (3*Km - 2*m)/(2.*(3*Km + m)))* (Kc + 1/(1 + 1j*tau*w)))/(aspectRatio*(1 + Kc)*np.pi)))/15.
    Klf = Kd + \
        fp*(Km*((1 + (3*Km)/(4.*m)) + \
        fc*(1 + (2*Km*(1 - (3*Km - 2*m)/(2.*(3*Km + m))))/(aspectRatio*m*np.pi)))*(1 + 3*(1 + Kc)*gam2))/((1 + Kc)*(1 + gam))
    Khf = fp*((gam - 3*gam2*(1 + Kc))*Km*(-((1 + (3*Km)/(4.*m))) + \
        fc*gam*(1 + (2*Km*(1 - (3*Km - 2*m)/(2.*(3*Km + m))))/(aspectRatio*m*np.pi))))/(gam*(1 +\
        gam)*(1 + Kc)*(1 - (1j*(1 + gam))/(gam*tau*w)))
    return [Klf+Khf, mu+mhf]

def AnisSquirtModel(lam, m, f, Kf, ee, eef, aspectRatio, tau, w, ratio):
    fc = (4.*aspectRatio*np.pi*ee)/3.
    ff = (4.*aspectRatio*np.pi*eef)/3.
    fp = f - fc - ff
    nu = lam/(2.*(lam + m))
    sc =(aspectRatio*m*np.pi)/(2.*(1 - nu))
    Kc = sc/Kf
    Kp = (4*m)/(3*Kf)
    gam = (3*fp*(1 + Kp)*sc)/(4.*fc*(1 + Kc)*m)
    gam2 = (gam*(1 - nu))/((1 + Kp)*(1 + nu))
    io = fc/(fc + aspectRatio*fp)
    bet = io*ff/fc
    tauf = ratio*tau
    D1 = (io*(-1j - bet*tau*w + tauf*w) + 3*gam2*(1 + Kc)*(io*(1j + bet*tau*w - tauf*w) + \
        (-1j + tau*w)*(1 + 1j*tauf*w)))/(3.*(1 + Kc)*(bet*(-1j + (1 + (-1 + gam)*io)*tau*w) - \
        (-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    D2 =  (bet*(-1j + tau*w))/((1 + Kc)*(bet*(-1j + \
        (1 + (-1 + gam)*io)*tau*w) -(-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    G1 = (tau*w)/((1 + Kc)*(-1j + tau*w))
    G2 = (1j/3*(1j + bet*tau*w - tauf*w)*(io + 1j*gam*io*tau*w - \
        3*1j*gam2*(-1 + io)*(1 + Kc)*(-1j + tau*w)))/((1 + Kc)*\
        (-1j + tau*w)*(bet*(-1j + (1 + (-1 + gam)*io)*tau*w) - \
        (-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    G3 = (bet*(-1j + gam*tau*w))/((1 + Kc)*(bet*(-1j + (1 + (-1 + gam)*io)*tau*w) - \
        (-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    F1 = (-3*gam2*(-1 + io)*(1 + Kc)*(-1j + tau*w) + \
        io*(-1j + gam*tau*w))/(3.*(1 + Kc)*(bet*(-1j + (1 + (-1 + gam)*io)*tau*w) - \
        (-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    F2 = (tauf*w*(gam + io - gam*io + 1j*gam*tau*w) + \
        bet*(-1j + (1 + (-1 + gam)*io)*tau*w))/((1 + Kc)*(bet*(-1j + (1 + (-1 + gam)*io)*tau*w) - \
        (-1j + tauf*w)*(-io + gam*(-1 + io - 1j*tau*w))))
    k = lam + (2*m)/3.
    L2 = k**2 + (16*m**2)/45.
    L4 = k**2 - (8*m**2)/45.
    c11 = lam + 2*m - \
        fp*(-((3*D1*k + D2*lam)*(1 + (3*k)/(4.*m))) + \
            (3*(1 - nu)*(3*lam**2 + 4*lam*m + (m**2*(36 + 20*nu))/(7 - 5*nu)))/(4.*m*(1 + nu))) - \
        fc*((32*m*(1 - nu))/(15.*aspectRatio*(2 - nu)*np.pi) - \
            G2*(3*k + (3*k**2)/sc) - G1*(k + L2/sc) - G3*(lam + (k*lam)/sc) + L2/sc) - \
        ff*(-3*F1*k*(1 + lam/sc) - F2*lam*(1 + lam/sc) + lam**2/sc)
    c33 = lam + 2*m - \
        fp*(-((1 + (3*k)/(4.*m))*(3*D1*k + D2*(lam + 2*m))) + \
            (3*(1 - nu)*(3*lam**2 + 4*lam*m + (m**2*(36 + 20*nu))/(7 - 5*nu)))/(4.*m*(1 + nu))) - \
        fc*((32*m*(1 - nu))/(15.*aspectRatio*(2 - nu)*np.pi) - \
            G2*(3*k + (3*k**2)/sc) - G1*(k + L2/sc) - G3*(lam + 2*m + (k*(lam + 2*m))/sc) + L2/sc) - \
        ff*(-3*F1*k*(1 + (lam + 2*m)/sc) -  F2*(lam + 2*m)*(1 + (lam + 2*m)/sc) + (lam + 2*m)**2/sc)
    c44 = m - \
        fp*(15*m*(1 - nu))/(7 - 5*nu) - \
        fc*((8*m*(1 - nu))/(5.*aspectRatio*(2 - nu)*np.pi) + (4*(1 - G1)*m**2)/(15.*sc)) - \
        ff*(4*m*(1 - nu))/(aspectRatio*(2 - nu)*np.pi)
    c12 = lam - \
        fp*(-((3*D1*k + D2*lam)*(1 + (3*k)/(4.*m))) + \
            (3*(1 - nu)*(3*lam**2 + 4*lam*m - (4*m**2*(1 + 5*nu))/(7 - 5*nu)))/(4.*m*(1 + nu))) - \
        fc*((-16*m*(1 - nu))/(15.*aspectRatio*(2 - nu)*np.pi) - \
            G2*(3*k + (3*k**2)/sc) - G1*(k + L4/sc) - G3*(lam + (k*lam)/sc) + L4/sc) - \
        ff*(-3*F1*k*(1 + lam/sc) - F2*lam*(1 + lam/sc) + lam**2/sc)
    c13 = lam - \
        fp*(-((1 + (3*k)/(4.*m))*(3*D1*k + D2*(lam + m))) + \
            (3*(1 - nu)*(3*lam**2 + 4*lam*m - (4*m**2*(1 + 5*nu))/(7 - 5*nu)))/(4.*m*(1 + nu))) - \
        fc*((-16*m*(1 - nu))/(15.*aspectRatio*(2 - nu)*np.pi) - \
            G2*(3*k + (3*k**2)/sc) - G1*(k + L4/sc) - G3*(lam + m + (k*(lam + m))/sc) + L4/sc) - \
        ff*(-3*F1*k*(1 + (lam + m)/sc) - F2*(lam + m + (lam*(lam + 2*m))/sc) + (lam*(lam + 2*m))/sc)
    return [c11, c13, c33, c44, c12]

def AnisHTIVelocities(c11, c33, c44, c12, c13, rho, theta):
    H = 4*(c13 + c44)**2*np.cos(theta)**2*np.sin(theta)**2 + (-((c33 - c44)*np.cos(theta)**2) + (c11 - c44)*np.sin(theta)**2)**2
    Vp = np.sqrt((c11*np.sin(theta)**2 + c33*np.cos(theta)**2 + c44 + np.sqrt(H))/(2*rho))
    Vsv = np.sqrt((c11*np.sin(theta)**2 + c33*np.cos(theta)**2 + c44 - np.sqrt(H))/(2*rho))
    Vsh = np.sqrt((1/2*(c11 - c12)*np.sin(theta)**2 + c44*np.cos(theta)**2)/rho)
    return [Vp, Vsv, Vsh]

