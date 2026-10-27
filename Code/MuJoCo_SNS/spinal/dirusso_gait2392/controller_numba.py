"""Preallocated compiled counterpart of controller.py for optimization.

Matches its corrected discrete equations; paper-to-model adaptations and
unresolved delay/balance details remain those of controller.py.
"""
import numpy as np
from numba import njit
import controller as C

GAIN_KEYS = ('ia_aff', 'ia_mono', 'ia_iain', 'ia_recip', 'ia_iain_mut',
             'ii_aff', 'ii_exc', 'ii_inh', 'ib_aff', 'ib_auto', 'ib_rev',
             'ibin_mut', 'rc_drive', 'rc_mn', 'rc_iain', 'rc_mut', 'cut_exc')


@njit(cache=True)
def _clip(x, lo=0., hi=1.):
    return min(max(x, lo), hi)


@njit(cache=True)
def _step(state, scratch, mnstate, mn, phi, W, g, sens, strikes,
          pitch, previous_pitch, kp, kd, target, omega, gamma, dt,
          ff, cut_aff, ant, ext):
    old0, old1 = phi[0], phi[1]
    phi[0] = (old0 + dt * (omega - gamma * np.sin(
        2*np.pi*(old0-old1)-np.pi)) / (2*np.pi)) % 1.
    phi[1] = (old1 + dt * (omega - gamma * np.sin(
        2*np.pi*(old1-old0)-np.pi)) / (2*np.pi)) % 1.
    for side in range(2):
        if strikes[side]:
            phi[side] = 0.
    derivative = 0. if np.isnan(previous_pitch[0]) else (pitch-previous_pitch[0])/dt
    previous_pitch[0] = pitch
    balance = -kp*(pitch-target)-kd*derivative
    k = dt/C.TAU
    for side in range(2):
        for m in range(9):
            ia, ii, ib = state[side,0,m], state[side,1,m], state[side,2,m]
            iain, rc, ibin = state[side,3,m], state[side,4,m], state[side,5,m]
            scratch[side,0,m] = ia+k*(g[0,m]*sens[side,0,m]-ia)
            scratch[side,1,m] = ii+k*(g[5,m]*sens[side,1,m]-ii)
            scratch[side,2,m] = ib+k*(g[8,m]*sens[side,2,m]-ib)
            drive_ia = g[2,m]*max(ia,0.)-g[14,m]*_clip(rc)
            drive_rc = g[12,m]*mn[side,m]
            drive_ib = max(ib,0.)
            for a in range(9):
                if ant[m,a]:
                    drive_ia -= g[4,a]*_clip(state[side,3,a])
                    drive_rc -= g[15,a]*_clip(state[side,4,a])
                    drive_ib -= g[11,a]*_clip(state[side,5,a])
            scratch[side,3,m] = iain+k*(drive_ia-iain)
            scratch[side,4,m] = rc+k*(drive_rc-rc)
            scratch[side,5,m] = ibin+k*(drive_ib-ibin)
            for channel, source in ((6,ib),(7,ii),(8,ii)):
                scratch[side,channel,m] = state[side,channel,m]+k*(
                    max(source,0.)-state[side,channel,m])
        for m in range(9):
            drive = 0.
            for p in range(5):
                dist = ((phi[side]-(0.1+0.2*p)+0.5)%1.)-0.5
                primitive = 0. if abs(dist)>0.2 else 0.5*(1.+np.cos(np.pi*dist/0.2))
                drive += W[side*9+m,p]*primitive
            drive *= ff
            drive += g[1,m]*max(scratch[side,0,m],0.) + g[6,m]*_clip(scratch[side,7,m])
            drive -= g[9,m]*_clip(scratch[side,5,m]) + g[13,m]*_clip(scratch[side,4,m])
            drive += g[16,m]*cut_aff*sens[side,3,0]
            if ext[m]:
                drive += g[10,m]*_clip(scratch[side,6,m])
            for a in range(9):
                if ant[m,a]:
                    drive -= g[3,a]*_clip(scratch[side,3,a]) + g[7,a]*_clip(scratch[side,8,a])
            mnstate[side,m] += (dt/C.TAU_A)*(drive-mnstate[side,m])
            value = _clip(mnstate[side,m])
            if m == 0 or m == 1 or m == 2:
                value += 0.1*balance[m]*(1. if sens[side,3,0]>0.2 else 0.)
            mn[side,m] = _clip(value)
        for channel in range(9):
            lo = 0. if channel == 4 or channel == 6 or channel == 7 else -50.
            for m in range(9):
                state[side,channel,m] = _clip(scratch[side,channel,m],lo,50.)


class CompiledDiRusso:
    def __init__(self, W, dt=.002, omega=5.11, gamma=1., gains=None,
                 balance_gains=None):
        defaults = C.default_gains()
        defaults.update(gains or {})
        self.g = np.stack([np.broadcast_to(defaults[k],(9,)) for k in GAIN_KEYS])
        self.ff, self.cut_aff = defaults['ff'], defaults['cut_aff']
        self.W = np.ascontiguousarray(W,dtype=float)
        self.state, self.scratch = np.zeros((2,9,9)), np.zeros((2,9,9))
        self.mnstate, self.mn = np.zeros((2,9)), np.zeros((2,9))
        self.phi = np.array([0., .5])
        self.previous_pitch = np.array([np.nan])
        self.dt, self.omega, self.gamma = dt, omega, gamma
        balances = balance_gains or dict(kp=8.,kd=1.5,target=0.)
        self.bal = {key:np.broadcast_to(value,(3,)).copy()
                    for key,value in balances.items()}
        self.ant = np.zeros((9,9),dtype=np.bool_)
        for m, names in enumerate(C.ANT_I):
            self.ant[m,names] = True
        self.ext = np.array([name in C.EXTENSOR for name in C.MUSCLES])

    def step_array(self, sens, strikes, pitch):
        _step(self.state, self.scratch, self.mnstate, self.mn, self.phi,
              self.W, self.g, sens, strikes, pitch, self.previous_pitch,
              self.bal['kp'],self.bal['kd'],self.bal['target'],self.omega,
              self.gamma,self.dt,self.ff,self.cut_aff,self.ant,self.ext)
        return self.mn
