"""Heel-loading transients armed by genuine preceding whole-foot unloading."""
import math
import numpy as np

class HeelStrikeEncoder:
    def __init__(self,dt=.002,min_swing_s=.08,decay_s=.15,capture_s=.05,require_whole=True):
        if min(dt,decay_s,capture_s)<=0 or min_swing_s<0:raise ValueError('Invalid contact-event times')
        self.dt=dt;self.min_swing=min_swing_s;self.capture=capture_s;self.require_whole=require_whole
        self.decay=math.exp(-dt/decay_s);self.loaded=np.zeros(2,dtype=bool);self.armed=np.zeros(2,dtype=bool)
        self.unloaded=np.zeros(2);self.contact_age=np.zeros(2);self.trace=np.zeros(2)
        self.last_onset=np.zeros(2,dtype=bool);self.count=np.zeros(2,dtype=int)
    def step(self,heel_bw,whole_bw):
        for side in (0,1):
            old=self.loaded[side];heel=float(heel_bw[side]);whole=float(whole_bw[side]) if self.require_whole else heel
            if heel>.05:self.loaded[side]=True
            elif heel<.025:self.loaded[side]=False
            if whole<.025:
                self.unloaded[side]+=self.dt;self.contact_age[side]=0.
                if self.unloaded[side]>=self.min_swing:self.armed[side]=True
            else:
                self.unloaded[side]=0.;self.contact_age[side]+=self.dt
                if self.contact_age[side]>self.capture:self.armed[side]=False
            onset=bool(self.loaded[side] and not old and self.armed[side])
            self.last_onset[side]=onset;self.trace[side]*=self.decay
            if onset:self.trace[side]=1.;self.armed[side]=False;self.count[side]+=1
        return self.trace
