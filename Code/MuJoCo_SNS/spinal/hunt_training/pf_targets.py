"""Monoarticular command primitives, separate from RF/hamstring commands.

These are normalized engineering targets from a feasible simulation teacher,
not measured EMG. Positive native hip coordinates correspond to extension.
"""
import numpy as np
from scipy.ndimage import gaussian_filter1d

def mono_command_targets(command,acts,pairs):
    acts=list(acts);out={}
    for side in ('r','l'):
        groups=[('hip_flexion',(('glut_max1','glut_max2','glut_max3'),('iliacus','psoas'))),
                ('knee_angle',(('vas_med','vas_int','vas_lat'),('bifemsh',)))]
        for joint,families in groups:
            for cell,group in zip(pairs[joint+'_'+side],families):
                indices=[acts.index(name+'_'+side) for name in group]
                wave=gaussian_filter1d(command[:,indices].mean(axis=1),10.,mode='wrap')
                out[cell]=2.75*np.clip(wave/max(float(np.percentile(wave,95)),1e-6),0.,1.)
    return out
