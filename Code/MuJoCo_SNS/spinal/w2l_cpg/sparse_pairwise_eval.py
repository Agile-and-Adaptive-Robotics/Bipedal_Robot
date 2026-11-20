"""Sparse SNS summation retaining NumPy's contiguous pairwise reduction tree.

No fastmath or reassociation. Empty columns are pruned within the same
eight-lane leaf accumulators, then the original binary tree is evaluated.
Only suitable for finite ordinary conductances, as checked by the caller.
"""
import numpy as np
from numba import njit
from preallocated_exact_eval import PreallocatedExactForward

def reduction_plan(n):
    leaves=[];joins=[]
    def branch(start,count):
        if count<=128:
            index=len(leaves);leaves.append((start,count));return index
        left=count//2;left-=left%8
        a=branch(start,left);b=branch(start+left,count-left)
        index=100000+len(joins);joins.append((a,b));return index
    root=branch(0,n)
    # Leaf IDs precede internal nodes; construction is postorder.
    count=len(leaves)
    remap=lambda i:i if i<100000 else count+i-100000
    return np.array(leaves,dtype=np.int64),np.array([(remap(a),remap(b)) for a,b in joins],dtype=np.int64).reshape(-1,2),remap(root)

@njit(cache=True)
def sparse_sums(v,post,pre,gmax,elo,erange,reversal,leaf,lane,tail,counts,joins,root,acc,values,out):
    acc[:]=0.;values[:]=0.
    for e in range(len(post)):
        p=post[e];f=(v[pre[e]]-elo[e])/erange[e]
        if f<0.:f=0.
        elif f>1.:f=1.
        g=gmax[e]*f;ge=g*reversal[e];b=leaf[e];k=lane[e]
        if tail[e]:
            # Tail is folded after the eight main lanes, not into a lane.
            acc[p,b,8+k,0]=g;acc[p,b,8+k,1]=ge
        else:
            acc[p,b,k,0]+=g;acc[p,b,k,1]+=ge
    for p in range(len(v)):
        for b in range(len(counts)):
            for kind in range(2):
                r=acc[p,b,:,kind]
                value=((r[0]+r[1])+(r[2]+r[3]))+((r[4]+r[5])+(r[6]+r[7]))
                for k in range(counts[b]%8):value+=r[8+k]
                values[p,b,kind]=value
        offset=len(counts)
        for j in range(len(joins)):
            for kind in range(2):values[p,offset+j,kind]=values[p,joins[j,0],kind]+values[p,joins[j,1],kind]
        out[p,0]=values[p,root,0];out[p,1]=values[p,root,1]

class SparsePairwiseForward(PreallocatedExactForward):
    def __init__(self,compiled):
        super().__init__(compiled)
        leaves,self.joins,self.root=reduction_plan(len(self.c.V));self.counts=leaves[:,1]
        self.leaf=np.zeros(len(self.pre),dtype=np.int64);self.lane=self.leaf.copy();self.tail=np.zeros(len(self.pre),dtype=np.bool_)
        for e,col in enumerate(self.pre):
            b=next(i for i,(start,count) in enumerate(leaves) if start<=col<start+count)
            local=col-leaves[b,0];body=self.counts[b]-self.counts[b]%8
            self.leaf[e]=b;self.tail[e]=local>=body;self.lane[e]=local-body if self.tail[e] else local%8
        n=len(self.c.V);self.acc=np.zeros((n,len(leaves),16,2));self.values=np.zeros((n,len(leaves)+len(self.joins),2));self.sums=np.zeros((n,2))
        assert np.all(np.isfinite(self.gmax)) and np.all(np.isfinite(self.reversal))

    def forward(self,x=None):
        c=self.c
        np.copyto(self.last_voltage,c.V);c.V_last=self.last_voltage
        if x is None:i_app=0
        elif self.input_onehot:
            np.take(x,self.iq,out=self.input_values);self.i_app[self.ip]=self.input_values;i_app=self.i_app
        else:i_app=np.matmul(c.input_connectivity,x)
        sparse_sums(c.V_last,self.post,self.pre,self.gmax,self.elo,self.erange,self.reversal,self.leaf,self.lane,self.tail,
                    self.counts,self.joins,self.root,self.acc,self.values,self.sums)
        i_syn=self.sums[:,1]-c.V_last*self.sums[:,0]
        if c.gated:
            a_inf=1/(1+c.k_a*np.exp(c.slope_a*(c.e_a-c.V_last)))
            b_inf=1/(1+c.k_b*np.exp(c.slope_b*(c.e_b-c.V_last)))
            c_inf=1/(1+c.k_c*np.exp(c.slope_c*(c.e_c-c.V_last)))
            tau_c=c.tau_max_c*c_inf*np.sqrt(c.k_c*np.exp(c.slope_c*(c.e_c-c.V_last)))
            c.b_gate_last,c.c_gate_last=np.copy(c.b_gate),np.copy(c.c_gate)
            c.b_gate=c.b_gate_last+c.dt*((b_inf-c.b_gate_last)/c.tau_max_b)
            c.c_gate=c.c_gate_last+c.dt*((c_inf-c.c_gate_last)/tau_c)
            i_ion=c.g_ion*(a_inf**c.pow_a)*(c.b_gate**c.pow_b)*(c.c_gate**c.pow_c)*(c.e_ion-c.V_last)
            i_gated=np.sum(i_ion,axis=0)
            c.V=c.V_last+c.time_factor_membrane*(-c.g_m*(c.V_last-c.V_rest)+c.i_b+i_syn+i_app+i_gated)
        else:c.V=c.V_last+c.time_factor_membrane*(-c.g_m*(c.V_last-c.V_rest)+c.i_b+i_syn+i_app)
        if self.output_onehot:np.take(c.V,self.oq,out=self.output_values);c.outputs=self.output_values
        else:c.outputs=np.matmul(c.output_voltage_connectivity,c.V)
        return c.outputs
