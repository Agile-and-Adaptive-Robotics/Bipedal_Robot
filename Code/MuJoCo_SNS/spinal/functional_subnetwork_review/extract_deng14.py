"""Extract literal wiring from ODE AST, without importing or running it."""
import ast
import json
from audit_saved_graphs import ROOT, OUT, digest, table

source=ROOT/'Code/MuJoCo_SNS/spinal/deng_cpg_ode.py'
tree=ast.parse(source.read_text())
env={}
edges=[]
def value(n):
    if isinstance(n,ast.Constant):return n.value
    if isinstance(n,ast.Name):return env[n.id]
    if isinstance(n,(ast.List,ast.Tuple,ast.Set)):return [value(x) for x in n.elts]
    if isinstance(n,ast.UnaryOp) and isinstance(n.op,ast.USub):return -value(n.operand)
    if isinstance(n,ast.JoinedStr):return ''.join(str(value(x.value)) if isinstance(x,ast.FormattedValue) else x.value for x in n.values)
    raise ValueError(ast.dump(n))
def visit(n):
    if isinstance(n,ast.Assign):
        try:v=value(n.value)
        except (ValueError,KeyError):return
        for target in n.targets:
            if isinstance(target,ast.Name):env[target.id]=v
            elif isinstance(target,ast.Tuple):
                for name,x in zip(target.elts,v):env[name.id]=x
    elif isinstance(n,ast.Expr) and isinstance(n.value,ast.Call) and isinstance(n.value.func,ast.Name) and n.value.func.id=='add':
        a,b,g,e,w=[value(x) for x in n.value.args]
        edges.append({'id':f'deng14:edge:{len(edges)}','source':a,'target':b,'gmax_uS':g,'Erev_mV':e,'threshold_mV':w[0],'saturation_mV':w[1],'source_line':n.lineno,'delay_ms':0,'kind':'graded'})
    elif isinstance(n,ast.For) and isinstance(n.target,ast.Name):
        for tag in value(n.iter):
            env[n.target.id]=tag
            for x in n.body:visit(x)
for statement in tree.body:
    if isinstance(statement,(ast.FunctionDef,ast.If)):continue
    visit(statement)
assert len(env['NAMES'])==14 and len(edges)==18
folder=OUT/'reduced_Deng14';folder.mkdir(exist_ok=True)
nodes=[]
for i,name in enumerate(env['NAMES']):
    nodes.append({'id':f'deng14:{i}','name':name,'type':'NaP half-center' if i in env['NA'] else 'graded leaky neuron','rest_mV':env['VREST_MN'] if i>=12 else env['VREST_HC'],'Cm_nF':env['CM'],'Gm_uS':env['GM'],'passive_tau_ms':env['CM']/env['GM'],'fixed_h_tau_ms':env['TAU_H'] if i in env['NA'] else None,'spike_threshold_reset':'not applicable','side':'single limb; no bilateral side assigned','source_sha256':digest(source)})
table(folder/'nodes.csv',nodes);table(folder/'edges.csv',edges)
result={'source_sha256':digest(source),'method':'AST literal extraction only; source not imported','neurons':14,'edges':18,'named_reference_call':"simulate('fixed', t_end=20000.0, dt=0.1, i_pulse=10.0, t_pulse=100.0, pulse_ms=1.0)",'time_units':'ms','physical_feedback_edges':0,'outputs':['MN_hip_ext','MN_hip_flx'],'PF_MN_ODE_uS':[2.565,3.632],'PF_MN_Simulink_uS':[2.6563,3.6325],'version_policy':'separate implementations; no merged parameter set','Table2':'graded model; NaP h is not GLIF adaptive threshold','simulation_executed':False}
(folder/'STATIC_VALIDATION.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
