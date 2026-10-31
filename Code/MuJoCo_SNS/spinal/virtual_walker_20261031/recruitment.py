"""Explicit opt-in mirrored synaptic calibration; no MN-current allocator."""
def apply_recruitment(controller,scale=1.,toe_mn_gmax=None,hip_scale=1.,hip_flexor_scale=1.):
 assert scale>0 and hip_scale>0 and hip_flexor_scale>0;edits=[];toe=[];hips=[];flexors=[]
 for k,e in enumerate(controller.edges):
  pre=controller.names[int(e[0])];post=controller.names[int(e[1])]
  if ' PF ' in pre and not pre.endswith(' IN') and ' MN ' in post and not post.endswith(' RE'):
   edits.append(dict(source=pre,target=post,source_gmax_uS=float(e[2]),effective_gmax_uS=float(e[2]*scale)));controller.edges[k,2]*=scale
   if ' Hip PF ' in pre and ' Hip MN ' in post:
    controller.edges[k,2]*=hip_scale
    if pre.endswith(' Hip PF flx') and post.endswith(' Hip MN flx'):controller.edges[k,2]*=hip_flexor_scale;flexors.append(dict(source=pre,target=post,extra_flexor_scale=hip_flexor_scale,effective_gmax_uS=float(controller.edges[k,2])))
    edits[-1]['effective_gmax_uS']=float(controller.edges[k,2]);hips.append(dict(source=pre,target=post,additional_scale=hip_scale,extra_flexor_scale=hip_flexor_scale if post.endswith(' flx') else 1.,effective_gmax_uS=float(controller.edges[k,2])))
  if pre in ('L toe contact','R toe contact') and ' Ank MN ' in post:
   if toe_mn_gmax is not None:
    assert toe_mn_gmax>=0;controller.edges[k,2]=toe_mn_gmax
    row=next(x for x in controller.contact_edge_edits if x['source']==pre and x['target']==post);row['effective_gmax_uS']=float(toe_mn_gmax);row['absolute_parameter_group']='MN_toe_ankle'
   toe.append(dict(source=pre,target=post,effective_gmax_uS=float(controller.edges[k,2])))
 if scale!=1.:assert len(edits)==12
 if hip_scale!=1.:assert len(hips)==4
 if hip_flexor_scale!=1.:assert len(flexors)==2
 if toe_mn_gmax is not None:assert len(toe)==2 and toe[0]['effective_gmax_uS']==toe[1]['effective_gmax_uS']
 for row in edits:
  other_source=('R' if row['source'].startswith('L') else 'L')+row['source'][1:];other_target=('R' if row['target'].startswith('L') else 'L')+row['target'][1:];other=next((x for x in edits if x['source']==other_source and x['target']==other_target),None)
  if other is not None:assert row['effective_gmax_uS']==other['effective_gmax_uS']
 return dict(pf_mn_recruitment_scale=scale,hip_pf_mn_scale=hip_scale,hip_flexor_pf_mn_scale=hip_flexor_scale,hip_pf_mn_edges=hips,hip_flexor_edges=flexors,pf_mn_edges=edits,toe_mn_edges=toe,scope='Actual conductance edits, bilateral absolute equality asserted; no neuron bias or force allocation.')
