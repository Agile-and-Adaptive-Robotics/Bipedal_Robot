"""Experimental spike-driven stance primitives, fitted offline, not allocated online."""
import json
import numpy as np
from numba import njit
from balance_pathway import BodyBalanceRelay

ROLES={'hip_flexor':'Hip MN flx','hip_extensor':'Hip MN ext',
       'knee_flexor':'Knee MN flx','knee_extensor':'Knee MN ext',
       'ankle_dorsiflexor':'Ank MN ext','ankle_plantarflexor':'Ank MN flx'}


@njit(cache=True)
def chemical_statistics(v,edges,indices,integral,peaks,dt):
    for k in range(len(indices)):
        e=edges[indices[k]];pre=int(e[0]);post=int(e[1]);a=min(1.,max(0.,(v[pre]-e[3])/(e[4]-e[3])));current=e[2]*a*(e[5]-v[post]);integral[k]+=current*dt;peaks[k]=max(peaks[k],abs(current))


@njit(cache=True)
def protect_transient(v,edges,indices,pf_bindings,gains,threshold):
    for k in range(len(indices)):
        pre,low,high=pf_bindings[k];activity=min(1.,max(0.,(v[int(pre)]-low)/(high-low)));quiet=min(1.,max(0.,(threshold-activity)/threshold));edges[indices[k],2]=gains[k]*quiet


class StanceRecruitmentRelay(BodyBalanceRelay):
    """Two new PF stance populations supply fixed, mirrored recruitment weights.

    Actual ipsilateral normal force and RG-E fraction gate the populations.
    Existing transient PFâ†’MN axons are unchanged. Neither desired tension nor
    inverse dynamics is evaluated during the coupled simulation.
    """
    learned=tuple(role+'_'+kind+'_uS' for role in ROLES for kind in ('exc','inh'))+('blend',)
    @classmethod
    def learned_keys(cls,configuration):
        return cls.learned+(('rg_ext_reference_mV',) if 'rg_ext_reference_mV' in configuration else ())
    def __init__(self,c,configuration):
        from spiking_sensory import encode_step
        self.c=c;self.cfg=configuration;assert configuration['pathway_kind']=='phase_stance_recruitment';self.encode=encode_step
        assert 0<=configuration['blend']<=1.02 and 0<configuration['off_load_BW']<configuration['full_load_BW']<=1
        assert 0<configuration['input_current_nA']<=4 and configuration['load_filter_tau_s']>=.01
        assert all(0<=configuration[k]<=15 for k in self.learned)
        if 'rg_ext_reference_mV' in configuration:assert 0<configuration['rg_ext_reference_mV']<=20
        self.indices=np.arange(len(c.v),len(c.v)+2);self.edge_specs=[]
        for key in ('v','newv','m','h','gs','gi','applied','tonic','last_test_impulse_nA_s'):
            old=getattr(c,key);setattr(c,key,np.r_[old,[-60.,-60.] if key in ('v','newv') else [0.,0.]])
        c.names.extend(['L stance-PF tonic relay','R stance-PF tonic relay']);c.params=np.vstack([c.params,np.tile([-60.,25.,1.,0.,0.,1.,1.,0.,1.,1.,0.],(2,1))]);edges=[]
        for pop,side in enumerate('LR'):
            for role,tail in ROLES.items():
                target=side+' '+tail;post=c.names.index(target);assert c.params[post,2]==1 and c.params[post,3]==0
                for kind,reversal in [('exc',-10.),('inh',-100.)]:
                    gain=configuration[role+'_'+kind+'_uS']*configuration['blend'];edges.append([self.indices[pop],post,gain,-60.,-40.,reversal]);self.edge_specs.append(dict(index=len(c.edges)+len(edges)-1,source=c.names[self.indices[pop]],target=target,role=role,kind=kind,reversal_mV=reversal,configured_gmax_uS=gain))
        c.edges=np.vstack([c.edges,edges]);self.edge_indices=np.array([e['index'] for e in self.edge_specs],dtype=np.int64)
        self.pf_bindings=[]
        for spec in self.edge_specs:
            post=c.names.index(spec['target']);original=[e for e in c.edges[:216] if int(e[1])==post and ' PF ' in c.names[int(e[0])] and not c.names[int(e[0])].endswith(' IN')];assert len(original)==1;edge=original[0];self.pf_bindings.append([edge[0],edge[3],edge[4]])
        self.pf_bindings=np.array(self.pf_bindings);self.configured_gains=np.array([s['configured_gmax_uS'] for s in self.edge_specs])
        if configuration.get('transient_quiet_PF_threshold',0):assert 0<configuration['transient_quiet_PF_threshold']<=.3
        self.v=np.tile(np.linspace(.01,.99,16),(2,1));self.ref=np.zeros((2,16));self.decoded=np.zeros(2);self.spikes=np.zeros((2,16),dtype=np.int8);self.counts=np.zeros(2,dtype=np.int64);self.total_counts=self.counts.copy();self.target=np.zeros(2);self.current=np.zeros(2);self.gates=np.zeros(2);self.tau=np.full(2,.025);self.delay=np.zeros(1);self.pointer=np.array(0,dtype=np.int64);self.measurement=np.zeros(8+len(edges));self.integral=np.zeros(len(edges));self.logs=[]
        self.phase_indices=np.array([[c.names.index(s+' RG flx'),c.names.index(s+' RG ext')] for s in 'LR'])
        self.encode(self.target,self.v.copy(),self.ref.copy(),self.decoded.copy(),self.spikes.copy(),self.counts.copy(),c.dt,200.,.005,self.tau,.001)
        chemical_statistics(c.v,c.edges,self.edge_indices,self.integral.copy(),self.measurement[8:].copy(),c.dt)

    def observe(self,model,data,previous,time,interval):
        from import_virtual_walker import mujoco
        root=mujoco.mj_name2id(model,mujoco.mjtObj.mjOBJ_BODY,'Root');bw=float(model.body_subtreemass[root])*9.81;forces=np.zeros(2) if previous is None else previous['foot_force_world'][:,1]/bw
        load=np.clip((forces-self.cfg['off_load_BW'])/(self.cfg['full_load_BW']-self.cfg['off_load_BW']),0,1);self.gates+=(load-self.gates)*(-np.expm1(-interval/self.cfg['load_filter_tau_s']));phase=np.empty(2)
        for k,(F,E) in enumerate(self.phase_indices):
            f=max(0.,self.c.v[F]-self.c.params[F,0]);e=max(0.,self.c.v[E]-self.c.params[E,0]);phase[k]=e/(f+e+1e-12)
            if 'rg_ext_reference_mV' in self.cfg:phase[k]*=min(1.,e/self.cfg['rg_ext_reference_mV'])
        self.current[:]=self.cfg['input_current_nA']*self.gates*phase;self.target[:]=self.current/20.;self.measurement[:8]=np.r_[forces,phase,self.target,self.current]

    def prepare(self):
        super().prepare();self.measurement[8:]=0.

    def advance(self):
        self.encode(self.target,self.v,self.ref,self.decoded,self.spikes,self.counts,self.c.dt,200.,.005,self.tau,.001)
        # 40Hz corresponds to full stance activation; only actual filtered
        # spikes populate this proxy, with the source synaptic clip preserved.
        self.proxy=-60.+100.*self.decoded;self.c.v[self.indices]=self.proxy
        if self.cfg.get('transient_quiet_PF_threshold',0):protect_transient(self.c.v,self.c.edges,self.edge_indices,self.pf_bindings,self.configured_gains,self.cfg['transient_quiet_PF_threshold'])
        chemical_statistics(self.c.v,self.c.edges,self.edge_indices,self.integral,self.measurement[8:],self.c.dt)

    def save(self,folder):
        keys=('t','counts','decoded','input_nA','load_phase_peak_current','load_gate','synaptic_integral_nA_s');np.savez_compressed(folder/'stance_recruitment_spikes.npz',**{k:np.array([row[i] for row in self.logs]) for i,k in enumerate(keys)})
        (folder/'STANCE_RECRUITMENT_ARCHITECTURE.json').write_text(json.dumps(dict(configuration=self.cfg,edges=self.edge_specs,actual_added_spiking_cells=32,spike_counts=self.total_counts.tolist(),scope='Experimental fixed neural stance primitive: two16-cell LIF populations,24chemical axons (zero gmax disconnected). Actual normal-load20ms filter times source RG-E fraction gates spike input;40Hz gives unit synaptic activation. Standardâˆ’10/âˆ’100mV reversals, mirrored weights fitted OFFLINE from source Hill/geometry target. Existing transient PF axons unchanged. No online desired force, inverse dynamics, pose clamp, or joint torque. Source108/216 unchanged; runtime110/240,736spiking+80graded cells. Peak-current fields are actual neural-substep maxima, not only2ms averages.'),indent=2),encoding='utf-8')


def probe(root,configuration_name='STANCE_RECRUITMENT_CONFIGURATION.json',output_name='ISOLATED_NEURAL_CHECK.json'):
    """Current and complete-state checks before a coupled trial."""
    import inspect
    from pathlib import Path
    from spiking_central_runtime import SpikingCentralController as C
    from recruitment import apply_recruitment
    from import_virtual_walker import mujoco
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());base=Path(fit['baseline']);record=json.loads((base/'SPIKING_COMPARISON.json').read_text());bank=json.loads((base.parent/'IMPORT_AUDIT.json').read_text())['muscles'];cfg=record['source_configuration'];settings=json.loads((root/configuration_name).read_text());model=mujoco.MjModel.from_xml_path(str(base/'model.xml'));data=mujoco.MjData(model);data.qpos[:]=record['result']['initial_qpos'];mujoco.mj_forward(model,data);allowed=set(inspect.signature(C.__bases__[0]).parameters)-{'source','bank'}
    def make(blend):
        C.configuration=record['configuration'];C.balance_configuration=dict(settings,blend=blend);C.sensory_spiking=True;C.clock_scale=1.;C.audit_state=False;C.reference_calculations=False;c=C(cfg['source'],bank,**{k:v for k,v in cfg.items() if k in allowed});apply_recruitment(c,cfg['pf_mn_recruitment_scale'],cfg['toe_mn_gmax'],cfg['hip_pf_mn_scale'],cfg['hip_flexor_pf_mn_scale']);c.apply_output_gain();return c
    zero=np.zeros(len(bank));c=make(1.);bw=record['result']['robot_mass_kg']*9.81;peak=np.zeros(24);continued=None
    def advance(c,k):
        previous=dict(foot_force_world=np.array([[0,.5*bw if k*.002>=1 else 0,0],[0,.5*bw if k*.002>=1 else 0,0]]));c.balance.observe(model,data,previous,k*.002,.002);c.step(k*.002,zero,zero,zero,{},.002,False)
    for k in range(1500):
        advance(c,k);peak=np.maximum(peak,c.balance.measurement[8:])
        if k==1000:
            initial=c.export_state();states=[]
            for repeat in range(2):
                c.import_state(initial)
                for j in range(1,33):advance(c,k+j)
                states.append(c.export_state())
            continued={name:bool(v.shape==states[1][name].shape and v.dtype==states[1][name].dtype and v.tobytes()==states[1][name].tobytes()) for name,v in states[0].items()};assert all(continued.values());c.import_state(initial)
    caps={row['role']:row['individual_axon_current_envelope_nA'] for row in fit['fits']};rows=[dict(**e,actual_substep_peak_nA=float(peak[k]),source_current_envelope_nA=caps[e['role']],ratio=float(peak[k]/max(1e-12,caps[e['role']]))) for k,e in enumerate(c.balance.edge_specs)];report=dict(all_states_byte_identical=all(continued.values()),state_fields=len(continued),added_spikes=c.balance.total_counts.tolist(),current_rows=rows,maximum_current_envelope_ratio=max(row['ratio'] for row in rows),scope='3s isolated actual736-spiking/80-graded controller: zero muscle afferents; normal-load step at1s gates new stance cells. Complete32-interval continuation replay includes physical input/phase observation. Does not establish coupled ground success.')
    report['configuration_source']=str(root/configuration_name);(root/output_name).write_text(json.dumps(report,indent=2),encoding='utf-8');return report


def write_checkpoint(root):
    """One compact ongoing ledger; avoid a report or script for each trial."""
    import hashlib
    from pathlib import Path
    root=Path(root);fit=json.loads((root/'OFFLINE_NEURAL_FIT.json').read_text());reference=Path(fit['baseline']);target=Path(fit['target']);full=json.loads(target.read_text());mechanical=json.loads((target.parent/'MECHANICAL_HOLD.json').read_text());cases=[]
    caps={r['role']:r['individual_axon_current_envelope_nA'] for r in fit['fits']}
    for folder in sorted((root/'W2L_VirtualWalker_Newer').glob('stance*')):
        if not (folder/'SPIKING_COMPARISON.json').exists():continue
        rec=json.loads((folder/'SPIKING_COMPARISON.json').read_text());r=rec['result'];a=rec['assessment'];z=np.load(folder/'trajectory.npz');bw=r['robot_mass_kg']*9.81;forces=np.stack([z[k][:,1] for k in ('rail_force_world','guide_force_world','harness_force')],axis=1);up=float(np.maximum(forces,0).sum(1).mean()/bw);down=float(np.maximum(-forces,0).sum(1).mean()/bw);row=dict(path=str(folder),mode=rec['mode'],feedback=rec['feedback'],air_gate=rec.get('air_gate'),gate=a['alternating_supported_stepping_gate'],failure=r['failure'],swings=[a['legs'][s]['unloaded_clear_swings'] for s in 'LR'],net_assistance_BW=r['mean_net_vertical_assistance_BW'],upward_channels_BW=up,downward_channels_BW=down,progress_m=r['root_forward_relative_to_belt_m'],knee_min_deg=r['anatomical_knee_min_deg'],near_arena_limit=any(v['near_limit'] for v in r['slide_travel_checks'].values()),warmed_wall_s=r['warmed_wall_s'],actual_spiking_cells=rec['spiking']['actual_spiking_cells'])
        s=np.load(folder/'stance_recruitment_spikes.npz');arch=json.loads((folder/'STANCE_RECRUITMENT_ARCHITECTURE.json').read_text());peaks=s['load_phase_peak_current'][:,8:].max(axis=0);row['actual_substep_current_peak_nA']=float(peaks.max());row['source_envelope_ratio']=float(max(peak/max(1e-12,caps[e['role']]) for peak,e in zip(peaks,arch['edges'])));row['added_spikes']=arch['spike_counts'];role_results=[]
        for role in ROLES:
            for command in full['shared_functional_torque_primitive']['motor_commands']:
                if command['role']!=role:continue
                side=0 if command['MN_name'].startswith('L') else 1;mask=(z['t']>=1)&(z['foot_vertical'][:,side]>.1*bw)
                if mask.any():role_results.append(dict(role=role,side='LR'[side],mean_loaded_tension_N=float(z['forces'][mask,command['bank_index']].mean()),mean_loaded_MN_mV=float(z['neuro'][mask,z['neuron_names'].tolist().index(command['MN_name'])].mean())))
        row['loaded_motor_outputs']=role_results;cases.append(row)
    report=dict(machine='easteregg2',source_sha256=fit['source_sha256'],preserved_baseline=str(reference),preserved_baseline_net_assistance_BW=.3647539631203412,source_voltage_equality_feasible=full['feasible'],shared_functional_primitive=full['shared_functional_torque_primitive'],mechanical_checks=mechanical,cases=cases,promoted=False,scope='New offline-derived neural stance primitive. No harness-free walking or reduced-support promotion. Existing704spiking+80graded baseline unchanged; addon736+80. Independent constant-voltage mechanics is not a neural pass.',implementation_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),Path(__file__).with_name('symmetric_stance.py'),Path(__file__).with_name('spiking_central_runtime.py')]})
    measured=root/'MEASURED_WRENCH_SUPPORT_FEASIBILITY.json'
    if measured.exists():
        audit=json.loads(measured.read_text());bad=[r for r in audit['samples'] if not r['foot_only_all_support_removed_feasible'] and 'loaded_foot_full_projection_m' in r];report['measured_support_audit']=dict(path=str(measured),valid_reconstructions=audit['valid_reconstruction_samples'],foot_only_feasible=audit['all_support_removed_feasible_on_valid'],loaded_obstructed_samples=len(bad),posterior_required_CoP_gap_m=[min(r['required_CoP_gap_m'] for r in bad),max(r['required_CoP_gap_m'] for r in bad)],scope=audit['scope'])
    equivalent=root/'AIR_DEFAULT_EQUIVALENCE.json'
    if equivalent.exists():report['air_baseline_array_equivalence']={k:v['all_values_equal'] for k,v in json.loads(equivalent.read_text()).items()}
    report['fresh_basin_status']='No material support/clearance improvement established. Tiny nominal support differences do not promote a candidate; no redundant fresh/basin campaign. Preserved baseline retains existing fresh and8/8 evidence.'
    for key,name in [('landing_phase_pair','LANDING_PHASE_PAIR_TARGETS_FORCE_FEASIBILITY.json'),('landing_primitive_gates','LANDING_PRIMITIVE_GATES.json'),('mechanical_drift_audit','MECHANICAL_HOLD_DRIFT_AUDIT.json')]:
        if (root/name).exists():report[key]=dict(path=str(root/name),data=json.loads((root/name).read_text()))
    interlock=root/'INTERLOCK_GATES.json'
    if interlock.exists():
        source=json.loads(interlock.read_text());report['hip_II_precontact']=dict(path=str(interlock),cases=[{k:v for k,v in c.items() if k in ('path','mode','air_gate','ground_gate','swings','net_support_BW','upward_BW','downward_BW','progress_m','failure','knee_min_deg','warmed_wall_s','spiking','graded','median_landing')} for c in source['cases']],promoted=False)
        verification=root/'LOAD_GATE_COUPLED_VERIFICATION.json'
        if verification.exists():report['hip_II_precontact']['load_gate_verification']=dict(path=str(verification),data=json.loads(verification.read_text()))
    (root/'CHECKPOINT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    lines=['# Symmetric stance primitive experiment','', '**Machine:** easteregg2. Existing704-spiking/80-graded baseline remains preferred.','',f'Baseline: `{reference}`','',f'Offline source-mechanics target: `{target}`','', 'Exact shared L/R motor voltages are infeasible with inherited force-law differences. One shared anatomical torque magnitude per functional pair is feasible, using each actual source moment arm and inverse sigmoid offline. Both feet carry205.613N (half411.225N body weight), equilibrium residual1.28eâˆ’13. This is not an online allocator.','', 'Independent source-Hill/MuJoCo check: equilibrium-preloaded muscles hold2s with all rail/guide contacts disabled, kneesâ‰¥2.96Â°, and root drift21.7mm; cold-muscle startup falls backward at1.634s. Neither condition is a neural standing/walking result.','', 'The new two16-cell spiking stance-PF populations use actual normal load (20ms filter) and source RG-E fraction. Fixed mirroredâˆ’10/âˆ’100mV chemical axons recruit MNs; inherited transient PFâ†’MN axons, force law, geometry, cadence, and supports are unchanged. Offline fitting uses actual source pre-attenuation current envelopes. Zero-force antagonist pathways are disconnected. Isolated72-field neural continuation is byte-identical; actual neural-substep current peaks are logged.','', 'Full tonic projection weakens propulsion and is rejected. Preserving transients by waiting for PF-quiet stance was rejected before simulation because the right hip-extensor had no qualifying samples. A calf-only projection is tested separately. No candidate is promoted merely for passing knee cycles.','', '| Case | True swings L/R | Net support %BW | Up/down %BW | Belt-relative progress m | Physical gate |','|---|---:|---:|---:|---:|---|']
    for row in cases:
        if row['mode']=='air':continue
        lines.append(f"| {Path(row['path']).name} | {row['swings'][0]}/{row['swings'][1]} | {100*row['net_assistance_BW']:.3f} | {100*row['upward_channels_BW']:.3f}/{100*row['downward_channels_BW']:.3f} | {row['progress_m']:+.4f} | {row['gate']} |")
    lines.extend(['','The swing gate requires1mm clearance for80ms,10mm peak clearance, load<.025BW, contralateral load>.05BW for at least half the interval, â‰¥3swings/leg, andâ‰¥80% alternation. Positive belt-relative progress, knee guards, current envelopes, and arena bounds remain separate acceptance requirements.','', 'Complete compact metrics, current peaks, actual loaded MN/tension outputs, and paths: `CHECKPOINT.json`. Fresh/basin checks are required before promotion. No manuscript files changed.'])
    if measured.exists():lines.extend(['','## Outcome and next mechanical target','','The full tonic candidate is rejected (negative progression, greater support, and arena-edge proximity). All three calf-only doses retain real alternating steps, but none lowers support below the preserved36.475%BW baseline. Therefore no new controller is promoted and no redundant finalist basin campaign is run. Saved air arrays match the original baseline exactly, including original neural and central/sensory spike arrays; added stance populations remain silent in air.','', 'The authoritative measured-wrench diagnostic reconstructs all24 inspected baseline phases. Foot-only support is compatible with only7of those24 same motions/wrenches. At14obstructed phases with a loaded foot, the required CoP lies28â€“163mm posterior to the actual contact and outside the full foot projection; heels already touch the belt. Three other obstructed samples are near-unloaded transitions. This is a trajectory/landing/balance target, not the old tiptoe or sensor-height issue. It does not imply that a changed gait/controller cannot walk.','', f'Validated measured-wrench audit: `{measured}`. Earlier differentiated-kinematics3D audit is preliminary because it omits parts of apparatus reconstruction; do not cite its feasibility count.','', 'Next: derive ROM-respecting landing/body-motion primitives and fit their actual neural output, rather than repeat global tonic or termination-gain sweeps.'])
    if (root/'LANDING_PRIMITIVE_GATES.json').exists():lines.extend(['','## Landing geometry and short mechanical witnesses','','Using actual transformed foot/toe corners, passive toe equilibrium, unchanged source ROM, and all14mechanical degrees of freedom, one right-stance landing target and its anatomical mirror are force-feasible. Static force witnesses require8.404/8.418%BW upward rail support. Independent constant-voltage0.5s source-Hill checks give preloaded mean upward support9.077/9.090%BW, no guide force or downward cancellation. These are short mechanical observations, not neural balance or sustained low-support walking. Initial development used the wrong local foot thickness axis; invalid results are preserved in IK_INVALID_LANDMARK_DEVELOPMENT.zip and excluded from accepted targets.','', 'Saved10ms velocity audit: final preloaded landing root speeds are approximately0.95mm/s, maximum joint speeds0.305deg/s, and last100ms root drift0.086mm. The older2s home witness ends at66.1mm/s forward speed. Neither is presented as a demonstrated steady equilibrium. See MECHANICAL_HOLD_DRIFT_AUDIT.json for estimated terminal velocities and exact scope.','', 'The first fixed load/RG-phase neural projection fitted to the landing pair passes air and yields10/9 actual ground swings, but increases net support to47.466%BW (47.472up/0.0055down), with+0.338m relative-belt progress. It is rejected. A post-contact recruitment fit does not solve pre-contact landing placement. No fresh/basin run is warranted for this inferior candidate. Next mechanism is audited against actual hip-extensor-II observability before adding a specific hip-retraction/knee-extension timing pathway.'])
    if interlock.exists():
        lines.extend(['','## Hip-II pre-contact recruitment and load gating','','Actual hip-extensor receptor/actuator length Jacobians have the same positive sign under hip flexion. Raw II-current/hip-angle correlations are0.91L/0.90R; decoded sensory spikes ripple, so static proxy separation alone is not a phase detector. Dedicated inhibitory INs were added instead of hijacking existing Ia cells, whose RG/PF collaterals would change central phase. Optional zero gain omits new topology and is byte-identical to the omitted option under the isolated input check. Original704-spiking/80-graded configuration remains preferred.','', 'Exact source-Hill replay at recorded pre-step poses matches every muscle tension within4eâˆ’13N. Inhibition cuts precontact knee-extensor tension349/340N to52/44N, yet knee-extension onset remains about298ms before touchdown: knee-flexor MNs remain nearâˆ’91mV and produce only6â€“7N, mostly passive. Direct heteronymous excitation changes MN voltage toâˆ’80mV but still gives only9N. This is an actuator recruitment-range problem in addition to timing, not passive knee-extensor tension.','', 'Two additional16-cell adapting LIF PF populations were calibrated with expected Hip-II inputs. Input0.2ÂµS fails recruitment;0.4ÂµS gives mean chemical activation0.983 flexed and effectively0 retracted. The source-current-derived stance projection1.493ÂµS has a134.4nA maximum atâˆ’100mV, below the178.2nA source-axon envelope. Air then has5/5 antiphase cycles and~47Â° knee excursions, but the coupled ground case worsens support and backward progression. This limited input/readout fit is not complete-cycle Hunt training.','', 'An explicit engineering normal-foot-loadâ†’graded inhibitory-INâ†’extraPF pathway uses previous-frame actual GRF,8*clip((F/BWâˆ’0.025)/0.075,0,1)nA, and2ÂµS inhibition. It is not labelled tendon Ib and applies no body force or joint torque. Original heel/toe and II routes remain intact. Recorded input matches the actual-force mapping exactly (0nA error,8nA maxima). After150ms continuous load>.1BW, addedPF activation is0.00035L/0.00441R. All air trajectory, central-spike, and sensory-spike arrays match before/after gating byte-for-byte, cropping only two inactive new IN channels. Complete52-field neural continuation is exact.','', '| Ground case | Spiking/graded cells | True swings L/R | Net %BW | Up/down %BW | Relative-belt progress m |','|---|---:|---:|---:|---:|---:|'])
        for row in source['cases']:
            if row['mode']=='air':continue
            lines.append(f"| {Path(row['path']).name} | {row['spiking']}/{row['graded']} | {row['swings'][0]}/{row['swings'][1]} | {100*row['net_support_BW']:.3f} | {100*row['upward_BW']:.3f}/{100*row['downward_BW']:.3f} | {row['progress_m']:+.4f} |")
        lines.extend(['','The phase-correct swing projection0.4874ÂµS with load gating preserves10/9 actual swings and positive progression, but its0.355percentage-point nominal support reduction is not a demonstrated improvement: median swing clearance falls44.4/47.0mmâ†’36.1/42.5mm, knee extension still begins298ms before touchdown, and CoP remains94.5mm anterior to CoM (baseline97.2mm). No candidate is promoted, and no fresh/basin campaign or further strength sweep is justified. The independent low-support static phase pair remains a partial mechanical reference, not a feasible complete gait or a neural success. Next work must resolve whole-cycle landing/body-motion feasibility rather than add generic recruitment gains.','', 'Reproduction for any candidate: run_spiking_central.py --mode supported --duration10 --configuration <full candidate JSON in this experiment folder> --sensory-spiking --initial-qpos D:\\Users\\Ben Bolen\\NeuroCampaign_20261031\\spiking\\stance_home_v7\\standing_home.npz --pool-inherited-hip-physics --belt-speed0.2 --out <this experiment folder> --tag <new tag>. Supply spaces before10 and0.2 in the actual command. Source guard and four original recruitment scalars load from GENTLE_HIP_0P6_FLEX_0P5.json. Full configurations and numerical outcomes are in INTERLOCK_GATES.json and per-case SPIKING_COMPARISON.json.','', 'All neural and Hill equations use the existing sparse Numba/preallocated path. Observed warmed10s coupled loops take5.25â€“6.65s in these conditions; no matched scalar-reference benchmark was run for the enlarged736-cell variants. Historical central-only speed ratios do not establish their current full-model speedup. No manuscript or canonical AnimatLab model changed.'])
    if interlock.exists():
        runner=Path(__file__).with_name('run_spiking_central.py');outer=root.parent
        lines=[line for line in lines if not line.startswith('Reproduction for any candidate:')]
        lines.extend(['','## Exact reproduction commands','','Run on easteregg2. First command reproduces the preferred preserved baseline. Second reproduces the unpromoted swing-primitive experiment. Choose a new tag if that output already exists.','', '```powershell',rf"& 'D:\Anaconda\envs\myo\python.exe' '{runner}' --mode supported --duration 10 --configuration '{outer/'CENTRAL_RESET_REST_CANDIDATE.json'}' --sensory-spiking --initial-qpos '{outer/'stance_home_v7/standing_home.npz'}' --pool-inherited-hip-physics --belt-speed 0.2 --out '{outer}' --tag baseline_reproduction_unique",'```','','```powershell',rf"& 'D:\Anaconda\envs\myo\python.exe' '{runner}' --mode supported --duration 10 --configuration '{root/'HIP_II_SWING_PRIMITIVE_CANDIDATE.json'}' --sensory-spiking --initial-qpos '{outer/'stance_home_v7/standing_home.npz'}' --pool-inherited-hip-physics --belt-speed 0.2 --out '{root}' --tag hipII_swing_reproduction_unique",'```'])
    (root/'EXPERIMENT_RECORD.md').write_text('\n'.join(lines)+'\n',encoding='utf-8');return report


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser();parser.add_argument('root');args=parser.parse_args();print(json.dumps(probe(args.root),indent=2))
