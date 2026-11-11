"""Exercise our generated Qt page against a temporary native model copy."""
import sys,json,tempfile,traceback
from pathlib import Path
from PySide6.QtCore import QTimer,Qt
from PySide6.QtWidgets import QApplication
BASE=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(BASE/'Code/MuJoCo_SNS/connectome_app'))
from app import Window
from native_w2l import read_model
from native_files import save_model
from project import flatten
folder=tempfile.TemporaryDirectory(prefix='connectome_ui_');p=Path(folder.name)/'W2L108.aproj'
source=BASE/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer.aproj'
save_model(read_model(source),p)
app=QApplication([]);w=Window(str(p));w.setAttribute(Qt.WidgetAttribute.WA_ShowWithoutActivating);w.setWindowFlag(Qt.WindowType.WindowDoesNotAcceptFocus);w.show()
errors=[]
def console(level,message,line,source):
 print('WEB',message,line,flush=True)
 if 'Error' in message or 'ReferenceError' in message:errors.append(message)
w.page.javaScriptConsoleMessage=console
def js(code,callback):w.page.runJavaScript(code,callback)
def check():
 js("JSON.stringify(typeof desktopBridge==='undefined'?{}:{ready:!!desktopBridge,path:models[activeTab].path,n:nodes.length,error:typeof circuitProjection})",loaded)
def loaded(value):
 try:
  state=json.loads(value or '{}')
  if not state.get('path'):QTimer.singleShot(500,check);return
  print('LOADED',state,flush=True)
  js(overview_check_script(),overview_checked)
 except Exception:traceback.print_exc();app.exit(1)
def overview_check_script():
 return """function smokeMeaning(){const f=flattenView(modelSpec());return JSON.stringify(f.es.map(e=>[e.id,e.fromId,e.toId,e.synapse]).sort((a,b)=>a[0].localeCompare(b[0])));}smokeOverviewNodes=nodes.map(n=>({id:n.id,scope:n.scope,label:n.label,y:n.y}));smokeBefore=smokeMeaning();smokeRG=nodes.find(n=>n.scope?.layer==='RG'&&n.scope.side==='Left');smokeRGId=smokeRG.id;smokeRGAfferents=edges.filter(e=>e.to===smokeRGId&&nodeById(e.from)?.scope?.layer==='Afferents').map(e=>[e.id,desktopEndpointNode(e,'to').id]);selectNode(smokeRGId);document.activeElement?.blur();window.dispatchEvent(new KeyboardEvent('keydown',{key:'G',ctrlKey:true,shiftKey:true,bubbles:true}));smokeUnpacked=!nodeById(smokeRGId);smokeSame=smokeBefore===smokeMeaning();smokeAfferentsSame=smokeRGAfferents.every(([id,target])=>desktopEndpointNode(edgeById(id),'to').id===target);undo();JSON.stringify({overview:smokeOverviewNodes,unpacked:smokeUnpacked,same:smokeSame,afferentsSame:smokeAfferentsSame,afferents:smokeRGAfferents.length,restored:!!nodeById(smokeRGId)})"""
def assert_overview(data):
 rg=[n for n in data['overview'] if n.get('scope',{}).get('layer')=='RG'];assert {n['scope']['side'] for n in rg}=={'Left','Right'}
 assert data['unpacked'] and data['same'] and data['afferentsSame'] and data['restored'],str(data)
def motor_script():
 return "enterSub(nodes.find(n=>n.scope?.layer==='MN'&&n.scope.side==='Left'));JSON.stringify({view:!!circuitContext,neurons:nodes.map(n=>({type:n.type,label:n.label,x:n.x,y:n.y})),types:nodes.map(n=>n.type),aliases:nodes.filter(n=>n.external_alias).length,edges:edges.length,document:modelSpec()})"
def overview_checked(value):
 try:
  assert_overview(json.loads(value));print('PASS W2L108 bilateral overview, Ctrl+Shift+G afferent bindings, and Undo.',flush=True)
  w.grab().save(str(BASE/'tmp/connectome_bilateral_overview.png'));js(motor_script(),motor)
 except Exception:traceback.print_exc();app.exit(1)
def motor(value):
 try:
  data=json.loads(value);ns,es,_=flatten(data['document']);assert(len(ns),len(es))==(108,216);assert data['view'];assert any(n['type']=='ADAPTER' for n in data['neurons']);assert any(n['type']=='RC' for n in data['neurons']);print('MOTOR',len(data['neurons']),data['edges'],flush=True)
  js("smokeEdge=edges.find(e=>!e.view_only);smokeEdge.synapse.gmax_uS=.456789;smokeEdgeId=smokeEdge.id;desktopSave('save');JSON.stringify({id:smokeEdgeId})",saved_request)
 except Exception:traceback.print_exc();app.exit(1)
def saved_request(value):
 global edited_id
 edited_id=json.loads(value)['id'];QTimer.singleShot(900,verify_save)
def verify_save():
 try:
  again=read_model(p);edge=next(e for e in flatten(again)[1] if e['id']==edited_id);assert edge['synapse']['gmax_uS']==.456789
  js("desktopSave('save');JSON.stringify({native:!!circuitContext.root.native,path:models[activeTab].path})",lambda _:QTimer.singleShot(700,capture))
 except Exception:traceback.print_exc();app.exit(1)
def capture():
 w.grab().save(str(BASE/'tmp/connectome_motor_preview.png'))
 js("leaveCircuit();JSON.stringify({context:!!circuitContext,layers:nodes.map(n=>n.label),document:modelSpec()})",finish)
def finish(value):
 try:
  data=json.loads(value);assert not data['context'];assert len(flatten(data['document'])[0])==108;assert not errors;print('PASS Qt page native import, paired MN, physical symbols, two native saves, and return to layers.',flush=True)
  next_library()
 except Exception:traceback.print_exc();app.exit(1)
families=iter([('syn6',794,3200),('s3k',410,1186),('w2lvar',888,7866)])
def next_library():
 global expected
 try:expected=next(families)
 except StopIteration:print('PASS all four desktop hierarchy and paired motor views.',flush=True);app.exit(0);return
 w.open_controller(expected[0]);QTimer.singleShot(700,library_opened)
def library_opened():
 js(overview_check_script(),library_overview)
def library_overview(value):
 try:
  assert_overview(json.loads(value));print('PASS',expected[0],'bilateral overview, keyboard unpack, afferents, and Undo.',flush=True);js(motor_script(),library_motor)
 except Exception:traceback.print_exc();app.exit(1)
def library_motor(value):
 try:
  data=json.loads(value);ns,es,_=flatten(data['document']);assert(len(ns),len(es))==expected[1:];assert data['view'];assert 'SN-Ia' in data['types'];assert 'SN-II' in data['types'];assert 'SN-Ib' in data['types'];assert data['aliases']>0
  if any(n['type']=='RC' for n in ns):assert 'RC' in data['types']
  print('PASS desktop',expected[0],'full graph',len(ns),len(es),'motor aliases',data['aliases'],flush=True)
  # A new tab opened from inside the circuit must not carry the old model
  # into that tab. Returning to the old tab must preserve the complete graph.
  js("smokeTab=activeTab;desktopNew();smokeBlank=modelSpec().nodes.length;switchTab(smokeTab);JSON.stringify({blank:smokeBlank,document:modelSpec()})",library_tabs)
 except Exception:traceback.print_exc();app.exit(1)
def library_tabs(value):
 try:
  data=json.loads(value);assert data['blank']==0;assert len(flatten(data['document'])[0])==expected[1];assert not errors;next_library()
 except Exception:traceback.print_exc();app.exit(1)
QTimer.singleShot(1000,check);QTimer.singleShot(90000,lambda:app.exit(2));code=app.exec();w.allow_close=True;w.close();folder.cleanup();sys.exit(code)
