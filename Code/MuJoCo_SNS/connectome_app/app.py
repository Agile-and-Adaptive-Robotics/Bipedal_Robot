"""Desktop host for the existing full connectome block editor."""
from __future__ import annotations
import argparse,json,os,sys
from pathlib import Path
from PySide6.QtCore import QObject,Signal,Slot,QUrl,QTimer
from PySide6.QtGui import QAction,QKeySequence,QDesktopServices
from PySide6.QtWidgets import QApplication,QMainWindow,QFileDialog,QMessageBox
from PySide6.QtWebEngineWidgets import QWebEngineView
from PySide6.QtWebEngineCore import QWebEnginePage,QWebEngineProfile,QWebEngineSettings
from PySide6.QtWebChannel import QWebChannel
from project import FileSession,ProjectError,normalize
from build_frontend import build

HERE=Path(__file__).parent
def workspace():
    base=Path('D:/Users')/Path.home().name if Path('D:/Users').exists() else Path.home()/'Documents'
    path=base/'Connectome_Editor_Projects';path.mkdir(parents=True,exist_ok=True);return path

class Bridge(QObject):
    opened=Signal(str);saved=Signal(str);report=Signal(str)
    def __init__(self,window):super().__init__();self.window=window;self.files=FileSession();self.states={}
    @Slot()
    def ready(self):self.window.initial_document(self.window.initial_path)
    @Slot(str)
    def closed(self,tab):self.states.pop(tab,None)
    @Slot(str)
    def changed(self,payload):
        try:
            data=json.loads(payload);self.states[data['tab']]=data;self.window.active=data['tab'];self.window.refresh_title()
        except Exception as e:self.report.emit('State error: '+str(e))
    def open_path(self,path,tab=None):
        import uuid
        tab=tab or uuid.uuid4().hex;doc=self.files.open(path,tab)
        self.opened.emit(json.dumps({'tab':tab,'path':str(Path(path).resolve()),'document':doc}))
    @Slot()
    def open(self):
        path,_=QFileDialog.getOpenFileName(self.window,'Open connectome',str(workspace()),'Connectomes (*.json *.aproj *.asim);;AnimatLab model (*.aproj *.asim);;Editor project (*.json);;All files (*)')
        if not path:return
        try:self.open_path(path)
        except Exception as e:self.report.emit('Open failed: '+str(e))
    @Slot(str)
    def save(self,payload):
        try:
            request=json.loads(payload);tab=request['tab'];kind=request.get('kind','save');path=None
            if kind!='save' or tab not in self.files.paths:
                native=request['document'].get('native') if kind!='subnetwork' else None
                current=self.files.paths.get(tab)
                suffix=native.get('primary_suffix','.aproj') if native else '.sns.json'
                filename=(current.stem+'_copy' if current and native else request.get('name','Connectome'))+suffix
                folder=current.parent if current else workspace()
                path,_=QFileDialog.getSaveFileName(self.window,'Save subnetwork copy' if kind=='subnetwork' else 'Save connectome as',str(folder/filename),'AnimatLab model (*'+suffix+')' if native else 'Connectome project (*.json)')
                if not path:return
                if not path.lower().endswith(suffix):path+=suffix
            target_tab=tab if kind!='subnetwork' else tab+':subnetwork-copy'
            destination,doc=self.files.save(request['document'],target_tab,path)
            self.saved.emit(json.dumps({'tab':tab,'kind':kind,'path':str(destination),'document':doc}))
        except Exception as e:self.report.emit('Save failed: '+str(e))
    @Slot(str)
    def build(self,payload):
        try:
            from runtime import build_project
            result=build_project(json.loads(payload));self.window.built=result
            info={k:v for k,v in result.items() if k not in ('network','compiled')}
            self.report.emit('Neural model built from the current editable document. No integration steps executed.\n'+json.dumps(info,indent=2))
        except Exception as e:self.report.emit('Build stopped: '+str(e))
    @Slot(str)
    def export_svg(self,payload):
        try:
            path,_=QFileDialog.getSaveFileName(self.window,'Export drawing',str(workspace()/'Connectome.svg'),'SVG (*.svg)')
            if path:Path(path).write_text(payload,encoding='utf-8');self.report.emit('SVG saved: '+path)
        except Exception as e:self.report.emit('SVG export failed: '+str(e))

class Page(QWebEnginePage):
    def javaScriptConsoleMessage(self,level,message,line,source):
        print(f'Web: {message} ({line})',flush=True)

class Window(QMainWindow):
    def __init__(self,path=None):
        super().__init__();self.initial_path=path;self.initial_loaded=False;self.active=None;self.built=None;self.allow_close=False;self.setWindowTitle('Connectome Editor');self.resize(1460,940)
        self.view=QWebEngineView(self);self.setCentralWidget(self.view)
        # Off-the-record browser profile: original editor localStorage is untouched.
        self.profile=QWebEngineProfile(self);self.page=Page(self.profile,self);self.view.setPage(self.page)
        self.page.settings().setAttribute(QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls,False)
        self.channel=QWebChannel(self.page);self.bridge=Bridge(self);self.channel.registerObject('desktop',self.bridge);self.page.setWebChannel(self.channel)
        cache=workspace()/'.app';page=build(cache/'Connectome_Editor.html');self.view.setUrl(QUrl.fromLocalFile(str(page.resolve())))
        menu=self.menuBar().addMenu('&File')
        for label,shortcut,js in [('New','Ctrl+N','desktopNew()'),('Open…','Ctrl+O','desktopBridge.open()'),('Save','Ctrl+S','desktopSave("save")'),('Save As…','Ctrl+Shift+S','desktopSave("save_as")'),('Save Subnetwork As…',None,'desktopSave("subnetwork")'),('Build Neural Model','F5','desktopBuild()')]:
            action=QAction(label,self)
            if shortcut:action.setShortcut(QKeySequence(shortcut))
            action.triggered.connect(lambda checked=False,script=js:self.page.runJavaScript(script));menu.addAction(action)
        library=menu.addMenu('Open controller')
        for family in ('syn6','s3k','w2lvar','W2L108'):
            description={'syn6':'syn6 — reviewed final','s3k':'s3k — saved standing','w2lvar':'w2lvar — saved standing','W2L108':'W2L108 — editable native model'}[family]
            action=QAction(description,self);action.triggered.connect(lambda checked=False,key=family:self.open_controller(key));library.addAction(action)
        folder_action=QAction('Open project folder',self);folder_action.triggered.connect(lambda:QDesktopServices.openUrl(QUrl.fromLocalFile(str(workspace()))));menu.addAction(folder_action)
        helpmenu=self.menuBar().addMenu('&Help');action=QAction('Commands',self);action.triggered.connect(lambda:QMessageBox.information(self,'Editor commands','Original editor commands are retained.\n\nDrag palette → canvas: place\nShift/Alt-drag or Alt-click two cells: connect\nRubber-band / Ctrl-click: select\nDrag selection: group move\nAlt-click wire: bend; drag pink handle\nCtrl+G: pack subnetwork; Ctrl+Shift+G: unpack\nDouble-click module: enter; F4/Esc: leave\nF2/F3: hide/show layers\nCtrl+C/V: copy/paste\nDelete: delete selected graph elements\nCtrl+Z/Y: undo/redo\nCtrl+S: save; Ctrl+Shift+S: save as\nF5: build neural model without running it.'));helpmenu.addAction(action)
    def open_controller(self,family):
        try:
            if family=='W2L108':
                path=workspace()/'Models/W2L108/W2L108.aproj'
                if not path.exists():
                    from native_w2l import read_model
                    from native_files import save_model
                    source=HERE.parents[2]/'tmp/functional_subnetwork_remote/W2L108_verified/editable/W2L_VirtualWalker_Newer.aproj'
                    save_model(read_model(source),path)
                self.bridge.open_path(path);return
            folder=workspace()/'Controllers';index=folder/'LIBRARY_INDEX.json'
            if index.exists():
                row=next(r for r in json.loads(index.read_text()) if r['family']==family);path=Path(row['hierarchy'])
            else:path=folder/(('syn6_reviewed' if family=='syn6' else family)+'_hierarchy.sns.json')
            if not path.exists():
                from controller_library import prepare
                prepare(path.parent)
            self.bridge.open_path(path)
        except Exception as e:self.bridge.report.emit('Could not open controller: '+str(e))
    def initial_document(self,path):
        if self.initial_loaded:return
        self.initial_loaded=True
        candidate=Path(path) if path else workspace()/'Paired_MN.sns.json'
        if not candidate.exists() and path is None:
            candidate.write_bytes((HERE/'examples/paired_mn.sns.json').read_bytes())
        if candidate.exists():
            try:self.bridge.open_path(candidate)
            except Exception as e:self.bridge.report.emit(str(e))
    def refresh_title(self):
        info=self.bridge.states.get(self.active,{})
        self.setWindowTitle(('• ' if info.get('dirty') else '')+'Connectome Editor — '+info.get('name','Untitled'))
    def closeEvent(self,event):
        dirty=[x for x in self.bridge.states.values() if x.get('dirty')]
        if self.allow_close or not dirty:event.accept();return
        choice=QMessageBox.question(self,'Unsaved connectome edits','Save unsaved projects before closing?',QMessageBox.StandardButton.Save|QMessageBox.StandardButton.Discard|QMessageBox.StandardButton.Cancel,QMessageBox.StandardButton.Save)
        if choice==QMessageBox.StandardButton.Cancel:event.ignore();return
        if choice==QMessageBox.StandardButton.Discard:event.accept();return
        event.ignore()
        def finish(items):
            try:
                for item in items:
                    if not item.get('dirty'):continue
                    tab=item['tab'];path=None
                    if tab not in self.bridge.files.paths:
                        path,_=QFileDialog.getSaveFileName(self,'Save '+item['name'],str(workspace()/(item['name']+'.sns.json')),'JSON (*.json)')
                        if not path:return
                    self.bridge.files.save(item['document'],tab,path)
                self.allow_close=True;self.close()
            except Exception as e:QMessageBox.warning(self,'Could not save',str(e))
        self.page.runJavaScript('desktopWorkspace()',finish)

def main():
    if hasattr(sys.stdout,'reconfigure'):sys.stdout.reconfigure(encoding='utf-8',errors='replace')
    parser=argparse.ArgumentParser();parser.add_argument('project',nargs='?');args=parser.parse_args()
    app=QApplication(sys.argv);app.setApplicationName('Connectome Editor');window=Window(args.project);window.show();sys.exit(app.exec())
if __name__=='__main__':main()
