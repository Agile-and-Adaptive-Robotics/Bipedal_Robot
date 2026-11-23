function record = test_sns_gui_voltage_edits(outputFolder)
% Edit the real property controls and check the saved JSON and generated masks.
if nargin < 1, outputFolder = tempname; end
if ~isfolder(outputFolder), mkdir(outputFolder); end
model = 'SNS_GuiCircuit';
assert(~bdIsLoaded(model), 'Close the open GUI-generated model before running this test.');
root = fileparts(fileparts(mfilename('fullpath')));
app = SNSGui();
app.f.Visible = 'off';
closeApp = onCleanup(@()delete(app.f)); %#ok<NASGU>
app.loadJson(fullfile(root,'results','reflex_validation','biological_knee_circuit.json'));

editField(app, 'block', 'Thr', -46.25);
editField(app, 'block', 'Slope', 11.25);
editField(app, 'edge', 'ThrPre', -47.75);
editField(app, 'edge', 'SlopePre', 7.25);
file = fullfile(outputFolder,'edited_circuit.json');
app.saveJson(file);
saved = jsondecode(fileread(file));
assert(saved.neurons(1).Thr_mV == -46.25 && saved.neurons(1).Slope_mV == 11.25);
assert(saved.synapses(1).ThrPre_mV == -47.75 && saved.synapses(1).SlopePre_mV == 7.25);
app.loadJson(file);
app.cbBuild([],[],outputFolder);
load_system(fullfile(outputFolder,[model '.slx']));
closeModel = onCleanup(@()close_system(model,0)); %#ok<NASGU>
cellPath = [model '/' matlab.lang.makeValidName(saved.neurons(1).name)];
synapses = find_system(model,'SearchDepth',1,'MaskType','SNS NonSpiking Synapse');
synapse = synapses(cellfun(@(s) str2double(get_param(s,'SlopePre')) == 7.25, synapses));
assert(isscalar(synapse));
assert(str2double(get_param(cellPath,'Thr')) == -46.25);
assert(str2double(get_param(cellPath,'Slope')) == 11.25);
assert(str2double(get_param(synapse{1},'ThrPre')) == -47.75);
record = struct('neuron',cellPath,'Thr_mV',-46.25,'Slope_mV',11.25, ...
    'synapse',synapse{1},'ThrPre_mV',-47.75,'SlopePre_mV',7.25);
end

function editField(app, selection, name, value)
app.selType = selection;
app.selIdx = 1;
app.fillPropertyPanel();
if strcmp(selection,'block'), params = app.blocks(1).params;
else, params = app.edges(1).params; end
fields = app.propChildren(cellfun(@(h) isa(h,'matlab.ui.control.NumericEditField'), app.propChildren));
field = fields{strcmp(params(:,1),name)};
field.Value = value;
field.ValueChangedFcn(field,[]);
end
