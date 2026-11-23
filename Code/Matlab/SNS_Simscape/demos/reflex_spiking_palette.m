function [name,cleanup]=reflex_spiking_palette(libraryFolder)
name='SNS_Library';cleanup=[];
if getSimulinkBlockHandle([name '/HybridSpikingSynapse'])>0
    return
end
% Create the optional spiking blocks in a private palette. They are copied
% into the experiment model, so subsequent simulations only need the base library.
file=[tempname '.slx']; [~,name]=fileparts(file);
new_system(name,'Library');save_system(name,file);
cleanup=onCleanup(@()removePalette(name,file));
snsDestinationLibrary=name; %#ok<NASGU> used by sns_build_spiking
run(fullfile(libraryFolder,'sns_build_spiking.m'));
end

function removePalette(name,file)
if bdIsLoaded(name),close_system(name,0);end
if isfile(file),delete(file);end
end
