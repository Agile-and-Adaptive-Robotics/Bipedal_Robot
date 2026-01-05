classdef HelloWorld < matlabshared.addon.LibraryBase
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access = private, Constant = true)
        READ_COMMAND = hex2dec('01')
    end  
    
    properties(Access = protected, Constant = true)
        LibraryName = 'ExampleAddon/HelloWorld'
        DependentLibraries = {}
        LibraryHeaderFiles = {}
        CppHeaderFile = fullfile(arduinoio.FilePath(mfilename('fullpath')), 'src', 'HelloWorld.h')
        CppClassName = 'HelloWorld'
    end    
    
    methods
        function obj = HelloWorld(parentObj)
            obj.Parent = parentObj;
            
        end 
        
        function out = read(obj)
            cmdID = obj.READ_COMMAND;
            inputs = [];
            output = sendCommand(obj, obj.LibraryName, cmdID, inputs);
            out = char(output');
        end
    end
end

