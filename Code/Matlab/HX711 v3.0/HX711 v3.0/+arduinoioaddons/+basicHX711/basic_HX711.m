% Copyright 2018, Nicholas Giacoboni
% 
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright 
%    notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the 
%    documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT 
% NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
% WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
% OF SUCH DAMAGE.
%**************************************************************************
classdef basic_HX711 < matlabshared.addon.LibraryBase
    properties(Access = protected)
        Pins %Data and Clock pin
    end
    properties (Access = private, Constant = true)
        READ_HX711 = hex2dec('01')
    end
    properties(Access = protected, Constant = true)
        LibraryName = 'basicHX711/basic_HX711'
        DependentLibraries = {}
        LibraryHeaderFiles = 'basic_HX711/basic_HX711.h'
        CppHeaderFile = fullfile(arduinoio.FilePath(mfilename('fullpath')), 'src', 'basic_HX711.h')
        CppClassName = 'basic_HX711'
    end
    
    methods(Hidden, Access = public)
        function obj = basic_HX711(parentObj,inputPins)
            obj.Parent = parentObj;
            obj.Pins = getTerminalsFromPins(obj.Parent,inputPins);
            configurePin(parentObj,inputPins{1},'DigitalInput'); % Data Pin
            configurePin(parentObj,inputPins{2},'DigitalOutput'); % Clock Pin
        end
    end
    
    methods(Access = public)
        function force = read_HX711(obj)
            cmdID = obj.READ_HX711;
            inputs = obj.Pins;
            value = sendCommand(obj, obj.LibraryName, cmdID, inputs);
            value(3)=bitshift(value(3),16);
            value(2)=bitshift(value(2),8);
            force = bitor(value(3),bitor(value(2),value(1)));
        end
    end
end

