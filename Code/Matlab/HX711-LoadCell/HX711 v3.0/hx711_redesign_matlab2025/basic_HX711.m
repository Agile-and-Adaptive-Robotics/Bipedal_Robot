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
            configurePin(parentObj,inputPins{1},'DigitalInput');  % Data Pin
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
