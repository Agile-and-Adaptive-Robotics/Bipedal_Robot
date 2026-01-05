%**************************************************************************
%     Copyright (C) 2019  Nicholas Giacoboni
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%**************************************************************************
classdef advanced_HX711 < matlabshared.addon.LibraryBase
    
    properties(Access = protected)
        Pins %Data and Clock pin
        Gain %Gain (32, 64, 128)
        Interrupt %If you want to use interrupt 
        %NOTE: comunication between Arduino and HX711's ADC is time
        %sensitive, I make sure you don't use interrupt during comunication
    end
    
    properties (Access = private, Constant = true)
        TEST = hex2dec('00')
        READ_HX711 = hex2dec('01')
        POWER_DOWN = hex2dec('02')
        POWER_UP = hex2dec('03')
    end
    
    properties(Access = protected, Constant = true)
        LibraryName = 'advancedHX711/advanced_HX711'
        DependentLibraries = {}
        LibraryHeaderFiles = 'advanced_HX711/advanced_HX711.h'
        CppHeaderFile = fullfile(arduinoio.FilePath(mfilename('fullpath')), 'src', 'advanced_HX711.h')
        CppClassName = 'advanced_HX711'
    end
    
    methods(Hidden, Access = public)
        % CONSTRUCTOR: LoadCell = addon(a,'ExampleAddon/advanced_HX711',...
        %   'Pins',{'D2','D3'},'Gain',128,'Interrupt',true);
        function obj = advanced_HX711(parentObj,varargin)
            % Check the number of input parameters
            if(nargin < 3)
                matlabshared.hwsdk.internal.localizedError(...
                    'MATLAB:narginchk:notEnoughInputs');
            elseif(nargin > 7)
                matlabshared.hwsdk.internal.localizedError(...
                    'MATLAB:narginchk:tooManyInputs');
            end
            % Default value
            defaultGain = 128;
            defaultInterrupt = false;
            % Error messages
            errorMsg1 = 'Selectable gain are 32, 64 and 128';
            errorMsg2 = 'Interrupt parameter must be true or false';
            %*******************Validation Functions*******************
            %Gain must be 32, 64 or 128
            validGain = @(x) assert(isnumeric(x) && isscalar(x) && ...
                ((x == 32) || (x == 64) || (x == 128)),errorMsg1);
            %You would or wouldn't like to use Interrupt (TRUE or FALSE)
            validInterrupt = @(x) assert(islogical(x),errorMsg2);
            % Input Parser scheme
            try
                p = inputParser; %Object
                addParameter(p, 'Pins',[]);
                addParameter(p, 'Gain', defaultGain,validGain);
                addParameter(p, 'Interrupt', defaultInterrupt,validInterrupt);
                parse(p, varargin{1:end});
            catch e
                throwAsCaller(e);
            end
            obj.Parent = parentObj; %Parent the add-on object to the Arduino object.
            obj.Pins = cellstr(p.Results.Pins); %Cell array
            % Input and gain selection is controlled by the number of the
            % input PD_SCK pulses (Gain 128 = 25 sck Pulses, Gain 32 = 26
            % sck Pulses, Gain 64 = 27 sck Pulses. 
            % Each PD_SCK pulse shifts out one bit, starting with the MSB 
            % bit first, until all 24 bits are shifted out.
            switch p.Results.Gain
                case 128
                    obj.Gain = 1;
                case 64
                    obj.Gain = 3;
                case 32
                    obj.Gain = 2;
            end
            obj.Interrupt = p.Results.Interrupt;
            %***********Pin config.**********************
            configurePin(parentObj,obj.Pins{1},'DigitalInput'); % Data Pin
            configurePin(parentObj,obj.Pins{2},'DigitalOutput'); % Clock Pin
            % "DigitalWrite" and "DigitalRead" in .h file need terminals
            % not physical pins like configurePin;
            obj.Pins = getTerminalsFromPins(obj.Parent,obj.Pins);
        end
    end
    
    methods(Access = public)
        
        % Read data
        function force = read_HX711(obj)
            cmdID = obj.READ_HX711;
            inputs = [obj.Pins(1) obj.Pins(2) obj.Gain obj.Interrupt];
            value = sendCommand(obj, obj.LibraryName, cmdID, inputs);
            value(3)=bitshift(value(3),16);
            value(2)=bitshift(value(2),8);
            force = bitor(value(3),bitor(value(2),value(1)));
        end
        
        % Power down
        function powerDown(obj)
            cmdID = obj.POWER_DOWN;
            inputs = [obj.Pins(1) obj.Pins(2) obj.Interrupt];
            msg = sendCommand(obj,obj.LibraryName, cmdID, inputs);
            disp(char(msg'));
        end   
        
        % Power up
        function powerUp(obj)
            cmdID = obj.POWER_UP;
            inputs = [obj.Pins(1) obj.Pins(2)]; 
            msg = sendCommand(obj,obj.LibraryName, cmdID, inputs);
            disp(char(msg'));
        end 
        
        % Comunication test (send back the "Hello World" string from
        % Arduino and show it in the command window;
        function comunication_test(obj)
            cmdID = obj.TEST;
            inputs = [];
            msg = sendCommand(obj, obj.LibraryName, cmdID, inputs);
            disp(char(msg'));
        end
        
    end
end

