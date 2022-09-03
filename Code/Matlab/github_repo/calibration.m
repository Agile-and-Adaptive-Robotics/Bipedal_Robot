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
%% CALIBRATE YOUR LOADCELL WITH ARDUINO
classdef calibration < dynamicprops
    properties(Access = public)
        n = 100               % Number of readings
        known_weight = 1      % Weight used for calibration
    end
    properties(Access = private)
        tare_weight = 0
        scale_factor = 1
    end
    methods(Access = public)
        %% CONSTRUCTOR
        function cal = calibration(varargin)
            %Check the number of inputs
            if nargin < 2
                error('Not enough input arguments.');
            elseif nargin > 2
                error('Too many input arguments.');
            end
            % Check the inputs
            if varargin{1} < 0 
                error('The number of readings must be positive');
            elseif ~isnumeric(varargin{1}) || ~isscalar(varargin{1})
                error('Incorrect type');
            end
            if varargin{2} == 0
                error('The known weight can not be zero');
            elseif ~isnumeric(varargin{2}) || ~isscalar(varargin{2})
                error('Incorrect type');
            end
            cal.n = varargin{1};
            cal.known_weight = varargin{2};
        end
        %% Tare
        function value = tare(varargin)
            if nargin < 2
                error('MyComponent:incorrectType',...
                    'Not enough input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order.');
            elseif nargin > 2
                error('MyComponent:incorrectType',...
                    'Too many input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order.');
            end
            obj = varargin{1};
            HX711_obj = varargin{2};
            value = 1:1:obj.n;
            for i=1:1:obj.n
                value(i) = read_HX711(HX711_obj);
                if check_conversion(obj,value(i))
                    value(i) = NaN;
                end
            end
            value = mean(value,'omitnan');
            obj.tare_weight = value;
        end
        %% Scale factor
        function value = scale(varargin)
            if nargin < 2
                error('MyComponent:incorrectType',...
                    'Not enough input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order.');
            elseif nargin > 2
                error('MyComponent:incorrectType',...
                    'Too many input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order.');
            end
            obj = varargin{1};
            HX711_obj = varargin{2};
            value = 1:1:obj.n;
            for i=1:1:obj.n
                value(i) = read_HX711(HX711_obj);
                if check_conversion(obj,value(i))
                    value(i) = NaN;
                end
            end
            ev_value = mean(value,'omitnan');
            value = (ev_value - obj.tare_weight)/obj.known_weight;
            obj.scale_factor = value;
        end
        %% Get weight
        function weight = get_weight(obj,HX711_obj,varargin)
            if nargin < 2
                error('MyComponent:incorrectType',...
                    'Not enough input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order.');
            elseif nargin > 3
                error('MyComponent:incorrectType',...
                    'Too many input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order. You can also provide the number of readings');
            end
            if ~isempty(varargin)
                if ~isnumeric(varargin{1}) || ~isscalar(varargin{1})
                    error('Incorrect type');
                end
            end
            if isempty(varargin)
                % Raw reading
                val = read_HX711(HX711_obj);
                if check_conversion(obj,val)
                    weight = NaN;
                else
                    weight = (val - obj.tare_weight)/obj.scale_factor;
                end
            else
                % Average of multiple readings
                k = varargin{1};
                value = 1:1:k;
                for i=1:1:k
                    value(i) = read_HX711(HX711_obj);
                    if check_conversion(obj,value(i))
                        value(i) = NaN;
                    end
                end
                ev_value = mean(value,'omitnan');
                weight = (ev_value - obj.tare_weight)/obj.scale_factor;
            end
        end
        %% Statistical function
        function [Av,Std] = stat(obj,HX711_obj,k)
            if nargin < 3
                error('MyComponent:incorrectType',...
                    'Not enough input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order then the number of readings.');
            elseif nargin > 3
                error('MyComponent:incorrectType',...
                    'Too many input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order then the number of readings');
            end
            if ~isnumeric(k) || ~isscalar(k)
                error('Incorrect type');
            end
            weight = 1:1:k;
            for i=1:1:k
                val = read_HX711(HX711_obj);
                weight(i) = (val - obj.tare_weight)/obj.scale_factor;
                if check_conversion(obj,val)
                    weight(i) = NaN;
                end
            end
            Av = mean(weight,'omitnan');
            Std = std(weight,'omitnan');
        end
        %% Plot Function
        function plot_data(obj,HX711_obj,k,varargin)
            if nargin < 3
                error('MyComponent:incorrectType',...
                    'Not enough input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order then the number of readings and the weight (if it is known).');
            elseif nargin > 5
                error('MyComponent:incorrectType',...
                    'Too many input arguments:\nprovide *calibration* and *HX711* %s',...
                    'class objects in this order then the number of readings and the weight (if it is known).');
            end
            if ~isempty(varargin)
                if ~isnumeric(varargin{1}) || ~isscalar(varargin{1})
                    error('Incorrect type');
                end
            end
            [average,std_dev] = stat(obj,HX711_obj,k);
            x = linspace(average-4*std_dev,average+4*std_dev,5000);
            y = (std_dev*sqrt(2*pi))^-1*exp(-(x-average).^2/(2*std_dev^2));
            hold on;
            plot(x,y,'b','LineWidth',2);
            xlim([average-5*std_dev average+8.5*std_dev]);
            ylim([0 max(y)*1.1]);
            if isempty(varargin) == 0
                plot([varargin{1} varargin{1}],ylim,'m');
            end
            plot([average-2*std_dev average-2*std_dev],ylim/2,'r--','LineWidth',1.2);
            plot([average+2*std_dev average+2*std_dev],ylim/2,'r--','LineWidth',1.2);
            plot([average average],ylim,'b--','LineWidth',1.2);
            lgd = legend('$\frac{1}{\sigma\sqrt{2\pi}}\,e^{-\frac{(x-M)^2}{2\sigma^2}}$',...
                'Known weight','$M-2\sigma$','$M+2\sigma$');
            lgd.Interpreter = 'latex';
            lgd.FontSize=16;
            xlabel('Weight');
            title('Normal Distribution');
            hold off;
        end
        %% Check HX711 data
        function incorrect = check_conversion(~,val)
            if (val == hex2dec('800000') || val == hex2dec('7FFFFF'))
                incorrect = true; 
            else
                incorrect = false;
            end
        end
    end
end