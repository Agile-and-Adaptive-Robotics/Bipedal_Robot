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
%% CALIBRATE YOUR LOADCELL WITH ARDUINO
classdef calibration 
    properties
        n                  % Number of readings
        known_weight       % Weight used for calibration
        tare_weight = 0
        scale_factor = 1
    end   
    methods(Access = public)
        %% CONSTRUCTOR
        function cal = calibration(varargin)
            if nargin < 2
                error('Not enough input arguments.');
            elseif nargin > 2
                error('Not enough input arguments.');
            end 
            % The number of readings must be positive
            if varargin{1} < 0
                error('The number of readings must be positive');
            end
            cal.n = varargin{1};
            cal.known_weight = varargin{2};
        end       
        %% Tare
        function tare_value = tare(varargin)
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
            end
            tare_value = mean(value);
        end
        %% Scale factor
        function scale_value = scale(varargin)
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
            end
            ev_value = mean(value);
            scale_value = (ev_value - obj.tare_weight)/obj.known_weight;
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
            if isempty(varargin)==1
                % Raw reading
                weight = (read_HX711(HX711_obj) - obj.tare_weight)/obj.scale_factor;
            else
                % Average of multiple readings
                k = varargin{1};
                value = 1:1:k;
                for i=1:1:k
                    value(i) = read_HX711(HX711_obj);
                end
                ev_value = mean(value);
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
            weight = 1:1:k;
            for i=1:1:k
                weight(i) = (read_HX711(HX711_obj) - obj.tare_weight)/obj.scale_factor;
            end
            Av = mean(weight);
            Std = std(weight);
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
            lgd = legend('$\frac{1}{\sigma\sqrt{2\pi}}\,e^{-\frac{(x-M)^2}{2\sigma^2}}$',...
                'Known weight','$M-2\sigma$','$M+2\sigma$');
            lgd.Interpreter = 'latex';
            lgd.FontSize=16;
            xlabel('Weight');
            title('Normal Distribution');
            hold off;
        end
    end
end