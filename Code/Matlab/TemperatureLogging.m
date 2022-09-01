%% Acquire and analyze data from an Arduino
%Original code for one temperature sensor
%Modified code for 5V Pressure sensor and HX711 Load Cell Amplifier

%% Connect to Arduino
% Use the arduino command to connect to an Arduino device.

a = arduino;

%% Take a single temperature measurement
% The datasheet for the TMP36 temperature sensor tells us that the voltage
% reading is directly proportional to temperature in Celsius with an 
% offset of 0.5V and a scale factor of 10 mV/°C (equivalent to 100 °C/V).
% Therefore the conversion can be represented as
%
% $T_C = (V-0.5)*100$
%
% We can read the output voltage, convert it to Celsius and convert the
% result to Farenheit as follows:
v = readVoltage(a,'A0');            %Pressure, kPa
kPa = v*158.05-129.72;
psi = kPa*0.1450377377;
fprintf('Temperature Reading:\n  %.1f kPA\n  %.1f psi\n',kPa,psi)

%% Record and plot 10 seconds of temperature data

ii = 0;
psi = zeros(1e4,1);
t = zeros(1e4,1);

tic
while toc < 10
    ii = ii + 1;
    % Read current voltage value
    v = readVoltage(a,'A0');
    % Calculate temperature from voltage (based on data sheet)
    kPa(ii) = v*158.05-129.72;
    psi(ii) = kPa*0.1450377377;
    % Get time since starting
    t(ii) = toc;
end

% Post-process and plot the data. First remove any excess zeros on the
% logging variables.
kPa = kPa(1:ii);
psi = psi(1:ii);
t = t(1:ii);
% Plot pressure and force versus time
figure
plot(t,psi,'-o')
xlabel('Elapsed time (sec)')
ylabel('Temperature (\kPa)')
title('Ten Seconds of Pressure Data')
set(gca,'xlim',[t(1) t(ii)])

%% Compute acquisition rate

timeBetweenDataPoints = diff(t);
averageTimePerDataPoint = mean(timeBetweenDataPoints);
dataRateHz = 1/averageTimePerDataPoint;
fprintf('Acquired one data point per %.3f seconds (%.f Hz)\n',...
    averageTimePerDataPoint,dataRateHz)

%% Why is my data so choppy?

measurableIncrementV = 5/1023;
measurableIncrementkPa = measurableIncrementV*158.05;
measurableIncrementpsi = measurableIncrementkPa*0.1450377377;
fprintf('The smallest measurable increment of this sensor by the Arduino is\n %-6.4f V\n %-6.2fkPa\n %-6.2fpsi\n',...
    measurableIncrementV,measurableIncrementC,measurableIncrementF);

%% Acquire and display live data

figure
h = animatedline;
ax = gca;
ax.YGrid = 'on';
ax.YLim = [0 600];

stop = false;
startTime = datetime('now');
while ~stop
    % Read current voltage value
    v = readVoltage(a,'A0');
    % Calculate temperature from voltage (based on data sheet)
    TempC = (v - 0.5)*100;
    psi = 9/5*TempC + 32;    
    % Get current time
    t =  datetime('now') - startTime;
    % Add points to animation
    addpoints(h,datenum(t),psi)
    % Update axes
    ax.XLim = datenum([t-seconds(15) t]);
    datetick('x','keeplimits')
    drawnow
    % Check stop condition
    stop = readDigitalPin(a,'D12');
end

%% Plot the recorded data

[timeLogs,tempLogs] = getpoints(h);
timeSecs = (timeLogs-timeLogs(1))*24*3600;
figure
plot(timeSecs,tempLogs)
xlabel('Elapsed time (sec)')
ylabel('Temperature (\circF)')

%% Smooth out readings with moving average filter

smoothedTemp = smooth(tempLogs,25);
tempMax = smoothedTemp + 2*9/5;
tempMin = smoothedTemp - 2*9/5;

figure
plot(timeSecs,tempLogs, timeSecs,tempMax,'r--',timeSecs,tempMin,'r--')
xlabel('Elapsed time (sec)')
ylabel('Temperature (\circF)')
hold on 

%%
% Plot the original and the smoothed temperature signal, and illustrate the
% uncertainty.

plot(timeSecs,smoothedTemp,'r')

%% Save results to a file

T = table(timeSecs',tempLogs','VariableNames',{'Time_sec','Temp_F'});
filename = 'Temperature_Data.xlsx';
% Write table to file 
writetable(T,filename)
% Print confirmation to command line
fprintf('Results table with %g temperature measurements saved to file %s\n',...
    length(timeSecs),filename)
