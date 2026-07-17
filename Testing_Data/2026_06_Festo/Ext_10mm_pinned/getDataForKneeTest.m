clc;
clear all;

%% Assuming readserialnumbers() returns a matrix with 4 columns: Time(ms), Position(degrees), Force (N), Pressure (kPa)
data = readserialnumbers(); 
%% save current readings into their own file
timestamp = string(datetime("now", "Format", "yyyy-MM-dd_HH-mm-ss"));

matFileName = "forceData_" + timestamp + ".mat";
csvFileName = "forceData_" + timestamp + ".csv";

% Save raw data variable as MATLAB file

save(matFileName, "data");

% Save data as CSV file too

writematrix(data, csvFileName);

fprintf("Data saved as:\n");
fprintf("%s\n", matFileName);
fprintf("%s\n", csvFileName);


%% Extract columns

timeData = data(:, 1); % given in ms
angleData = data(:,2); % given in degrees
forceData = data(:, 3); % given in N 
pressureData = data(:,4); % given in kPa
%averageForceData = data(:,3);
% converting data
timeDataSeconds=timeData/1000; %converting ms to s
% pressure1DataUpdated = pressure1Data*6.89476; %converting from psig to kpa
% pressure2DataUpdated = pressure2Data*6.89476; %converting from psig to kpa
% pressureDifference = pressure1DataUpdated - pressure2DataUpdated

plot(timeDataSeconds, forceData);
%%
% Sampling rate and cutoff frequency for Butterworth filter
Fs = 1000; % Example: 1000 Hz sample rate
Fc = 60;  % Example: 120 Hz cutoff frequency
order = 12;
%  (12th order)
[b, a] = butter(order, Fc / (Fs / 2), 'low');

% Pre-allocate arrays for storing filtered data
%filteredForceData = zeros(10000,1);
%filteredPressureData = zeros(length(validStartTimes), plot_length);
%
% Apply the filter to each valid segment of forceData and pressureData
%for i = 1:length(timeData)
%     startIdx = validStartTimes(i);
%     pressureStartIdx = startIdx -10; % Adjust the start index for pressure data
%     if pressureStartIdx < 1
%         pressureStartIdx = 1; % Ensure the start index is within bounds
%     end
    filteredForceData = filtfilt(b, a, forceData);
    %filteredPressureData(i, :) = filtfilt(b, a, pressureData(pressureStartIdx:pressureStartIdx+plot_length-1));
%end
%
figure;
plot(timeDataSeconds, filteredForceData);
figure;
plot(timeDataSeconds, forceData);


%% THIS SECTION OF CODE IS FROM OLD STUFF FOR REFERENCE %%  
% %% Assign the entire output data to a variable
% simulatedData_40ms_16Hz = out.statespaceresponse_40ms_16Hz;
% 
% %% Assuming 'Data' is stored as a matrix, you can extract columns like this:
% simulatedtime_40ms_16Hz = t+0.01 ;        % Extract time (assuming it's a separate property)
% %%
% Data1_40ms_16Hz = simulatedData_40ms_16Hz.Data(:,1);  % Extract the first column pressure before
% Data2_40ms_16Hz = simulatedData_40ms_16Hz.Data(:,2);  % Extract the second column pressure after
% Data3_40ms_16Hz = simulatedData_40ms_16Hz.Data(:,3);  % Extract the third column pressure in system
% %% Plot Force vs Time
% figure;
% plot(timeDataSeconds, forceData, timeDataSeconds, valveontimeData);
% title('Force vs Time');
% xlabel('Time (seconds)');
% ylabel('Force (N)');
% grid on;
% 
% %% Plot Pressure1 vs Time
% % 
% % figure;
% % plot(timeDataSeconds, pressure1DataUpdated,timeDataSeconds, pressure2DataUpdated, timeDataSeconds, valveontimeData,simulatedtime_10Hz,Data2_10Hz);
% % title('Pressure (before) and Pressure (after) vs Time');
% % xlabel('Time (seconds)');
% % ylabel('Pressure (kPa)');
% % legend('before BPA','after BPA','pulse timer','simulation');
% % grid on;
% pressureDifference = pressure2Data40ms16Hz-Data1_40ms_16Hz;
% figure;
% plot(timer, pressure2Data40ms16Hz,simulatedtime_40ms_16Hz,Data1_40ms_16Hz, timer, pressureDifference);
% title('Pressure (before) and Pressure (after) vs Time');
% xlabel('Time (seconds)');
% ylabel('Pressure (kPa)');
% legend('after BPA','simulation','difference');
% grid on;
% 
% %% Plot Pressure vs Time
% figure;
% plot(timeDataSeconds, pressure2DataUpdated, timeDataSeconds, valveontimeData);
% title('Pressure2(after) vs Time');
% xlabel('Time (seconds)');
% ylabel('Pressure (kPa)');
% grid on;
% 
% %% plot pulse timer
% figure;
% plot(timeDataSeconds, valveontimeData);
% title('valve time vs Time');
% xlabel('Time (seconds)');
% ylabel('valve open or closed');
% grid on;
% %% to correct for the timer differences
% %% Define a common time vector
% timer= 0:0.001:3.05;
% %
% %% Find the common range between the two time sets and create a new time vector
% minTime = max(min(timer), min(simulatedtime_40ms_16Hz));
% maxTime = min(max(timer), max(simulatedtime_40ms_16Hz));
% commonTimeVector = linspace(minTime, maxTime, length(timer));  % Create a common time vector
% commonTimeVector=commonTimeVector(:);
% % Resample both data sets onto the common time vector
% % the 40 ms pulses
% pressure1Data40ms16Hz_resampled = interp1(timer, pressure1Data40ms16Hz, commonTimeVector, 'linear');
% pressure2Data40ms16Hz_resampled = interp1(timer, pressure2Data40ms16Hz, commonTimeVector, 'linear');
% simulatedPressure40ms16Hz_resampled = interp1(simulatedtime_40ms_16Hz, output_40ms(:,1), commonTimeVector, 'linear');
% pressure1Data40ms16Hz_resampled = pressure1Data40ms16Hz_resampled(:);
% pressure2Data40ms16Hz_resampled = pressure2Data40ms16Hz_resampled(:);
% simulatedPressure40ms16Hz_resampled = simulatedPressure40ms16Hz_resampled(:);
% 
% % Calculate the pressure difference
% pressure2Difference40ms16Hz = abs(pressure2Data40ms16Hz_resampled - output_40ms(:,1));%simulatedPressure12Hz_resampled
% pressure1Difference40ms16Hz = abs(pressure1Data40ms16Hz_resampled - output_40ms(:,2));%simulatedPressure12Hz_resampled
% percentage_Difference_Pre_40ms_16Hz= (pressure1Difference40ms16Hz ./ pressure1Data40ms16Hz_resampled)*100;
% percentage_Difference_Post_40ms_16Hz= (pressure2Difference40ms16Hz ./ pressure2Data40ms16Hz_resampled)*100;
% % Plot the results
% figure;
% % plot(commonTimeVector, pressure1Data12Hz_resampled, commonTimeVector, simulatedPressure12Hz_resampled, commonTimeVector, pressureDifference12Hz);
% plot(commonTimeVector-0.01, pressure2Data40ms16Hz_resampled,commonTimeVector-0.01,output_40ms(:,1)/1.32, commonTimeVector-0.01, pressure2Difference40ms16Hz ,commonTimeVector-0.01, percentage_Difference_Post_40ms_16Hz)
% title('Pressure (after) vs Time (Resampled)');
% xlabel('Time (seconds)');
% ylabel('Pressure (kPa)');
% legend('after BPA', 'simulation (resampled)', 'difference','percentage difference');
% grid on;
% %
% figure;
% plot(commonTimeVector, pressure1Data40ms16Hz_resampled,commonTimeVector,output_40ms(:,2), commonTimeVector, pressure1Difference40ms16Hz,commonTimeVector, percentage_Difference_Pre_40ms_16Hz);
% title('Pressure (before)  vs Time (Resampled)');
% xlabel('Time (seconds)');
% ylabel('Pressure (kPa)');
% legend('before BPA', 'simulation (resampled)', 'difference', 'percentage difference');
% %
% % figure;
% % plot(commonTimeVector, percentage_Difference_Post);
% %%
% figure;
% plot(commonTimeVector, pressure2Data40ms16Hz_resampled,commonTimeVector,output(:,1))
