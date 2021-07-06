%This function is written to be used with 2 arduino sketches
% (1) SLoadCell_ZeroFactorSketch, when protocol_id == '1'
% (2) Knee TorqueTest, when protocol_id == '2'
% Minor Adjustments are needed to this function depending on which Arduino 
% sketch it is being used with

function [Data, Stats] = KneeTest(protocol_id, totals) 
total = totals;
    %for protocol_id == '1', total == 150
    %for protocol_id == '2', total needs to be at least 6000 for the system
    %to reach equilibrium


format short g      %remove scientific notation

%Initialize serial port 
s = serialport('COM5',9600);

flush(s,'input');

while s.NumBytesAvailable() < 1
    writeline(s,num2str(protocol_id));
    disp(s.NumBytesAvailable());
end

flush(s,'input');
pause(0.1);

while s.NumBytesAvailable() > 1 
    disp('wait');
end

writeline(s,num2str(total));
disp(total);


svalues = zeros(total,3);

yyaxis left     %graph force on left axis in blue
Force = animatedline('color','blue');
ylim([-50,300]);
ylabel('Force (N)');

yyaxis right    %graph pressure on right in red
Pressure = animatedline('color','red');
ylim([15,700]);
ylabel('Pressure (kPa)');

%the load cell data will sometimes spike unexpectedly, these variables and
% the if statement ignore data that is greater than 100 N from the previous value

prev = 1;   %since the spikes can span over multiple data readings, 
%this variable keeps track of the last data point not in the spike and 
%compares new data against it until it finds one that is within 100 N,
%which indicates that the spike has ended

for i = 1:total
        svalues(i,1) = str2double(readline(s))*4.4482216; %read info from load cell into column 1 and convert from lbs to N                                   
        svalues(i,2) = str2double(readline(s))*395/512-115; %read info from pressure censor into column 2 and convert from analog input (0-1023) to kPa
        svalues(i,3) = str2double(readline(s))/1000; %put a timestamp in column 3 for the information
        % read into columns 1 and 2


    if i > 1
        if abs(svalues(i-prev,1)-svalues(i,1))<100      %check if next value is 
        %more than 100 N away from the last known value
            prev = 1;
            addpoints(Force,svalues(i,3),svalues(i,1)); %if so, add the force data
        else
            svalues(i,1) = svalues(i-prev,1);   %otherwise, replace it
            disp('spike');                      %with the last known data point
            disp(svalues(i,3));
            prev=prev+1;                        %increase the amount of readings
        end                                     %since the last value outside of the spike
                    
        
        addpoints(Pressure,svalues(i,3),svalues(i,2));  %add pressure and time data
        drawnow                                         %update graph
    end
end

addpoints(Force,svalues(i,3),svalues(i,1));     %add data to graph
addpoints(Pressure,svalues(i,3),svalues(i,2));
drawnow  

%save data as function output Data
Data = svalues;

%Use the following bit of code to find some basic statistics about the data
%that has been collected.  The operating assumption here is that collecting 
%6000 data points will be enough for the system to have reached equilibriam,
%and that the last 500 data points will exists outside of the transient
%system response.

% **NOTE: If protocol_id == '1', either comment out the following code,
% or change the range of data on which the basic statistics are calculated.

if (protocol_id == 2)
    
forceData = svalues(1:total,1);
pressureData = svalues(1:total,2);

%Force stats
stats = zeros(6,2);
stats(1,1) = mean(forceData);
stats(2,1) = median(forceData);
stats(3,1) = mode(forceData);
stats(4,1) = min(forceData);
stats(5,1) = max(forceData);
stats(6,1) = std(forceData);

%Pressure stats
stats(1,2) = mean(pressureData);
stats(2,2) = median(pressureData);
stats(3,2) = mode(pressureData);
stats(4,2) = min(pressureData);
stats(5,2) = max(pressureData);
stats(6,2) = std(pressureData);

rows = {'Mean','Median','mode','min','max','Standard Deviation'};
columns = {'Force','Pressure'};

stats = array2table(stats,'RowNames',rows,'VariableNames',columns);
end
Stats = stats;
% disp(Stats);


