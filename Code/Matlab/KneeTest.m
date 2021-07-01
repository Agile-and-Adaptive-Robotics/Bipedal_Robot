%This function is written to be used with 2 arduino sketches
% (1) SLoadCell_ZeroFactorSketch, when protocol_id == '1'
% (2) Knee TorqueTest, when protocol_id == '2'
% Minor Adjustments are needed to this function depending on which Arduino 
% sketch it is being used with

function [Data, Stats] = KneeTest(protocol_id)
format short g
%Initialize serial port 
s = serialport('COM5',9600);

%Open serial port s

%Send 'protocol_id' over serial until there is data on the 
%incoming buffer.  It was expected that only the fwrite(s, protocol_id)
%command should be needed, but when this is the case, there seems to a 
%timing issue with Arduino. The While statement was constructed as a work 
% around

while s.NumBytesAvailable < 1
    write(s, num2str(protocol_id),'string');
end



%prepare a cell array to receive ASCii data from the incoming buffer
total = 2000;  
  %for protocol_id == '1', total == 150
  %for protocol_id == '2', total needs to be at least 6000 for the system
  %to reach equilibrium
  
svalues = zeros(total,3);

yyaxis left     %graph force on left axis in blue
Force = animatedline('color','blue');
ylim([-2,300]);

yyaxis right    %graph pressure on right in red
Pressure = animatedline('color','red');
ylim([15,700]);

for i = 1:total

       svalues(i,1) = str2double(readline(s))*4.4482216; %read info from load cell into column 1 and convert from lbs to N
       svalues(i,2) = str2double(readline(s))*395/512-115; %read info from pressure censor into column 2 and convert from analog input (0-1023) to kPa
       svalues(i,3) = str2double(readline(s))/1000; %put a timestamp in column 3 for the information
       % read into columns 1 and 2
       
       addpoints(Force,svalues(i,3),svalues(i,1));     %add data to graph
       addpoints(Pressure,svalues(i,3),svalues(i,2));
       drawnow                                              %update graph
       
       
end


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
%Force stats
Stats = zeros(4,2);
Stats(1,1) = mean(svalues(1:total,1));
Stats(2,1) = median(svalues(1:total,1));
Stats(3,1) = mode(svalues(1:total,1));
Stats(4,1) = std(svalues(1:total,1));

%Pressure stats
Stats(1,2) = mean(svalues(1:total,2));
Stats(2,2) = median(svalues(1:total,2));
Stats(3,2) = mode(svalues(1:total,2));
Stats(4,2) = std(svalues(1:total,2));
end
Stats;


