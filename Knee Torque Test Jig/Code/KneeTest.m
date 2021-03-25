%This function is written to be used with 2 arduino sketches
% (1) SLoadCell_ZeroFactorSketch, when protocol_id == '1'
% (2) Knee TorqueTest, when protocol_id == '2'
% Minor Adjustments are needed to this function depending on which Arduino 
% sketch it is being used with

function [Data, Stats] = KneeTest(protocol_id)

%Initialize serial port 
s = serial('COM3', 'Baudrate', 9600);

%Open serial port s
fopen(s);

%Send 'protocol_id' over serial until there is data on the 
%incoming buffer.  It was expected that only the fwrite(s, protocol_id)
%command should be needed, but when this is the case, there seems to a 
%timing issue with Arduino. The While statement was constructed as a work 
% around
while (s.BytesAvailable == 0)
    fwrite(s, protocol_id)
end


%prepare a cell array to receive ASCii data from the incoming buffer
total = 6000;  
  %for protocol_id == '1', total == 150
  %for protocol_id == '3', total needs to be at least 6000 for the system
  %to reach equilibrium
svalues = cell (total,2);

tic %start timer
for i = 1:total

       svalues{i,1} = fgetl(s); %read information from buffer into column 1
       svalues{i,2} = toc; %put a timestamp in column 2 for the information
       % read into column 1 
end


%the following while loop converts the ascii values from column 1 of
%svalues into a numeric array
a = 1;
while ~isempty(svalues{a}) && a<total
        data(a,1) = str2num(svalues{a})*4.4482216; %%Convert lb to N
        data(a,2) = svalues{a,2};
        a = a+1;
end


%save data as function output Data
Data = data;

%Create a Force vs. time plot with the information in Data
Y = data(:,1);
X = data(:,2);
plot(X,Y,'.')
xlabel('Time (sec)')
ylabel('Force(N)')

%Use the following bit of code to find some basic statistics about the data
%that has been collected.  The operating assumption here is that collecting 
%6000 data points will be enough for the system to have reached equilibriam,
%and that the last 500 data points will exists outside of the transient
%system response.

% **NOTE: If protocol_id == '1', either comment out the following code,
% or change the range of data on which the basic statistics are calculated.
Stats = zeros(4,1);
Stats(1,1) = mean(data(5500:5999,1));
Stats(2,1) = median(data(5500:5999,1));
Stats(3,1) = mode(data(5500:5999,1));
Stats(4,1) = std(data(5500:5999,1));
Stats;

fclose(s)
