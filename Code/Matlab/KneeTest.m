%This function is written to be used with 2 arduino sketches
% (1) SLoadCell_ZeroFactorSketch, when protocol_id == '1'
% (2) Knee TorqueTest, when protocol_id == '2'
% Minor Adjustments are needed to this function depending on which Arduino 
% sketch it is being used with

function [Data, Stats] = KneeTest(protocol_id)

%Initialize serial port 
s = serial('COM4','Baudrate', 9600);

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
  %for protocol_id == '2', total needs to be at least 6000 for the system
  %to reach equilibrium
svalues = cell (total,3);

tic %start timer
for i = 1:total

       svalues{i,1} = fgetl(s); %read info from load cell into column 1
       svalues{i,2} = fget1(2); %read info from pressure censor into column 2
       svalues{i,3} = toc; %put a timestamp in column 3 for the information
       % read into columns 1 and 2
end


%the following while loop converts the ascii values from column 1 of
%svalues into a numeric array
a = 1;
while ~isempty(svalues{a}) && a<total
        data(a,1) = str2num(svalues{a,1})*4.4482216; %%Convert lb to N
        data(a,2) = str2num(svalues{a,2});
        data(a,3) = svalues{a,3};
        a = a+1;
end


%save data as function output Data
Data = data;

%Create a Force vs. time plot with the information in Data
Force = data(:,1);
X = data(:,3);
yyaxis left
xlabel('Time (sec)');
ylabel('Force(N)');
plot(X,Force,'.');

%Graph pressure against right y axis
Pressure = data(:,2);
yyaxis right
ylabel('Pressure( )');
plot(X,Pressure,'.');

%Use the following bit of code to find some basic statistics about the data
%that has been collected.  The operating assumption here is that collecting 
%6000 data points will be enough for the system to have reached equilibriam,
%and that the last 500 data points will exists outside of the transient
%system response.

% **NOTE: If protocol_id == '1', either comment out the following code,
% or change the range of data on which the basic statistics are calculated.

%Force stats
Stats = zeros(4,2);
Stats(1,1) = mean(data(5500:5999,1));
Stats(2,1) = median(data(5500:5999,1));
Stats(3,1) = mode(data(5500:5999,1));
Stats(4,1) = std(data(5500:5999,1));

%Pressure stats
Stats(1,2) = mean(data(5500:5999,2));
Stats(2,2) = median(data(5500:5999,2));
Stats(3,2) = mode(data(5500:5999,2));
Stats(4,2) = std(data(5500:5999,2));
Stats;

fclose(s)
