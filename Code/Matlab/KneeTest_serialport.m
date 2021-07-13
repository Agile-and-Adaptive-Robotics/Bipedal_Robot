%This function is written to be used with 2 arduino sketches
% (1) SLoadCell_ZeroFactorSketch, when protocol_id == '1'
% (2) Knee TorqueTest, when protocol_id == '2'
% Minor Adjustments are needed to this function depending on which Arduino 
% sketch it is being used with

function [Data, Stats] = KneeTest_serialport(protocol_id)


%Initialize serial port 
s = serialport("COM8",9600);
flush(s);
s.UserData = struct("Force",[],"Time",[],"Count",1);
total = 1500; 
% svalues = cell (total, 2);

% s = serial('COM4','Baudrate', 9600);
%Open serial port s
%fopen(s);

%Send 'protocol_id' over serial until there is data on the 
%incoming buffer.  It was expected that only the fwrite(s, protocol_id)
%command should be needed, but when this is the case, there seems to a 
%timing issue with Arduino. The While statement was constructed as a work 
% around
while (s.NumBytesAvailable == 0)
    write(s, protocol_id, "char")
end

%prepare a cell array to receive ASCii data from the incoming buffer
 
  %for protocol_id == '1', total == 150
  %for protocol_id == '2', total needs to be at least 6000 for the system
  %to reach equilibrium

tic %start timer

configureCallback(s,"terminator",@readDueData);

Force = src.UserData.Force*4.4482216; %%Convert lb to N
Time = src.UserData.Time;

%save data as function output Data
Data = [Force, Time];


% for i = 1:total
% 
%        svalues{i,1} = readDueData(s); %read information from buffer into column 1
%        svalues{i,2} = toc; %put a timestamp in column 2 for the information
%                             %  read into column 1 
% end
% 
% 
% % the following while loop converts the ascii values from column 1 of
% % svalues into a numeric array
% data = cell (total,2);
% a = 1;
% 
% while ~isempty(svalues{a}) && a<total
%         data(a,1) = src.UserData.Force(end+1)*4.4482216; %%Convert lb to N
%         data(a,2) = toc;
%         a = a+1;
% end

% %Create a Force vs. time plot with the information in Data
% Y = data(:,1);
% X = data(:,2);
% plot(X,Y,'.')
% xlabel('Time (sec)')
% ylabel('Force(N)')

%Use the following bit of code to find some basic statistics about the data
%that has been collected.  The operating assumption here is that collecting 
%6000 data points will be enough for the system to have reached equilibriam,
%and that the last 500 data points will exists outside of the transient
%system response.

% **NOTE: If protocol_id == '1', either comment out the following code,
% or change the range of data on which the basic statistics are calculated.
beg = total-500;
fin = total -1;
Stats = zeros(4,1);
Stats(1,1) = mean(Force(beg:fin,1));
Stats(2,1) = median(Force(beg:fin,1));
Stats(3,1) = mode(Force(beg:fin,1));
Stats(4,1) = std(Force(beg:fin,1));
% Stats;

function readDueData(src, ~)

% Read the ASCII data from the serialport object.
force = readline(src);

% Convert the string data to numeric type and save it in the UserData
% property of the serialport object.
src.UserData.Force(end+1) = str2double(force);

% % Update Time on clock
 src.UserData.Time(end+1) = toc;

% Update the Count value of the serialport object.
src.UserData.Count = src.UserData.Count + 1;

% If 1001 data points have been collected from the Arduino, switch off the
% callbacks and plot the data.
if src.UserData.Count > Total+1
    configureCallback(src, "off");
    plot(src.UserData.Time(2:end),src.UserData.Force(2:end),'.')
    xlabel('Time (sec)')
    ylabel('Force(N)');
end
end


end
