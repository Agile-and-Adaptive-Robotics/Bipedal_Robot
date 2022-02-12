
%This function is written to be used with 2 arduino sketches
%(1) SLoadCell_ZeroFactorSketch, when protocol_id == 1
%(2) Knee TorqueTest, when protocol_id == 2
%if protocol_id == 2, then the 'total' variable is how many data points
%to collect (should be at least 10000 to get reliable data)

function [Data, Stats] = KneeTest(protocol_id,port,varargin) 

    %Allows for the "total" argument be left blank if you are not using
    %KneeTorqueTest.ino, if left blank it will default to 10000
    p = inputParser;
    addRequired(p,'protocol_id');
    addRequired(p,'port');
    addOptional(p,'total',1000);
    parse(p,protocol_id,port,varargin{:});
    
    %Creates variables equal to the argument values
    protocol_id = num2str(p.Results.protocol_id);
    total = p.Results.total;
    
    %Creates a string in the format "COMx" where x is the port number to be
    %able to connect to the arduino
    port = strjoin({'COM',num2str(p.Results.port)},'');

    
    %remove scientific notation
    format short g

    %Initialize serial port 
    s = serialport(port,57600);

    while s.NumBytesAvailable < 1               
        write(s,string(protocol_id),'string');
    end
    %the arduino sends sends the string "running" when it receives the
    %protocol_id, which increases the value of s.NumBytesAvailable and 
    %breaks the loop

    write(s,string(total),'string')
    readline(s)
    %writes the amount of data points to collect, reads once to clear
    %the "running" string from the buffer before reading the data

    svalues = zeros(total,3);   %creates array for the data
    clf;                        %clears graph from any previous tests

    yyaxis left     %graph force on left axis in blue
    Force = animatedline('color','blue');
    ylim([0,50]);
    ylabel('Force (N)');

    yyaxis right    %graph pressure on right in red
    Pressure = animatedline('color','red');
    ylim([0,620]);
    ylabel('Pressure (kPa)');

%     the load cell data will sometimes spike unexpectedly which
%     causes problems in data collection

    prev = 1;   
%     since the spikes can span over multiple data readings, 
%     this variable keeps track of the last data point not in the spike 
%     and compares new data against it until it finds one that is within
%     100 N, which indicates that the spike has ended

    for i = 1:total
        
        svalues(i,1) = str2double(readline(s))*4.4482216;                                
        svalues(i,2) = str2double(readline(s))*153.62*(5/1023)-133.21;
        svalues(i,3) = str2double(readline(s))/1000;

%         read data to each column and convert units when needed
%         column 1 is force, converting lbs to N
%         column 2 is pressure, converting analog value (1-1023) to
%         kPa
%         column 3 is the time that the data was collected,
%         converting milliseconds to seconds

        if i > 1
            
            %i > 1 because the spike detection references a previous i
            %value, which doesn't exist if the loop started at the minimum
            %i value
            
                if abs(svalues(i-prev,1)-svalues(i,1))<200     
                    
                    %check if next value is more than
                    %100 N away from the last known value
                    prev = 1;
                    addpoints(Force,svalues(i,3),svalues(i,1)); 
                    %if so, add the force data
                    
                else
                    
                    svalues(i,1) = svalues(i-prev,1);   
                    %otherwise, replace it with the last known data point
                    
                    disp('spike');                      
                    disp(svalues(i,3));
                    prev = prev+1;                        
                    %increase the amount of points since 
                    %the last value outside of the spike
                    
                end                                     

            addpoints(Pressure,svalues(i,3),svalues(i,2));  
            %add pressure and time data
            
        end
    end

    addpoints(Force,svalues(i,3),svalues(i,1));     %add data to graph
    addpoints(Pressure,svalues(i,3),svalues(i,2));
    drawnow  

    %Use the following bit of code to find some basic statistics about the data
    %that has been collected.  The operating assumption here is that collecting 
    %6000 data points will be enough for the system to have reached equilibriam,
    %and that the last 500 data points will exists outside of the transient
    %system response.

    forceData = svalues(500:total,1);
    pressureData = svalues(500:total,2);

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

    %Makes table for stats
    rows = {'Mean','Median','mode','min','max','Standard Deviation'};
    columns = {'Force','Pressure'};
    stats = array2table(stats,'RowNames',rows,'VariableNames',columns);

    Data = svalues;
    Stats = stats;

end

