%This is an introductory function for reading in serial data

% function data = readserialnumbers()
% 
% %Initialize the serial port on the correct port, with a baud rate
% s = serialport('COM6', 500000);
% 
% %Determine how many lines of data to collect, and initialize a cell array
% %to store in the data from the serial port
% % when sending in arduino data, seperate variables using a serial print comma and
% % restart the line with final variable using serial println.  
% % i.e. Serial.print(timer);
% %    Serial.print(",");
% %    Serial.print(torqueNm, 3);
% %    Serial.print(",");
% %    Serial.println(torqueAvg, 3);
% length_to_collect = 285000; %number of points to collect
% stext = cell(length_to_collect,1);
% 
% %% Read the data and store in a matrix
% %Read all the data (ASCII) into a cell array
% for i = 1:length_to_collect+1  %First read almost always starts in the middle of a number (Or is 'arduino is ready'), so it is likely not useable, so add 1
%     stext{i} = readline(s);
% end
% 
% %Convert the data from ASCII to numbers
% a = 2; %initialize counting variable (first read is always 'junk')
% b = 1; %initialize storing variable
% while a<=length_to_collect+1 %While a is less than or equal to the amount of data we need to collect
%     if(~isempty(str2num(stext{a})))  %Make sure the data we are reading is actually a number
%         data(b,:) = str2num(stext{a});  %If it is a number, store it
%         b = b+1;  %Increment storage variable
%     end
%     a = a+1;   %Increment read variable
% end
% 
% %%
% %close the port before ending the function
% clear s
% This function reads serial data from Arduino and stores it in a numeric matrix

function data = readserialnumbers()

    % Serial settings
    port = 'COM12';
    baudRate = 115200;

    % Number of valid data rows to collect
    length_to_collect = 200;

    % Expected number of comma-separated values per line
    expectedNumCols = 4;

    % Open serial port
    s = serialport(port, baudRate);
    configureTerminator(s, "LF");
    s.Timeout = 10;

    % Clear old serial data from the buffer
    flush(s);

    % Preallocate data matrix
    data = NaN(length_to_collect, expectedNumCols);

    % Read and discard first line because it is often incomplete junk
    readline(s);

    % Counter for valid stored rows
    b = 1;

    while b <= length_to_collect

        % Read one line from serial
        lineText = strtrim(readline(s));

        % Split line by commas
        parts = split(lineText, ",");

        % Only process lines with the expected number of values
        if numel(parts) == expectedNumCols

            % Convert text values to numbers
            numericValues = str2double(parts).';

            % Only store if all values converted correctly
            if all(~isnan(numericValues))
                data(b, :) = numericValues;
                b = b + 1;
            else
                disp("Skipped bad numeric line: " + lineText);
            end

        else
            disp("Skipped wrong-size line: " + lineText);
        end
    end

    % Close serial connection
    clear s

end