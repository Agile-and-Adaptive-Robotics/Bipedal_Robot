function fr_glog2(logFile, line)
%FR_GLOG2 Append one line to a log file (single writer per file).
fid = fopen(logFile, 'a');
if fid == -1
    return
end
fprintf(fid, '%s\n', line);
fclose(fid);
end
