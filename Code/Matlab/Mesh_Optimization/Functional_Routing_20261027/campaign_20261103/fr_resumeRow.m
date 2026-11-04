function ra = fr_resumeRow(resArr, name)
%FR_RESUMEROW The resume-scan row of one group (default row if absent).
idx = find(strcmp({resArr.name}, name), 1);
if isempty(idx)
    ra = fr_resRow();
else
    ra = resArr(idx);
end
end
