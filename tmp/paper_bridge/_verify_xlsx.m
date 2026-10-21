files = {'PinnedFlexor_Preopt_GoF.xlsx','PinnedFlexor_GoF.xlsx','PinnedExtensor_GoF.xlsx','AdoptedCorrectionTerms.xlsx','Route_PinnedBrackets.xlsx','Route_Biomimetic10mm.xlsx','Route_Redesign20mm.xlsx'};
base = 'D:/GitHub/Bipedal_Robot/Documentation/Reports and Papers/Knee_Torque_Test/Tables/';
for f = 1:numel(files)
  fprintf('==================== %s\n', files{f});
  try
    raw = readcell([base files{f}], 'TextType', 'string');
    [r, c] = size(raw);
    fprintf('-- %dx%d\n', r, c);
    for i = 1:min(r, 13)
      row = '';
      for j = 1:min(c, 8)
        v = raw{i,j};
        if ismissing(string(v)), v = ''; end
        row = [row char(string(v)) ' | '];
      end
      if strtrim(row) ~= '', fprintf('%s\n', row); end
    end
  catch ME
    fprintf('ERROR: %s\n', ME.message);
  end
end
