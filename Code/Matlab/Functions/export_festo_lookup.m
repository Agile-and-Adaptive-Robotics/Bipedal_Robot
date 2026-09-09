% Export the Festo normalized-force lookup tables (f_10, f20, f40) from
% FestoLookup.mat to CSV grids for the Python/MuJoCo BPA muscle port.
%
% festo4.m loads these and calls them as f(rel, P) where rel = relative
% strain (contraction / KMAX) and P = pressure normalized by 620 kPa.
% This script inspects what the objects are and exports dense evaluation
% grids so Python can interpolate identically.

cd('C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/Functions');
S = load('FestoLookup.mat');
disp(fieldnames(S));

outdir = 'C:/Users/Ben/Documents/GitHub/Bipedal_Robot/Code/Matlab/Functions';

names = fieldnames(S);
json = '{';
for k = 1:numel(names)
    obj = S.(names{k});
    fprintf('== %s: class=%s\n', names{k}, class(obj));
    switch class(obj)
        case 'sfit'
            % Closed-form surface: port the formula + coefficients
            cn = coeffnames(obj);
            cv = coeffvalues(obj);
            fprintf('   formula: %s\n', formula(obj));
            coeffs = '';
            for j = 1:numel(cn)
                fprintf('   %s = %.10g\n', cn{j}, cv(j));
                if j > 1, coeffs = [coeffs ', ']; end %#ok<AGROW>
                coeffs = [coeffs sprintf('"%s": %.12g', cn{j}, cv(j))]; %#ok<AGROW>
            end
            json = [json sprintf('%s"%s": {"formula": "%s", "coeffs": {%s}}', ...
                mergeStr(k), names{k}, formula(obj), coeffs)]; %#ok<AGROW>
        case 'griddedInterpolant'
            gv = obj.GridVectors;
            [REL, P] = ndgrid(gv{1}, gv{2});
            M = [REL(:), P(:), reshape(obj.Values, [], 1)];
            fn = fullfile(outdir, sprintf('festo_lookup_%s.csv', names{k}));
            writematrix(M, fn);
            fprintf('   wrote %s (%d rows)\n', fn, size(M, 1));
        otherwise
            fprintf('   unhandled class, skipping\n');
    end
end
json = [json '}'];
fid = fopen(fullfile(outdir, 'festo_lookup_coeffs.json'), 'w');
fwrite(fid, json); fclose(fid);
fprintf('wrote festo_lookup_coeffs.json\n');
disp('EXPORT DONE');

function s = mergeStr(k)
    if k > 1, s = ', '; else, s = ''; end
end

