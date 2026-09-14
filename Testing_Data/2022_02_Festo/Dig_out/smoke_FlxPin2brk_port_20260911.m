%% smoke_FlxPin2brk_port_20260911.m
% Equivalence check for the minimizeFlxPin two-bracket port (2026-09-11).
% minimizeFlxPin now carries the minimizeFlxPin2brk mechanics with the
% two-rotation frames fixed (K = [X1,X2,X1], K2 = [X1,X1,X2]) and the +5.3 deg
% encoder offset on test 4 (40cm-tendon) instead of test 3 (47cm).
%
% Expected:
%   - prediction curves M_p identical on ALL tests (port vs 2brk, same Xi,
%     both finite-Xi and the (0,Inf,Inf) escape)
%   - GoF identical on tests 1, 2, 5 (48/46/41 cm)
%   - GoF DIFFERENT on tests 3 and 4 (offset attribution moved)
%   - Aexp identical on tests 1, 2, 5; differs on 3 and 4

here = fileparts(mfilename('fullpath'));              % ...\Dig_out
root = fileparts(fileparts(fileparts(here)));         % ...\Bipedal_Robot
addpath(genpath(fullfile(root, 'Code', 'Matlab')));
cd(fullfile(root, 'Testing_Data', '2022_02_Festo'));
addpath(pwd);   % so auto-created pool workers find the data mats via path

Xi = [0.0089, 5.62e4, 1.85e4];   % settled 2trans_noT3 pick 107 (dissertation)

[fA, bA]  = minimizeFlxPin(Xi(1), Xi(2), Xi(3));
[fB, bB]  = minimizeFlxPin2brk(Xi(1), Xi(2), Xi(3), [], true, '2trans');
[fA0, bA0] = minimizeFlxPin(0, Inf, Inf);
[fB0, bB0] = minimizeFlxPin2brk(0, Inf, Inf, [], true, '2trans');

lbl = {'48cm','46cm','47cm','40cm-t','41cm'};
fprintf('\n=== finite Xi [%.4f, %.2e, %.2e]: port vs 2brk GoF ===\n', Xi);
for j = 1:5
    fprintf('%-6s port [%6.3f %6.4f %6.3f]  2brk [%6.3f %6.4f %6.3f]  dMax %.3g\n', ...
        lbl{j}, fA(j,:), fB(j,:), max(abs(fA(j,:) - fB(j,:))));
end
fprintf('\n=== escape (0, Inf, Inf): port vs 2brk GoF ===\n');
for j = 1:5
    fprintf('%-6s port [%6.3f %6.4f %6.3f]  2brk [%6.3f %6.4f %6.3f]  dMax %.3g\n', ...
        lbl{j}, fA0(j,:), fB0(j,:), max(abs(fA0(j,:) - fB0(j,:))));
end

dP = zeros(5,1); dP0 = zeros(5,1);
for j = 1:5
    dP(j)  = max(abs(bA(j).M_p(:)  - bB(j).M_p(:)),  [], 'omitnan');
    dP0(j) = max(abs(bA0(j).M_p(:) - bB0(j).M_p(:)), [], 'omitnan');
end
fprintf('\nmax |M_p port - M_p 2brk| per test (finite Xi): %s\n', sprintf('%.3g ', dP));
fprintf('max |M_p port - M_p 2brk| per test (escape):    %s\n', sprintf('%.3g ', dP0));

fprintf('\nAexp attribution (max |port - 2brk|, deg; nonzero = offset moved):\n');
for j = 1:5
    fprintf('  test %d (%s): %.4g\n', j, lbl{j}, max(abs(bA(j).Aexp - bB(j).Aexp)));
end

ok = all(dP < 1e-12) && all(dP0 < 1e-12) ...
  && max(abs(fA([1 2 5],:) - fB([1 2 5],:)), [], 'all') < 1e-12 ...
  && max(abs(fA0([1 2 5],:) - fB0([1 2 5],:)), [], 'all') < 1e-12 ...
  && max(abs(fA(3,:) - fB(3,:)), [], 'all') > 1e-9 ...
  && max(abs(fA(4,:) - fB(4,:)), [], 'all') > 1e-9;
fprintf('\nSMOKE RESULT: %s\n', string(ok));
