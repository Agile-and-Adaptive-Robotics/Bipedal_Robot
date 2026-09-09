%Dig_FlxBio_handtune.m — hand-tuned Xi guesses evaluated on the biomimetic flexor
%(straight calc), scored identically to Dig_FlxBio_dubfilt (325kPa FVU excluded from
%the score because that baseline is degenerate). Prints the hand table and the
%combined top 8 (hand guesses + the 329 GA-mined unique candidates).
cd('D:/GitHub/Bipedal_Robot/Testing_Data/2022_02_Festo');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Functions/ModernRobotics');
addpath('D:/GitHub/Bipedal_Robot/Code/Matlab/Robot_Data');

x0s  = [0.004, 0.005, 0.0065, 0.0075, 0.009, 0.010];          %Xi0: zero-torque crossing sits near here
x12  = [1e5, 5e5, 1e6, 2e6, 5e6, 2e7; ...                     %Xi1 guesses
        1e4, 1.5e4, 2e4, 1.6e4, 1.5e4, 1.6e4];                %Xi2 guesses (Xi1 > 10*Xi2 zone covered)
npair = size(x12, 2);
hand = zeros(numel(x0s)*npair, 3); hIdx = 0;
for i0 = 1:numel(x0s)
    for p = 1:npair
        hIdx = hIdx + 1;
        hand(hIdx,:) = [x0s(i0), x12(1,p), x12(2,p)];
    end
end

[b10, b620, b325] = minimizeFlx(0, Inf, Inf);
bR = [b10(1), b620(1), b325(1)]; bF = [b10(2), b620(2), b325(2)];

nh = size(hand,1);
hscore = Inf(nh,1); hmet = nan(nh,9);
for i = 1:nh
    try
        [u, v620, v325] = minimizeFlx(hand(i,1), hand(i,2), hand(i,3));
        hmet(i,:) = [u(1:3), v620(1:3), v325(1:3)];
        R = hmet(i,[1 4 7]) ./ bR; V = hmet(i,[2 5 8]) ./ bF;
        hscore(i) = mean([R, V(1:2)]);
    catch
        hscore(i) = Inf;
    end
end

[~, hord] = sort(hscore, 'ascend');
fprintf('===== HAND-TUNED GUESSES (%d), best first — all 3 GoF per case =====\n', nh);
fprintf('%-5s %8s %9s %9s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %6s\n', ...
    'rank','Xi0(m)','Xi1','Xi2','R10','F10','M10','R620','F620','M620','R325','F325','M325','score');
for q = 1:nh
    i = hord(q);
    fprintf('%-5d %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f\n', ...
        q, hand(i,1), hand(i,2), hand(i,3), hmet(i,1), hmet(i,2), hmet(i,3), ...
        hmet(i,4), hmet(i,5), hmet(i,6), hmet(i,7), hmet(i,8), hmet(i,9), hscore(i));
end

%% Combined top 8: hand guesses + GA-mined unique candidates
Dg = load('Dig_out/Dig_FlxBio_dubfilt_results.mat');
d  = Dg.dub_filt_results; met = Dg.met; n = size(d,1);
score = Inf(n,1);
for i = 1:n
    R = met(i,[1 4 7]) ./ bR; V = met(i,[2 5 8]) ./ bF;
    score(i) = mean([R, V(1:2)]);
end
allXi  = [d; hand];
allSc  = [score; hscore];
allMet = [met; hmet];
allSrc = [repmat({"GA"},n,1); repmat({"hand"},nh,1)];
[~, ord] = sort(allSc, 'ascend');
fprintf('\n===== COMBINED TOP 8 (GA + hand) =====\n');
fprintf('%-5s %5s %8s %9s %9s | %5s %5s %5s | %5s %5s %5s | %5s %5s %5s | %6s\n', ...
    'rank','src','Xi0(m)','Xi1','Xi2','R10','F10','M10','R620','F620','M620','R325','F325','M325','score');
shown = 0;
for q = 1:numel(allSc)
    i = ord(q);
    if ~isfinite(allSc(i)), continue; end
    shown = shown + 1;
    fprintf('%-5d %-5s %8.4f %9.3g %9.3g | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %5.2f %5.2f %5.2f | %6.3f\n', ...
        q, allSrc{i}, allXi(i,1), allXi(i,2), allXi(i,3), ...
        allMet(i,1), allMet(i,2), allMet(i,3), allMet(i,4), allMet(i,5), allMet(i,6), ...
        allMet(i,7), allMet(i,8), allMet(i,9), allSc(i));
    if shown == 8, break; end
end
save('Dig_out/Dig_FlxBio_handtune_results.mat', 'hand', 'hscore', 'hmet');
fprintf('\nSaved hand results to Dig_out\\Dig_FlxBio_handtune_results.mat\n');
