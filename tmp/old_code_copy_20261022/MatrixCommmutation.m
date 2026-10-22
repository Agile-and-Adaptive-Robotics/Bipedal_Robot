%% Matrix commutation
% Set up rotation matrices Rx, Ry, Rz. Change order of multiplication.
% compare results.

%clean up workspace
clc; clear; close all;

% Create angles alpha, beta, gamma, and distances x, y, z. Make them real
% positive valued.
syms alpha beta gamma x y z real positive

% Call helper functions at the bottom of script
R_z = Rz(gamma);
R_y = Ry(beta);
R_x = Rx(alpha);
p = [x y z]';

R{1} = R_z*R_y*R_x;
R{2} = R_x*R_y*R_z;
R{3} = R_y*R_x*R_z;

pnew = cell(size(R));

for i = 1:length(R)
    pnew{i}=R{i}*p;
    pnew{i}
end

% Helper functions and nested functions are at the bottom of scripts
% These ones are for rotation matrices
function RotZ = Rz(gamm)
    RotZ = [cos(gamm) -sin(gamm) 0; sin(gamm) cos(gamm) 0; 0 0 1];
end

function RotY = Ry(bet)
    RotY = [cos(bet) 0 sin(bet); 0 1 0; -sin(bet) 0 cos(bet)];
end

function RotX = Rx(alph)
    RotX = [1 0 0; 0  cos(alph) -sin(alph); 0 sin(alph) cos(alph)];
end
    