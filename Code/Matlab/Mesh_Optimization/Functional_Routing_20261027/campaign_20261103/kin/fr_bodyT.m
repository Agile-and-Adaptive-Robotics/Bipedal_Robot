function T = fr_bodyT(kin, q, body)
%FR_BODYT Body-in-ground transform at q (4x4, column convention):
% p_ground = T(1:3,1:3)*p_body + T(1:3,4).
% Composition convention: see fr_kinematics header (validated against the
% geometry extraction's bodyTransformsInGround).
if strcmp(body, 'ground')
    T = eye(4);
    return
end
S = kin.chain.(body);
J = kin.joints.(S.jointName);
Tj = fr_jointT(kin, q, J);                                    % off_P <- off_C
% the joint's OWN scoped offset frames (names collide across joints)
Xpo = fr_offInBody(J.frames.(J.parentFrame));                 % bodyP <- off_P
Xco = fr_offInBody(J.frames.(J.childFrame));                  % bodyC <- off_C
Tjfull = Xpo * Tj / Xco;                                 % bodyP <- bodyC
T = fr_bodyT(kin, q, S.parentBody) * Tjfull;
end

function T = fr_jointT(~, q, J)
R = eye(3);
t = [0; 0; 0];
for k = 1:6
    A = J.axes{k};
    v = fr_axisValue(A, q);
    if k <= 3
        R = R * fr_rotAxis(A.axis(:), v);
    else
        t = t + A.axis(:) * v;
    end
end
T = eye(4);
T(1:3, 1:3) = R;
T(1:3, 4) = t;
end

function v = fr_axisValue(A, q)
switch A.type
    case 'linear'
        qv = 0;
        if ~isempty(A.coords); qv = q.(A.coords{1}); end
        v = A.coeffs(1) * qv + A.coeffs(2);
    case 'constant'
        v = A.value;
    case 'spline'
        qv = 0;
        if ~isempty(A.coords); qv = q.(A.coords{1}); end
        v = fr_splinemm(A.x, A.y, qv);
    otherwise
        error('fr_bodyT:AxisType', 'unknown axis type %s', A.type);
end
end

function X = fr_offInBody(F)
X = eye(4);
X(1:3, 1:3) = F.R;
X(1:3, 4) = F.t(:);
end
