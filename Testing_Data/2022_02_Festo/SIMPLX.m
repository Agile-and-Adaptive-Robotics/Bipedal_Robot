clear; clc; close all;
%Create a Nelder-Mead simplex algorith to search for 3-Dimensional search
load KneeFlx_10mm_42cm.mat TorqueR phiD Location Name CrossPoint Dia T_Pam fitting
Theoretical = TorqueR(:,3)';   %Load theoretical torque
load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque modp mdl1 gofp2 TorqueMean pres 
pres = mean(pres);                  %Make pressure a scalar value
y1 = Torque;                        %Make it just the data
y = feval(mdl1, X1);               %Make y1 the curve fit
y3 = diff(y);
c = [1 1];                          %Weights c1 & c2 for functions 1 & 2, respectively

n = 3;              %number of dimensions
a = 0.001;          %size of tetrahedron
p = a/(n*sqrt(2))*(sqrt(n+1)+n-1);
q = a/(n*sqrt(2))*(sqrt(n+1)-1);
pq = p*eye(n)+q*(ones(n)-eye(n));
X = zeros((n+1),n);
rest = 0.415;
tendon = 0.012;
kmax = 0.350;
X(1,:) = [rest, tendon, kmax];

for k = 1:n
    dX = pq(k,:);
    X(k + 1, : ) = X(1, :) + dX;   
end

% figure('Color', 'w'), hold on, grid on, rotate3d on
% plot3( X( 1, 1 ), X( 1, 2 ), X( 1, 3 ), '.b', 'Markersize', 20 )
% plot3( X( 2, 1 ), X( 2, 2 ), X( 2, 3 ), '.r', 'Markersize', 20 )
% plot3( X( 3, 1 ), X( 3, 2 ), X( 3, 3 ), '.k', 'Markersize', 20 )
% plot3( X( 4, 1 ), X( 4, 2 ), X( 4, 3 ), '.m', 'Markersize', 20 )

gamma = 1.5;       %Expansion value
beta = 0.5;        %Contraction value
epsilon = 0.001;   %Convergence value

f = zeros((n+1),1);
k = 1;
while k < 101                                 %Set max iterations to 100
   for i = 1:n+1
       f(i,1) = minimizeFlx10mm_nest(X(i,:),X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c) 
%        f(i,1) = (1-x(i,1))^2+(2-x(i,2))^2;
   end    
    A = [X, f];
    B = sortrows(A,(n+1));
    x = [B(1,1), B(1,2); B(2,1), B(2,2); B(3,1), B(3,2)];
    xc = 0.5*(x(1,:)+x(2,:));  %Centroid
    xr = 2*xc-x(3,:); %Reflection point
    fr = (1-xr(1,1))^2+(2-xr(1,2))^2;
        if fr < f(1)  %expansion
            xnew = (1+gamma)*xc-gamma*x(3,:);
        elseif fr >= f(3) %contraction
            xnew = (1-beta)*xc+beta*x(3,:);
        elseif fr < f(3) && fr > f(2) %contraction
            xnew = (1+beta)*xc-beta*x(3,:);
        end
    fnew = (1-xnew(1,1))^2+(2-xnew(1,2))^2;
    fc = (1-xc(1,1))^2+(2-xc(1,2))^2;
    Q = (((f(1)-fc)^2+(f(2)-fc)^2+(fnew-fc)^2)/3)^(1/2);
       if Q <= epsilon
           disp(xnew)
           disp(fnew)
           disp(Q)
           break
       elseif Q > epsilon
           k = k+1;
           x(3,:) = xnew;    
       end 
end    