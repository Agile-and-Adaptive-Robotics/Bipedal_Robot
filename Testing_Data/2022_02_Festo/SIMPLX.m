clear; close all;
%Create a Nelder-Mead simplex algorith to search for 3-Dimensional search
load KneeFlx_10mm_42cm.mat TorqueR phiD Location Name CrossPoint Dia T_Pam fitting
Theoretical = TorqueR(:,3)';   %Load theoretical torque
load Plot_KneeFlx_10mm_42cm.mat X1 Angle Torque modp mdl1 gofp2 TorqueMean pres 
pres = mean(pres);                  %Make pressure a scalar value
y1 = Torque';                        %Make it just the data
y = feval(mdl1, X1);               %Make y1 the curve fit
y3 = diff(y);
weight = [0.75 0.25];               %Assign relative weight to  functions 1 and 2
weight_norm = norm(weight);         %Find magnitude of vector "weight"
c = weight/weight_norm;                          %Weights c1 & c2 for functions 1 & 2, respectively, in unit vector form

n = 3;              %number of dimensions
a = 0.003;          %size of tetrahedron
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

alpha = 1;         %Reflection value
gamma = 2;         %Expansion value
beta = 0.5;        %Contraction value
sigma = 0.5;       %Shrink value
epsilon = 0.001;   %Convergence value

f = zeros((n+1),1);
k = 1;
iter = 250;                                 %number of iterations
while k < iter+1                               %Set max iterations to 100
   for i = 1:n+1
       f(i,1) = minimizeFlx10mm_nest(X(i,:),X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c); 

   end    
    A = [X, f];
    B = sortrows(A,(n+1));                  %sort A from low to high
    x = B(:,1:n);
    xc = (1/n)*sum(x((1:n),:));             %Centroid
    xr = (1+alpha)*xc-alpha*x((n+1),:);     %Reflection point
    fr = minimizeFlx10mm_nest(xr,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);

        if fr < f(1)                        %Expansion
            xnew = (1-gamma)*xc+gamma*xr;
        elseif fr >= f(n+1)                     %Contraction, fr>f(n+1)
            xnew = (1-beta)*xc+beta*x(n+1,:);
        elseif fr < f(n+1) && fr >= f(n)        %Contraction, fr=<f(n+1)
            xnew = (1+beta)*xc-beta*x(n+1,:);
        else
            x = x(1,:)+sigma*(x-x(1,:));        %Shrink
        end
    fnew = minimizeFlx10mm_nest(xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
    fc = minimizeFlx10mm_nest(xc,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
    Q = sqrt((1/n+1)*(sum(f(1:n)-fc)^2+(fnew-fc)^2));
       if Q <= epsilon
           xnew
           fnew
           Q
           k
           break
       elseif Q > epsilon
           k = k+1;
           x((n+1),:) = xnew;
           if k == iter+1
                xnew
                fnew
                Q
                k
           end
       end 
end    