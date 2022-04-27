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

alpha = 1;         %Reflection value
gamma = 1.5;         %Expansion value
beta = 0.5;        %Contraction value
sigma = 0.5;       %Shrink value
epsilon = 0.001;   %Convergence value

f = zeros((n+1),1); 
k = 1;
iter = 150;                                   %number of iterations
while k < iter+1                              %Run loop for maximum number of iterations if needed
    
    for i = 1:n+1
       f(i,1) = minimizeFlx10mm_nest(X(i,:),X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c); 

    end
    
    A = [X, f];
    B = sortrows(A,(n+1));                  %sort A from low to high
    X = B(:,1:n);
    Xc = (1/n)*sum(X((1:n),:));             %Centroid
    Xr = (1+alpha)*Xc-alpha*X((n+1),:);     %Reflection point
    fr = minimizeFlx10mm_nest(Xr,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);

    
        if fr < f(n) && fr > f(1)              %Reflection
            Xnew = Xr;
            fnew = fr;
            
        elseif fr < f(1)                        %Expansion
            Xnew = (1+gamma)*Xc-gamma*X(n+1,:);
            fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);

        elseif fr >= f(n+1)                     %Contraction, fr>f(n+1)
            Xnew = (1-beta)*Xc+beta*X(n+1,:);
            fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
                
        elseif fr < f(n+1) && fr > f(n)        %Contraction, fr=<f(n+1)
            Xnew = (1+beta)*Xc-beta*X(n+1,:);
            fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
        
        else                                    %Shrink
            X = X(1,:)+sigma*(X-X(1,:));
            Xnew = X(1,:);
            fnew = f(1);
        end
        
%         if fr < f(n) && fr >= f(1)              %Reflection
% %             Xnew = Xr;
% %             fnew = fr;
%             X(n+1,:) = Xr;
%             f(n+1) = fr;
%             
%         elseif fr < f(1)                        %Expansion
%             Xnew = (1-gamma)*Xc+gamma*Xr;
%             fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
%             if fnew > f(1)
%                 X(n+1,:) = Xr;
%                 f(n+1) = fr;
%             else
%                 X(n+1,:) = Xnew;
%                 f(n+1)=fnew;
%             end
%         elseif fr >= f(n+1)                     %Contraction, fr>f(n+1)
%             Xnew = (1-beta)*Xc+beta*X(n+1,:);
%             fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
%             if fnew >= f(n+1)                       %Shrink
%                 X = X(1,:)+sigma*(X-X(1,:));
%                 Xnew = X(1,:);
%                 fnew = f(1);
%                 for i = 2:n+1
%                     f(i,1) = minimizeFlx10mm_nest(X(i,:),X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
%                 end
%             else
%                 X(n+1,:) = Xnew;
%                 f(n+1) = fnew;
%             end
%                 
%         elseif fr < f(n+1) && fr >= f(n)        %Contraction, fr=<f(n+1)
%             Xnew = (1+beta)*Xc-beta*X(n+1,:);
%             fnew = minimizeFlx10mm_nest(Xnew,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
%             if fnew >= fr                       %Shrink
%                 X = X(1,:)+sigma*(X-X(1,:));
%                 Xnew = X(1,:);
%                 fnew = f(1);
%                 for i = 2:n+1
%                     f(i,1) = minimizeFlx10mm_nest(X(i,:),X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
%                 end
%             else
%                 X(n+1,:) = Xnew;
%                 f(n+1) = fnew;
%             end
% 
%         end


    fc = minimizeFlx10mm_nest(Xc,X1,Name, Location, CrossPoint, Dia, T_Pam, fitting, pres, phiD, y1, y3, c);
    Q = sqrt(1/(n+1)*sum((f(1:n+1)-fc).^2));
       if Q <= epsilon
           X(1,:)
           f(1)
           X(n+1,:)
           f(n+1)
           Q
           k
           break
       elseif Q > epsilon
           k = k+1;
           X(n+1,:) = Xnew
           if k == iter+1
                X(1,:)
                f(1)
                X(n+1,:)
                f(n+1)
                Q
                k
           end
       end 
end    