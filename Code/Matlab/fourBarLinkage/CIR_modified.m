% ----- Rotation Center Determination ----------
% ----- Neider Nadid® 11/03/2020 ---------------
clc; close all; clear;
%% ------------ Dimensions of mechanism ----------------------------------
a = 51.4; 
b = 39.76;
c = 46.86; 
x_A = 4.5;
y_A = -395.83;
x_B = 22.28; 
y_B = -411.04;
A = [x_A;y_A];
B = [x_B;y_B];
%% ------------------------------------------------------------------------
phi = 180:1:270+45;
t = 1; % Configuration may be either -1 or 1
for i=1:size(phi,2)
    
    u_AC = [cos(phi(i)*pi/180);sin(phi(i)*pi/180)];
    C = A+a*u_AC;
    s = ((B-C)'*(B-C))^0.5;
    l = (b^2+s^2-c^2)/(2*s);
    h = t*(b^2-l^2)^0.5;
    N = (1/s)*[l -h;h l];
    D = C+N*(B-C);
    u=(1/a)*(C-A); 
    v=(1/c)*(D-B);
    D1=[(B(1)-A(1)) v(1);(B(2)-A(2)) v(2)]; 
    DD=[u(1) v(1);u(2) v(2)];
    landa=det(D1)/det(DD);
    I=A+landa*u_AC; % Instant centre of rotation 
    x_I(i) = I(1);
    y_I(i) = I(2);
    x = [A(1) C(1) D(1) B(1)];
    y = [A(2) C(2) D(2) B(2)];
    plot(x,y,'g','LineWidth',2); grid on; hold on; axis equal; axis([-60 120 -490 -320]) 
    plot(x,y,'o r','LineWidth',2)
    plot(x_I,y_I,'m','LineWidth',1)
    plot([C(1) I(1) D(1)],[C(2) I(2) D(2)],'-- b','LineWidth',1)
    pause (0.01)
    drawnow
    hold off
end

