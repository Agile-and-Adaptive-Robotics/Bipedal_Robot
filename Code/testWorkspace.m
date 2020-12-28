% clc                 %Command Line Clear
% clear               %Clear the workspace of stored variables
% close all           %Close all open figures
% 
% gamma = 0.5;                    %Shape factor for force length curve of contractile element
% epsilon = 0.6;                  %Passive muscle strain due to maximum isometric force
% kPE = 4;                        %Exponential shape facor for passive muscle element
% 
% tsL = 0.12;                     %Tendon slack length
% alphaP = 0.05235988;            %Optimal pennation angle
% ofL = 0.121;                    %Optimal muscle fiber length
% 
% mtL = Add_Mag2.MuscleLength(1)                   %Muscle-tendon length
% 
% mif = 343;                      %maximum isometric force
% 
% %Estimate the nominal muscle length, without the tendon
% if mtL == tsL
%     mtL0 = (mtL*1.001 - tsL)/ofL;
% else
%     mtL0 = (mtL - tsL)/ofL;
% end
% 
% %Set function solver parameters
% options = optimoptions('fsolve','Display','none','FunctionTolerance',0.001);
% 
% mL = fsolve(@fBalance, mtL0, options);
% 
% 
% 
% %Force of the contracile muscle element
% fL = exp(-(Lmn - 1).^2/gamma);
% 
% %Force of the passive elastic muscle element
% fPE = (exp(kPE*(Lmn - 1)/epsilon - 1)/(exp(kPE) - 1));
% 
% %Force of the tendon elastic element
% fT = 37.5/Lts*(Lmt - Lmn.*sqrt(1 - (sin(alphap)./Lmn).^2) - Lts);
% 
% cosAlpha = sqrt(1 - (sin(alphaP)/x0)^2);
% 
% 
%     
%     
% figure
% hold on
% plot(Lmn, fL)
% plot(Lmn, fPE)
% hold off


% addpath('Previous Optimization Code\Functions')
% ftest = forz(mtL, mif, ofL, tsL, alphaP)
% 
% 
Lmn = linspace(0, 2);           %Nominal muscle length
gamma = 0.5;                    %Shape factor for force length curve of contractile element
epsilon = 0.6;                  %Passive muscle strain due to maximum isometric force
kPE = 4;                        %Exponential shape facor for passive muscle element

Lts = 0.12;                     %Tendon slack length
alphap = 0.05235988;            %Optimal pennation angle
Lofl = 0.121;                   %Optimal muscle fiber length


mif = 343;                      %maximum isometric force

%Force of the contracile muscle element
fL = exp(-(Lmn - 1).^2/gamma);

%Force of the passive elastic muscle element
fPE = (exp(kPE*(Lmn - 1)/epsilon) - 1)/(exp(kPE) - 1);

%Force of the tendon elastic element
% fT = 37.5/Lts*(Lmt - Lmn.*sqrt(1 - (sin(alphap)./Lmn).^2) - Lts);

%
% cosalpha = sqrt(((1 - (sin(alphap)/Lmn)).^2));

figure
hold on
plot(Lmn, fL)
plot(Lmn, fPE)
% plot(Lmn, fT)
hold off
% function muscleF = fBalance(nmL)
%     fPE = (exp(kPE*(nmL - 1)/epsilon - 1)/(exp(kPE) - 1));
%     
%     fL = exp(-(nmL - 1).^2/gamma);
% 
%     cosAlpha = sqrt(1 - (sin(alphaP)/nmL)^2);
%     
%     fT = 37.5/tsL*(mtL - nmL*cosAlpha - tsL);
%     
%     muscleF = (fL + fPE)*cosAlpha - fT;
% end