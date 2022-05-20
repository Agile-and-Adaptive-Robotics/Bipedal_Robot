function F = festo3(Lmt, rest, dia, pres, kmax, ten, Fitting)
%Inputs:
%Lmt == muscle-tendon length, scalar
%rest == resting length of artificial muscle, "size" from Size function
%dia == diameter of Festo tube, from Size function
%long == longest musculotendon length
%pres == pressure in kPa
%kmax == maximum measured contraction length. Input as length, will
%convert to percent in the code.
%ten == tendon length
%Fitting == fitting length
%Outputs:
%F == Force, N
if nargin == 5              %Use if comparing measured muscle length to resting
    tendon = 0;
    fitting = 0;
elseif nargin == 6          %Use if measured muscle length and "tendon" where tendon can be a stand in for everything else in the Lmt that isn't active muscle
    tendon = ten;
    fitting = 0;
elseif nargin == 7          %Use if muscle length, tendon, and fitting sizes are known
    tendon = ten;
    fitting = Fitting;
else
    fprintf('Invalid number of arguments\n')
end

load ForceStrainTable.mat ForceStrain
            
kmax = (rest-kmax)/rest; %Convert maximum contraction from length into percent

X = linspace(0,620,19); %Pressure for interpolation
Y = linspace(0,1,30);   %Relative strain range for interpolation

k = zeros(size(Lmt,1),1);
rel = zeros(size(Lmt,1),1);
F = zeros(size(Lmt,1),1);
            for i = 1:size(Lmt, 1)
                k(i,1) = (rest-(Lmt(i,1)-tendon-2*fitting))/rest; %current strain
                rel(i,1) = k(i,1)/kmax; %relative strain               
                if rel(i,1) >= 0 && rel(i,1) <=1
                    F(i,1) = interp2(X, Y, ForceStrain(:,2:20), pres, rel(i), 'linear');
                elseif k(i,1) < 0 && k(i,1) >= -0.03
                    x1 = [ -.03   -.02  -.01      0]';
                    z1 = [741.9    613   523  458.2]';
                    z2 = [759.2  629.3 539.3  473.3]';
                    z3 = [785.1  653.5 562.6  495.9]';
                    x = [x1; x1; x1];
                    y = [580*ones(length(x1),1);
                         600*ones(length(x1),1);
                         630*ones(length(x1),1)];
                    z = [z1; z2; z3]; 
                    BPAFit = fit([x, y],z,'linearinterp','Normalize','on');
                    F(i,1) = BPAFit(contract,pres);
                elseif rel(i,1) > 1
                    F(i,1) = 0;
                else
                    F(i,1) = NaN;
                end
            end

%If diameter is not 10 mm, then upscale force (but probably better to use
%festo4.m function
if dia == 20
    F = (1500/630)*F;
end

if dia == 40
    F = (6000/630)*F;
end

end

%Reference:
%Hunt, Alexander J., Alexander Graber-Tilton, and Roger D. Quinn. "Modeling length effects of braided pneumatic actuators."
%In ASME 2017 International Design Engineering Technical Conferences and Computers and Information in Engineering Conference, pp. V05AT08A008-V05AT08A008.
%American Society of Mechanical Engineers, 2017.