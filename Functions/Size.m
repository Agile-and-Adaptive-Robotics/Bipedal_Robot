function [ID, size, longest]= Size(length, mif)
%function Size for sizing of Festo artificial muscles
%Inputs:
%Length is the musculotendon length, either a column vector or a matrix
%mif is Maximum Isometric Force
%Outputs:
%ID is inside diameter of the specified festo muscle
%size is the length of the festo muscle
%long is the longest artificial musculo-tendon length
%Variables:
%Fc is a percentage of max theoretical festo muscle force
%x is length of fitting from end to attachment point.

%Notes:
%1) If error calculating ID based on force, then two muscles in parallel
%should be considered.
%2) "Size Calc error1" if no change between shortest and largest lengths
%for give dof(s). 
%3) "Size Calc error2" if error calculating size because delta is too large 
%compared to smallest length, then consider disregarding muscle or changing 
%attachment points(e.g. short interior groin and hip muscles).

long = max(max(length));
short = min(min(length));
delta = long - short;  
F1 = 630;             %Maximum theoretical force for 10mm I.D. Festo muscle
F2 = 1500;            %Maximum theoretical force for 20mm I.D. Festo muscle
F3 = 6000;            %Maximum theoretical force for 40mm I.D. Festo muscle
Fc = 0.9;          
x = 0.0125;           %Length of air fittings

if mif < Fc*F1
    ID = 10;
elseif mif >= Fc*F1 && mif < Fc*F2
    ID = 20;
elseif mif >= Fc*F2 && mif <= Fc*F3
    ID = 40;
else
    ID = 'Use additional muscles';
end

opt_max = 1.09;  %Maximum allowable length change for optimal performance (Festo)
%opt_min = 1;
max_def = 1.25;  %Maximum allowable length change (Festo)


size = delta/(1 - 1/opt_max);
if delta == 0
    size = 'Size calc error1';
    return
end    

if size+2*x > long
    for def = opt_max:0.001:max_def
        if size+2*x > long && def <= max_def
            size = delta/(1 - 1/def);
        elseif size+2*x <= long
            def;
            break
%         elseif def == max_def
%             def
%             size = 'Size calc error2';
        end
    end
end

if nargout > 2
    longest = long;
end

end