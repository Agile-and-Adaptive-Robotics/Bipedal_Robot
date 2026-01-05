function [ T ] = scale( x )
%SCALE Scale the domain by building a matrix T, which when multiplied by
%the domain gives values that are roughly the same order of magnitude.
%   Detailed explanation goes here
    n = size(x,1);
    T = zeros(n);
    for i=1:n
        val = abs(x(i));
        if val > 0
            oom = floor(log10(val));
        else
            oom = 0;
        end
        T(i,i) = 1/(10^oom);
    end
end

