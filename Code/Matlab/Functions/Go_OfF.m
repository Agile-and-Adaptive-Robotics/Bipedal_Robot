function [RMSE, fvu, maxResid] = Go_OfF(y,val)
%(Go)odness (Of) (F)it
%The lower the better for all values, therefore this function is good to
%use with optimization algorithms.

%Inputs:
%   y == Data being fit
%   val == Expected values at each data point
%Outputs
%  RMSE == Root Mean Square Error
%  fvu == Fraction of variance unexplained, adjusted
%  maxResid == Maximum residual

%If it's a matrix, make sure it is sorted. Then create interpolant models
%to find values at the 
if size(val,2)>1 && size(y,2)>1
    Yq = y(:,1);        % query points from y matrix
    y = y(:,2);
    if ~issortedrows(val,1)
        val = sortedrows(val,1);
    else
    end
    [a,~,c] = size(val);
    xData = NaN(a,c);
    yData = NaN(a,c);
    mdl = cell(c);
    value = cell(c);
        for i = 1:c
        xData(:,i) = val(:,1,i);
        yData(:,i) = val(:,2,i);
        mdl{i} = griddedInterpolant(xData(:,i),yData(:,i));
        value{i} = mdl{i}(Yq(:,1));
        end
    s = [mdl; value];
    rowNames = {'models'; 'values'};
    r = cell2struct(s(2,:), rowNames(2),1);
        for i = 1:c
        val(:,i) = r(i).values;
        end
elseif size(val,2)>1 && size(y,2)==1
    val = val(:,2);
else
end

yresid = y-val;                     %residual error
maxResid = max(abs(yresid));             %largest residual
SSresid = sum(yresid.^2,'omitnan'); %Sum of squares of the residual
SStot = sum((y-mean(y,1)).^2);        %total sum of squares
n = (sum(~isnan(val)));             %number of data points
RMSE = sqrt(SSresid./n);             % RMSE for function 1
fvu = SSresid*(n-1)./(SStot.*n);      %Fraction of Variance Unexplained (FVU),adjusted, Rsq = 1 - FVU 

end

