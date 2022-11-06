function [Y_of_X, r_sq, C_final] = linear_fit_and_r_sq(X,Y,bnds)
    
%     plane = @(C,X,Y) C(1) + C(2).*X + C(3).*Y; %X, Y must be columns
%     plane_diff = @(C,X,Y,Z) plane(C,X(:),Y(:)) - Z(:);
    
    line = @(C,X) C(1) + C(2).*X; %X must be a column
    line_diff = @(C,X,Y) line(C,X(:)) - Y(:); %Y must be a column
    
%     f = @(C) 0.5*sum(plane_diff(C,X,Y,Z).^2);
    f = @(C) 0.5*sum(line_diff(C,X,Y).^2);
    if isempty(bnds)
        bnds = [-1000,1000;-1000,1000];
    end
    
    pop_size = 50000;
%     C_GA = GA(f,bnds,false,false,'min',[50,.95,1e-6],pop_size,0.5,'continuous',.001,[],1,false,0,1);
    C_GA = GA(f,bnds,{[],[]},[1000,.95,1e-6,1],5000,.5,'single',.001,1,true,1,1);
%     C_final = opti(f,[mean(Z);.001;.001],{[],[]},{bnds,0,0},[],[],[1e6,1e-6,1e-6],'min','bfgs',[],'linesearch',{'fit',[1e-4,0.8,0.5,20],[],[]},[],1);
%     C_final = opti(f,C_GA,{[],[]},{bnds,0,0},[],[],[1e6,1e-6,1e-6],'min','bfgs',[],'linesearch',{'fit',[1e-4,0.8,0.5,20],[],[]},[],1);
    C_final = bnd_bfgs(f,C_GA,bnds,{[],[]},[],[],[],1,1);
    
    Y_of_X = @(X) line(C_final,X);
    
    Yfit = Y_of_X(X); %Model's prediction of Z
    Y_resid = Y - Yfit; %Error from the actual data
    SSresid = sum(Y_resid.^2); %Sum of the squared residuals
    SStotal = sum((Y-mean(Y)).^2); %Sum of the data squared
    r_sq = 1 - SSresid/SStotal; %Definition of r_sq   
    
end