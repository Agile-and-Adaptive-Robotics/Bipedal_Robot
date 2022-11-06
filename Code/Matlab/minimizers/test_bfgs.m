
clear all


% f = @(x) [0.5*x'*x;(1+x(1)-x(2)).^2];
% f = @(x) 0.5*x'*x;
% x0 = [-1;-1];
% bnds = [-20,20;-10,10;-5,5];
% x0 = [-10;-10;-20];
% bnds = [-20,-5;-20,-5;-30,-10];
% cons = {[1,-1,0],[0]};
% bnds = [];
% cons = [];
% x0 = -4;
% bnds = [-5,5];
% cons = {[],[]};
% pop_size = 101;

% f = @(x) [100*(x(2)-x(1).^2).^2 + (1-x(1)).^2;(x(1)-x(2)).^2];
f = @(x) 100*(x(2)-x(1).^2).^2 + (1-x(1)).^2;
% x0 = [-1;1];
x0 = [0;0]
% bnds = [-1,1;-1,1];
% bnds = [-1,1;-1,1];
% cons = [];
% pop_size = 1001;

% f = @(x) sum(x);
% bnds = [-1,1;-1,1];
% cons = [];
% pop_size = 1001;


% f = @(x) (x(1)-1e6).^4 + (x(1)-1e6).^2 + (x(2)-1e-6).^4 + (x(2)-1e-6).^2;

% f = @(x) 1e-3*(x'*x) + (1-cos(x(1))) + (1-cos(x(2)));

% x0 = [-10;1];
% x0 = [-1;-1;-2];
% bnds = [];
% bnds = [-20,-5;-20,20];
% bnds = [-20,1;-20,20];
% bnds = [-1.1,1;-1,1.1];
% bnds = [-20,20;-10,10;-5,5];
% bnds = [-1.01,1;-1.01,1;-2.01,1];
% cons = {[1,1,-1],[0]};
% cons = {[1,-1,0],[0]};
% conv_crit = [1e3,1e-14,1e-14,0];
direction_options = [false];
search_options = {'backtrack',1e-4,0.5,15,0.8,1};
plot_options = [1,2,3];
print_options = 1;
parallelize = false;

% x_final = bfgs(f,x0,bnds,cons,conv_crit,direction_options,search_options,plot_options,print_options)
% x_final = bnd_bfgs(f,x0,bnds,cons,conv_crit,direction_options,search_options,plot_options,print_options)

% x_ga = GA(f,bnds,cons,[1000,.99,1e-9,0],pop_size,[],'single',.001,randi(100,1,1),parallelize,1,2)
% x_final = bnd_bfgs(f,x0,bnds,[],conv_crit,direction_options,search_options,plot_options,print_options)
x_final = bfgs(f,x0,[],[],0,[2000],direction_options,search_options,plot_options,print_options);
% x_final = nelder_mead(f,x0,{bnds,1e-6},cons,[],[],true,1)
% x_final = nelder_mead(f,x0,[],[],[1000,1e-15],[],true,1)
% x_final = bnd_nelder_mead(f,x0,bnds,cons,[],[],[1,2,3,4,5,6],1)












