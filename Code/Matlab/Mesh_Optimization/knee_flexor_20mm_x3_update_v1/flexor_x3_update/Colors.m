%% Colors that are accessible
c = cell(7,1);
c{1} = '#FFD700'; % gold
c{2} = '#FFB14E'; % orange
c{3} = '#FA8775'; % light orange
c{4} = '#EA5F94'; % pink
c{5} = '#CD34B5'; % magenta
c{6} = '#9D02D7'; % magenta 2
c{7} = '#0000FF'; % indigo

d = hex2rgb(c); 


% Output will be 
% d =[1.0000    0.8431         0;
%     1.0000    0.6941    0.3059;
%     0.9804    0.5294    0.4588;
%     0.9176    0.3725    0.5804;
%     0.8039    0.2039    0.7098;
%     0.6157    0.0078    0.8431;
%          0         0    1.0000];