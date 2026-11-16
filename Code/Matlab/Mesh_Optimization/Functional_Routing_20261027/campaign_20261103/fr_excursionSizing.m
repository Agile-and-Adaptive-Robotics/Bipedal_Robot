function T=fr_excursionSizing(C,path_m,bay_m,humanFmax_N,forceScale)
% First-pass stroke and force screen for a route sampled at common poses.
% Rest at the longest path gives zero nominal contraction at that pose.
% Candidate-route travel and sleeve bay must use the same pose ordering.
if nargin<5; forceScale=1; end
path_m=path_m(:); bay_m=bay_m(:);
assert(numel(bay_m)==1 || numel(bay_m)==numel(path_m));
assert(all(isfinite(path_m)) && all(isfinite(bay_m)) && all(bay_m>0));
assert(isscalar(humanFmax_N) && humanFmax_N>0);
takeup=max(path_m)-path_m;
halfFmax10=.5*maxBPAforce(max(path_m),'10',C.pressure);
rows=cell(6,1);
k=0;
for dia=[10 20]
    key=num2str(dia);
    limit=1-C.LminFraction.(['d' key]);
    for gain=1:3
        k=k+1;
        reeving=fr_travelReeving(gain);
        stroke=takeup/gain;
        minimum=max(stroke)/limit;
        maximum=min(bay_m+stroke);
        feasible=minimum<=maximum && maximum>0;
        restForce=NaN; worstForce=NaN; count=NaN;
        if feasible
            % The largest fitting rest length provides the best-case force
            % screen. The exact solver can select a shorter viable sleeve.
            perBPA=maxBPAforce(maximum,key,C.pressure);
            f=perBPA*festo4(dia,stroke/(maximum*limit),C.pressure)/gain;
            restForce=perBPA*festo4(dia,0,C.pressure)/gain;
            worstForce=min(f);
            if worstForce>0
                count=ceil(forceScale*humanFmax_N/worstForce);
            end
        end
        rows{k}=struct('diameter_mm',dia,'gain',gain,'ideal_tackle_sheaves',reeving.sheaveCount, ...
            'required_stroke_mm',1000*max(stroke),'minimum_rest_mm',1000*minimum, ...
            'maximum_proxy_rest_mm',1000*maximum,'axial_screen_pass',feasible, ...
            'rest_cable_force_per_BPA_N',restForce,'minimum_cable_force_per_BPA_N',worstForce, ...
            'parallel_count_screen',count,'initial_10mm_half_Fmax_N',halfFmax10, ...
            'human_Fmax_with_scale_N',forceScale*humanFmax_N);
    end
end
T=struct2table(vertcat(rows{:}));
end
