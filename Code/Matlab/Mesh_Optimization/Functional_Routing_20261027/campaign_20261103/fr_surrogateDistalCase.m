function [chosen,info] = fr_surrogateDistalCase(C,kin,chosen,branch,budget,mountSearch,retained,allowDistalSizing)
% Optimize continuous geometry and lengths on a coarse joint grid.
if nargin<6; mountSearch=false; end
if nargin<7; retained=[]; end
if nargin<8; allowDistalSizing=false; end
if mountSearch && ~isfield(chosen.caseDef,'distalSleeveGuideRow')
    chosen=fr_prepareAnkleMountSearch(C,kin,chosen);
end
if ~isempty(retained)
    chosen=fr_prepareAnkleWrap(C,chosen,kin);
end
cs=chosen.caseDef;
x0=fr_packX(cs,chosen.design);
x0=min(cs.ub,max(cs.lb,x0));
free=find(cs.ub-cs.lb>1e-10);
if mountSearch; free=free(free<=3*numel(cs.movable)); end
assert(~isempty(free),'No continuous decision variables');
lower=cs.lb(free); upper=cs.ub(free);
start=x0(free);
varySpread=mountSearch && cs.seed.nBPA>1;
if varySpread
    spreadIndex=numel(start)+1;
    lower(end+1)=1;
    upper(end+1)=max(3.2,cs.bundleSpread);
    start(end+1)=max(1,cs.bundleSpread);
end
if allowDistalSizing
    assert(mountSearch && any(strcmp(cs.rowOwning,'toes_r')));
    toe=fr_toeShiftPatches(C,kin,cs);
    restIndex=numel(start)+1;
    restRow=3*numel(cs.movable)+1;
    lower(end+1:end+2)=[cs.lb(restRow)/10;0];
    upper(end+1:end+2)=[cs.ub(restRow)/10;15];
    start(end+1:end+2)=[min(upper(restIndex),max(lower(restIndex),chosen.design.L0));0];
end
coupled=~strcmp(cs.name,'soleus');
opts=struct('quiet',true,'envelopes',{{}});
G=[];
if coupled
    N=4;
    G=fr_ankleCoupledSetup(C,kin,cs,N);
    source=fullfile(C.camp,'run_20261105_distal_stages','ankle_coupled', ...
        [cs.name '_' branch '_coupled_torque.csv']);
    D=readtable(source);
    a=unique(D.ankle_deg); b=unique(D{:,2});
    F=griddedInterpolant({a,b},reshape(D.human_ankle_Nm,numel(a),numel(b)));
    H1=F(G.ankle_deg,G.second_deg);
    F.Values=reshape(D.human_second_Nm,numel(a),numel(b));
    H2=F(G.ankle_deg,G.second_deg);
    if ~isempty(retained)
        [inner,innerDetail]=fr_evalCoupledAnkle(C,kin,retained.caseDef,retained.design,branch,N);
        assert(isequal(inner.ankle_deg,G.ankle_deg) && isequal(inner.second_deg,G.second_deg));
    end
end
spec=jsondecode(fileread(fullfile(C.camp,'objective_spec.json')));
terms=spec.terms;
if isstruct(terms); terms=num2cell(terms); end
shape=weight(terms,'torque_shape'); lengthWeight=weight(terms,'actuator_length');
attachmentWeight=weight(terms,'attachment_deviation');
capacityWeight=weight(terms,'torque_capacity_deficit');
romWeight=weight(terms,'usable_rom_shortfall');
overlapWeight=weight(terms,'retained_bpa_corridor_soft');
if allowDistalSizing; restWeight=weight(terms,'bpa_rest_length_preference'); end
evals=0;
firstError='';
bestScore=inf;
bestX=x0;
bestCase=cs;
initialScore=objective(start);
starts=start';
screen=[];
if ~isempty(retained)
    [seedX,home]=fr_seedAnkleShaftWrap(C,kin,cs,chosen.design);
    candidates=seedX(:,free);
    candidates=[candidates,repmat(start(numel(free)+1:end)',size(candidates,1),1)];
    % Each screened candidate uses the same coupled length refit and
    % constraints as the search, and counts against its evaluation budget.
    count=min(size(candidates,1),max(0,budget-3));
    scores=zeros(count,1);
    for k=1:count; scores(k)=objective(candidates(k,:)'); end
    [~,order]=sort([initialScore;scores]);
    proposals=[start';candidates(1:count,:)];
    starts=proposals(order(1:min(6,numel(order))),:);
    screen=struct('homeMargin_m',home(1:count,2), ...
        'homeLargestTurn_deg',rad2deg(home(1:count,3)), ...
        'coupledObjective',scores,'seedObjective',initialScore);
    fprintf('%s %s: screened %d shaft routes, seed objective %.4g -> %.4g\n', ...
        cs.name,branch,count,initialScore,min([initialScore;scores]));
end
o=optimoptions('surrogateopt','Display','off','UseParallel',false, ...
    'MaxFunctionEvaluations',budget-evals,'InitialPoints',starts);
[~,~,exitFlag,out]=surrogateopt(@objective,lower,upper,o);
assert(isfinite(bestScore) && bestScore<1e10,'No finite evaluations: %s',firstError);
chosen.caseDef=bestCase;
chosen.design=fr_unpackX(bestCase,bestX);
[chosen.f,chosen.rep]=fr_objective(C,bestCase,kin,bestX,branch,opts);
info=struct('evaluations',evals,'initialScore',initialScore,'score',bestScore, ...
    'solverOutput',out,'exitFlag',exitFlag,'firstEvaluationError',firstError, ...
    'gridSize',1+coupled*3,'seedScreen',screen);

    function f=objective(z)
        evals=evals+1;
        x=x0; x(free)=z(1:numel(free));
        candidate=cs;
        grid=G;
        if varySpread
            candidate.bundleSpread=z(spreadIndex);
            if coupled; grid.cs.bundleSpread=z(spreadIndex); end
        end
        try
            d=fr_unpackX(candidate,x);
            movement=mean(vecnorm(d.uv-chosen.design.uv,2,2))/25;
            if allowDistalSizing
                [candidate,shift]=fr_shiftToeRoute(candidate,toe,z(restIndex+1));
                movement=movement+mean(shift)/.025;
                grid=fr_ankleCoupledSetup(C,kin,candidate,N);
            end
            if ~isempty(retained)
                [candidate,d]=fr_materializeAnkleWrap(C,kin,candidate,d);
                grid=fr_ankleCoupledSetup(C,kin,candidate,N);
            end
            if mountSearch
                if allowDistalSizing
                    d=fr_refitMountLengths(C,kin,candidate,d,branch,grid,z(restIndex));
                else
                    d=fr_refitMountLengths(C,kin,candidate,d,branch,grid);
                end
                x=fr_packX(candidate,d);
            end
            if coupled
                layerCost=0;
                if isempty(retained)
                    e=fr_evalCoupledAnkle(C,kin,candidate,d,branch,N,grid);
                else
                    [e,outerDetail]=fr_evalCoupledAnkle(C,kin,candidate,d,branch,N,grid);
                    layer=fr_compareAnkleLayers(C,kin,e.poses,outerDetail.sleeves, ...
                        innerDetail.sleeves,candidate.fittingLength,retained.caseDef.fittingLength);
                    layerCost=1e5*sum(max(0,-[layer.orderingMargin_m,layer.rigidMargin_m]-1e-4)) ...
                        +overlapWeight*layer.sleeveOverlap_m/.020;
                end
                bad=e.noForce;
                hard=mean(bad)+sum(max(0,-[e.boneMargin_m,e.leadMargin_m,e.footMargin_m]+(-1e-4)))+ ...
                    max(0,-e.bundleMargin_m-1e-9)+max(0,e.tendonRequired_m-d.tendon-1e-6);
                hard=hard+sum(max(0,-[e.distalGuideMargin_m,e.wrapMargin_m,e.shaftCableMargin_m, ...
                    e.bundleHardwareMargin_m]-1e-4));
                f=1e5*hard+signedCost(e.ankle_Nm,H1)+signedCost(e.second_Nm,H2)+layerCost;
                if isfield(e,'stopTurn_deg')
                    f=f+12*((e.guideTurn_deg/90)^2+(e.stopTurn_deg/90)^2);
                end
                f=f+lengthWeight*max(0,d.L0+d.tendon-cs.L0Seed-cs.tendonSeed)/ ...
                    max(cs.L0Seed+cs.tendonSeed,eps);
                f=f+attachmentWeight*movement;
                if allowDistalSizing
                    f=f+restWeight*d.L0/max(cs.L0Seed,eps);
                end
            else
                f=fr_objective(C,candidate,kin,x,branch,opts);
            end
            if ~isfinite(f); f=1e10; end
        catch ME
            if isempty(firstError); firstError=ME.message; end
            f=1e10;
        end
        if f<bestScore; bestScore=f; bestX=x; bestCase=candidate; end
    end

    function f=signedCost(robot,human)
        deficit=max(0,abs(human)-robot.*sign(human));
        f=capacityWeight*sum(deficit)/max(sum(abs(human)),eps)+ ...
          romWeight*mean(deficit>1e-6)+shape*mean((robot-human).^2)/max(mean(human.^2),eps);
    end
end

function w=weight(terms,name)
for k=1:numel(terms)
    if strcmp(terms{k}.name,name); w=terms{k}.weight; return; end
end
error('Missing objective weight %s',name);
end
