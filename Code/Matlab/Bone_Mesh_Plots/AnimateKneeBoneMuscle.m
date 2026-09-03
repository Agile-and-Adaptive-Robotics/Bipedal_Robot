function fig = AnimateKneeBoneMuscle(T, T_ICR_t1, phi, pos, p1, p2, Bifemsh, varargin)
% AnimateKneeBoneMuscle
% Human bones/paths use T. A supplied robot Location uses T_Pam.
% Location is M-by-3-by-N with rows before CrossPoint in the femur frame
% and rows CrossPoint:end already in ICR coordinates, exactly as passed to
% MonoPamDataExplicit_balanceX3. Do NOT apply T_ICR_t1 to those rows again.
%
% New route-aware inputs (old seven-argument calls still work):
%   'Location', Location, 'CrossPoint', CrossPoint, 'T_Pam', T_Pam
%   'Static', true          render only pos; no animation or export
%   'FrameIndices', frames  explicit playback order, overriding FrameStep
% Bifemsh can be one muscle or a cell array of muscles with Location/Cross.
% Empty XLim/YLim/ZLim fit all displayed bone and muscle coordinates.
% For routed Y-up XYZ animations only, AnimationPaddingX extends +X;
% AnimationPaddingZ extends both ends of Z. Both are in metres and default
% to 0.10 and 0.05. Padding is added after the existing limits are computed.
% FullSkeleton=true/false selects full skeleton / femur+tibia in either mode.
% If omitted, static defaults to StaticFullSkeleton (true); animation to false.
% Full skeleton adds pelvis, sacrum, lumbar spine and talus/calcaneus/toes.
% In animation the foot follows the human tibia pose; proximal bones stay fixed.
% StaticPadding adds 0.15 m on each side of the static scene by default;
% pass a scalar or [X Y Z] in metres. Static limits always contain the bones.
% RobotColor / HumanColor / EndColor default to c{5} / c{3} / c{6}
% from the supplied Colors.m. Pass actual c entries to override these.
% DisplayAxisMap remaps the displayed coordinates AFTER DisplayRotation.
% BOTH display matrices default to eye(3): native X/Y/Z are preserved.
% DisplayAxisMap=[] also means identity. There is no implicit axis swap.
% For Y-up native coordinates, rotate the CAMERA with CameraOrbitDeg;
% do not swap point coordinates just to obtain a different viewpoint.
% Axis limits and AnimationPaddingX/Z refer to the final displayed axes.

    p = inputParser;
    addParameter(p,'FemurFile','Femur_Mesh_Points.xlsx');
    addParameter(p,'TibiaFile','Tibia_Mesh_Points.xlsx');
    addParameter(p,'SpineFile','Spine_Mesh_Points.xlsx');
    addParameter(p,'SacrumFile','Sacrum_Mesh_Points.xlsx');
    addParameter(p,'PelvisFile','Pelvis_R_Mesh_Points.xlsx');
    addParameter(p,'TalusFile','Talus_Mesh_Points.xlsx');
    addParameter(p,'CalcaneusFile','Calcaneus_Mesh_Points.xlsx');
    addParameter(p,'ToesFile','Toes_Mesh_Points.xlsx');
    addParameter(p,'UnitScale',1);
    addParameter(p,'Hip',[-0.0707,-0.0661,0.0835]);
    % Preserve the identity DisplayRotation in your uploaded file.
    addParameter(p,'DisplayRotation',eye(3));
    addParameter(p,'DisplayAxisMap',eye(3));
    addParameter(p,'PauseTime',0.25);
    addParameter(p,'FrameStep',1);
    addParameter(p,'Loop',true);
    addParameter(p,'ExportGif',false);
    addParameter(p,'GifFile','KneeBoneMuscle.gif');
    addParameter(p,'ExportVideo',false);
    addParameter(p,'VideoFile','KneeBoneMuscle.mp4');
    addParameter(p,'FrameRate',20);
    addParameter(p,'PauseAtFrames',[]);
    addParameter(p,'XLim',[-0.5,0.05]);
    addParameter(p,'YLim',[-0.225,0.0348]);
    addParameter(p,'ZLim',[-1,0.04]);
    addParameter(p,'Location',[]);
    addParameter(p,'CrossPoint',[]);
    addParameter(p,'T_Pam',[]);
    addParameter(p,'HumanLabels',{});
    addParameter(p,'Static',false);
    addParameter(p,'StaticFullSkeleton',true);
    addParameter(p,'FullSkeleton',[]);
    addParameter(p,'StaticPadding',0.15);
    addParameter(p,'FrameIndices',[]);
    addParameter(p,'CameraOrbitDeg',0);
    addParameter(p,'AnimationPaddingX',0.10);
    addParameter(p,'AnimationPaddingZ',0.05);
    addParameter(p,'RobotColor','#CD34B5'); % c{5}: magenta PAM
    addParameter(p,'HumanColor','#FA8775'); % c{3}: light orange muscle
    addParameter(p,'EndColor','#9D02D7');   % c{6}
    parse(p,varargin{:});
    opt = p.Results;

    phi = phi(:);
    N = numel(phi);
    validateattributes(phi,{'numeric'},{'real','finite','nonempty'});
    validateattributes(pos,{'numeric'},{'scalar','integer','>=',1,'<=',N});
    validateattributes(p1,{'numeric'},{'real','finite','numel',3});
    validateattributes(p2,{'numeric'},{'real','finite','numel',3});
    validateattributes(opt.FrameStep,{'numeric'},{'scalar','integer','positive'});
    validateattributes(opt.PauseTime,{'numeric'},{'scalar','finite','nonnegative'});
    validateattributes(opt.FrameRate,{'numeric'},{'scalar','finite','positive'});
    validateattributes(opt.Hip,{'numeric'},{'real','finite','numel',3});
    validateattributes(opt.DisplayRotation,{'numeric'},{'real','finite','size',[3,3]});
    if isempty(opt.DisplayAxisMap)
        opt.DisplayAxisMap = eye(3);
    end
    validateattributes(opt.DisplayAxisMap,{'numeric'},{'real','finite','size',[3,3]});
    if norm(opt.DisplayAxisMap.'*opt.DisplayAxisMap-eye(3),'fro') > 1e-10
        error('AnimateKneeBoneMuscle:DisplayAxisMap', ...
            'DisplayAxisMap must preserve lengths (an orthogonal 3-by-3 matrix).')
    end
    validateattributes(opt.AnimationPaddingX,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(opt.AnimationPaddingZ,{'numeric'},{'real','finite','scalar','nonnegative'});
    validateattributes(opt.StaticFullSkeleton,{'logical','numeric'},{'scalar','binary'});
    showFullSkeleton = opt.Static && opt.StaticFullSkeleton;
    if ~isempty(opt.FullSkeleton)
        validateattributes(opt.FullSkeleton,{'logical','numeric'},{'scalar','binary'});
        showFullSkeleton = logical(opt.FullSkeleton);
    end
    validateattributes(opt.StaticPadding,{'numeric'},{'real','finite','vector','nonnegative'});
    if ~ismember(numel(opt.StaticPadding),[1,3])
        error('AnimateKneeBoneMuscle:StaticPadding', ...
            'StaticPadding must be a scalar or an [X Y Z] vector in metres.')
    end
    staticPadding = opt.StaticPadding(:).';
    if isscalar(staticPadding)
        staticPadding = repmat(staticPadding,1,3);
    end
    validateTransforms(T,N,'T');
    opt.Hip = reshape(opt.Hip,1,3);
    p1 = reshape(p1,1,3);
    p2 = reshape(p2,1,3);
    Rplot = opt.DisplayRotation*opt.DisplayAxisMap;
    hasRoute = ~isempty(opt.Location);

    if hasRoute
        validateRoute(opt.Location,opt.CrossPoint,N,'Location');
        validateTransforms(opt.T_Pam,N,'T_Pam');
        endLabel = 'pEnd';
    else
        % Legacy two-point call: keep its original home-pose conversion.
        validateTransforms(T_ICR_t1,N,'T_ICR_t1');
        v2 = transformPoints(T_ICR_t1(:,:,pos),p2);
        endLabel = 'p2';
    end

    if iscell(Bifemsh)
        humans = Bifemsh(:);
    else
        humans = {Bifemsh};
    end
    nHuman = numel(humans);
    for h = 1:nHuman
        validateRoute(humans{h}.Location,humans{h}.Cross,N,'Human Location');
    end
    humanLabels = opt.HumanLabels;
    if isempty(humanLabels)
        if nHuman == 1
            humanLabels = {'Human muscle path'};
        else
            humanLabels = cellstr(compose('Human muscle %d',1:nHuman));
        end
    end
    if numel(humanLabels) ~= nHuman
        error('AnimateKneeBoneMuscle:HumanLabels', ...
            'Provide one HumanLabels entry per human muscle.')
    end
    humanLabels = cellstr(string(humanLabels));

    if opt.Static
        idxList = pos;
    elseif ~isempty(opt.FrameIndices)
        validateattributes(opt.FrameIndices,{'numeric'}, ...
            {'vector','integer','>=',1,'<=',N});
        idxList = reshape(opt.FrameIndices,1,[]);
    else
        idxFlex = pos:-opt.FrameStep:1;
        if idxFlex(end) ~= 1
            idxFlex(end+1) = 1;
        end
        idxReturn = 1:opt.FrameStep:pos;
        if idxReturn(end) ~= pos
            idxReturn(end+1) = pos;
        end
        idxList = [idxFlex,idxReturn(2:end)];
    end
    if (opt.ExportGif || opt.ExportVideo) && opt.Loop
        opt.Loop = false;
    end

    % Load extra meshes once only when full skeleton was selected.
    femur0 = loadCloud(opt.FemurFile,opt.UnitScale);
    tibia0 = loadCloud(opt.TibiaFile,opt.UnitScale);
    femurPlot = (femur0 + opt.Hip)*Rplot;
    extraBones = struct('Points',{},'Label',{},'Tag',{});
    footTibia = [];
    if showFullSkeleton
        [extraBones,footTibia] = skeletonMeshes(opt,T(:,:,pos),Rplot);
    end

    % Cache the small route arrays, not the large transformed bone meshes.
    % Static and moving plots use these exact same per-frame coordinates.
    robotPaths = cell(1,N);
    humanPaths = cell(nHuman,N);
    limits = [min(femurPlot,[],1);max(femurPlot,[],1)];
    for b = 1:numel(extraBones)
        limits = includePoints(limits,extraBones(b).Points);
    end
    poseIndices = unique([pos,idxList]);
    for i = poseIndices
        if hasRoute
            robotPaths{i} = routeFrame(opt.Location,opt.CrossPoint,opt.T_Pam(:,:,i),i);
        else
            robotPaths{i} = [p1;transformPoints(T(:,:,i),v2)];
        end
        limits = includePoints(limits,(robotPaths{i}+opt.Hip)*Rplot);
        for h = 1:nHuman
            humanPaths{h,i} = routeFrame(humans{h}.Location,humans{h}.Cross,T(:,:,i),i);
            limits = includePoints(limits,(humanPaths{h,i}+opt.Hip)*Rplot);
        end
        if opt.Static || showFullSkeleton || isempty(opt.XLim) || isempty(opt.YLim) || isempty(opt.ZLim)
            tibiaPlot = (transformPoints(T(:,:,i),tibia0)+opt.Hip)*Rplot;
            limits = includePoints(limits,tibiaPlot);
        end
        if showFullSkeleton
            footPlot = (transformPoints(T(:,:,i),footTibia)+opt.Hip)*Rplot;
            limits = includePoints(limits,footPlot);
        end
    end

    figureName = 'Knee bones and muscle routes - animation';
    if opt.Static
        figureName = 'Knee bones and muscle routes - static';
    end
    fig = figure('Name',figureName,'Color','w');
    setappdata(fig,'KneeBonePlotterSource',mfilename('fullpath'));
    setappdata(fig,'KneeBoneDisplayMap',Rplot);
    setappdata(fig,'KneeBoneFullSkeleton',showFullSkeleton);
    setappdata(fig,'isPaused',false);
    setappdata(fig,'stopAnimation',false);
    set(fig,'KeyPressFcn',@localKeyPress);
    ax = axes('Parent',fig);
    hold(ax,'on');
    axis(ax,'equal');
    view(ax,3);
    set(ax,'FontName','Arial','FontSize',10,'FontWeight','bold', ...
        'LineWidth',2,'Box','off','XMinorTick','on','YMinorTick','on', ...
        'ZMinorTick','on','TickLength',[0.025,0.05]);
    grid(ax,'off');

    % Filled markers stay visible even if a previous script set the root
    % MarkerEdgeColor to none. No global graphics defaults are changed here.
    hFemur = plot3(ax,femurPlot(:,1),femurPlot(:,2),femurPlot(:,3), ...
        'o','LineStyle','none','Color',[0.2,0.2,1.0], ...
        'MarkerSize',2,'MarkerFaceColor',[0.2,0.2,1.0], ...
        'MarkerEdgeColor','none','Tag','FemurMesh');
    hTibia = plot3(ax,NaN,NaN,NaN,'o','LineStyle','none', ...
        'Color',[0.1,0.1,0.1],'MarkerSize',2, ...
        'MarkerFaceColor',[0.1,0.1,0.1],'MarkerEdgeColor','none','Tag','TibiaMesh');
    hExtraLegend = gobjects(numel(extraBones),1);
    hFoot = gobjects(0);
    for b = 1:numel(extraBones)
        points = extraBones(b).Points;
        hBone = plot3(ax,points(:,1),points(:,2),points(:,3),'o', ...
            'LineStyle','none','Color','b','MarkerSize',2, ...
            'MarkerFaceColor','b','MarkerEdgeColor','none', ...
            'Tag',extraBones(b).Tag,'HandleVisibility','off');
        if strcmp(extraBones(b).Tag,'FootMesh')
            hFoot = hBone;
        end
        hExtraLegend(b) = plot3(ax,NaN,NaN,NaN,'o', ...
            'LineStyle','none','MarkerSize',8,'MarkerFaceColor','b', ...
            'MarkerEdgeColor','none');
    end
    hBPA = plot3(ax,NaN,NaN,NaN,'o-','Color',opt.RobotColor, ...
        'LineWidth',2,'MarkerSize',5,'MarkerFaceColor',opt.RobotColor, ...
        'MarkerEdgeColor','none','Tag','RobotRoute');
    hHuman = gobjects(nHuman,1);
    for h = 1:nHuman
        hHuman(h) = plot3(ax,NaN,NaN,NaN,'o-','Color',opt.HumanColor, ...
            'LineWidth',2,'MarkerSize',4,'MarkerFaceColor',opt.HumanColor, ...
            'MarkerEdgeColor','none','Tag',sprintf('HumanRoute%d',h));
    end
    hP1 = plot3(ax,NaN,NaN,NaN,'o','Color','r','MarkerFaceColor','r', ...
        'MarkerSize',8,'MarkerEdgeColor','none','Tag','RobotOrigin');
    hP2 = plot3(ax,NaN,NaN,NaN,'o','Color',opt.EndColor,'MarkerFaceColor',opt.EndColor, ...
        'MarkerSize',8,'MarkerEdgeColor','none','Tag','RobotEnd');
    xlabel(ax,'X, m');
    ylabel(ax,'Y, m');
    zlabel(ax,'Z, m');
    % Legend-only markers; leave the actual bone-cloud markers small.
    hFemurLegend = plot3(ax,NaN,NaN,NaN,'o', ...
        'LineStyle','none','MarkerSize',8, ...
        'MarkerFaceColor',hFemur.MarkerFaceColor, ...
        'MarkerEdgeColor','none');
    
    hTibiaLegend = plot3(ax,NaN,NaN,NaN,'o', ...
        'LineStyle','none','MarkerSize',8, ...
        'MarkerFaceColor',hTibia.MarkerFaceColor, ...
        'MarkerEdgeColor','none');
    lgd = legend(ax,[hFemurLegend;hTibiaLegend;hExtraLegend;hBPA;hHuman;hP1;hP2], ...
        [{'Femur','Tibia'}, {extraBones.Label}, {'Optimized BPA path'}, ...
        humanLabels(:).',{'p1',endLabel}], ...
        'Location','best','AutoUpdate','off');
    xlim(ax,plotLimits(opt.XLim,limits(:,1)));
    ylim(ax,plotLimits(opt.YLim,limits(:,2)));
    zlim(ax,plotLimits(opt.ZLim,limits(:,3)));
    if opt.Static || showFullSkeleton
        % Include every selected bone, even with old knee-only limits.
        % Animation fits the WHOLE sweep once, not separately per frame.
        scenePadding = zeros(1,3);
        if opt.Static
            scenePadding = staticPadding;
        end
        xlim(ax,[min(ax.XLim(1),limits(1,1))-scenePadding(1), ...
                 max(ax.XLim(2),limits(2,1))+scenePadding(1)]);
        ylim(ax,[min(ax.YLim(1),limits(1,2))-scenePadding(2), ...
                 max(ax.YLim(2),limits(2,2))+scenePadding(2)]);
        zlim(ax,[min(ax.ZLim(1),limits(1,3))-scenePadding(3), ...
                 max(ax.ZLim(2),limits(2,3))+scenePadding(3)]);
    end
    isYUpDisplay = norm([0,1,0]*Rplot-[0,1,0]) < 1e-10;
    if ~opt.Static && hasRoute && isYUpDisplay
        % Fixed limits for the entire animation; no camera tracking or zoom.
        % Leave -X and Y alone. Move +X outward and add lateral Z room.
        animationXLim = xlim(ax);
        animationZLim = zlim(ax);
        xlim(ax,animationXLim + [0,opt.AnimationPaddingX]);
        zlim(ax,animationZLim + [-opt.AnimationPaddingZ,opt.AnimationPaddingZ]);
    end
    
    xlabel(ax,'X, m');

    % if isequal(Rplot,eye(3))
    %     % Actual anatomical coordinates: no permutation or sign inversion.
    %     ylabel(ax,'Y, m');
    %     zlabel(ax,'Z, m');
    % 
    %     % Preserve the previous starting viewing direction, expressed
    %     % in anatomical coordinates instead of the old display coordinates.
    %     cameraDirection = [ ...
    %         cosd(30)*sind(-37.5), ...
    %         sind(30), ...
    %         cosd(30)*cosd(-37.5)];
    % 
    %     view(ax,cameraDirection);
    % else
    %     % Existing display convention for unchanged flexor calls.
    %     ylabel(ax,'-Z, m');
    %     zlabel(ax,'Y, m');
    % end
    
    anatomicalUp = [0,1,0]*Rplot;
    camup(ax,anatomicalUp);
    
    if opt.CameraOrbitDeg ~= 0
        camorbit(ax,opt.CameraOrbitDeg,0,'data',anatomicalUp);
    end
    % Place the three rulers at one lower, rear corner.
    % Applies to the anatomical XYZ display, with Y pointing up.
    % if isequal(Rplot,eye(3))
    %     cameraDirection = campos(ax) - camtarget(ax);
    % 
    %     % Choose the X/Z boundaries farther from the camera.
    %     xCorner = ax.XLim(1 + (cameraDirection(1) < 0));
    %     yCorner = ax.YLim(1);    % Bottom of the Y-up scene
    %     zCorner = ax.ZLim(1 + (cameraDirection(3) < 0));
    % 
    %     % X ruler: constant Y and Z.
    %     ax.XAxis.FirstCrossoverValue  = yCorner;
    %     ax.XAxis.SecondCrossoverValue = zCorner;
    % 
    %     % Y ruler: constant X and Z.
    %     ax.YAxis.FirstCrossoverValue  = xCorner;
    %     ax.YAxis.SecondCrossoverValue = zCorner;
    % 
    %     % Z ruler: constant X and Y.
    %     ax.ZAxis.FirstCrossoverValue  = xCorner;
    %     ax.ZAxis.SecondCrossoverValue = yCorner;
    % end
    drawPose(pos);
    drawnow;
    if opt.Static
        return
    end

    firstGifFrame = true;
    videoObj = [];
    if opt.ExportVideo
        videoObj = VideoWriter(opt.VideoFile,'MPEG-4');
        videoObj.FrameRate = opt.FrameRate;
        open(videoObj);
        videoCleanup = onCleanup(@() close(videoObj));
    end
    keepGoing = true;
    while keepGoing && isgraphics(fig)
        for i = idxList
            while isgraphics(fig) && getappdata(fig,'isPaused') && ...
                    ~getappdata(fig,'stopAnimation')
                drawnow;
                pause(0.05);
            end
            if ~isgraphics(fig) || getappdata(fig,'stopAnimation')
                keepGoing = false;
                break
            end
            drawPose(i);
            drawnow;
            if ~isgraphics(fig)
                keepGoing = false;
                break
            end
            if any(i == opt.PauseAtFrames)
                pause;
            end
            if ~isgraphics(fig)
                keepGoing = false;
                break
            end
            if opt.ExportGif || opt.ExportVideo
                frame = getframe(fig);
                if opt.ExportGif
                    [A,map] = rgb2ind(frame2im(frame),256);
                    if firstGifFrame
                        imwrite(A,map,opt.GifFile,'gif','LoopCount',Inf, ...
                            'DelayTime',1/opt.FrameRate);
                        firstGifFrame = false;
                    else
                        imwrite(A,map,opt.GifFile,'gif','WriteMode','append', ...
                            'DelayTime',1/opt.FrameRate);
                    end
                end
                if opt.ExportVideo
                    writeVideo(videoObj,frame);
                end
            end
            pause(opt.PauseTime);
        end
        keepGoing = keepGoing && opt.Loop;
    end

    function drawPose(i)
        tibiaPlot = (transformPoints(T(:,:,i),tibia0)+opt.Hip)*Rplot;
        robotPlot = (robotPaths{i}+opt.Hip)*Rplot;
        setPath(hTibia,tibiaPlot);
        if showFullSkeleton
            setPath(hFoot,(transformPoints(T(:,:,i),footTibia)+opt.Hip)*Rplot);
        end
        setPath(hBPA,robotPlot);
        for h = 1:nHuman
            setPath(hHuman(h),(humanPaths{h,i}+opt.Hip)*Rplot);
        end
        setPath(hP1,robotPlot(1,:));
        setPath(hP2,robotPlot(end,:));
        % Sum every segment; a wrapped route is not its endpoint chord.
        pathLength = sum(sqrt(sum(diff(robotPaths{i},1,1).^2,2)));
        title(ax,sprintf('Knee angle = %.1f deg, BPA path length = %.3f m', ...
            phi(i)*180/pi,pathLength),'FontName','Arial','FontSize',12, ...
            'FontWeight','bold');
    end
end

function [bones,footTibia] = skeletonMeshes(opt,Tpose,Rplot)
    % Same mesh origins and metre offsets as MuscleBonePlotting. Sacrum and
    % pelvis are already in the world/pelvis frame: do not add Hip to them.
    back = [-0.1007,0.0815,0];
    ankle = [0,-0.43,0];
    subtalar = [-0.04877,-0.04195,0.00792];
    mtp = [0.1788,-0.002,0.00108];
    pelvis = loadCloud(opt.PelvisFile,opt.UnitScale)*Rplot;
    sacrum = loadCloud(opt.SacrumFile,opt.UnitScale)*Rplot;
    spine = (loadCloud(opt.SpineFile,opt.UnitScale)+back)*Rplot;

    % Foot meshes first enter the tibia frame. Apply the SAME selected-pose
    % transform as the tibia, then Hip; do not add a second Knee offset.
    footTibia = [loadCloud(opt.TalusFile,opt.UnitScale)+ankle; ...
        loadCloud(opt.CalcaneusFile,opt.UnitScale)+ankle+subtalar; ...
        loadCloud(opt.ToesFile,opt.UnitScale)+ankle+subtalar+mtp];
    foot = (transformPoints(Tpose,footTibia)+opt.Hip)*Rplot;
    bones = struct('Points',{pelvis,sacrum,spine,foot}, ...
        'Label',{'Pelvis / hip','Sacrum','Lumbar spine','Foot'}, ...
        'Tag',{'PelvisMesh','SacrumMesh','SpineMesh','FootMesh'});
end

function validateTransforms(T,N,label)
    if ~isnumeric(T) || size(T,1) ~= 4 || size(T,2) ~= 4 || ...
            size(T,3) ~= N || ndims(T) > 3 || any(~isfinite(T(:)))
        error('AnimateKneeBoneMuscle:Transforms', ...
            '%s must be a finite 4-by-4-by-numel(phi) array.',label)
    end
end

function validateRoute(Location,Cross,N,label)
    if ~isnumeric(Location) || isempty(Location) || size(Location,2) ~= 3 || ...
            ndims(Location) > 3 || ...
            ~(size(Location,3) == N || size(Location,3) == 1) || ...
            any(~isfinite(Location(:)))
        error('AnimateKneeBoneMuscle:Location', ...
            '%s must be finite M-by-3 or M-by-3-by-numel(phi).',label)
    end
    validateattributes(Cross,{'numeric'}, ...
        {'scalar','integer','>=',1,'<=',size(Location,1)+1});
end

function points = loadCloud(file,scale)
    points = readmatrix(file)*scale;
    if size(points,2) < 3
        error('AnimateKneeBoneMuscle:BoneMesh','%s needs three coordinate columns.',file)
    end
    points = points(:,1:3);
    points = points(all(isfinite(points),2),:);
    if isempty(points)
        error('AnimateKneeBoneMuscle:BoneMesh','%s contains no finite points.',file)
    end
end

function points = transformPoints(T,points)
    % Batched equivalent of RowVecTrans(T,points(j,:)); T is one 4x4 matrix.
    points = points*T(1:3,1:3).' + T(1:3,4).';
end

function points = routeFrame(Location,Cross,T,i)
    frameIndex = min(i,size(Location,3));
    points = Location(:,:,frameIndex);
    points(Cross:end,:) = transformPoints(T,points(Cross:end,:));
end

function limits = includePoints(limits,points)
    limits = [min(limits(1,:),min(points,[],1)); ...
              max(limits(2,:),max(points,[],1))];
end

function limits = plotLimits(requested,dataLimits)
    if isempty(requested)
        padding = max(0.005,0.05*(dataLimits(2)-dataLimits(1)));
        limits = [dataLimits(1)-padding,dataLimits(2)+padding];
    else
        validateattributes(requested,{'numeric'},{'real','finite','numel',2});
        if requested(2) <= requested(1)
            error('AnimateKneeBoneMuscle:Limits','Axis limits must be increasing.')
        end
        limits = reshape(requested,1,2);
    end
end

function setPath(h,points)
    set(h,'XData',points(:,1),'YData',points(:,2),'ZData',points(:,3));
end

function localKeyPress(src,event)
    switch event.Key
        case 'space'
            setappdata(src,'isPaused',~getappdata(src,'isPaused'));
        case 'escape'
            setappdata(src,'stopAnimation',true);
    end
end
