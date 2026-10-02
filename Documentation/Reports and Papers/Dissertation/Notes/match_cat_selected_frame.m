scriptDir = fileparts(mfilename('fullpath'));
dissertationDir = fileparts(scriptDir);
videoPath = fullfile(dissertationDir, 'CHATGPT_staging_20260925', ...
    'Decerebrate Cat walks and exhibits multiple gait patterns.mp4');
shot = imread('C:\Users\BENBOL~1\AppData\Local\Temp\codex-clipboard-54f0c016-cc73-4ccc-8dd3-8221d9f4d3bb.png');
fprintf('Screenshot size: %d x %d\n', size(shot,2), size(shot,1));
% Displayed video rectangle, excluding pillarbox bars and transport controls.
target = mean(double(shot(1:910,337:1568,1:3)),3);
[xq,yq] = meshgrid(linspace(1,size(target,2),320),linspace(1,size(target,1),240));
target = interp2(target,xq,yq,'linear');
roiRows = 45:205; roiCols = 40:280;
targetROI = target(roiRows,roiCols); targetROI = targetROI(:);
targetROI = (targetROI-mean(targetROI))/std(targetROI);
vr = VideoReader(videoPath);
scores = []; times = []; frames = {}; index = 0;
while hasFrame(vr)
    time = vr.CurrentTime; frame = readFrame(vr); index = index+1;
    gray = mean(double(frame),3); roi = gray(roiRows,roiCols); roi = roi(:);
    roi = (roi-mean(roi))/std(roi);
    scores(index) = mean((roi-targetROI).^2);
    times(index) = time; frames{index} = frame;
end
[~,order] = sort(scores);
best = order(1);
fprintf('Best match: frame %d, %.6f seconds, score %.8f\n',best,times(best),scores(best));
disp([order(1:10)',times(order(1:10))',scores(order(1:10))']);
outputDir = fullfile(scriptDir,'cat_selected_match');
if ~exist(outputDir,'dir'), mkdir(outputDir); end
imwrite(frames{best},fullfile(outputDir,'selected_original.png'));
save(fullfile(outputDir,'match_metadata.mat'),'scores','times','best');
fig = figure('Visible','off','Color','w','Position',[100 100 1400 1000]);
tl = tiledlayout(fig,3,3,'TileSpacing','compact','Padding','compact');
for i = 1:9
    k=order(i); ax=nexttile(tl); image(ax,frames{k}); axis(ax,'image'); axis(ax,'off');
    title(ax,sprintf('Frame %d | %.3f s | score %.4f',k,times(k),scores(k)));
end
exportgraphics(fig,fullfile(outputDir,'match_contact_sheet.png'),'Resolution',130); close(fig);
