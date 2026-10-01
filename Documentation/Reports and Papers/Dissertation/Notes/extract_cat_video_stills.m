scriptDir = fileparts(mfilename('fullpath'));
dissertationDir = fileparts(scriptDir);
videoPath = fullfile(dissertationDir, 'CHATGPT_staging_20260925', ...
    'Decerebrate Cat walks and exhibits multiple gait patterns.mp4');
outputDir = fullfile(scriptDir, 'cat_video_stills');
if ~exist(outputDir, 'dir'), mkdir(outputDir); end
vr = VideoReader(videoPath);
sampleTimes = linspace(0, max(0, vr.Duration - 1 / vr.FrameRate), 20);
fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1800 1200]);
tl = tiledlayout(fig, 4, 5, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('Supplied cat footage | %.2f s | %.1f fps', vr.Duration, vr.FrameRate));
for j = 1:numel(sampleTimes)
    vr.CurrentTime = sampleTimes(j);
    frame = readFrame(vr);
    imwrite(frame, fullfile(outputDir, sprintf('cat_frame_%02d.png', j)));
    ax = nexttile(tl); image(ax, frame); axis(ax, 'image'); axis(ax, 'off');
    title(ax, sprintf('%02d | %.2f s', j, sampleTimes(j)), 'FontSize', 10);
end
exportgraphics(fig, fullfile(outputDir, 'contact_sheet.png'), 'Resolution', 140);
close(fig);
fprintf('Duration %.3f s; rate %.3f fps; frame %d x %d\n', vr.Duration, vr.FrameRate, vr.Width, vr.Height);
disp(sampleTimes);
