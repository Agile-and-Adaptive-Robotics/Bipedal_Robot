scriptDir = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(fileparts(fileparts(scriptDir))));
videoDir = fullfile(repoRoot, 'Pictures', 'PictureAnalysis');
outputDir = fullfile(scriptDir, 'video_contact_sheets');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

videos = {
    '20230510_161942.mp4'
    '20230419_230001_1.mp4'
};

for i = 1:numel(videos)
    videoPath = fullfile(videoDir, videos{i});
    vr = VideoReader(videoPath);
    sampleTimes = linspace(0, max(0, vr.Duration - 1 / vr.FrameRate), 12);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 1600 1000]);
    tl = tiledlayout(fig, 3, 4, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s | %.2f s | %.1f fps', videos{i}, vr.Duration, vr.FrameRate), ...
        'Interpreter', 'none', 'FontSize', 14, 'FontWeight', 'bold');

    for j = 1:numel(sampleTimes)
        vr.CurrentTime = sampleTimes(j);
        frame = readFrame(vr);
        ax = nexttile(tl);
        image(ax, frame);
        axis(ax, 'image');
        axis(ax, 'off');
        title(ax, sprintf('%.2f s', sampleTimes(j)), 'FontSize', 10);
    end

    [~, stem] = fileparts(videos{i});
    outPath = fullfile(outputDir, [stem '_contact.png']);
    exportgraphics(fig, outPath, 'Resolution', 160);
    close(fig);
    fprintf('%s\n', outPath);
end
