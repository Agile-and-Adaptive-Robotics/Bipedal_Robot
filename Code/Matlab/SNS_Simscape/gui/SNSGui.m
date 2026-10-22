classdef SNSGui < handle
%SNSGUI  Block-library GUI for building SNS Simulink networks.
%
% Ben's ask (2026-10-21): an AnimatLab-style network builder that is
% easier and faster — palette of the SNS_Library blocks, place-and-wire
% canvas in the library's own icon language, every mask value editable in
% a property panel, model generation through the EXISTING
% sns_build_from_json machinery, and loading of the existing network
% JSONs (the tuned production network and the winner slice among them).
%
% Connections follow the 2026-09-22 architecture: a connection IS a
% one-input synapse onto the target's syn-port stack, so drawing a
% connection creates the synapse, and its canvas marker (white triangle
% = excitatory, black dot = inhibitory) derives from the Esyn sign,
% exactly like the library icons.
%
% Programmatic uifigure class on purpose: .m code only, version
% controllable, no binary .mlapp. Run:  app = SNSGui()
% then try  app.demoKneeReflex()   (builds a demo circuit) or
%          app.loadJson('<path to a network json>')
%
% Metrics instrumentation: every user-equivalent click increments
% app.nClicks; load/build/save timings land in app.timings. The demo
% builder reports clicks-to-circuit so "easier and faster" stays a
% measured claim (README_SNS_Simscape.md, GUI section).

properties
    f                    % uifigure
    ax                   % canvas uiaxes
    lstPalette           % uilistbox
    propPanel            % uigridlayout for the property editor
    propChildren = {}    % components of the current property panel
    btnWiring            % connect-mode toggle
    lblStatus            % status line
    selType = 'none'     % 'block' | 'edge' | 'none'
    selIdx = 0
    mode = 'idle'        % 'idle' | 'place' | 'wire'
    wireSrc = 0          % source block index while wiring
    % ---- circuit model ----
    blocks = struct('name', {}, 'type', {}, 'x', {}, 'y', {}, ...
                    'params', {})          % params: {name; value} cell col
    edges = struct('src', {}, 'dst', {}, 'x', {}, 'y', {}, ...
                   'params', {})           % synapse connections
    % ---- metrics ----
    nClicks = 0
    timings = struct('loadJson', [], 'build', [], 'saveJson', [])
    nextId = 1
end

properties (Constant)
    % palette: type, display label, params {name default}, glyph kind
    Catalog = { ...
      'NonSpikingNeuron',  'Graded neuron',        {'Vrest' -52; 'Gm' 0.1; 'Cm' 5; 'Thr' -45; 'Slope' 5}, 'neuron'
      'NonSpikingSynapse', 'Graded synapse',       {'gmax' 1; 'Esyn' -72; 'ThrPre' -55; 'SlopePre' 1},    'synapse'
      'SpikingLIFNeuron',  'Spiking LIF neuron',   {'Vrest' -70; 'Vth' 30; 'Rm' 10; 'tau' 10},            'lif'
      'IaMuscleSpindle',   'Ia spindle encoder',   {'Imax' 10; 'Wl' 1; 'Wv' 0.5},                         'ia'
      'IbGolgiTendon',     'Ib tendon encoder',    {'Imax' 10; 'Kf' 1},                                   'ib'
      'MuscleActivation',  'Activation dynamics',  {'tauAct' 50; 'A0' 0},                                 'pentagon'
      'SynSum',            'Synapse-sum junction', {},                                                    'junction'
      'BPA_10mm',          'BPA 10 mm',            {'Rest' 0.20; 'Kmax' 0.165; 'TendonL' 0; 'FittingL' 0}, 'mckibben'
      'BPA_20mm',          'BPA 20 mm',            {'Rest' 0.20; 'Kmax' 0.165; 'TendonL' 0; 'FittingL' 0}, 'mckibben'
      'BPA_40mm',          'BPA 40 mm',            {'Rest' 0.20; 'Kmax' 0.165; 'TendonL' 0; 'FittingL' 0}, 'mckibben'
      'BioMuscle',         'Biological muscle',    {'Fmax' 1000; 'l_opt' 0.10; 'l_slack' 0.05; 'alpha' 0; 'v_max' 10; 'A_hill' 0.25; 'kpe' 4; 'e_pas' 0.6}, 'pennate'
      'BPAForce',          'Generic BPA force',    {'Fmax' 500; 'epsMax' 0.25},                           'mckibben'
      'SpikingNeuron',     'Spiking neuron',       {'Vrest' -70; 'Vth' 30; 'Vreset' -75; 'Gm' 0.1; 'Cm' 5}, 'lif'
      'SpikingSynapse',    'Spiking synapse',      {'gmax' 1; 'ginc' 0.1; 'tau_syn' 20; 'Esyn' 0},        'synapse'
      'HybridSpikingSynapse', 'Hybrid synapse',    {'gmax' 1; 'ginc' 0.1; 'tau_syn' 20; 'Esyn' 0; 'ThrPre' -55}, 'synapse'}
    SynapseDefaults = {'gmax' 1; 'Esyn' -72; 'ThrPre' -55; 'SlopePre' 1}
end

methods
    function app = SNSGui()
        % SNS_Library + builders must be on the path for Build/Load
        here = fileparts(mfilename('fullpath'));
        addpath(fileparts(here));          % SNS_Simscape root
        addpath(fullfile(fileparts(here), 'demos'));
        app.f = uifigure('Name', 'SNS Circuit Builder (SNS_Library)', ...
            'Position', [60 60 1280 760], 'Color', [0.97 0.97 0.97]);
        grid = uigridlayout(app.f, [2 3], ...
            'ColumnWidth', {230, '1x', 300}, 'RowHeight', {36, '1x'});
        % ---- toolbar ----
        uibutton(grid, 'push', 'Text', 'New', 'FontWeight', 'bold', ...
            'ButtonPushedFcn', @(s, e) app.newCircuit());
        uibutton(grid, 'push', 'Text', 'Load JSON...', ...
            'ButtonPushedFcn', @app.cbLoad);
        uibutton(grid, 'push', 'Text', 'Build Simulink model', ...
            'FontWeight', 'bold', 'ButtonPushedFcn', @app.cbBuild);
        palWrap = uigridlayout(grid);
        palWrap.Layout.Row = 2; palWrap.Layout.Column = 1;
        palWrap.RowHeight = {20, '1x'}; palWrap.Padding = [2 2 2 2];
        uilabel(palWrap, 'Text', 'SNS_Library palette (pick, then click canvas)', ...
            'FontName', 'Arial', 'FontWeight', 'bold');
        app.lstPalette = uilistbox(palWrap, ...
            'Items', app.Catalog(:, 2)', 'Value', 'Graded neuron', ...
            'FontSize', 11);
        app.lstPalette.ValueChangedFcn = @app.cbPalettePick;
        app.btnWiring = uibutton(grid, 'state', 'Text', 'Connect mode', ...
            'ValueChangedFcn', @app.cbWireMode);
        app.btnWiring.Layout.Row = 2; app.btnWiring.Layout.Column = 2;
        app.lblStatus = uilabel(grid, 'Text', 'Ready.', ...
            'FontName', 'Arial', 'FontColor', [0.25 0.25 0.25]);
        app.lblStatus.Layout.Row = 2; app.lblStatus.Layout.Column = 3;
        % ---- canvas ----
        app.ax = uiaxes(grid, 'Box', 'off');
        app.ax.Layout.Row = 3; app.ax.Layout.Column = 2;
        hold(app.ax, 'on');
        app.ax.XLim = [0 100]; app.ax.YLim = [0 100];
        app.ax.XGrid = 'on'; app.ax.YGrid = 'on'; app.ax.GridAlpha = 0.10;
        axis(app.ax, 'equal');
        app.ax.ButtonDownFcn = @app.cbCanvasClick;
        % ---- property panel ----
        app.propPanel = uigridlayout(grid);
        app.propPanel.Layout.Row = 3; app.propPanel.Layout.Column = 3;
        app.propPanel.RowHeight = {22, 22, 26, 26, 26, 26, 26, 26, ...
            26, 26, 26, 26};
        app.propPanel.ColumnWidth = {105, '1x'};
        h = uilabel(app.propPanel, 'Text', 'Properties', ...
            'FontWeight', 'bold', 'FontSize', 12);
        h.Layout.Row = 1; h.Layout.Column = 1;
        app.draw();
    end

    % ------------------------------------------------ palette / placement
    function cbPalettePick(app, ~, ~)
        app.nClicks = app.nClicks + 1;
        app.mode = 'place';
        app.setStatus(sprintf('Click the canvas to place %s.', ...
            app.lstPalette.Value));
    end

    function cbCanvasClick(app, ~, ~)
        app.nClicks = app.nClicks + 1;
        cp = app.ax.CurrentPoint;
        p = cp(1, 1:2);
        switch app.mode
            case 'place'
                lbl = app.lstPalette.Value;
                t = app.Catalog{strcmp(app.Catalog(:, 2), lbl), 1};
                app.addBlock(t, p(1), p(2));
                app.mode = 'idle';
                app.setStatus(sprintf('Placed %s.', t));
            case 'wire'
                hit = app.hitBlock(p);
                if hit == 0, return; end
                if app.wireSrc == 0
                    app.wireSrc = hit;
                    app.setStatus('Connect: now click the TARGET block.');
                elseif hit ~= app.wireSrc
                    app.addEdgeFromBlocks(app.wireSrc, hit);
                    app.wireSrc = 0;
                end
            otherwise
                app.selType = 'none'; app.selIdx = 0;
                app.draw(); app.fillPropertyPanel();
        end
    end

    function cbWireMode(app, btn, ~)
        app.nClicks = app.nClicks + 1;
        if btn.Value
            app.mode = 'wire'; app.wireSrc = 0;
            app.setStatus('Connect: click SOURCE block, then TARGET.');
        else
            app.mode = 'idle'; app.wireSrc = 0;
            app.setStatus('Connect mode off.');
        end
    end

    function cbBlockClick(app, idx)
        app.nClicks = app.nClicks + 1;
        if app.mode == 'wire'
            if app.wireSrc == 0
                app.wireSrc = idx;
                app.setStatus('Connect: now click the TARGET block.');
            elseif idx ~= app.wireSrc
                app.addEdgeFromBlocks(app.wireSrc, idx);
                app.wireSrc = 0;
            end
            return;
        end
        app.selType = 'block'; app.selIdx = idx;
        app.draw(); app.fillPropertyPanel();
    end

    function cbEdgeClick(app, idx)
        app.nClicks = app.nClicks + 1;
        app.selType = 'edge'; app.selIdx = idx;
        app.draw(); app.fillPropertyPanel();
    end

    % ------------------------------------------------------- circuit data
    function addBlock(app, type, x, y)
        c = app.Catalog(strcmp(app.Catalog(:, 1), type), :);
        if isempty(c)
            app.setStatus(sprintf('Unknown block type %s.', type)); return;
        end
        nm = sprintf('%s_%d', type, app.nextId); app.nextId = app.nextId + 1;
        b = struct('name', nm, 'type', type, 'x', x, 'y', y, ...
                   'params', {c{3}});
        app.blocks(end+1) = b; %#ok<AGROW>
        app.selType = 'block'; app.selIdx = numel(app.blocks);
        app.draw(); app.fillPropertyPanel();
    end

    function addEdgeFromBlocks(app, src, dst)
        a = app.blocks(src); b = app.blocks(dst);
        e = struct('src', src, 'dst', dst, ...
            'x', (a.x + b.x) / 2 - 5, 'y', (a.y + b.y) / 2, ...
            'params', {app.SynapseDefaults});
        app.edges(end+1) = e; %#ok<AGROW>
        app.selType = 'edge'; app.selIdx = numel(app.edges);
        app.draw(); app.fillPropertyPanel();
        app.setStatus(sprintf('%s -> %s: synapse added. Esyn < 0 is inhibitory.', ...
            a.name, b.name));
    end

    function newCircuit(app)
        app.blocks = struct('name', {}, 'type', {}, 'x', {}, 'y', {}, ...
                            'params', {});
        app.edges = struct('src', {}, 'dst', {}, 'x', {}, 'y', {}, ...
                           'params', {});
        app.selType = 'none'; app.selIdx = 0; app.nextId = 1;
        app.draw(); app.fillPropertyPanel();
        app.setStatus('New circuit.');
    end

    % --------------------------------------------------------- rendering
    function draw(app)
        cla(app.ax); hold(app.ax, 'on');
        for k = 1:numel(app.edges)
            e = app.edges(k);
            src = app.blocks(e.src); dst = app.blocks(e.dst);
            sel = strcmp(app.selType, 'edge') && app.selIdx == k;
            app.drawEdge(k, src, dst, e, sel);
        end
        for k = 1:numel(app.blocks)
            b = app.blocks(k);
            sel = strcmp(app.selType, 'block') && app.selIdx == k;
            app.drawBlockGlyph(k, b, sel);
        end
        if isempty(app.blocks) && isempty(app.edges)
            text(app.ax, 50, 50, 'Pick a palette block, then click here.', ...
                'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                'FontSize', 11, 'Color', [0.55 0.55 0.55]);
        end
    end

    function drawEdge(app, k, src, dst, e, sel)
        % Orthogonal route with per-edge lanes so wires read cleanly:
        % leave the source's right edge, run right of both blocks, drop
        % to the target row, and enter the target's left edge.
        lw = 1.4; if sel, lw = 2.6; end
        ys = src.y; ye = dst.y;
        xs = src.x + 2.8;
        xe = dst.x - 6.0;
        lane = dst.x - 9.0 - 2.5 * mod(k - 1, 4);
        lane = min(lane, xe - 1.5);
        poly = [xs ys; lane ys; lane ye; xe ye];
        plot(app.ax, poly(:, 1), poly(:, 2), '-', ...
            'Color', [0.25 0.25 0.25], 'LineWidth', lw, 'HitTest', 'off');
        mx = xe + 0.6; my = ye;   % marker just before the target
        esyn = app.paramValue(e.params, 'Esyn', -72);
        isCurrent = ~strcmp(src.type, 'NonSpikingNeuron');
        if isCurrent
            % afferent-current edge: plain line, small injection dot
            plot(app.ax, mx, my, 'o', 'MarkerSize', 3, ...
                'MarkerEdgeColor', [0.25 0.25 0.25], 'MarkerFaceColor', ...
                [0.25 0.25 0.25], 'HitTest', 'off');
        elseif esyn >= 0    % excitatory: white triangle, black edge
            patch(app.ax, 'Faces', 1:3, 'Vertices', ...
                [mx-1.6 my-1.4; mx+1.2 my; mx-1.6 my+1.4], ...
                'FaceColor', 'w', 'EdgeColor', 'k', 'LineWidth', 1.0, ...
                'HitTest', 'off');
        else            % inhibitory: solid black dot
            th = linspace(0, 2*pi, 17);
            patch(app.ax, mx + 1.1 * cos(th), my + 1.1 * sin(th), ...
                [0.05 0.05 0.05], 'EdgeColor', 'none', 'HitTest', 'off');
        end
        % click target: near-transparent fat dot over the marker
        % (skipped on large networks to keep loads fast)
        if numel(app.edges) <= 200
            th2 = linspace(0, 2*pi, 13);
            patch(app.ax, mx + 2.4 * cos(th2), my + 2.4 * sin(th2), ...
                [1 1 1], 'FaceAlpha', 0.02, 'EdgeColor', 'none', ...
                'ButtonDownFcn', @(s, ev) app.cbEdgeClick(k));
        end
    end

    function drawBlockGlyph(app, k, b, sel)
        kind = app.Catalog{strcmp(app.Catalog(:, 1), b.type), 4};
        x = b.x; y = b.y;
        showLabel = numel(app.blocks) <= 40;
        switch kind
            case 'neuron'
                th = linspace(0, 2*pi, 33);
                app.gPatch(k, x + 2.4*cos(th), y + 2.4*sin(th), ...
                    [0.99 0.96 0.81]);
                xw = linspace(x-1.5, x+1.5, 25);
                yw = y - 0.8 + 1.4*exp(-((xw - x + 0.25)/0.7).^2);
                plot(app.ax, xw, yw, 'k-', 'LineWidth', 0.8, 'HitTest', 'off');
                plot(app.ax, [x-1.5 x+1.5], [y-0.8 y-0.8], 'k-', ...
                    'LineWidth', 0.8, 'HitTest', 'off');
            case 'lif'
                th = linspace(0, 2*pi, 33);
                app.gPatch(k, x + 2.4*cos(th), y + 2.4*sin(th), ...
                    [0.99 0.96 0.81]);
                xs = x - 1.7 + 0.34*(0:9);
                ys = y - 0.8 * ones(size(xs));
                ys(2:2:end) = y + 1.4;
                plot(app.ax, xs, ys, 'k-', 'LineWidth', 0.8, 'HitTest', 'off');
            case 'ia'
                app.gCapsule(k, x, y, 'Ia', [0.90 0.93 0.87], 3.0, 1.6);
            case 'ib'
                app.gCapsule(k, x, y, 'Ib', [0.93 0.93 0.90], 3.0, 1.6);
            case 'pentagon'
                th = linspace(pi/2, pi/2 + 2*pi, 6);
                app.gPatch(k, x + 2.3*cos(th), y + 2.3*sin(th), ...
                    [0.96 0.78 0.80]);
            case 'junction'
                app.gPatch(k, x + [-2 -2 2 2], y + [-1.5 1.5 1.5 -1.5], ...
                    [0.93 0.93 0.93]);
            case 'mckibben'
                app.gMuscle(k, x, y, 6.0, 2.5, sns_icon_geom('BPA'));
            case 'pennate'
                app.gMuscle(k, x, y, 6.0, 2.7, sns_icon_geom('BIO'));
            otherwise  % standalone synapse block
                plot(app.ax, [x-1.6 x+0.6], [y y], 'k-', ...
                    'LineWidth', 1.6, 'HitTest', 'off');
                patch(app.ax, 'Faces', 1:3, 'Vertices', ...
                    [x+0.6 y-1.0; x+2.2 y; x+0.6 y+1.0], ...
                    'FaceColor', 'w', 'EdgeColor', 'k', 'HitTest', 'off');
        end
        if sel
            plot(app.ax, x + [-4.4 4.4 4.4 -4.4 -4.4], ...
                 y + [3.6 3.6 -3.6 -3.6 3.6], '-', ...
                 'Color', [0.00 0.45 0.74], 'LineWidth', 1.6, 'HitTest', 'off');
        end
        if showLabel
            text(app.ax, x, y - 4.0, strrep(b.type, '_', '\_'), ...
                'HorizontalAlignment', 'center', 'FontName', 'Arial', ...
                'FontSize', 8, 'HitTest', 'off');
        end
    end

    function gPatch(app, k, px, py, fc)
        patch(app.ax, px, py, fc, 'EdgeColor', 'k', 'LineWidth', 1.4, ...
            'ButtonDownFcn', @(s, ev) app.cbBlockClick(k));
    end

    function gCapsule(app, k, x, y, lbl, fc, hw, hh)
        t = linspace(0, 1, 24);
        xw = [linspace(x-hw, x+hw, 24), fliplr(linspace(x-hw, x+hw, 24))];
        yy = [hh*(1-(2*t-1).^2), -fliplr(hh*(1-(2*t-1).^2))];
        app.gPatch(k, xw, y + yy, fc);
        text(app.ax, x, y, lbl, 'HorizontalAlignment', 'center', ...
            'FontName', 'Arial', 'FontSize', 8, 'HitTest', 'off');
    end

    function gMuscle(app, k, x, y, w, h, g) %#ok<INUSD>
        % side-view McKibben / pennate muscle from the shared geometry
        sx = w / (2 * 0.95); sy = h / (2 * 0.42);
        mx = @(xq) x + xq * sx; my = @(yq) y + yq * sy;
        patch(app.ax, mx(g.bodyX), my(g.bodyY), g.fill, ...
            'EdgeColor', 'k', 'LineWidth', 1.2, ...
            'ButtonDownFcn', @(s, ev) app.cbBlockClick(k));
        segs = [];
        if isfield(g, 'braid'), segs = g.braid; end
        if isfield(g, 'fibers'), segs = g.fibers; end
        for i = 1:size(segs, 1)
            plot(app.ax, mx(segs(i, [1 3])), my(segs(i, [2 4])), ...
                'k-', 'LineWidth', 0.5, 'HitTest', 'off');
        end
        if isfield(g, 'capR')
            patch(app.ax, mx(g.capR(:, 1)), my(g.capR(:, 2)), g.cap, ...
                'EdgeColor', 'k', 'LineWidth', 1.0, 'HitTest', 'off');
            patch(app.ax, mx(g.capL(:, 1)), my(g.capL(:, 2)), g.cap, ...
                'EdgeColor', 'k', 'LineWidth', 1.0, 'HitTest', 'off');
        end
        if isfield(g, 'tendonR')
            patch(app.ax, mx(g.tendonR(:, 1)), my(g.tendonR(:, 2)), ...
                g.tendon, 'EdgeColor', 'k', 'LineWidth', 1.0, 'HitTest', 'off');
            patch(app.ax, mx(g.tendonL(:, 1)), my(g.tendonL(:, 2)), ...
                g.tendon, 'EdgeColor', 'k', 'LineWidth', 1.0, 'HitTest', 'off');
        end
    end

    % ---------------------------------------------------- property editor
    function fillPropertyPanel(app)
        for k = 1:numel(app.propChildren)
            if isvalid(app.propChildren{k}), delete(app.propChildren{k}); end
        end
        app.propChildren = {};
        if strcmp(app.selType, 'none')
            h = uilabel(app.propPanel, 'Text', ...
                'Select a block or a connection marker.', ...
                'FontName', 'Arial');
            app.place(h, 3, 1); app.track(h);
            return;
        end
        if strcmp(app.selType, 'block')
            b = app.blocks(app.selIdx);
            ttl = sprintf('%s   (%s)', b.name, b.type);
        else
            e = app.edges(app.selIdx);
            ttl = sprintf('%s -> %s   (synapse)', ...
                app.blocks(e.src).name, app.blocks(e.dst).name);
        end
        h = uilabel(app.propPanel, 'Text', ttl, 'FontName', 'Arial', ...
            'FontWeight', 'bold');
        app.place(h, 2, 1); app.track(h);
        if strcmp(app.selType, 'block'), P = app.blocks(app.selIdx).params;
        else, P = app.edges(app.selIdx).params; end
        for r = 1:size(P, 1)
            h = uilabel(app.propPanel, 'Text', P{r, 1}, ...
                'FontName', 'Arial');
            app.place(h, r + 2, 1); app.track(h);
            h = uieditfield(app.propPanel, 'numeric', ...
                'Value', P{r, 2}, 'FontName', 'Arial', ...
                'ValueChangedFcn', @(s, ev) app.cbParamEdit(P{r, 1}));
            app.place(h, r + 2, 2); app.track(h);
        end
    end

    function place(app, h, r, c) %#ok<INUSD>
        h.Layout.Row = r; h.Layout.Column = c;
    end

    function h = track(app, h)
        app.propChildren{end+1} = h; %#ok<AGROW>
    end

    function cbParamEdit(app, pname)
        app.nClicks = app.nClicks + 1;
        row = find(strcmp(app.propChildren2Names(), pname), 1);
        if isempty(row), return; end
        v = app.propFieldsByOrder{row}.Value;
        if strcmp(app.selType, 'block')
            P = app.blocks(app.selIdx).params;
            P{strcmp(P(:, 1), pname), 2} = v;
            app.blocks(app.selIdx).params = P;
        else
            P = app.edges(app.selIdx).params;
            P{strcmp(P(:, 1), pname), 2} = v;
            app.edges(app.selIdx).params = P;
        end
        app.draw();   % Esyn sign flips the connection marker
        app.setStatus(sprintf('%s = %g', pname, v));
    end

    % ------------------------------------------------------ JSON + build
    function cbLoad(app, ~, ~)
        app.nClicks = app.nClicks + 1;
        [fn, pn] = uigetfile('*.json', 'Pick a network JSON');
        if isequal(fn, 0), return; end
        app.loadJson(fullfile(pn, fn));
    end

    function loadJson(app, path)
        t0 = tic;
        J = jsondecode(fileread(path));
        N = asStructArray(J.neurons);
        SY = asStructArray(J.synapses);
        nN = numel(N); nS = numel(SY);
        app.newCircuit();
        for k = 1:nN
            nm = N(k).name;
            if contains(nm, 'MN_'), c = 5;
            elseif startsWith(nm, {'Ia_', 'II_', 'Ib_'}), c = 4;
            elseif startsWith(nm, 'IBEXC'), c = 3;
            elseif startsWith(nm, 'PF'), c = 2;
            else, c = 1; end
            P = {'Vrest' -52; 'Gm' 0.1; 'Cm' 5; 'Thr' -45; 'Slope' 5};
            if isfield(N, 'Vrest_mV'), P{1, 2} = N(k).Vrest_mV; end
            if isfield(N, 'Gm_uS'),    P{2, 2} = N(k).Gm_uS; end
            if isfield(N, 'Cm_nF'),    P{3, 2} = N(k).Cm_nF; end
            app.blocks(end+1) = struct('name', nm, ...
                'type', 'NonSpikingNeuron', ...
                'x', 8 + 17 * c, 'y', 90 - 3.4 * mod(k, 24), ...
                'params', {P}); %#ok<AGROW>
        end
        pos = containers.Map;
        for k = 1:numel(app.blocks), pos(app.blocks(k).name) = k; end
        for k = 1:nS
            s = SY(k);
            if isKey(pos, s.src) && isKey(pos, s.dst)
                a = app.blocks(pos(s.src)); b = app.blocks(pos(s.dst));
                app.edges(end+1) = struct('src', pos(s.src), ...
                    'dst', pos(s.dst), 'x', (a.x + b.x)/2 - 5, ...
                    'y', (a.y + b.y)/2, 'params', ...
                    {{'gmax', s.g_uS; 'Esyn', s.Esyn_mV; ...
                      'ThrPre', s.ThrPre_mV; ...
                      'SlopePre', s.SlopePre_mV}}); %#ok<AGROW>
            end
        end
        app.timings.loadJson = toc(t0);
        app.selType = 'none';
        app.draw(); app.fillPropertyPanel();
        app.setStatus(sprintf(['Loaded: %d neurons, %d synapses from %s ' ...
            '(%.2f s). Block labels hide above 40 cells.'], nN, nS, ...
            path, app.timings.loadJson));
    end

    function cbSave(app, ~, ~)
        app.nClicks = app.nClicks + 1;
        [fn, pn] = uiputfile('circuit.json', 'Save network JSON');
        if isequal(fn, 0), return; end
        t0 = tic;
        app.saveJson(fullfile(pn, fn));
        app.timings.saveJson = toc(t0);
        app.setStatus('JSON saved.');
    end

    function saveJson(app, path)
        S = struct();
        S.meta = struct('source', 'SNSGui canvas circuit', 'app', 'SNSGui');
        neurons = struct('name', {}, 'Vrest_mV', {}, 'Gm_uS', {}, ...
                         'Cm_nF', {});
        for k = 1:numel(app.blocks)
            b = app.blocks(k);
            if ~strcmp(b.type, 'NonSpikingNeuron'), continue; end
            neurons(end+1) = struct('name', b.name, ...
                'Vrest_mV', app.paramValue(b.params, 'Vrest', -52), ...
                'Gm_uS', app.paramValue(b.params, 'Gm', 0.1), ...
                'Cm_nF', app.paramValue(b.params, 'Cm', 5)); %#ok<AGROW>
        end
        S.neurons = neurons;
        syn = struct('src', {}, 'dst', {}, 'g_uS', {}, 'Esyn_mV', {}, ...
                     'ThrPre_mV', {}, 'SlopePre_mV', {});
        for k = 1:numel(app.edges)
            e = app.edges(k);
            if ~strcmp(app.blocks(e.src).type, 'NonSpikingNeuron'), continue; end
            if ~strcmp(app.blocks(e.dst).type, 'NonSpikingNeuron'), continue; end
            syn(end+1) = struct('src', app.blocks(e.src).name, ...
                'dst', app.blocks(e.dst).name, ...
                'g_uS', app.paramValue(e.params, 'gmax', 1), ...
                'Esyn_mV', app.paramValue(e.params, 'Esyn', -72), ...
                'ThrPre_mV', app.paramValue(e.params, 'ThrPre', -55), ...
                'SlopePre_mV', app.paramValue(e.params, 'SlopePre', 1)); %#ok<AGROW>
        end
        S.synapses = syn;
        % Afferent-sourced edges (Ia/Ib/spindle/tendon blocks emit
        % currents) export as input currents onto the target neuron's
        % Iapp port — the same treatment the production network gives
        % its Ia_/II_/Ib_ inputs.
        ins = struct('dst', {});
        for k = 1:numel(app.edges)
            e = app.edges(k);
            if strcmp(app.blocks(e.src).type, 'NonSpikingNeuron'), continue; end
            if ~strcmp(app.blocks(e.dst).type, 'NonSpikingNeuron'), continue; end
            ins(end+1) = struct('dst', app.blocks(e.dst).name); %#ok<AGROW>
        end
        S.inputs = ins;
        fid = fopen(path, 'w');
        fwrite(fid, jsonencode(S, 'PrettyPrint', true)); fclose(fid);
    end

    function cbBuild(app, ~, ~)
        app.nClicks = app.nClicks + 1;
        here = fileparts(mfilename('fullpath'));
        tmpJ = fullfile(here, 'circuit_last.json');
        app.saveJson(tmpJ);
        mdl = 'SNS_GuiCircuit';
        slx = fullfile(here, [mdl '.slx']);
        t0 = tic;
        sns_build_from_json(tmpJ, mdl, slx);
        app.timings.build = toc(t0);
        nNeural = app.countType('NonSpikingNeuron');
        nOther = numel(app.blocks) - nNeural;
        app.setStatus(sprintf(['Built %s (%.2f s; %d neurons, %d ' ...
            'synapses). Canvas sensory/muscle blocks outside the neural ' ...
            'layer are kept in the JSON, not carried into the model.'], ...
            slx, app.timings.build, nNeural, numel(app.edges)));
    end

    % -------------------------------------------------------- demo helper
    function demoKneeReflex(app)
        % Reduced knee-reflex circuit for demos and click accounting,
        % faithful to the demo architecture: afferent encoders deliver
        % currents to sensory interneurons, which synapse onto the MNs
        % (reciprocal and autogenic inhibition included).
        app.newCircuit();
        app.addBlock('NonSpikingNeuron', 46, 72);   % 1 MN_ext
        app.addBlock('NonSpikingNeuron', 46, 26);   % 2 MN_flex
        app.addBlock('NonSpikingNeuron', 26, 84);   % 3 SN_E (Ia ext)
        app.addBlock('NonSpikingNeuron', 26, 60);   % 4 SN_I (Ib ext)
        app.addBlock('NonSpikingNeuron', 26, 30);   % 5 SN_F (Ia flex)
        app.addBlock('IaMuscleSpindle', 10, 92);    % 6 Ia ext
        app.addBlock('IbGolgiTendon', 10, 66);      % 7 Ib ext
        app.addBlock('IaMuscleSpindle', 10, 22);    % 8 Ia flex
        app.addEdgeFromBlocks(6, 3);   % Ia ext current -> SN_E
        app.addEdgeFromBlocks(7, 4);   % Ib ext current -> SN_I
        app.addEdgeFromBlocks(8, 5);   % Ia flex current -> SN_F
        app.addEdgeFromBlocks(3, 1);   % SN_E -> MN_ext
        app.edges(end).params{2, 2} = 0;   % Esyn 0: excitatory
        app.addEdgeFromBlocks(4, 1);   % SN_I -> MN_ext (autogenic inh)
        app.addEdgeFromBlocks(5, 1);   % SN_F -> MN_ext (reciprocal inh)
        app.addEdgeFromBlocks(3, 2);   % SN_E -> MN_flex (reciprocal inh)
        app.selType = 'edge'; app.selIdx = 4;
        app.draw(); app.fillPropertyPanel();
        app.setStatus('Demo knee-reflex circuit (reduced). Select a connection.');
    end

    % ------------------------------------------------------------- utils
    function names = propChildren2Names(app)
        % param-name list in the order the edit fields were created
        if strcmp(app.selType, 'block'), P = app.blocks(app.selIdx).params;
        elseif strcmp(app.selType, 'edge'), P = app.edges(app.selIdx).params;
        else, names = {}; return; end
        names = P(:, 1);
    end

    function fields = propFieldsByOrder(app)
        fields = app.propChildren(arrayfun(@(h) isa(h, 'matlab.ui.control.NumericEditField'), ...
            app.propChildren));
    end

    function n = countType(app, type)
        n = sum(strcmp({app.blocks.type}, type));
    end

    function v = paramValue(~, P, nm, dflt)
        v = dflt;
        row = strcmp(P(:, 1), nm);
        if any(row), v = P{row, 2}; end
    end

    function hit = hitBlock(app, p)
        hit = 0; best = 4.5;
        for k = 1:numel(app.blocks)
            d = norm([app.blocks(k).x - p(1), app.blocks(k).y - p(2)]);
            if d < best, best = d; hit = k; end
        end
    end

    function setStatus(app, s)
        app.lblStatus.Text = s;
        drawnow;
    end
end
end

function A = asStructArray(v)
% jsondecode returns a CELL of structs when array elements carry
% different field sets (the winner-slice JSONs do); normalize to the
% common-field struct array the loader expects.
if isstruct(v), A = v; return; end
A = struct([]);
for k = 1:numel(v)
    e = v{k};
    if isempty(A)
        A = e;
    else
        fn = intersect(fieldnames(A), fieldnames(e));
        A = rmfield(A, setdiff(fieldnames(A), fn));
        e = rmfield(e, setdiff(fieldnames(e), fn));
        A(end+1) = e;
    end
end
end
