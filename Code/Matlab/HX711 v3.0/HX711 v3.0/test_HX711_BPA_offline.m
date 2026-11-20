function test_HX711_BPA_offline()
    %TEST_HX711_BPA_OFFLINE Hardware-free smoke test for HX711_BPA.
    %
    % Run on any machine after pulling the repo to confirm the app builds
    % and its non-hardware logic works (guards, known-factor entry,
    % conversions, save format, txt sidecar, run auto-increment, save-dir
    % resolution, Clean). No Arduino needed. Prints PASS/FAIL summary.
    %
    %   >> test_HX711_BPA_offline

    here = fileparts(mfilename('fullpath'));
    addpath(here);
    fprintf('Constructing HX711_BPA (hidden)...\n');
    app = HX711_BPA(here, true);  % keepHidden = true: no window flash
    res = app.offlineSelfTest(fullfile(tempdir, 'hx711_bpa_selftest'));
    delete(app);
    if res.allPass
        fprintf('HX711_BPA OFFLINE SELF TEST: PASS\n');
    else
        warning('HX711_BPA OFFLINE SELF TEST: FAIL - see per-check lines above.');
    end
end
