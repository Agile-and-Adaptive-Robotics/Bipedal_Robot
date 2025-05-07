function ExtPin10mm_compare()
    % Load BPA dataset
    load ExtPinBPASet.mat ke

    % === Build training data ===
    rel_all = []; theta_all = []; rest_all = []; ten_all = []; Lm_all = []; Fm_all = []; D_all = []; P_all = []; Fexp_all = [];

    for i = 1:numel(ke)
        k = ke(i);
        Lm = k.Lm_h;
        rel = computeRel(k, Lm);  % original rel from geometry
        rel_all = [rel_all; rel];
        theta_all = [theta_all; k.A_h];
        rest_all = [rest_all; repmat(k.rest, numel(rel), 1)];
        ten_all = [ten_all; repmat(k.ten, numel(rel), 1)];
        Lm_all = [Lm_all; Lm];
        Fm_all = [Fm_all; repmat(k.Fm, numel(rel), 1)];
        D_all  = [D_all; repmat(k.dBPA, numel(rel), 1)];
        P_all  = [P_all; repmat(k.P, numel(rel), 1)];
        Fexp_all = [Fexp_all; k.Mexp ./ k.mA_h];
    end

    % === Standardize predictors ===
    X_raw = [rest_all, theta_all, ten_all];
    mu = mean(X_raw);  sigma = std(X_raw);
    X_std = (X_raw - mu) ./ sigma;
    X_std_aug = [X_std, ones(size(X_std,1),1)];

    % === Regression model ===
    coeff_std = (X_std_aug' * X_std_aug) \ (X_std_aug' * (Fexp_all - rel_all));

    % === Unstandardized form (optional display only) ===
    coeff_unstd = coeff_std(1:3) ./ sigma';
    offset = coeff_std(4) - dot(coeff_std(1:3) ./ sigma', mu);
    fprintf(['\n✅ Correction model for rel:\nrel_corr = rel_raw ' ...
             '+ %.4f * L_rest + %.4f * theta + %.4f * tendon + %.4f\n'], ...
            coeff_unstd(1), coeff_unstd(2), coeff_unstd(3), offset);

    % === Apply correction and recompute torques ===
    for i = 1:numel(ke)
        k = ke(i);
        N = numel(k.Lm_h);
        rel_raw = computeRel(k, k.Lm_h);

        % Predict correction
        X = [repmat(k.rest,N,1), k.A_h, repmat(k.ten,N,1)];
        X_std = (X - mu) ./ sigma;
        rel_corr = rel_raw + X_std * coeff_std(1:3) + coeff_std(4);

        % New force prediction from rel_corr
        P_norm = k.P / 620;
        F_corr = f_festo(rel_corr, P_norm, k.dBPA) .* k.Fm;

        % New torque
        M_corr = F_corr .* k.mA_h;

        % Store
        ke(i).M_h_corr = M_corr;
        ke(i).F_h_corr = F_corr;
    end

    % === Plot comparison ===
    figure('Color','w'); tiledlayout(2,2);
    for i = 1:4
        nexttile; hold on;
        scatter(ke(i).A_h, ke(i).M_h, 25, [0.4 0.6 1], 'DisplayName', 'Original');
        scatter(ke(i).A_h, ke(i).M_h_corr, 25, [1 0.3 0.3], 'DisplayName', 'Corrected');
        scatter(ke(i).Aexp, ke(i).Mexp, 25, 'k', 'DisplayName', 'Experimental');
        xlabel('\theta_k (deg)'); ylabel('Torque (Nm)');
        title(sprintf('BPA #%d', i)); legend(); grid on;
    end
end

%% --- Utility: Compute relative contraction
function rel = computeRel(k, Lm)
    rest = k.rest;
    ten  = k.ten;
    fitn = k.fitn;
    KMAX = (rest - k.Kmax) / rest;
    contraction = (rest - (Lm - ten - 2*fitn)) / rest;
    rel = contraction / KMAX;
end

