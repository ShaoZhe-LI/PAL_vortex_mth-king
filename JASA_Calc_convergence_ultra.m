%% ============================================================
% Ultrasound convergence analysis: N_FHT sweep only
%
% PURPOSE:
% 1) Keep only the ultrasound calculation.
% 2) Analyze multiple target points at the same time.
% 3) Sweep N_FHT:
%    - Compute the ultrasound field using the transform method.
%    - Snap each requested target point to the nearest transform-grid point.
%    - Compute the direct-integral reference value at the actual grid point.
%    - Compute the signed SPL difference:
%        DeltaSPL = SPL_transform - SPL_direct
% 4) Plot clean convergence curves in one figure.
% 5) Save figures, results, replotting data, and parameter txt files.
%
% DEPENDENCIES:
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%   - make_source_velocity.m
%% ============================================================

clear; clc; close all;

for i = 1:2
    close all;

    %% ===================== save and display options =====================
    save_figures      = true;
    save_calc_results = true;
    save_params_txt   = true;
    show_figures      = true;

    %% ===================== physical parameters =====================
    a     = 0.05;
    v0    = 0.172;
    c     = 343;
    rho0  = 1.21;
    pref  = 2e-5;

    fu = 40e3;
    fa = 0.5e3;
    f1 = fu;
    f2 = fu + fa;

    if i == 1
        m1 = 0;
    elseif i == 2
        m1 = 3;
    end

    m2 = m1;

    %% ===================== target points =====================
    if m1 == 0
        targets = [ ...
            0.0, 0.3;
            0.0, 1.0;
            0.3, 0.3];

        target_labels = { ...
            '(0, 0, 0.3) m', ...
            '(0, 0, 1.0) m', ...
            '(0.3, 0, 0.3) m'};
    end

    if m1 == 1
        targets = [ ...
            0.02, 0.3;
            0.06, 1.0;
            0.30, 0.3];

        target_labels = { ...
            '(0.02, 0, 0.3) m', ...
            '(0.06, 0, 1.0) m', ...
            '(0.3, 0, 0.3) m'};
    end

    if m1 == 3
        targets = [ ...
            0.05, 0.3;
            0.15, 1.0;
            0.30, 0.3];

        target_labels = { ...
            '(0.05, 0, 0.3) m', ...
            '(0.15, 0, 1.0) m', ...
            '(0.3, 0, 0.3) m'};
    end

    nTarget = size(targets, 1);

    %% ===================== sweep parameters =====================
    delta_fixed  = 0.001;
    N_FHT_fixed  = 16384;

    N_FHT_list = unique(round(logspace(log10(128), log10(32768), 42)));
    N_FHT_list = sort([N_FHT_list, N_FHT_fixed]);

    %% ===================== other calculation parameters =====================
    rho_max = 0.5;
    zu_max  = 3.0;

    green_R_min = 1e-12;

    %% ===================== parallel settings =====================
    use_parallel = true;
    num_workers  = 20;

    %% ===================== direct-reference calculation parameters =====================
    direct_ref_cfg = struct();
    direct_ref_cfg.src_block_size = 120000;
    direct_ref_cfg.N_FHT          = 65536;
    direct_ref_cfg.delta          = c / f2 / 24;
    direct_ref_cfg.dis_coe        = 64;
    direct_ref_cfg.num_workers    = num_workers;

    %% ===================== save path =====================
    time_tag = datestr(now, 'mmdd_HHMMSS');

    save_root_parent = 'result_convergence';
    save_root = fullfile(save_root_parent, sprintf('NFHTconv_%s_m%d', time_tag, m1));

    if (save_figures || save_calc_results || save_params_txt) && ~exist(save_root, 'dir')
        mkdir(save_root);
    end

    %% ===================== parallel pool initialization =====================
    if use_parallel
        p = gcp('nocreate');

        if isempty(p)
            parpool('local', num_workers);
        else
            if p.NumWorkers ~= num_workers
                delete(p);
                parpool('local', num_workers);
            end
        end
    end

    %% ===================== save parameters =====================
    if save_params_txt
        meta = struct();

        meta.a = a;
        meta.v0 = v0;
        meta.c = c;
        meta.rho0 = rho0;
        meta.pref = pref;

        meta.fu = fu;
        meta.fa = fa;
        meta.f1 = f1;
        meta.f2 = f2;

        meta.m1 = m1;
        meta.m2 = m2;

        meta.rho_max = rho_max;
        meta.zu_max = zu_max;
        meta.green_R_min = green_R_min;

        meta.targets = targets;
        meta.target_labels = target_labels;

        meta.N_FHT_list = N_FHT_list;
        meta.N_FHT_fixed = N_FHT_fixed;
        meta.delta_fixed = delta_fixed;

        meta.direct_ref_cfg = direct_ref_cfg;

        meta.use_parallel = use_parallel;
        meta.num_workers = num_workers;

        meta.save_root = save_root;
        meta.time_tag = time_tag;

        meta.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';

        local_write_all_params_txt(fullfile(save_root, 'all_parameters.txt'), meta);
    end

    %% ===================== build source and medium configuration =====================
    source_cfg = build_source_cfg(a, v0, m1, f1, fa, f2);
    medium_cfg = build_medium_cfg(c, rho0, pref);

    %% ===================== prepare direct-reference source discretization =====================
    fprintf('Preparing direct-reference source discretization...\n');

    calc_cfg_ref = build_calc_cfg_for_direct_reference( ...
        direct_ref_cfg.N_FHT, rho_max, zu_max, direct_ref_cfg.delta, ...
        direct_ref_cfg.dis_coe, direct_ref_cfg.num_workers);

    [source_prep_ref, ~, ~] = make_source_velocity(source_cfg, medium_cfg, calc_cfg_ref);

    Xs_ref = source_prep_ref.dim.Xs(:);
    Ys_ref = source_prep_ref.dim.Ys(:);
    Zs_ref = source_prep_ref.dim.Zs(:);

    vn1_ref = source_prep_ref.dim.Vn_pts_f1(:);
    vn2_ref = source_prep_ref.dim.Vn_pts_f2(:);

    if isfield(source_prep_ref.dim, 'dA_pts') && ~isempty(source_prep_ref.dim.dA_pts)
        dA_ref = source_prep_ref.dim.dA_pts(:);
    else
        dA_ref = source_prep_ref.dim.dA;
    end

    if isscalar(dA_ref)
        q1_ref = vn1_ref * dA_ref;
        q2_ref = vn2_ref * dA_ref;
    else
        q1_ref = vn1_ref .* dA_ref;
        q2_ref = vn2_ref .* dA_ref;
    end

    %% ============================================================
    % Sweep N_FHT
    %% ============================================================
    fprintf('\n============================================================\n');
    fprintf('N_FHT sweep only with fixed delta\n');
    fprintf('============================================================\n');

    res_NFHT = struct();

    res_NFHT.N_FHT_list = N_FHT_list(:);
    res_NFHT.delta_fixed = delta_fixed;
    res_NFHT.targets_req = targets;
    res_NFHT.target_labels = target_labels;
    res_NFHT.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';

    nCase = numel(N_FHT_list);

    res_NFHT.rho_actual = zeros(nCase, nTarget);
    res_NFHT.z_actual   = zeros(nCase, nTarget);

    res_NFHT.p1_transform = complex(zeros(nCase, nTarget));
    res_NFHT.p2_transform = complex(zeros(nCase, nTarget));
    res_NFHT.p1_direct    = complex(zeros(nCase, nTarget));
    res_NFHT.p2_direct    = complex(zeros(nCase, nTarget));

    res_NFHT.spl1_transform = zeros(nCase, nTarget);
    res_NFHT.spl2_transform = zeros(nCase, nTarget);
    res_NFHT.spl1_direct    = zeros(nCase, nTarget);
    res_NFHT.spl2_direct    = zeros(nCase, nTarget);

    res_NFHT.spl_diff1 = zeros(nCase, nTarget);
    res_NFHT.spl_diff2 = zeros(nCase, nTarget);

    res_NFHT.rel_err_amp1 = zeros(nCase, nTarget);
    res_NFHT.rel_err_amp2 = zeros(nCase, nTarget);

    for ii = 1:nCase
        N_FHT_now = N_FHT_list(ii);

        fprintf('\n[N_FHT sweep] case %d / %d : N_FHT = %d, delta = %.6e m\n', ...
            ii, nCase, N_FHT_now, delta_fixed);

        out = evaluate_one_transform_case_and_direct_reference_multi_targets( ...
            N_FHT_now, delta_fixed, ...
            targets, ...
            source_cfg, ...
            c, rho0, pref, f1, f2, ...
            rho_max, zu_max, ...
            green_R_min, ...
            Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref, ...
            source_prep_ref.k1, source_prep_ref.k2, ...
            source_prep_ref.w1, source_prep_ref.w2, ...
            source_prep_ref.medium.rho0, ...
            direct_ref_cfg.src_block_size);

        res_NFHT.rho_actual(ii, :) = out.rho_actual(:).';
        res_NFHT.z_actual(ii, :)   = out.z_actual(:).';

        res_NFHT.p1_transform(ii, :) = out.p1_transform(:).';
        res_NFHT.p2_transform(ii, :) = out.p2_transform(:).';
        res_NFHT.p1_direct(ii, :)    = out.p1_direct(:).';
        res_NFHT.p2_direct(ii, :)    = out.p2_direct(:).';

        res_NFHT.spl1_transform(ii, :) = out.spl1_transform(:).';
        res_NFHT.spl2_transform(ii, :) = out.spl2_transform(:).';
        res_NFHT.spl1_direct(ii, :)    = out.spl1_direct(:).';
        res_NFHT.spl2_direct(ii, :)    = out.spl2_direct(:).';

        res_NFHT.spl_diff1(ii, :) = out.spl_diff1(:).';
        res_NFHT.spl_diff2(ii, :) = out.spl_diff2(:).';

        res_NFHT.rel_err_amp1(ii, :) = out.rel_err_amp1(:).';
        res_NFHT.rel_err_amp2(ii, :) = out.rel_err_amp2(:).';
    end

    %% ============================================================
    % Prepare data for replotting
    %% ============================================================
    plot_data = struct();

    plot_data.N_FHT = res_NFHT.N_FHT_list;
    plot_data.delta_fixed = delta_fixed;
    plot_data.targets_req = targets;
    plot_data.target_labels = target_labels;
    plot_data.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';

    plot_data.spl_diff_f1 = res_NFHT.spl_diff1;
    plot_data.spl_diff_f2 = res_NFHT.spl_diff2;

    plot_data.spl_transform_f1 = res_NFHT.spl1_transform;
    plot_data.spl_transform_f2 = res_NFHT.spl2_transform;
    plot_data.spl_direct_f1    = res_NFHT.spl1_direct;
    plot_data.spl_direct_f2    = res_NFHT.spl2_direct;

    plot_data.rho_actual = res_NFHT.rho_actual;
    plot_data.z_actual   = res_NFHT.z_actual;

    plot_data.rel_err_amp1 = res_NFHT.rel_err_amp1;
    plot_data.rel_err_amp2 = res_NFHT.rel_err_amp2;

    %% ============================================================
    % Clean convergence plot
    %% ============================================================
    N_FHT = plot_data.N_FHT(:);

    if isfield(plot_data, 'spl_diff_f2')
        spl_diff = plot_data.spl_diff_f2;
    elseif isfield(plot_data, 'spl_diff_f1')
        spl_diff = plot_data.spl_diff_f1;
    else
        error('No SPL difference data was found in plot_data.');
    end

    if isfield(plot_data, 'target_labels')
        target_labels_plot = plot_data.target_labels;
    else
        target_labels_plot = arrayfun(@(k) sprintf('Target %d', k), ...
            1:size(spl_diff, 2), 'UniformOutput', false);
    end

    if isstring(target_labels_plot)
        target_labels_plot = cellstr(target_labels_plot);
    elseif ischar(target_labels_plot)
        target_labels_plot = {target_labels_plot};
    end

    target_labels_plot = cellfun(@(s) regexprep(s, '\s*m\)', ') m'), ...
        target_labels_plot, 'UniformOutput', false);

    nTarget_plot = size(spl_diff, 2);

    N_FHT_ref = 16384;

    style_spec = { ...
        {'Color',[1 0 0],          'LineStyle','-',  'LineWidth',2.1}, ...
        {'Color',[0 0.45 0.74],    'LineStyle','--', 'LineWidth',2.1}, ...
        {'Color',[0 0.6 0],        'LineStyle',':',  'LineWidth',2.3}, ...
        {'Color',[0.49 0.18 0.56], 'LineStyle','-.', 'LineWidth',2.1}, ...
        {'Color',[0.85 0.33 0.10], 'LineStyle','-',  'LineWidth',2.1}, ...
        {'Color',[0.30 0.30 0.30], 'LineStyle','--', 'LineWidth',2.1}};

    legend_fontsize = 19;
    legend_position = [0.45 0.72 0.37 0.20];

    fig1 = local_create_figure(show_figures, ...
        'Color', 'w', ...
        'Position', [100 100 680 540]);

    ax = axes('Parent', fig1);
    hold(ax, 'on');

    for it = 1:nTarget_plot
        style_idx = mod(it - 1, numel(style_spec)) + 1;
        semilogx(ax, N_FHT, spl_diff(:, it), style_spec{style_idx}{:});
    end

    yline(ax, 0, 'k:', 'LineWidth', 1.2);

    if N_FHT_ref >= min(N_FHT) && N_FHT_ref <= max(N_FHT)
        xline(ax, N_FHT_ref, 'k--', 'LineWidth', 1.4);
    end

    ax.Box = 'on';
    ax.Layer = 'top';
    ax.XScale = 'log';

    ax.FontName = 'Times New Roman';
    ax.FontSize = 22;
    ax.LineWidth = 1.2;

    ax.TickDir = 'in';
    ax.TickLength = [0.018 0.018];
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'off';

    ax.XColor = [0 0 0];
    ax.YColor = [0 0 0];

    grid(ax, 'off');

    x_min = min(N_FHT);
    x_max = max(N_FHT);

    xtick_candidates = [1e2 1e3 1e4];
    xticks_use = xtick_candidates(xtick_candidates >= x_min & xtick_candidates <= x_max);

    if isempty(xticks_use)
        xticks_use = unique([x_min; x_max]);
    elseif numel(xticks_use) == 1
        xticks_use = unique([x_min; xticks_use(:); x_max]);
    end

    ax.XTick = xticks_use(:).';

    if x_min == x_max
        xlim(ax, [0.9 * x_min, 1.1 * x_max]);
    else
        xlim(ax, [x_min, x_max]);
    end

    y_all = spl_diff(:);
    y_all = y_all(isfinite(y_all));

    if isempty(y_all)
        y_lim = [-5 5];
    else
        y_min = min(y_all);
        y_max = max(y_all);

        if y_min == y_max
            y_min = y_min - 1;
            y_max = y_max + 1;
        end

        pad = 0.08 * (y_max - y_min);
        y_lim = [y_min - pad, y_max + pad];

        y_lim(1) = 5 * floor(y_lim(1) / 5);
        y_lim(2) = 5 * ceil(y_lim(2) / 5);
    end

    ylim(ax, y_lim);

    xlabel(ax, '$N_{\mathrm{FHT}}$', ...
        'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 30);

    ylabel(ax, '$\Delta$ SPL (dB)', ...
        'Interpreter', 'latex', ...
        'FontName', 'Times New Roman', ...
        'FontSize', 30);

    lgd = legend(ax, target_labels_plot, 'Location', 'northeast');
    lgd.Box = 'off';
    lgd.FontName = 'Times New Roman';
    lgd.FontSize = legend_fontsize;

    lgd.Units = 'normalized';
    lgd.Position = legend_position;

    ax.Units = 'normalized';
    ax.Position = [0.17 0.23 0.77 0.70];

    %% ============================================================
    % Save figures and results
    %% ============================================================
    if save_figures
        png_path = fullfile(save_root, 'SPLdiff_vs_NFHT_logx_clean.png');
        fig_path = fullfile(save_root, 'SPLdiff_vs_NFHT_logx_clean.fig');
        pdf_path = fullfile(save_root, 'SPLdiff_vs_NFHT_logx_clean.pdf');

        saveas(fig1, png_path);
        savefig(fig1, fig_path);

        try
            exportgraphics(fig1, pdf_path, 'ContentType', 'vector');
        catch
            set(fig1, 'PaperPositionMode', 'auto');
            print(fig1, pdf_path, '-dpdf', '-painters');
        end
    end

    if save_calc_results
        save(fullfile(save_root, 'convergence_results.mat'), ...
            'res_NFHT', 'plot_data', 'targets', 'target_labels', ...
            'N_FHT_list', 'delta_fixed', '-v7.3');

        save(fullfile(save_root, 'plot_data_only.mat'), ...
            'plot_data', '-v7.3');
    end

    if ~show_figures
        close(fig1);
    end

    fprintf('\nAll done.\n');
    fprintf('Results folder: %s\n', save_root);

    if save_figures
        fprintf('Saved files:\n');
        fprintf('  %s\n', png_path);
        fprintf('  %s\n', fig_path);
        fprintf('  %s\n', pdf_path);
    end
end

%% ============================================================
% Create figure with optional display
%% ============================================================
function fig_handle = local_create_figure(show_figures, varargin)

if show_figures
    fig_handle = figure(varargin{:});
else
    fig_handle = figure(varargin{:}, 'Visible', 'off');
end

end

%% ============================================================
% One transform case plus direct reference at multiple targets
%% ============================================================
function out = evaluate_one_transform_case_and_direct_reference_multi_targets( ...
    N_FHT, delta, ...
    targets, ...
    source_cfg, ...
    c, rho0, pref, f1, f2, ...
    rho_max, zu_max, ...
    green_R_min, ...
    Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref, ...
    k1_ref, k2_ref, w1_ref, w2_ref, rho0_ref, ...
    src_block_size_ref)

w1 = 2 * pi * f1;
w2 = 2 * pi * f2;

k1 = w1 / c + 1j * AbsorpAttenCoef(f1);
k2 = w2 / c + 1j * AbsorpAttenCoef(f2);

Nh = 1.2 * rho_max;
NH = 1.2 * w2 / c;

n_FHT = 0:N_FHT - 1;
[a_solve, k0, x1, x0] = solve_kappa0(N_FHT, n_FHT);

xh = (x1 * Nh).';
z  = 0:delta:zu_max;
Nz = numel(z);

syms rho_v

a = source_cfg.a;
v0 = source_cfg.v0;
m1 = source_cfg.m;
m2 = source_cfg.m;

switch source_cfg.profile
    case 'Uniform'
        vs_sym1 = v0 * (heaviside(rho_v) - heaviside(rho_v - a));
        vs_sym2 = vs_sym1;

    case 'Focus'
        F = source_cfg.F;

        vs_sym1 = v0 * exp(-1j * real(k1) * sqrt(rho_v.^2 + F^2)) ...
            .* (heaviside(rho_v) - heaviside(rho_v - a));

        vs_sym2 = v0 * exp(-1j * real(k2) * sqrt(rho_v.^2 + F^2)) ...
            .* (heaviside(rho_v) - heaviside(rho_v - a));

    case 'Vortex-m'
        vs_sym1 = v0 * (heaviside(rho_v) - heaviside(rho_v - a));
        vs_sym2 = vs_sym1;

    otherwise
        error('Unknown source profile.');
end

Nh_v = 1.1 * a;
NH_v = NH;
xh_v = (x1 * Nh_v).';

vs_f1 = matlabFunction(vs_sym1);
vs_f2 = matlabFunction(vs_sym2);

vs1 = vs_f1(xh_v);
vs2 = vs_f2(xh_v);

Vs1 = m_FHT(vs1, N_FHT, 1, Nh_v, NH_v, a_solve, x0, x1, k0, m1);
Vs2 = m_FHT(vs2, N_FHT, 1, Nh_v, NH_v, a_solve, x0, x1, k0, m2);

G1_raw = build_green_space_g_transform_raw( ...
    xh, z, k1, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G2_raw = build_green_space_g_transform_raw( ...
    xh, z, k2, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G1 = (-4 * pi * 1j) * G1_raw;
G2 = (-4 * pi * 1j) * G2_raw;

p1_transform_full = compute_ultrasound_pressure_from_G( ...
    G1, Vs1, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m1, rho0, c, k1);

p2_transform_full = compute_ultrasound_pressure_from_G( ...
    G2, Vs2, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m2, rho0, c, k2);

nTarget = size(targets, 1);

rho_actual = zeros(nTarget, 1);
z_actual   = zeros(nTarget, 1);

p1_transform = complex(zeros(nTarget, 1));
p2_transform = complex(zeros(nTarget, 1));
p1_direct    = complex(zeros(nTarget, 1));
p2_direct    = complex(zeros(nTarget, 1));

spl1_transform = zeros(nTarget, 1);
spl2_transform = zeros(nTarget, 1);
spl1_direct    = zeros(nTarget, 1);
spl2_direct    = zeros(nTarget, 1);

spl_diff1 = zeros(nTarget, 1);
spl_diff2 = zeros(nTarget, 1);

rel_err_amp1 = zeros(nTarget, 1);
rel_err_amp2 = zeros(nTarget, 1);

for it = 1:nTarget
    rho_t = targets(it, 1);
    z_t   = targets(it, 2);

    [~, id_rho] = min(abs(xh - rho_t));
    [~, id_z]   = min(abs(z  - z_t));

    rho_actual(it) = xh(id_rho);
    z_actual(it)   = z(id_z);

    p1_transform(it) = p1_transform_full(id_rho, id_z);
    p2_transform(it) = p2_transform_full(id_rho, id_z);

    [p1_direct(it), p2_direct(it)] = calc_ultra_direct_single_point_blocked( ...
        Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref, ...
        k1_ref, k2_ref, w1_ref, w2_ref, rho0_ref, ...
        [rho_actual(it), 0, z_actual(it)], src_block_size_ref);

    spl1_transform(it) = local_pressure_to_spl(abs(p1_transform(it)), pref);
    spl2_transform(it) = local_pressure_to_spl(abs(p2_transform(it)), pref);

    spl1_direct(it) = local_pressure_to_spl(abs(p1_direct(it)), pref);
    spl2_direct(it) = local_pressure_to_spl(abs(p2_direct(it)), pref);

    spl_diff1(it) = spl1_transform(it) - spl1_direct(it);
    spl_diff2(it) = spl2_transform(it) - spl2_direct(it);

    rel_err_amp1(it) = abs(p1_transform(it) - p1_direct(it)) ...
        / max(abs(p1_direct(it)), 1e-16);

    rel_err_amp2(it) = abs(p2_transform(it) - p2_direct(it)) ...
        / max(abs(p2_direct(it)), 1e-16);

    fprintf('    target %d requested : rho = %.6f m, z = %.6f m\n', ...
        it, rho_t, z_t);

    fprintf('    target %d actual    : rho = %.6f m, z = %.6f m\n', ...
        it, rho_actual(it), z_actual(it));

    fprintf('    target %d f2 SPL(transform/direct) = %.6f / %.6f dB, DeltaSPL = %.6e dB\n', ...
        it, spl2_transform(it), spl2_direct(it), spl_diff2(it));
end

out = struct();

out.rho_actual = rho_actual;
out.z_actual   = z_actual;

out.p1_transform = p1_transform;
out.p2_transform = p2_transform;
out.p1_direct    = p1_direct;
out.p2_direct    = p2_direct;

out.spl1_transform = spl1_transform;
out.spl2_transform = spl2_transform;
out.spl1_direct    = spl1_direct;
out.spl2_direct    = spl2_direct;

out.spl_diff1 = spl_diff1;
out.spl_diff2 = spl_diff2;

out.rel_err_amp1 = rel_err_amp1;
out.rel_err_amp2 = rel_err_amp2;

end

%% ============================================================
% Build source configuration
%% ============================================================
function source_cfg = build_source_cfg(a, v0, m1, f1, fa, f2)

source_cfg = struct();

source_cfg.profile = 'Vortex-m';

source_cfg.a = a;
source_cfg.v0 = v0;
source_cfg.v_ratio = 1;
source_cfg.m = m1;
source_cfg.F = 0.2;

source_cfg.f1 = f1;
source_cfg.fa = fa;
source_cfg.f2 = f2;

source_cfg.internal = struct();

end

%% ============================================================
% Build medium configuration
%% ============================================================
function medium_cfg = build_medium_cfg(c, rho0, pref)

medium_cfg = struct();

medium_cfg.c0 = c;
medium_cfg.rho0 = rho0;
medium_cfg.pref = pref;

medium_cfg.use_absorp = true;
medium_cfg.atten_handle = @(f) AbsorpAttenCoef(f);

medium_cfg.internal = struct();

end

%% ============================================================
% Build calculation configuration for the direct reference
%% ============================================================
function calc_cfg = build_calc_cfg_for_direct_reference( ...
    N_FHT, rho_max, zu_max, delta, dis_coe, num_workers)

calc_cfg = struct();

calc_cfg.fht = struct();
calc_cfg.fht.N_FHT = N_FHT;
calc_cfg.fht.rho_max = rho_max;
calc_cfg.fht.Nh_scale = 1.2;
calc_cfg.fht.NH_scale = 1.2;
calc_cfg.fht.Nh_v_scale = 1.1;
calc_cfg.fht.zu_max = zu_max;
calc_cfg.fht.za_max = 0.1;
calc_cfg.fht.delta = delta;

calc_cfg.dim = struct();
calc_cfg.dim.method = 'rayleigh';
calc_cfg.dim.use_parallel = true;
calc_cfg.dim.num_workers = num_workers;
calc_cfg.dim.use_freq = 'f2';
calc_cfg.dim.dis_coe = dis_coe;
calc_cfg.dim.margin = 1;
calc_cfg.dim.src_discretization = 'polar';

calc_cfg.dim.grid_mode = 'uniform';
calc_cfg.dim.ds_rho = 32;
calc_cfg.dim.ds_rho_src = 16;
calc_cfg.dim.ds_rho_obs = 32;
calc_cfg.dim.uniform_dx = delta / 2;
calc_cfg.dim.uniform_dz = delta / 1;

calc_cfg.internal = struct();

end

%% ============================================================
% Direct ultrasound calculation at a single observation point
%% ============================================================
function [p1, p2] = calc_ultra_direct_single_point_blocked( ...
    Xs, Ys, Zs, q1, q2, k1, k2, w1, w2, rho0, obs_point, src_block_size)

x_obs = obs_point(1);
y_obs = obs_point(2);
z_obs = obs_point(3);

phi1 = 0;
phi2 = 0;

Nsrc = numel(Xs);

for ib = 1:src_block_size:Nsrc
    ie = min(ib + src_block_size - 1, Nsrc);
    idx = ib:ie;

    R = sqrt((x_obs - Xs(idx)).^2 + ...
             (y_obs - Ys(idx)).^2 + ...
             (z_obs - Zs(idx)).^2);

    R(R < 1e-12) = 1e-12;

    G1 = exp(1j * k1 .* R) ./ (4 * pi * R);
    G2 = exp(1j * k2 .* R) ./ (4 * pi * R);

    phi1 = phi1 + sum(G1 .* q1(idx));
    phi2 = phi2 + sum(G2 .* q2(idx));
end

p1 = -2j * rho0 * w1 .* phi1;
p2 = -2j * rho0 * w2 .* phi2;

end

%% ============================================================
% Compute ultrasound pressure from spectral-domain Green function
%% ============================================================
function p_out = compute_ultrasound_pressure_from_G( ...
    G_spec, Vs, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use, rho0, c0, k_use)

F = G_spec .* Vs;

phi = -1j * m_FHT(F, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use);

p_out = 1j * rho0 * c0 * real(k_use) .* phi;

end

%% ============================================================
% Zeroth-order Hankel transform of the spatial-domain Green function
%% ============================================================
function G_raw = build_green_space_g_transform_raw( ...
    rho_vec, z_vec, k_use, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0)

[RHO, Z] = ndgrid(rho_vec, z_vec);

RR = sqrt(RHO.^2 + Z.^2);
RR_use = max(RR, green_R_min);

g_space = exp(1j * k_use * RR_use) ./ (4 * pi * RR_use);

G_raw = m_FHT(g_space, N_FHT, numel(z_vec), Nh, NH, a_solve, x0, x1, k0, 0);

end

%% ============================================================
% Convert pressure amplitude to SPL
%% ============================================================
function spl = local_pressure_to_spl(p_amp, pref)

spl = 20 * log10(p_amp / max(pref, eps) / sqrt(2) + eps);

end

%% ============================================================
% Save all parameters to a txt file
%% ============================================================
function local_write_all_params_txt(txt_path, S)

fid = fopen(txt_path, 'w');

if fid < 0
    error('Cannot open txt file for writing: %s', txt_path);
end

fprintf(fid, '============================================================\n');
fprintf(fid, 'ALL PARAMETERS\n');
fprintf(fid, 'Generated time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '============================================================\n\n');

local_dump_any(fid, 'params', S, 0);

fclose(fid);

end

%% ============================================================
% Recursively write arbitrary variables
%% ============================================================
function local_dump_any(fid, name, val, indent)

sp = repmat(' ', 1, indent);

if isstruct(val)
    fprintf(fid, '%s%s = struct\n', sp, name);

    fn = fieldnames(val);

    for i = 1:numel(fn)
        local_dump_any(fid, sprintf('%s.%s', name, fn{i}), val.(fn{i}), indent + 2);
    end

    return;
end

if islogical(val) && isscalar(val)
    fprintf(fid, '%s%s = %s\n', sp, name, mat2str(val));
    return;
end

if isnumeric(val)
    if isscalar(val)
        if isreal(val)
            fprintf(fid, '%s%s = %.16g\n', sp, name, val);
        else
            fprintf(fid, '%s%s = %.16g %+.16gi\n', sp, name, real(val), imag(val));
        end
    else
        sz = size(val);
        fprintf(fid, '%s%s : numeric array, size = [%s]\n', sp, name, num2str(sz));

        vec = val(:);
        nshow = min(numel(vec), 20);

        if isreal(vec)
            fprintf(fid, '%s  first %d values = [', sp, nshow);
            fprintf(fid, ' %.8g', vec(1:nshow));
            fprintf(fid, ' ]\n');
        else
            fprintf(fid, '%s  first %d values = [', sp, nshow);

            for k = 1:nshow
                fprintf(fid, ' %.8g%+.8gi', real(vec(k)), imag(vec(k)));
            end

            fprintf(fid, ' ]\n');
        end
    end

    return;
end

if ischar(val) || (isstring(val) && isscalar(val))
    fprintf(fid, '%s%s = %s\n', sp, name, char(string(val)));
    return;
end

if isa(val, 'function_handle')
    fprintf(fid, '%s%s = %s\n', sp, name, func2str(val));
    return;
end

if iscell(val)
    sz = size(val);
    fprintf(fid, '%s%s : cell, size = [%s]\n', sp, name, num2str(sz));

    nshow = min(numel(val), 10);

    for i = 1:nshow
        local_dump_any(fid, sprintf('%s{%d}', name, i), val{i}, indent + 2);
    end

    return;
end

try
    fprintf(fid, '%s%s = %s\n', sp, name, evalc('disp(val)'));
catch
    fprintf(fid, '%s%s = <unprintable type: %s>\n', sp, name, class(val));
end

end