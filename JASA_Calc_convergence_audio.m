%% ============================================================
% Audio convergence analysis: sweep N_FHT and delta (3 target points)
%
% REVISED VERSION:
% 1) Explicitly distinguish m1 and m2 in the transform method.
% 2) When m1 ~= m2, the direct ultrasound field does not use a mixed full cache.
% 3) The direct ultrasound field is assembled from pure-mode caches:
%      p1_direct_rz <- p1_direct_rz from pure cache (m1,m1)
%      p2_direct_rz <- p2_direct_rz from pure cache (m2,m2)
% 4) If a pure-mode cache does not exist, it is recomputed and cached.
% 5) make_source_velocity.m is not modified.
%
% DEPENDENCIES:
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%   - make_source_velocity.m
%% ============================================================

clear; clc; close all;

%% ===================== save and display options =====================
save_figures      = true;
save_calc_results = true;
save_params_txt   = true;
show_figures      = true;

for i = 1:2

close all;

%% ===================== basic physical parameters =====================
a     = 0.05;
v0    = 0.172;
c     = 343;
rho0  = 1.21;
beta  = 1.2;
pref  = 2e-5;

fu = 40e3;
fa = 0.5e3;
f1 = fu;
f2 = fu + fa;

if i == 1
    m1 = 0;
    m2 = 0;
end
if i == 2
    m1 = 3;
    m2 = 3;
end
if i == 3
    m1 = 0;
    m2 = 3;
end

ma = m2 - m1;

%% ===================== target points: 3 points =====================
if ma == 0
    targets = [ ...
        0.0, 0.3; ...
        0.0, 1.0; ...
        0.3, 0.3];

    target_labels = { ...
        '(0, 0, 0.3) m', ...
        '(0, 0, 1.0) m', ...
        '(0.3, 0, 0.3) m'};
end

if ma == 1
    targets = [ ...
        0.02, 0.3; ...
        0.06, 1.0; ...
        0.30, 0.30];

    target_labels = { ...
        '(0.02, 0, 0.3) m', ...
        '(0.06, 0, 1.0) m', ...
        '(0.3, 0, 0.3) m'};
end

if ma == 3
    targets = [ ...
        0.05, 0.3; ...
        0.15, 1.0; ...
        0.30, 0.30];

    target_labels = { ...
        '(0.05, 0, 0.3 m)', ...
        '(0.15, 0, 1.0 m)', ...
        '(0.3, 0, 0.3 m)'};
end

nTarget = size(targets,1);

%% ===================== sweep parameters =====================
delta_fixed  = 0.001;
N_FHT_fixed  = 16384;

N_FHT_list = unique(round(logspace(log10(512), log10(32768), 42)));
N_FHT_list = sort([N_FHT_list, N_FHT_fixed]);

delta_list = logspace(log10(0.0005), log10(0.5), 42);
delta_list = sort([delta_list, delta_fixed]);

%% ===================== other calculation parameters =====================
rho_max = 0.5;
zu_max  = 15.0;
za_max_default = 1.1 + delta_fixed;

green_R_min = 1e-12;

%% ===================== parallel settings =====================
use_parallel = true;
num_workers  = 20;

%% ===================== direct-reference calculation parameters =====================
direct_ref_cfg = struct();
direct_ref_cfg.src_block_size_ultra = 120000;
direct_ref_cfg.obs_block_size_rz    = 2048;
direct_ref_cfg.audio_q_block_size   = 200000;

direct_ref_cfg.N_FHT       = 65536;
direct_ref_cfg.delta       = c / f2 / 8;
direct_ref_cfg.dis_coe     = 32;
direct_ref_cfg.num_workers = num_workers;

%% ===================== audio direct-rz-phi grid parameters =====================
direct2d = struct();
direct2d.dr   = direct_ref_cfg.delta;
direct2d.dz   = direct_ref_cfg.delta;
direct2d.Nphi = 180;

audio_direct_use_z_mirror = true;

%% ===================== save path =====================
time_tag  = datestr(now, 'mmdd_HHMMSS');
save_root_parent = 'result_convergence';
save_root = fullfile(save_root_parent, sprintf('AudioConv_%s_m1_%d__m2_%d', time_tag, m1, m2));

if (save_figures || save_calc_results || save_params_txt) && ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% ===================== direct ultrasound cache path =====================
cache_root = 'data_pu';
if ~exist(cache_root, 'dir')
    mkdir(cache_root);
end

% The mixed folder is kept only for bookkeeping.
% The actual direct field uses pure-mode caches.
pu_cache_folder_name = build_pu_cache_folder_name( ...
    m1, m2, rho_max, zu_max, ...
    direct_ref_cfg.N_FHT, direct_ref_cfg.delta, direct_ref_cfg.dis_coe, ...
    direct2d.Nphi);

pu_cache_dir = fullfile(cache_root, pu_cache_folder_name);
pu_cache_mat = fullfile(pu_cache_dir, 'pu_direct_rz_cache.mat');
pu_cache_txt = fullfile(pu_cache_dir, 'pu_cache_parameters.txt');

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

%% ===================== build source / medium / calc cfg =====================
source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2);
medium_cfg = build_medium_cfg(c, rho0, beta, pref);

%% ===================== direct-rz grid =====================
rho_rz = build_uniform_grid_0_to_max(rho_max, direct2d.dr);
z_rz   = build_uniform_grid_0_to_max(zu_max,  direct2d.dz);

%% ============================================================
% Direct ultrasound field:
%   m1 == m2:
%       Use p1 and p2 from the pure cache (m1,m1).
%
%   m1 ~= m2:
%       Use pure-mode splitting:
%       p1_direct_rz <- p1_direct_rz from pure cache (m1,m1)
%       p2_direct_rz <- p2_direct_rz from pure cache (m2,m2)
%
%   This avoids the risk that a mixed cache internally uses only one source.m.
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Preparing direct ultrasound field on uniform (rho,z) grid ...\n');
fprintf('Target modal pair: m1 = %d, m2 = %d, ma = %d\n', m1, m2, ma);
fprintf('Pure-mode cache root: %s\n', cache_root);
fprintf('============================================================\n');

t_direct_all = tic;

[p1_pure_m1, p2_pure_m1, cache_info_m1] = get_or_compute_pure_mode_ultrasound_cache( ...
    m1, ...
    rho_rz, z_rz, ...
    cache_root, ...
    rho_max, zu_max, za_max_default, ...
    direct_ref_cfg, direct2d, ...
    a, v0, f1, fa, f2, ...
    medium_cfg, beta, rho0, c, ...
    use_parallel);

if m2 == m1
    p1_direct_rz = p1_pure_m1;
    p2_direct_rz = p2_pure_m1;

    cache_info_m2 = cache_info_m1;

else
    [~, p2_pure_m2, cache_info_m2] = get_or_compute_pure_mode_ultrasound_cache( ...
        m2, ...
        rho_rz, z_rz, ...
        cache_root, ...
        rho_max, zu_max, za_max_default, ...
        direct_ref_cfg, direct2d, ...
        a, v0, f1, fa, f2, ...
        medium_cfg, beta, rho0, c, ...
        use_parallel);

    p1_direct_rz = p1_pure_m1;
    p2_direct_rz = p2_pure_m2;
end

q_direct_rz = conj(p1_direct_rz) .* p2_direct_rz ...
    * beta*(2*pi*fa)/(1j*rho0^2*c^4);

time_direct_rz_once = toc(t_direct_all);

loaded_from_cache = cache_info_m1.loaded_from_cache && cache_info_m2.loaded_from_cache;
loaded_from_split_cache = (m1 ~= m2);

fprintf('Direct ultrasound prepared. elapsed = %.2f s\n', time_direct_rz_once);
fprintf('    p1 source cache: %s\n', cache_info_m1.cache_dir);
fprintf('    p2 source cache: %s\n', cache_info_m2.cache_dir);

%% ===================== save parameters =====================
if save_params_txt
    meta = struct();
    meta.a = a;
    meta.v0 = v0;
    meta.c = c;
    meta.rho0 = rho0;
    meta.beta = beta;
    meta.pref = pref;
    meta.fu = fu;
    meta.fa = fa;
    meta.f1 = f1;
    meta.f2 = f2;
    meta.m1 = m1;
    meta.m2 = m2;
    meta.ma = ma;
    meta.rho_max = rho_max;
    meta.zu_max = zu_max;
    meta.za_max_default = za_max_default;
    meta.green_R_min = green_R_min;
    meta.targets = targets;
    meta.target_labels = target_labels;
    meta.N_FHT_list = N_FHT_list;
    meta.delta_fixed = delta_fixed;
    meta.N_FHT_fixed = N_FHT_fixed;
    meta.delta_list = delta_list;
    meta.direct_ref_cfg = direct_ref_cfg;
    meta.direct2d = direct2d;
    meta.audio_direct_use_z_mirror = audio_direct_use_z_mirror;
    meta.use_parallel = use_parallel;
    meta.num_workers = num_workers;
    meta.save_root = save_root;
    meta.time_tag = time_tag;
    meta.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';

    meta.pu_cache_dir_bookkeeping = pu_cache_dir;
    meta.pu_cache_mat_bookkeeping = pu_cache_mat;
    meta.direct_cache_strategy = 'pure-mode split: p1 from (m1,m1), p2 from (m2,m2)';
    meta.p1_cache_info = cache_info_m1;
    meta.p2_cache_info = cache_info_m2;

    local_write_all_params_txt(fullfile(save_root, 'all_parameters.txt'), meta);
end

%% ===================== representative 2D case cache =====================
plot2d_case = struct();
plot2d_case.from_sweep = '';
plot2d_case.case_index = [];
plot2d_case.N_FHT = [];
plot2d_case.delta = [];
plot2d_case.rho_grid_audio = [];
plot2d_case.z_grid_audio = [];
plot2d_case.pa_transform_full = [];

%% ============================================================
% Sweep 1: N_FHT
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Sweep 1/2: N_FHT (delta fixed)\n');
fprintf('============================================================\n');

res_NFHT = init_audio_result_struct(numel(N_FHT_list), nTarget);
res_NFHT.N_FHT_list     = N_FHT_list(:);
res_NFHT.delta_fixed    = delta_fixed;
res_NFHT.targets_req    = targets;
res_NFHT.target_labels  = target_labels;
res_NFHT.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';
res_NFHT.case_time_sec  = zeros(numel(N_FHT_list),1);

for ii = 1:numel(N_FHT_list)
    t_case = tic;
    N_FHT_now = N_FHT_list(ii);

    fprintf('\n[N_FHT sweep] case %d / %d : N_FHT = %d, delta = %.6e m\n', ...
        ii, numel(N_FHT_list), N_FHT_now, delta_fixed);

    out = evaluate_one_audio_transform_case_with_fixed_direct_q( ...
        N_FHT_now, delta_fixed, 2.0 + delta_fixed, ...
        targets, ...
        source_cfg, ...
        c, rho0, beta, pref, f1, f2, fa, ...
        rho_max, zu_max, ...
        green_R_min, ma, ...
        rho_rz, z_rz, q_direct_rz, direct2d, ...
        direct_ref_cfg.audio_q_block_size, ...
        use_parallel, audio_direct_use_z_mirror);

    res_NFHT = assign_audio_result_struct(res_NFHT, ii, out);
    res_NFHT.case_time_sec(ii) = toc(t_case);

    if ii == numel(N_FHT_list)
        plot2d_case.from_sweep = 'NFHT';
        plot2d_case.case_index = ii;
        plot2d_case.N_FHT = N_FHT_now;
        plot2d_case.delta = delta_fixed;
        plot2d_case.rho_grid_audio = out.rho_grid_audio;
        plot2d_case.z_grid_audio   = out.z_grid_audio;
        plot2d_case.pa_transform_full = out.pa_transform_full;
    end

    fprintf('[N_FHT sweep] case done. elapsed = %.2f s\n', res_NFHT.case_time_sec(ii));
end

%% ============================================================
% Sweep 2: delta
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Sweep 2/2: delta (N_FHT fixed)\n');
fprintf('============================================================\n');

res_delta = init_audio_result_struct(numel(delta_list), nTarget);
res_delta.delta_list     = delta_list(:);
res_delta.N_FHT_fixed    = N_FHT_fixed;
res_delta.targets_req    = targets;
res_delta.target_labels  = target_labels;
res_delta.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';
res_delta.case_time_sec  = zeros(numel(delta_list),1);

for ii = 1:numel(delta_list)
    t_case = tic;
    delta_now = delta_list(ii);
    za_max_now = 2.0 + delta_now;

    fprintf('\n[delta sweep] case %d / %d : N_FHT = %d, delta = %.6e m\n', ...
        ii, numel(delta_list), N_FHT_fixed, delta_now);

    out = evaluate_one_audio_transform_case_with_fixed_direct_q( ...
        N_FHT_fixed, delta_now, za_max_now, ...
        targets, ...
        source_cfg, ...
        c, rho0, beta, pref, f1, f2, fa, ...
        rho_max, zu_max, ...
        green_R_min, ma, ...
        rho_rz, z_rz, q_direct_rz, direct2d, ...
        direct_ref_cfg.audio_q_block_size, ...
        use_parallel, audio_direct_use_z_mirror);

    res_delta = assign_audio_result_struct(res_delta, ii, out);
    res_delta.case_time_sec(ii) = toc(t_case);

    need_update_plot2d = false;
    if isempty(plot2d_case.delta)
        need_update_plot2d = true;
    elseif delta_now < plot2d_case.delta - 1e-15
        need_update_plot2d = true;
    elseif abs(delta_now - plot2d_case.delta) <= 1e-15 && N_FHT_fixed > plot2d_case.N_FHT
        need_update_plot2d = true;
    end

    if need_update_plot2d
        plot2d_case.from_sweep = 'delta';
        plot2d_case.case_index = ii;
        plot2d_case.N_FHT = N_FHT_fixed;
        plot2d_case.delta = delta_now;
        plot2d_case.rho_grid_audio = out.rho_grid_audio;
        plot2d_case.z_grid_audio   = out.z_grid_audio;
        plot2d_case.pa_transform_full = out.pa_transform_full;
    end

    fprintf('[delta sweep] case done. elapsed = %.2f s\n', res_delta.case_time_sec(ii));
end

%% ============================================================
% Prepare data for replotting
%% ============================================================
plot_data = struct();
plot_data.N_FHT = N_FHT_list(:);
plot_data.delta_fixed = delta_fixed;
plot_data.delta_list = delta_list(:);
plot_data.N_FHT_fixed = N_FHT_fixed;
plot_data.targets_req = targets;
plot_data.target_labels = target_labels;
plot_data.deltaSPL_definition = 'DeltaSPL = SPL_transform - SPL_direct';
plot_data.time_direct_rz_once = time_direct_rz_once;
plot_data.loaded_from_cache = loaded_from_cache;
plot_data.loaded_from_split_cache = loaded_from_split_cache;
plot_data.pu_cache_dir = pu_cache_dir;
plot_data.p1_cache_info = cache_info_m1;
plot_data.p2_cache_info = cache_info_m2;

plot_data.res_NFHT = res_NFHT;
plot_data.res_delta = res_delta;

%% ============================================================
% Prepare 2D plotting data, downsampled to several hundred points
%% ============================================================
plot2d_view = struct();
plot2d_view.meta = plot2d_case;

if ~isempty(plot2d_case.pa_transform_full)
    rho_full = plot2d_case.rho_grid_audio(:);
    z_full   = plot2d_case.z_grid_audio(:);
    pa_full  = plot2d_case.pa_transform_full;

    n_rho_keep = min(400, numel(rho_full));
    n_z_keep   = min(400, numel(z_full));

    idx_rho_ds = unique(round(linspace(1, numel(rho_full), n_rho_keep)));
    idx_z_ds   = unique(round(linspace(1, numel(z_full),   n_z_keep)));

    rho_ds = rho_full(idx_rho_ds);
    z_ds   = z_full(idx_z_ds);
    pa_ds  = pa_full(idx_rho_ds, idx_z_ds);

    plot2d_view.rho_grid = rho_ds;
    plot2d_view.z_grid   = z_ds;
    plot2d_view.pa_2d    = pa_ds;

    p_abs_ds = abs(pa_ds);
    p_ref_ds = max(p_abs_ds(:));
    plot2d_view.spl_2d_norm = 20*log10(p_abs_ds / max(p_ref_ds, eps) + eps);

    [~, idx_rho0_ds] = min(abs(rho_ds));
    plot2d_view.rho_axis_value = rho_ds(idx_rho0_ds);
    plot2d_view.pa_axis = pa_ds(idx_rho0_ds, :);
    plot2d_view.spl_axis_norm = 20*log10(abs(plot2d_view.pa_axis) / max(abs(plot2d_view.pa_axis), eps) + eps);
else
    plot2d_view.rho_grid = [];
    plot2d_view.z_grid   = [];
    plot2d_view.pa_2d    = [];
    plot2d_view.spl_2d_norm = [];
    plot2d_view.rho_axis_value = [];
    plot2d_view.pa_axis = [];
    plot2d_view.spl_axis_norm = [];
end

%% ============================================================
% Clean plotting
%% ============================================================
if isstring(target_labels)
    target_labels_plot = cellstr(target_labels);
elseif ischar(target_labels)
    target_labels_plot = {target_labels};
else
    target_labels_plot = target_labels;
end

target_labels_plot = cellfun(@(s) regexprep(s, '\s*m\)', ') m'), ...
    target_labels_plot, 'UniformOutput', false);

N_FHT_ref = N_FHT_fixed;
delta_ref = delta_fixed;

style_spec = { ...
    {'Color',[1 0 0],          'LineStyle','-',  'LineWidth',2.1}, ...
    {'Color',[0 0.45 0.74],    'LineStyle','--', 'LineWidth',2.1}, ...
    {'Color',[0 0.6 0],        'LineStyle',':',  'LineWidth',2.3}, ...
    {'Color',[0.49 0.18 0.56], 'LineStyle','-.', 'LineWidth',2.1}, ...
    {'Color',[0.85 0.33 0.10], 'LineStyle','-',  'LineWidth',2.1}, ...
    {'Color',[0.30 0.30 0.30], 'LineStyle','--', 'LineWidth',2.1}};

legend_fontsize = 19;

legend_pos_fig1 = [0.45 0.72 0.37 0.20];
legend_pos_fig2 = [0.25 0.72 0.37 0.20];

%% ============================================================
% Figure 1: DeltaSPL vs N_FHT
%% ============================================================
fig1 = local_create_figure(show_figures, ...
    'Color', 'w', ...
    'Position', [100 100 680 540]);

ax1 = axes('Parent', fig1);
hold(ax1, 'on');

nTarget1 = size(res_NFHT.spl_diff, 2);

for it = 1:nTarget1
    style_idx = mod(it - 1, numel(style_spec)) + 1;
    semilogx(ax1, res_NFHT.N_FHT_list, res_NFHT.spl_diff(:, it), ...
        style_spec{style_idx}{:});
end

yline(ax1, 0, 'k:', 'LineWidth', 1.2);

if N_FHT_ref >= min(res_NFHT.N_FHT_list) && N_FHT_ref <= max(res_NFHT.N_FHT_list)
    xline(ax1, N_FHT_ref, 'k--', 'LineWidth', 1.4);
end

local_set_axis_style_box(ax1);
local_set_log_xticks_nfht_sparse(ax1, res_NFHT.N_FHT_list);
local_set_y_margin(ax1, res_NFHT.spl_diff(:));

xlabel(ax1, '$N_{\mathrm{FHT}}$', ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 30);

ylabel(ax1, '$\Delta$ SPL (dB)', ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 30);

lgd1 = legend(ax1, target_labels_plot, 'Location', 'northeast');
lgd1.Box = 'off';
lgd1.FontName = 'Times New Roman';
lgd1.FontSize = legend_fontsize;
lgd1.Units = 'normalized';
lgd1.Position = legend_pos_fig1;

ax1.Units = 'normalized';
ax1.Position = [0.17 0.23 0.77 0.70];

%% ============================================================
% Figure 2: DeltaSPL vs delta_z
%% ============================================================
fig2 = local_create_figure(show_figures, ...
    'Color', 'w', ...
    'Position', [140 120 680 540]);

ax2 = axes('Parent', fig2);
hold(ax2, 'on');

nTarget2 = size(res_delta.spl_diff, 2);

for it = 1:nTarget2
    style_idx = mod(it - 1, numel(style_spec)) + 1;
    semilogx(ax2, res_delta.delta_list, res_delta.spl_diff(:, it), ...
        style_spec{style_idx}{:});
end

yline(ax2, 0, 'k:', 'LineWidth', 1.2);

if delta_ref >= min(res_delta.delta_list) && delta_ref <= max(res_delta.delta_list)
    xline(ax2, delta_ref, 'k--', 'LineWidth', 1.4);
end

local_set_axis_style_box(ax2);
local_set_log_xticks_delta(ax2, res_delta.delta_list);
local_set_y_margin(ax2, res_delta.spl_diff(:));

xlabel(ax2, '$\Delta_z$ (m)', ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 30);

ylabel(ax2, '$\Delta$ SPL (dB)', ...
    'Interpreter', 'latex', ...
    'FontName', 'Times New Roman', ...
    'FontSize', 30);

lgd2 = legend(ax2, target_labels_plot, 'Location', 'northwest');
lgd2.Box = 'off';
lgd2.FontName = 'Times New Roman';
lgd2.FontSize = legend_fontsize;
lgd2.Units = 'normalized';
lgd2.Position = legend_pos_fig2;

ax2.Units = 'normalized';
ax2.Position = [0.17 0.23 0.77 0.70];

%% ============================================================
% Save clean figures
%% ============================================================
if save_figures
    png1_path = fullfile(save_root, 'Audio_DeltaSPL_vs_NFHT_replot.png');
    fig1_path = fullfile(save_root, 'Audio_DeltaSPL_vs_NFHT_replot.fig');
    pdf1_path = fullfile(save_root, 'Audio_DeltaSPL_vs_NFHT_replot.pdf');

    png2_path = fullfile(save_root, 'Audio_DeltaSPL_vs_deltaz_replot.png');
    fig2_path = fullfile(save_root, 'Audio_DeltaSPL_vs_deltaz_replot.fig');
    pdf2_path = fullfile(save_root, 'Audio_DeltaSPL_vs_deltaz_replot.pdf');

    local_save_figure_png_fig_pdf(fig1, png1_path, fig1_path, pdf1_path);
    local_save_figure_png_fig_pdf(fig2, png2_path, fig2_path, pdf2_path);

    fprintf('Saved figure files:\n');
    fprintf('  %s\n', png1_path);
    fprintf('  %s\n', fig1_path);
    fprintf('  %s\n', pdf1_path);
    fprintf('  %s\n', png2_path);
    fprintf('  %s\n', fig2_path);
    fprintf('  %s\n', pdf2_path);
end

if ~show_figures
    close(fig1);
    close(fig2);
end

%% ============================================================
% Save results
%% ============================================================
if save_calc_results
    save(fullfile(save_root, 'audio_convergence_results.mat'), ...
        'res_NFHT', 'res_delta', 'plot_data', 'plot2d_case', 'plot2d_view', ...
        'targets', 'target_labels', ...
        'N_FHT_list', 'delta_fixed', 'delta_list', 'N_FHT_fixed', ...
        'time_direct_rz_once', 'loaded_from_cache', 'loaded_from_split_cache', ...
        'pu_cache_dir', 'cache_info_m1', 'cache_info_m2', ...
        '-v7.3');

    save(fullfile(save_root, 'audio_plot_data_only.mat'), ...
        'plot_data', 'plot2d_case', 'plot2d_view', '-v7.3');
end

fprintf('\nAll done.\n');
fprintf('Direct ultrasound precompute/load time = %.2f s\n', time_direct_rz_once);
fprintf('p1 cache dir = %s\n', cache_info_m1.cache_dir);
fprintf('p2 cache dir = %s\n', cache_info_m2.cache_dir);
fprintf('Results folder: %s\n', save_root);

end

%% ============================================================
% Initialize result structure
%% ============================================================
function S = init_audio_result_struct(nCase, nTarget)
S = struct();
S.rho_actual     = zeros(nCase, nTarget);
S.z_actual       = zeros(nCase, nTarget);

S.pa_transform   = complex(zeros(nCase, nTarget));
S.pa_direct      = complex(zeros(nCase, nTarget));

S.spl_transform  = zeros(nCase, nTarget);
S.spl_direct     = zeros(nCase, nTarget);
S.spl_diff       = zeros(nCase, nTarget);

S.rel_err_amp    = zeros(nCase, nTarget);
end

%% ============================================================
% Assign one case result
%% ============================================================
function S = assign_audio_result_struct(S, idx, out)
S.rho_actual(idx, :)    = out.rho_actual(:).';
S.z_actual(idx, :)      = out.z_actual(:).';

S.pa_transform(idx, :)  = out.pa_transform(:).';
S.pa_direct(idx, :)     = out.pa_direct(:).';

S.spl_transform(idx, :) = out.spl_transform(:).';
S.spl_direct(idx, :)    = out.spl_direct(:).';
S.spl_diff(idx, :)      = out.spl_diff(:).';

S.rel_err_amp(idx, :)   = out.rel_err_amp(:).';
end

%% ============================================================
% One transform case plus fixed-direct-q reference for multiple target points
%% ============================================================
function out = evaluate_one_audio_transform_case_with_fixed_direct_q( ...
    N_FHT, delta, za_max, ...
    targets, ...
    source_cfg, ...
    c, rho0, beta, pref, f1, f2, fa, ...
    rho_max, zu_max, ...
    green_R_min, ma, ...
    rho_rz, z_rz, q_direct_rz, direct2d, ...
    audio_q_block_size, ...
    use_parallel, use_z_mirror)

w1 = 2*pi*f1;
w2 = 2*pi*f2;
wa = 2*pi*fa;

k1 = w1/c + 1j*AbsorpAttenCoef(f1);
k2 = w2/c + 1j*AbsorpAttenCoef(f2);
ka = wa/c + 1j*AbsorpAttenCoef(fa);

Nh = 2 * rho_max;
NH = 4 * w2 / c;

n_FHT = 0:N_FHT-1;
[a_solve, k0, x1, x0] = solve_kappa0(N_FHT, n_FHT);

xh = (x1 * Nh).';
z  = 0:delta:zu_max;
z_audio = 0:delta:za_max;

Nz  = numel(z);

%% ---- source velocity ----
syms rho_v
a = source_cfg.a;
v0 = source_cfg.v0;

% Key fix: do not set m1 = source_cfg.m and m2 = source_cfg.m here.
m1 = source_cfg.m1;
m2 = source_cfg.m2;

switch source_cfg.profile
    case 'Uniform'
        vs_sym1 = v0 * (heaviside(rho_v) - heaviside(rho_v-a));
        vs_sym2 = vs_sym1;

    case 'Focus'
        F = source_cfg.F;
        vs_sym1 = v0 * exp(-1j*real(k1)*sqrt(rho_v.^2+F^2)) ...
            .* (heaviside(rho_v) - heaviside(rho_v-a));
        vs_sym2 = v0 * exp(-1j*real(k2)*sqrt(rho_v.^2+F^2)) ...
            .* (heaviside(rho_v) - heaviside(rho_v-a));

    case 'Vortex-m'
        vs_sym1 = v0 * (heaviside(rho_v) - heaviside(rho_v-a));
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

%% ---- ultrasound transform ----
G1_raw = build_green_space_g_transform_raw( ...
    xh, z, k1, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G2_raw = build_green_space_g_transform_raw( ...
    xh, z, k2, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G1 = (-4*pi*1j) * G1_raw;
G2 = (-4*pi*1j) * G2_raw;

p1_transform = compute_ultrasound_pressure_from_G( ...
    G1, Vs1, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m1, rho0, c, k1);

p2_transform = compute_ultrasound_pressure_from_G( ...
    G2, Vs2, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m2, rho0, c, k2);

q_transform = conj(p1_transform) .* p2_transform * beta*wa/(1j*rho0^2*c^4);

mask_rho = double(xh(:) <= rho_max);
q_transform = q_transform .* repmat(mask_rho, 1, size(q_transform,2));

%% ---- audio transform ----
Ga_raw = build_green_space_g_transform_raw( ...
    xh, z, ka, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);
Ga = (-4*pi*1j) * Ga_raw;

[pa_transform_full, ~] = compute_paW_from_Gr00( ...
    q_transform, z, z_audio, Ga, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa);

%% ---- target points ----
nTarget = size(targets,1);

rho_actual    = zeros(nTarget,1);
z_actual      = zeros(nTarget,1);
pa_transform  = complex(zeros(nTarget,1));
pa_direct     = complex(zeros(nTarget,1));
spl_transform = zeros(nTarget,1);
spl_direct    = zeros(nTarget,1);
spl_diff      = zeros(nTarget,1);
rel_err_amp   = zeros(nTarget,1);

rho_eval = zeros(nTarget,1);
z_eval   = zeros(nTarget,1);

for it = 1:nTarget
    rho_t = targets(it,1);
    z_t   = targets(it,2);

    [~, id_rho] = min(abs(xh - rho_t));
    [~, id_z]   = min(abs(z_audio - z_t));

    rho_actual(it) = xh(id_rho);
    z_actual(it)   = z_audio(id_z);

    rho_eval(it) = rho_actual(it);
    z_eval(it)   = z_actual(it);

    pa_transform(it) = pa_transform_full(id_rho, id_z);
end

%% ---- audio direct at target points using fixed q_direct_rz ----
pa_direct = compute_audio_direct_points_from_qrzphi_parallel( ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d.dr, direct2d.dz, direct2d.Nphi, ma, ...
    ka, rho0, wa, ...
    rho_eval, z_eval, ...
    audio_q_block_size, ...
    use_parallel, use_z_mirror);

%% ---- SPL and error ----
for it = 1:nTarget
    spl_transform(it) = local_pressure_to_spl(abs(pa_transform(it)), pref);
    spl_direct(it)    = local_pressure_to_spl(abs(pa_direct(it)), pref);
    spl_diff(it)      = spl_transform(it) - spl_direct(it);

    rel_err_amp(it)   = abs(pa_transform(it) - pa_direct(it)) / max(abs(pa_direct(it)), 1e-16);

    fprintf('    target %d requested : rho = %.6f m, z = %.6f m\n', it, targets(it,1), targets(it,2));
    fprintf('    target %d actual    : rho = %.6f m, z = %.6f m\n', it, rho_actual(it), z_actual(it));
    fprintf('    target %d audio SPL(transform/direct) = %.6f / %.6f dB, DeltaSPL = %.6e dB\n', ...
        it, spl_transform(it), spl_direct(it), spl_diff(it));
end

out = struct();
out.rho_actual    = rho_actual;
out.z_actual      = z_actual;
out.pa_transform  = pa_transform;
out.pa_direct     = pa_direct;
out.spl_transform = spl_transform;
out.spl_direct    = spl_direct;
out.spl_diff      = spl_diff;
out.rel_err_amp   = rel_err_amp;

out.rho_grid_audio = xh(:);
out.z_grid_audio   = z_audio(:);
out.pa_transform_full = pa_transform_full;
end

%% ============================================================
% Get or compute pure-mode direct ultrasound cache
%% ============================================================
function [p1_direct_rz, p2_direct_rz, cache_info] = get_or_compute_pure_mode_ultrasound_cache( ...
    m_use, ...
    rho_rz, z_rz, ...
    cache_root, ...
    rho_max, zu_max, za_max, ...
    direct_ref_cfg, direct2d, ...
    a, v0, f1, fa, f2, ...
    medium_cfg, beta, rho0, c, ...
    use_parallel)

cache_info = struct();
cache_info.m = m_use;
cache_info.loaded_from_cache = false;
cache_info.cache_dir = '';
cache_info.cache_mat = '';
cache_info.elapsed_sec = NaN;

pu_cache_folder_name = build_pu_cache_folder_name( ...
    m_use, m_use, rho_max, zu_max, ...
    direct_ref_cfg.N_FHT, direct_ref_cfg.delta, direct_ref_cfg.dis_coe, ...
    direct2d.Nphi);

pu_cache_dir = fullfile(cache_root, pu_cache_folder_name);
pu_cache_mat = fullfile(pu_cache_dir, 'pu_direct_rz_cache.mat');
pu_cache_txt = fullfile(pu_cache_dir, 'pu_cache_parameters.txt');

cache_info.cache_dir = pu_cache_dir;
cache_info.cache_mat = pu_cache_mat;

fprintf('\n------------------------------------------------------------\n');
fprintf('Preparing pure-mode ultrasound cache: m = %d\n', m_use);
fprintf('Cache folder: %s\n', pu_cache_dir);
fprintf('------------------------------------------------------------\n');

t0 = tic;

if exist(pu_cache_mat, 'file')
    fprintf('Found pure-mode cache. Loading...\n');

    S_cache = load(pu_cache_mat, ...
        'p1_direct_rz', 'p2_direct_rz', 'rho_rz', 'z_rz');

    if ~isequal(rho_rz(:), S_cache.rho_rz(:)) || ~isequal(z_rz(:), S_cache.z_rz(:))
        error(['Pure-mode cached rho/z grid does not match current setup. ' ...
               'Please delete cache folder: %s'], pu_cache_dir);
    end

    p1_direct_rz = S_cache.p1_direct_rz;
    p2_direct_rz = S_cache.p2_direct_rz;

    cache_info.loaded_from_cache = true;
    cache_info.elapsed_sec = toc(t0);

    fprintf('Pure-mode cache loaded. m = %d, elapsed = %.2f s\n', ...
        m_use, cache_info.elapsed_sec);
    return;
end

fprintf('Pure-mode cache not found. Computing m = %d direct ultrasound field...\n', m_use);

if ~exist(pu_cache_dir, 'dir')
    mkdir(pu_cache_dir);
end

source_cfg_mode = build_source_cfg(a, v0, m_use, m_use, f1, fa, f2);

calc_cfg_ref = build_calc_cfg_for_direct_reference( ...
    direct_ref_cfg.N_FHT, rho_max, zu_max, za_max, direct_ref_cfg.delta, ...
    direct_ref_cfg.dis_coe, direct_ref_cfg.num_workers);

[source_prep_ref, ~, ~] = make_source_velocity(source_cfg_mode, medium_cfg, calc_cfg_ref);

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

[p1_direct_rz, p2_direct_rz] = compute_ultrasound_direct_rz_grid_parallel_from_ref_sources( ...
    rho_rz, z_rz, ...
    Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref, ...
    source_prep_ref.k1, source_prep_ref.k2, ...
    source_prep_ref.w1, source_prep_ref.w2, ...
    source_prep_ref.medium.rho0, ...
    direct_ref_cfg.src_block_size_ultra, ...
    direct_ref_cfg.obs_block_size_rz, ...
    use_parallel);

q_direct_rz = conj(p1_direct_rz) .* p2_direct_rz ...
    * beta*(2*pi*fa)/(1j*rho0^2*c^4);

cache_meta = struct();
cache_meta.cache_type = 'pure-mode';
cache_meta.m1 = m_use;
cache_meta.m2 = m_use;
cache_meta.ma = 0;
cache_meta.rho_max = rho_max;
cache_meta.zu_max = zu_max;
cache_meta.N_FHT = direct_ref_cfg.N_FHT;
cache_meta.delta = direct_ref_cfg.delta;
cache_meta.dis_coe = direct_ref_cfg.dis_coe;
cache_meta.Nphi = direct2d.Nphi;
cache_meta.a = a;
cache_meta.v0 = v0;
cache_meta.f1 = f1;
cache_meta.f2 = f2;
cache_meta.fa = fa;
cache_meta.c = c;
cache_meta.rho0 = rho0;
cache_meta.beta = beta;
cache_meta.time_tag = datestr(now, 'yyyy-mm-dd HH:MM:SS');

save(pu_cache_mat, ...
    'p1_direct_rz', 'p2_direct_rz', 'q_direct_rz', ...
    'rho_rz', 'z_rz', 'cache_meta', '-v7.3');

local_write_all_params_txt(pu_cache_txt, cache_meta);

cache_info.loaded_from_cache = false;
cache_info.elapsed_sec = toc(t0);

fprintf('Pure-mode direct ultrasound cache saved: %s\n', pu_cache_mat);
fprintf('Pure-mode computation finished. m = %d, elapsed = %.2f s\n', ...
    m_use, cache_info.elapsed_sec);
end

%% ============================================================
% Direct ultrasound calculation on the (rho,z) grid using reference source discretization
%% ============================================================
function [p1_rz, p2_rz] = compute_ultrasound_direct_rz_grid_parallel_from_ref_sources( ...
    rho_grid, z_grid, ...
    Xs, Ys, Zs, q1, q2, ...
    k1, k2, w1, w2, rho0, ...
    src_block_size, obs_block_size, use_parallel)

[RHO, Z] = ndgrid(rho_grid, z_grid);
xo = RHO(:);
yo = zeros(size(xo));
zo = Z(:);

Nobs = numel(xo);
nBlocks = ceil(Nobs / obs_block_size);

p1_cell = cell(nBlocks,1);
p2_cell = cell(nBlocks,1);

fprintf('    [ultra-direct-rz] Nsrc = %d, Nobs = %d\n', numel(Xs), Nobs);
fprintf('    [ultra-direct-rz] src_block_size = %d, obs_block_size = %d, nBlocks = %d\n', ...
    src_block_size, obs_block_size, nBlocks);

t_progress = tic;

if use_parallel
    D = parallel.pool.DataQueue;
    progress_count = 0;
    afterEach(D, @local_update_progress);

    parfor iblk = 1:nBlocks
        [p1_cell{iblk}, p2_cell{iblk}] = calc_ultra_direct_obs_block_from_xyz( ...
            iblk, obs_block_size, xo, yo, zo, ...
            Xs, Ys, Zs, q1, q2, ...
            k1, k2, w1, w2, rho0, src_block_size);
        send(D, iblk);
    end
else
    for iblk = 1:nBlocks
        [p1_cell{iblk}, p2_cell{iblk}] = calc_ultra_direct_obs_block_from_xyz( ...
            iblk, obs_block_size, xo, yo, zo, ...
            Xs, Ys, Zs, q1, q2, ...
            k1, k2, w1, w2, rho0, src_block_size);
        fprintf('        rz obs progress: %5d / %5d blocks (%.1f%%), elapsed = %.1f s\n', ...
            iblk, nBlocks, 100*iblk/nBlocks, toc(t_progress));
    end
end

p1_vec = complex(zeros(Nobs,1));
p2_vec = complex(zeros(Nobs,1));

for iblk = 1:nBlocks
    ib = (iblk-1)*obs_block_size + 1;
    ie = min(iblk*obs_block_size, Nobs);
    idx_obs = ib:ie;

    p1_vec(idx_obs) = p1_cell{iblk};
    p2_vec(idx_obs) = p2_cell{iblk};
end

p1_rz = reshape(p1_vec, size(RHO));
p2_rz = reshape(p2_vec, size(RHO));

    function local_update_progress(~)
        progress_count = progress_count + 1;
        fprintf('        rz obs progress: %5d / %5d blocks (%.1f%%), elapsed = %.1f s\n', ...
            progress_count, nBlocks, 100*progress_count/nBlocks, toc(t_progress));
    end
end

%% ============================================================
% Direct ultrasound calculation for one observation block
%% ============================================================
function [p1_blk, p2_blk] = calc_ultra_direct_obs_block_from_xyz( ...
    iblk, obs_block_size, xo, yo, zo, ...
    Xs, Ys, Zs, q1, q2, k1, k2, w1, w2, rho0, src_block_size)

Nobs = numel(xo);
ib = (iblk-1)*obs_block_size + 1;
ie = min(iblk*obs_block_size, Nobs);
idx_obs = ib:ie;

p1_blk = complex(zeros(numel(idx_obs),1));
p2_blk = complex(zeros(numel(idx_obs),1));

for is = 1:src_block_size:numel(Xs)
    je = min(is + src_block_size - 1, numel(Xs));
    idx_src = is:je;

    dx = xo(idx_obs) - Xs(idx_src).';
    dy = yo(idx_obs) - Ys(idx_src).';
    dz = zo(idx_obs) - Zs(idx_src).';

    R = sqrt(dx.^2 + dy.^2 + dz.^2);
    R(R < 1e-12) = 1e-12;

    G1 = exp(1j * k1 .* R) ./ (4*pi*R);
    G2 = exp(1j * k2 .* R) ./ (4*pi*R);

    p1_blk = p1_blk + (-2j * rho0 * w1) * (G1 * q1(idx_src));
    p2_blk = p2_blk + (-2j * rho0 * w2) * (G2 * q2(idx_src));
end
end

%% ============================================================
% Compute ultrasound pressure using the given spectral-domain G
%% ============================================================
function p_out = compute_ultrasound_pressure_from_G( ...
    G_spec, Vs, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use, rho0, c0, k_use)

F = G_spec .* Vs;
phi = -1j * m_FHT(F, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use);
p_out = 1j * rho0 * c0 * real(k_use) .* phi;
end

%% ============================================================
% Zeroth-order Hankel transform of spatial-domain g, raw version
%% ============================================================
function G_raw = build_green_space_g_transform_raw( ...
    rho_vec, z_vec, k_use, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0)

[RHO, Z] = ndgrid(rho_vec, z_vec);
RR = sqrt(RHO.^2 + Z.^2);
RR_use = max(RR, green_R_min);

g_space = exp(1j * k_use * RR_use) ./ (4*pi * RR_use);
G_raw = m_FHT(g_space, N_FHT, numel(z_vec), Nh, NH, a_solve, x0, x1, k0, 0);
end

%% ============================================================
% Compute audio pa_W from the given Gr00
%% ============================================================
function [pa_W, phia_W] = compute_paW_from_Gr00( ...
    q_full, z, z_audio, Gr00, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa)

Nz  = numel(z);
Nza = numel(z_audio);

absz = [-fliplr(z(2:end)) z];
Nz1  = numel(absz);
N_conv = Nz1 + Nza - 1;

Qr00 = m_FHT(q_full, N_FHT, Nz, Nh, NH, a_solve, x0, x1, k0, ma);
Qr0  = [fliplr(Qr00(:,2:end)) Qr00];

[M_0, N_0] = size(Qr0);
Qr = [Qr0, zeros(M_0, N_conv-N_0)];
Q  = (fft(Qr.')).';

Gr0 = [fliplr(Gr00(:,2:end)) Gr00];
Gr  = [Gr0, zeros(M_0, N_conv-N_0)];
G   = (fft(Gr.')).';

Pa   = Q .* G;
par0 = (ifft(Pa.')).';
par  = par0(:, N_conv-Nza+1:N_conv);

phia0 = m_FHT(par, N_FHT, Nza, NH, Nh, a_solve, x0, x1, k0, ma);
phia_W = -phia0 * delta * 1j / 2;
pa_W   = 1j * rho0 * wa * phia_W;
end

%% ============================================================
% Direct audio calculation at multiple target points from q(rho,z,phi)
%% ============================================================
function pa_pts = compute_audio_direct_points_from_qrzphi_parallel( ...
    rho_grid, z_grid_pos, q_rz_pos, ...
    dr, dz, Nphi, ma, ...
    ka, rho0, wa, ...
    rho_eval, z_eval, ...
    q_block_size, ...
    use_parallel, use_z_mirror)

if use_z_mirror
    [z_grid_use, q_rz_use] = build_mirrored_qrz_about_z0(z_grid_pos, q_rz_pos);
    fprintf('    [audio-direct-rzphi] z-mirror enabled.\n');
else
    z_grid_use = z_grid_pos(:);
    q_rz_use = q_rz_pos;
    fprintf('    [audio-direct-rzphi] z-mirror disabled.\n');
end

rho_grid   = rho_grid(:);
z_grid_use = z_grid_use(:);

Npt = numel(rho_eval);
pa_pts = complex(zeros(Npt,1));

fprintf('    [audio-direct-rzphi] size(q_rz_use) = [%d, %d], eval pts = %d\n', ...
    numel(rho_grid), numel(z_grid_use), Npt);

z_block_size   = min(max(64, floor(q_block_size / max(numel(rho_grid),1))), 512);
phi_block_size = 8;

for ip = 1:Npt
    pa_pts(ip) = calc_audio_direct_single_point_from_qrz_blockwise( ...
        rho_grid, z_grid_use, q_rz_use, ...
        dr, dz, Nphi, ma, ...
        ka, rho0, wa, ...
        rho_eval(ip), z_eval(ip), ...
        z_block_size, phi_block_size);

    fprintf('        audio point progress: %4d / %4d (%.1f%%)\n', ...
        ip, Npt, 100*ip/Npt);
end
end

%% ============================================================
% Mirror q(rho,z) about z = 0
%% ============================================================
function [z_out, q_out] = build_mirrored_qrz_about_z0(z_in, q_in)

z_in = z_in(:);
mask_pos = z_in > 1e-12;

z_neg = -flipud(z_in(mask_pos));
q_neg = fliplr(q_in(:, mask_pos));

z_out = [z_neg; z_in];
q_out = [q_neg, q_in];
end

%% ============================================================
% Direct audio calculation at one point, blockwise over z and phi
%% ============================================================
function p_obs = calc_audio_direct_single_point_from_qrz_blockwise( ...
    rho_grid, z_grid, q_rz, ...
    dr, dz, Nphi, ma, ...
    ka, rho0, wa, ...
    rho_obs, z_obs, ...
    z_block_size, phi_block_size)

phi_all = (0:Nphi-1).' * (2*pi/Nphi);
dphi = 2*pi / Nphi;

Nr = numel(rho_grid);
Nz = numel(z_grid);

I_total = 0;

for iz = 1:z_block_size:Nz
    iz_end = min(iz + z_block_size - 1, Nz);
    z_blk = z_grid(iz:iz_end);
    q_blk = q_rz(:, iz:iz_end);

    rho_row = rho_grid(:);
    z_row   = z_blk(:).';

    for iphi = 1:phi_block_size:Nphi
        iphi_end = min(iphi + phi_block_size - 1, Nphi);
        phi_blk = phi_all(iphi:iphi_end);

        cos_phi   = cos(phi_blk).';
        phase_phi = exp(1j * ma * phi_blk).';

        Rxy2 = rho_obs^2 + rho_row.^2 - 2*(rho_row * (rho_obs * cos_phi));
        Dz2  = (z_obs - z_row).^2;

        R2 = reshape(Rxy2, [Nr, 1, numel(phi_blk)]) + reshape(Dz2, [1, numel(z_blk), 1]);
        R  = sqrt(R2);
        R(R < 1e-12) = 1e-12;

        G = exp(1j * ka .* R) ./ (4*pi*R);

        Qbase = q_blk .* (rho_row * ones(1, numel(z_blk))) * (dr * dz);
        Qfull = reshape(Qbase, [Nr, numel(z_blk), 1]) ...
              .* reshape(phase_phi * dphi, [1, 1, numel(phi_blk)]);

        I_total = I_total + sum(Qfull(:) .* G(:));
    end
end

p_obs = -1j * rho0 * wa * I_total;
end

%% ============================================================
% Build source cfg
%% ============================================================
function source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2)
source_cfg = struct();

source_cfg.profile = 'Vortex-m';
source_cfg.a = a;
source_cfg.v0 = v0;
source_cfg.v_ratio = 1;

source_cfg.m1 = m1;
source_cfg.m2 = m2;

% Kept for compatibility with make_source_velocity.
% For pure-mode cache calls, m1 = m2 = m_use, so this is correct.
% The main transform calculation no longer relies on source_cfg.m
% to distinguish the two modal orders.
source_cfg.m  = m1;

source_cfg.F = 0.2;

source_cfg.f1 = f1;
source_cfg.fa = fa;
source_cfg.f2 = f2;

source_cfg.internal = struct();
end

%% ============================================================
% Build medium cfg
%% ============================================================
function medium_cfg = build_medium_cfg(c, rho0, beta, pref)
medium_cfg = struct();

medium_cfg.c0 = c;
medium_cfg.rho0 = rho0;
medium_cfg.beta = beta;
medium_cfg.pref = pref;

medium_cfg.use_absorp = true;
medium_cfg.atten_handle = @(f) AbsorpAttenCoef(f);

medium_cfg.internal = struct();
end

%% ============================================================
% Build calc cfg for the direct reference
%% ============================================================
function calc_cfg = build_calc_cfg_for_direct_reference( ...
    N_FHT, rho_max, zu_max, za_max, delta, dis_coe, num_workers)

calc_cfg = struct();

calc_cfg.fht = struct();
calc_cfg.fht.N_FHT = N_FHT;
calc_cfg.fht.rho_max = rho_max;
calc_cfg.fht.Nh_scale = 1.2;
calc_cfg.fht.NH_scale = 4;
calc_cfg.fht.Nh_v_scale = 1.1;
calc_cfg.fht.zu_max = zu_max;
calc_cfg.fht.za_max = za_max;
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
% Build a uniform grid over [0, maxv]
%% ============================================================
function g = build_uniform_grid_0_to_max(maxv, d)
n = floor(maxv / d);
g = (0:n) * d;
if isempty(g) || abs(g(end) - maxv) > 1e-12
    g = [g, maxv];
end
g = unique(g, 'stable');
end

%% ============================================================
% Convert pressure amplitude to SPL
%% ============================================================
function spl = local_pressure_to_spl(p_amp, pref)
spl = 20 * log10(p_amp / max(pref, eps) / sqrt(2) + eps);
end

%% ============================================================
% Build pu cache subfolder name
%% ============================================================
function folder_name = build_pu_cache_folder_name(m1, m2, rho_max, zu_max, N_FHT, delta, dis_coe, Nphi)

s_rho   = local_num_to_tag(rho_max);
s_zu    = local_num_to_tag(zu_max);
s_delta = local_num_to_tag(delta);

folder_name = sprintf('m1%dm2%drm%s_zu%s_nf%d_dt%s_dc%d_np%d', ...
    m1, m2, s_rho, s_zu, N_FHT, s_delta, dis_coe, Nphi);
end

%% ============================================================
% Convert a numeric value to a compact tag
%% ============================================================
function s = local_num_to_tag(x)

if x == 0
    s = '0';
    return;
end

ax = abs(x);

if ax >= 1e-2 && ax < 1e3
    s = sprintf('%.6f', x);
else
    s = sprintf('%.6e', x);
end

s = strrep(s, '-', 'm');
s = strrep(s, '+', '');
s = strrep(s, '.', 'p');

while ~isempty(s) && s(end) == '0'
    s(end) = [];
end
if ~isempty(s) && s(end) == 'p'
    s(end) = [];
end

s = strrep(s, 'em0', 'em');
s = strrep(s, 'ep0', 'ep');
end

%% ============================================================
% Save all parameters to txt
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

%% ============================================================
% Set the common boxed axis style
%% ============================================================
function local_set_axis_style_box(ax)

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

end

%% ============================================================
% Set sparse logarithmic ticks for N_FHT
%% ============================================================
function local_set_log_xticks_nfht_sparse(ax, x)

x = x(:);
x = x(isfinite(x) & x > 0);

if isempty(x)
    return;
end

xmin = min(x);
xmax = max(x);

tick_candidates = [1e2 1e3 1e4];
tick_use = tick_candidates(tick_candidates >= xmin & tick_candidates <= xmax);

if isempty(tick_use)
    tick_use = unique([xmin; xmax]);
elseif numel(tick_use) == 1
    tick_use = unique([xmin; tick_use(:); xmax]);
end

ax.XTick = tick_use(:).';

if xmin == xmax
    xlim(ax, [0.9 * xmin, 1.1 * xmax]);
else
    xlim(ax, [xmin, xmax]);
end

end

%% ============================================================
% Set logarithmic ticks for delta_z
%% ============================================================
function local_set_log_xticks_delta(ax, x)

x = x(:);
x = x(isfinite(x) & x > 0);

if isempty(x)
    return;
end

xmin = min(x);
xmax = max(x);

pmin = floor(log10(xmin));
pmax = ceil(log10(xmax));

major_ticks = 10.^(pmin:pmax);
major_ticks = major_ticks(major_ticks >= xmin & major_ticks <= xmax);

if isempty(major_ticks)
    major_ticks = unique([xmin; xmax]);
elseif numel(major_ticks) == 1
    major_ticks = unique([xmin; major_ticks(:); xmax]);
end

ax.XTick = major_ticks(:).';

if xmin == xmax
    xlim(ax, [0.9 * xmin, 1.1 * xmax]);
else
    xlim(ax, [xmin, xmax]);
end

end

%% ============================================================
% Set y-axis margins and align limits to 5 dB
%% ============================================================
function local_set_y_margin(ax, y)

y = y(:);
y = y(isfinite(y));

if isempty(y)
    ylim(ax, [-5 5]);
    return;
end

ymin = min(y);
ymax = max(y);

if ymin == ymax
    ymin = ymin - 1;
    ymax = ymax + 1;
end

pad = 0.08 * (ymax - ymin);
yl = [ymin - pad, ymax + pad];

yl(1) = 5 * floor(yl(1) / 5);
yl(2) = 5 * ceil(yl(2) / 5);

ylim(ax, yl);

end

%% ============================================================
% Save one figure as png, fig, and vector pdf
%% ============================================================
function local_save_figure_png_fig_pdf(fig_handle, png_path, fig_path, pdf_path)

saveas(fig_handle, png_path);
savefig(fig_handle, fig_path);

try
    exportgraphics(fig_handle, pdf_path, 'ContentType', 'vector');
catch
    set(fig_handle, 'PaperPositionMode', 'auto');
    print(fig_handle, pdf_path, '-dpdf', '-painters');
end

end

%% ============================================================
% Safely set x-axis limits and avoid errors when there is only one point
%% ============================================================
function local_set_x_limit(ax, x)
x = x(:);
x = x(isfinite(x));

if isempty(x)
    return;
end

xmin = min(x);
xmax = max(x);

if xmin == xmax
    if xmin == 0
        xlim(ax, [-1 1]);
    else
        span = 0.1 * abs(xmin);
        if span == 0
            span = 1;
        end
        xlim(ax, [xmin - span, xmax + span]);
    end
else
    xlim(ax, [xmin xmax]);
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