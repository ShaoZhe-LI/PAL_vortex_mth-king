%% ============================================================
% Fixed-case comparison of audio field:
%   king method vs direct method
%   with local-effect correction
%
% Fixed parameters:
%   N_FHT = 16384
%   delta = 0.001
%
% Other parameters are kept as consistent as possible with the original Audio convergence analysis code.
%
% Calculation contents:
%   1) King-method audio field, without local effect.
%   2) Direct-integration audio field, without local effect.
%   3) Call the ultrasound pressure / velocity solvers at preset target points, using both methods.
%   4) Add local-effect correction for both methods.
%   5) Compare axial / radial results after including local effect.
%
% External dependencies:
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%   - make_source_velocity.m
%   - calc_ultrasound_field.m
%   - calc_ultrasound_velocity_field.m
% Updated.
%% ============================================================
%% ============================================================
% NOTE ON CACHE REUSE
%
% New caches store the physical and source parameters in cache_meta,
% including v0, a, f1, f2, fa, c, rho0, beta, source.profile,
% m1, m2, rho_max, zu_max, N_FHT, delta, dis_coe, and Nphi.
%
% When a cache is loaded, all parameters except v0 must match the current
% setup. If only v0 differs, the cached pressure fields are rescaled by
%   alpha_v0 = v0_current / v0_cache.
%
% The virtual-source term q is either rescaled by alpha_v0^2 for a full
% cache, or recomputed from the rescaled p1 and p2 for split caches.
%
% Old caches without cache_meta.v0 are not reused, because their velocity
% normalization cannot be verified safely.
%
% For mixed mode m1 ~= m2, this script still avoids mixed full-cache reuse
% and uses pure-mode split caches:
%   p1 <- cache with m1 = m2 = m1
%   p2 <- cache with m1 = m2 = m2
%% ============================================================

clear; clc; close all;

for i = 1:3

clearvars -except i
close all;

%% ===================== Save and display options =====================
save_figures      = true;
save_calc_results = true;
save_params_txt   = true;
show_figures      = true;

%% ===================== Basic physical parameters =====================
a     = 0.05;
v0    = 0.108;
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
    m2 = m1;
end
if i == 2
    m1 = 3;
    m2 = m1;
end
if i == 3
    m1 = 0;
    m2 = 3;
end

ma = m2 - m1;

% Cache scaling is determined from cache_meta.v0.
% No fixed reference velocity is assumed.

%% ===================== Fixed FHT parameters =====================
N_FHT = 16384 * 1;
delta = 0.001;

%% ===================== Other calculation parameters =====================
rho_max = 0.5;
zu_max  = 15.0;
za_max  = 5 + delta;

green_R_min = 1e-12;

%% ===================== Parallel settings =====================
use_parallel = true;
num_workers  = 20;

%% ===================== Direct-reference calculation settings =====================
direct_ref_cfg = struct();
direct_ref_cfg.src_block_size_ultra = 120000;
direct_ref_cfg.obs_block_size_rz    = 2048;
direct_ref_cfg.audio_q_block_size   = 200000;
direct_ref_cfg.audio_point_batch   = [];

direct_ref_cfg.N_FHT       = 65536;
direct_ref_cfg.delta       = c / f2 / 8;
direct_ref_cfg.dis_coe     = 32;
direct_ref_cfg.num_workers = num_workers;

%% ===================== Audio direct-rz-phi grid settings =====================
direct2d = struct();
direct2d.dr   = direct_ref_cfg.delta;
direct2d.dz   = direct_ref_cfg.delta;
direct2d.Nphi = 180;

audio_direct_use_z_mirror = true;

%% ===================== Line-result settings =====================
nLinePts = 100;

rho_axis_fixed = [];
theta_line_deg = [];

fprintf('case %d: m1 = %d, m2 = %d\n', i, m1, m2);

if m1 == m2
    disp('m1 = m2, axial and radial');

    z_axis_line_min = 1e-2;
    z_axis_line_max = 5.0;
    rho_axis_fixed  = 0.0;

    z_axis_req   = logspace(log10(z_axis_line_min), log10(z_axis_line_max), nLinePts).';
    rho_axis_req = rho_axis_fixed * ones(nLinePts,1);

    z_radial_fixed = 0.2;
    rho_radial_min = 1e-3;
    rho_radial_max = 0.5;

    rho_radial_req = logspace(log10(rho_radial_min), log10(rho_radial_max), nLinePts).';
    z_radial_req   = z_radial_fixed * ones(nLinePts,1);

else
    disp('m1 ~= m2, oblique and radial');

    theta_line_deg = 7.5;

    z_axis_line_min = 1e-2;
    z_axis_line_max_by_rho = rho_max / tand(theta_line_deg);
    z_axis_line_max = min(za_max, z_axis_line_max_by_rho);

    z_axis_req   = logspace(log10(z_axis_line_min), log10(z_axis_line_max), nLinePts).';
    rho_axis_req = z_axis_req * tand(theta_line_deg);

    z_radial_fixed = 0.2;
    rho_radial_min = 1e-3;
    rho_radial_max = 0.5;

    rho_radial_req = logspace(log10(rho_radial_min), log10(rho_radial_max), nLinePts).';
    z_radial_req   = z_radial_fixed * ones(nLinePts,1);
end

%% ===================== Save path =====================
time_tag  = datestr(now, 'mmdd_HHMMSS');
save_root_parent = 'result_fixed_compare_final';
save_root = fullfile(save_root_parent, sprintf('AudioFixed_LocalCorr_%s_m1_%d_m2_%d', time_tag, m1, m2));

if (save_figures || save_calc_results || save_params_txt) && ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% ===================== Direct-ultrasound cache path =====================
cache_root = 'data_pu';
if ~exist(cache_root, 'dir')
    mkdir(cache_root);
end

pu_cache_folder_name = build_pu_cache_folder_name( ...
    m1, m2, rho_max, zu_max, ...
    direct_ref_cfg.N_FHT, direct_ref_cfg.delta, direct_ref_cfg.dis_coe, ...
    direct2d.Nphi);

pu_cache_dir = fullfile(cache_root, pu_cache_folder_name);
pu_cache_mat = fullfile(pu_cache_dir, 'pu_direct_rz_cache.mat');
pu_cache_txt = fullfile(pu_cache_dir, 'pu_cache_parameters.txt');

%% ===================== Parallel pool initialization =====================
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

%% ===================== Save parameters =====================
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
    meta.N_FHT = N_FHT;
    meta.delta = delta;
    meta.rho_max = rho_max;
    meta.zu_max = zu_max;
    meta.za_max = za_max;
    meta.green_R_min = green_R_min;
    meta.direct_ref_cfg = direct_ref_cfg;
    meta.direct2d = direct2d;
    meta.audio_direct_use_z_mirror = audio_direct_use_z_mirror;
    meta.use_parallel = use_parallel;
    meta.num_workers = num_workers;
    meta.save_root = save_root;
    meta.time_tag = time_tag;
    meta.pu_cache_dir = pu_cache_dir;
    meta.pu_cache_mat = pu_cache_mat;
    meta.local_effect_formula = 'p_corr = p_wo - [rho0/2 * conj(v1)·v2 - (w1/w2 + w2/w1 - 1) * conj(p1)p2 / (2 rho0 c^2)]';

    meta.axis_line = struct();
    meta.axis_line.nLinePts = nLinePts;
    meta.axis_line.z_min = z_axis_line_min;
    meta.axis_line.z_max = z_axis_line_max;
    meta.axis_line.rho_req = rho_axis_req;
    meta.axis_line.z_req = z_axis_req;

    if m1 ~= m2
        meta.axis_line.type = 'oblique';
        meta.axis_line.theta_deg = theta_line_deg;
        meta.axis_line.rho_fixed = [];
    else
        meta.axis_line.type = 'axial';
        meta.axis_line.theta_deg = [];
        meta.axis_line.rho_fixed = rho_axis_fixed;
    end

    meta.radial_line = struct();
    meta.radial_line.nLinePts = nLinePts;
    meta.radial_line.z_fixed = z_radial_fixed;
    meta.radial_line.rho_min = rho_radial_min;
    meta.radial_line.rho_max = rho_radial_max;
    meta.radial_line.rho_req = rho_radial_req;
    meta.radial_line.z_req = z_radial_req;

    meta.cache_reuse_note = 'Cache scaling is determined from cache_meta.v0; no fixed v0 reference is assumed.';

    local_write_all_params_txt(fullfile(save_root, 'all_parameters.txt'), meta);
end

%% ===================== Build source / medium / calc cfg =====================
source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2);
medium_cfg = build_medium_cfg(c, rho0, beta, pref);

%% ===================== Prepare direct-reference source discretization =====================
fprintf('Preparing direct-reference source discretization...\n');

calc_cfg_ref = build_calc_cfg_for_direct_reference( ...
    direct_ref_cfg.N_FHT, rho_max, zu_max, za_max, direct_ref_cfg.delta, ...
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

%% ===================== Direct-rz grid =====================
rho_rz = build_uniform_grid_0_to_max(rho_max, direct2d.dr);
z_rz   = build_uniform_grid_0_to_max(zu_max,  direct2d.dz);

%% ============================================================
% Direct-ultrasound cache loading logic
%
% m1 == m2:
%   A full cache can be used.
%
% m1 ~= m2:
%   A mixed full cache is not used.
%   Pure-mode split caches are enforced:
%       p1 <- cache with m1 = m2 = m1
%       p2 <- cache with m1 = m2 = m2
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Preparing direct ultrasound field on uniform (rho,z) grid ...\n');
fprintf('Cache folder: %s\n', pu_cache_dir);
fprintf('============================================================\n');

time_direct_rz_once = NaN;
loaded_from_cache = false;
loaded_from_split_cache = false;

use_full_cache = (m1 == m2) && exist(pu_cache_mat, 'file');

if use_full_cache
    fprintf('Found full cached direct-ultrasound data. Loading...\n');
    t_load = tic;

    S_cache = load(pu_cache_mat, ...
        'p1_direct_rz', 'p2_direct_rz', 'q_direct_rz', ...
        'rho_rz', 'z_rz', 'cache_meta');

    if ~isfield(S_cache, 'cache_meta')
        error(['Cached data does not contain cache_meta. ' ...
               'Old caches without metadata cannot be safely reused. ' ...
               'Please delete cache folder and regenerate: %s'], pu_cache_dir);
    end

    if ~isequal(rho_rz(:), S_cache.rho_rz(:)) || ~isequal(z_rz(:), S_cache.z_rz(:))
        error(['Cached rho/z grid does not match current setup. ' ...
               'Please delete cache folder: %s'], pu_cache_dir);
    end

    alpha_v0 = local_validate_cache_meta_and_get_v0_scale( ...
        S_cache.cache_meta, ...
        a, v0, c, rho0, beta, f1, f2, fa, ...
        source_cfg.profile, ...
        m1, m2, rho_max, zu_max, ...
        direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
        direct_ref_cfg.dis_coe, direct2d.Nphi);

    p1_direct_rz = S_cache.p1_direct_rz * alpha_v0;
    p2_direct_rz = S_cache.p2_direct_rz * alpha_v0;

    if isfield(S_cache, 'q_direct_rz') && ~isempty(S_cache.q_direct_rz)
        q_direct_rz = S_cache.q_direct_rz * (alpha_v0^2);
    else
        q_direct_rz = conj(p1_direct_rz) .* p2_direct_rz ...
            * beta*(2*pi*fa)/(1j*rho0^2*c^4);
    end

    time_direct_rz_once = toc(t_load);
    loaded_from_cache = true;
    loaded_from_split_cache = false;

    fprintf('Full cache loaded. elapsed = %.2f s\n', time_direct_rz_once);
    fprintf('Cached pu scaled by v0 ratio: current/cache = %.8f / %.8f = %.8f\n', ...
        v0, S_cache.cache_meta.v0, alpha_v0);

else
    if m1 ~= m2
        fprintf('Mixed-mode case detected: m1 = %d, m2 = %d.\n', m1, m2);
        fprintf('Full mixed cache will NOT be used. Trying pure-mode split-cache reuse...\n');
    else
        fprintf('Full cache not found. Trying pure-mode split-cache reuse...\n');
    end

    t_load = tic;

    info_p1 = find_reusable_single_field_cache( ...
        cache_root, 'p1', m1, ...
        rho_max, zu_max, ...
        direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
        direct_ref_cfg.dis_coe, direct2d.Nphi);

    info_p2 = find_reusable_single_field_cache( ...
        cache_root, 'p2', m2, ...
        rho_max, zu_max, ...
        direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
        direct_ref_cfg.dis_coe, direct2d.Nphi);

    if info_p1.found && info_p2.found
        fprintf('Pure-mode split caches found.\n');
        fprintf('  p1 from: %s\n', info_p1.cache_dir);
        fprintf('  p2 from: %s\n', info_p2.cache_dir);

        S1 = load(info_p1.cache_mat, 'p1_direct_rz', 'rho_rz', 'z_rz', 'cache_meta');
        S2 = load(info_p2.cache_mat, 'p2_direct_rz', 'rho_rz', 'z_rz', 'cache_meta');

        if ~isfield(S1, 'cache_meta')
            error(['Split cache for p1 does not contain cache_meta. ' ...
                   'Old caches without metadata cannot be safely reused. ' ...
                   'Please regenerate cache: %s'], info_p1.cache_mat);
        end
        if ~isfield(S2, 'cache_meta')
            error(['Split cache for p2 does not contain cache_meta. ' ...
                   'Old caches without metadata cannot be safely reused. ' ...
                   'Please regenerate cache: %s'], info_p2.cache_mat);
        end

        if ~isequal(rho_rz(:), S1.rho_rz(:)) || ~isequal(z_rz(:), S1.z_rz(:))
            error('Split cache for p1 has mismatched rho/z grid.');
        end
        if ~isequal(rho_rz(:), S2.rho_rz(:)) || ~isequal(z_rz(:), S2.z_rz(:))
            error('Split cache for p2 has mismatched rho/z grid.');
        end

        alpha_v0_p1 = local_validate_cache_meta_and_get_v0_scale( ...
            S1.cache_meta, ...
            a, v0, c, rho0, beta, f1, f2, fa, ...
            source_cfg.profile, ...
            m1, m1, rho_max, zu_max, ...
            direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
            direct_ref_cfg.dis_coe, direct2d.Nphi);

        alpha_v0_p2 = local_validate_cache_meta_and_get_v0_scale( ...
            S2.cache_meta, ...
            a, v0, c, rho0, beta, f1, f2, fa, ...
            source_cfg.profile, ...
            m2, m2, rho_max, zu_max, ...
            direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
            direct_ref_cfg.dis_coe, direct2d.Nphi);

        p1_direct_rz = S1.p1_direct_rz * alpha_v0_p1;
        p2_direct_rz = S2.p2_direct_rz * alpha_v0_p2;

        q_direct_rz = conj(p1_direct_rz) .* p2_direct_rz ...
            * beta*(2*pi*fa)/(1j*rho0^2*c^4);

        time_direct_rz_once = toc(t_load);
        loaded_from_cache = true;
        loaded_from_split_cache = true;

        fprintf('Pure-mode split cache loaded. elapsed = %.2f s\n', time_direct_rz_once);
        fprintf('p1 cache v0 scaling: current/cache = %.8f / %.8f = %.8f\n', ...
            v0, S1.cache_meta.v0, alpha_v0_p1);
        fprintf('p2 cache v0 scaling: current/cache = %.8f / %.8f = %.8f\n', ...
            v0, S2.cache_meta.v0, alpha_v0_p2);

    else
        if m1 ~= m2
            error(['Pure-mode split cache not found for mixed-mode case m1=%d, m2=%d.\n' ...
                   'Because make_source_velocity currently only supports a single m, ' ...
                   'this script stops here to avoid generating an incorrect mixed-mode cache.\n' ...
                   'Please first generate pure-mode caches m1=m2=%d and m1=m2=%d.'], ...
                   m1, m2, m1, m2);
        end

        fprintf('Split-cache reuse failed. Recomputing full direct ultrasound field for single-mode case ...\n');

        if ~exist(pu_cache_dir, 'dir')
            mkdir(pu_cache_dir);
        end

        t_direct_rz = tic;

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

        time_direct_rz_once = toc(t_direct_rz);
        fprintf('Direct ultrasound rz-grid finished. elapsed = %.2f s\n', time_direct_rz_once);

        cache_meta = local_build_direct_cache_meta( ...
            source_cfg, medium_cfg, ...
            fu, fa, f1, f2, ...
            m1, m2, ma, ...
            rho_max, zu_max, ...
            direct_ref_cfg.N_FHT, direct_ref_cfg.delta, ...
            direct_ref_cfg.dis_coe, direct2d.Nphi);
        cache_meta.loaded_from_cache = false;
        cache_meta.loaded_from_split_cache = false;

        save(pu_cache_mat, ...
            'p1_direct_rz', 'p2_direct_rz', 'q_direct_rz', ...
            'rho_rz', 'z_rz', 'cache_meta', '-v7.3');

        local_write_all_params_txt(pu_cache_txt, cache_meta);
        fprintf('Direct ultrasound cache saved: %s\n', pu_cache_mat);
    end
end

%% ============================================================
% Single fixed case: compute transform + direct + local correction
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Fixed case calculation WITH local-effect correction: N_FHT = %d, delta = %.6f m\n', N_FHT, delta);
fprintf('============================================================\n');

t_case = tic;

out_fixed = evaluate_fixed_case_transform_and_direct_with_local_effect( ...
    N_FHT, delta, za_max, ...
    source_cfg, medium_cfg, ...
    c, rho0, beta, pref, f1, f2, fa, ...
    rho_max, zu_max, ...
    green_R_min, ma, ...
    rho_rz, z_rz, q_direct_rz, direct2d, ...
    direct_ref_cfg.audio_q_block_size, ...
    direct_ref_cfg.audio_point_batch, ...
    use_parallel, audio_direct_use_z_mirror, ...
    rho_axis_req, z_axis_req, rho_radial_req, z_radial_req);

time_fixed_case = toc(t_case);
fprintf('Fixed case done. elapsed = %.2f s\n', time_fixed_case);

%% ============================================================
% Organize data
%% ============================================================
axis_res = out_fixed.axis_res;
radial_res = out_fixed.radial_res;

%% ============================================================
% Clean plotting and saving
%% ============================================================
p_ref = 20e-6;
p_floor = 1e-16;

rho_show_min = 1e-2;
z_bad = 0.072;

if m1 ~= m2
    use_oblique_line = true;
    line_xlab_text = sprintf(['$\\mathrm{Axial\\ coordinate,}\\ z\\ (\\mathrm{m}),\\ ', ...
                              '\\rho=z\\tan %.1f^\\circ$'], theta_line_deg);
    line_file_prefix = 'Oblique';
else
    use_oblique_line = false;
    line_xlab_text = '$\mathrm{Axial\ distance,}\ z\ (\mathrm{m})$';
    line_file_prefix = 'Axis';
end

line_xlab_interpreter = 'latex';

radial_xlab_text = '$\mathrm{Radial\ distance,}\ \rho\ (\mathrm{m})$';
spl_ylab_text    = '$\mathrm{SPL}\ (\mathrm{dB})$';
diff_ylab_text   = '$|\Delta \mathrm{SPL}|\ (\mathrm{dB})$';

sty_dim_wo  = {'Color',[0.90 0.35 0.95], 'LineStyle','-',  'LineWidth',3.0};
sty_dim_w   = {'Color',[0.25 0.70 0.90], 'LineStyle','-',  'LineWidth',3.0};
sty_king_wo = {'Color',[0.85 0.00 0.00], 'LineStyle','--', 'LineWidth',3.0};
sty_king_w  = {'Color',[0.10 0.15 0.95], 'LineStyle','--', 'LineWidth',3.0};

sty_err_w   = {'Color',[0.10 0.15 0.95], 'LineStyle','-',  'LineWidth',2.8};
sty_err_wo  = {'Color',[0.95 0.00 0.00], 'LineStyle','-',  'LineWidth',2.8};

idx_axis_plot = axis_res.z_actual > 0;

if exist('z_bad', 'var') && ~isempty(z_bad)
    [~, idx_bad_local] = min(abs(axis_res.z_actual - z_bad));
    idx_axis_plot(idx_bad_local) = false;
end

z_plot_axis = axis_res.z_actual(idx_axis_plot);

spl_t_axis_wo = 20*log10(max(abs(axis_res.pa_transform_wo(idx_axis_plot)), p_floor) / p_ref);
spl_d_axis_wo = 20*log10(max(abs(axis_res.pa_direct_wo(idx_axis_plot)),    p_floor) / p_ref);
spl_t_axis    = 20*log10(max(abs(axis_res.pa_transform(idx_axis_plot)),    p_floor) / p_ref);
spl_d_axis    = 20*log10(max(abs(axis_res.pa_direct(idx_axis_plot)),       p_floor) / p_ref);

diff_spl_axis_wo = abs(spl_t_axis_wo - spl_d_axis_wo);
diff_spl_axis    = abs(spl_t_axis    - spl_d_axis);

rho_plot_radial = radial_res.rho_actual(:);

spl_t_radial_wo = 20*log10(max(abs(radial_res.pa_transform_wo(:)), p_floor) / p_ref);
spl_d_radial_wo = 20*log10(max(abs(radial_res.pa_direct_wo(:)),    p_floor) / p_ref);
spl_t_radial    = 20*log10(max(abs(radial_res.pa_transform(:)),    p_floor) / p_ref);
spl_d_radial    = 20*log10(max(abs(radial_res.pa_direct(:)),       p_floor) / p_ref);

diff_spl_radial_wo = abs(spl_t_radial_wo - spl_d_radial_wo);
diff_spl_radial    = abs(spl_t_radial    - spl_d_radial);

rho_show_max = max(rho_plot_radial(isfinite(rho_plot_radial)));
idx_radial_show = rho_plot_radial >= rho_show_min & ...
                  rho_plot_radial <= rho_show_max & ...
                  isfinite(rho_plot_radial);

fig_pos = [100 100 980 520];
ax_pos  = [0.17 0.32 0.78 0.46];

fs_ax   = 18;
fs_lab  = 24;
fs_leg  = 19;

fs_x_ax  = fs_ax  * 1.35;
fs_y_ax  = fs_ax  * 1.05;
fs_x_lab = fs_lab * 1.35;
fs_y_lab = fs_lab * 1.20;

fig1 = local_create_figure(show_figures, 'Color','w', 'Position', fig_pos);
ax1 = axes('Parent', fig1);
hold(ax1, 'on');

p1 = plot(ax1, z_plot_axis, spl_d_axis_wo, sty_dim_wo{:});
p2 = plot(ax1, z_plot_axis, spl_d_axis,    sty_dim_w{:});
p3 = plot(ax1, z_plot_axis, spl_t_axis_wo, sty_king_wo{:});
p4 = plot(ax1, z_plot_axis, spl_t_axis,    sty_king_w{:});

set(ax1, 'XScale', 'log');
xlim(ax1, [min(z_plot_axis), max(z_plot_axis)]);

local_set_clean_axis_style(ax1);
ax1.XAxis.FontSize = fs_x_ax;
ax1.YAxis.FontSize = fs_y_ax;
local_set_log_xticks_line(ax1, z_plot_axis);
local_set_y_margin_general(ax1, [spl_d_axis_wo; spl_d_axis; spl_t_axis_wo; spl_t_axis]);

xlabel(ax1, line_xlab_text, ...
    'FontSize', fs_x_lab, ...
    'Interpreter', line_xlab_interpreter);

ylabel(ax1, spl_ylab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_y_lab);

ax1.ActivePositionProperty = 'position';
ax1.Position = ax_pos;

lgd1 = legend(ax1, [p1 p2 p3 p4], ...
    {'DIM (W/o local)', 'DIM (W/ local)', 'King (W/o local)', 'King (W/ local)'}, ...
    'NumColumns', 2, ...
    'Location', 'northoutside');
lgd1.Box = 'on';
lgd1.FontSize = fs_leg;
lgd1.Interpreter = 'tex';
lgd1.AutoUpdate = 'off';

ax1.Position = ax_pos;

fig2 = local_create_figure(show_figures, 'Color','w', 'Position', fig_pos + [20 20 0 0]);
ax2 = axes('Parent', fig2);
hold(ax2, 'on');

plot(ax2, z_plot_axis, diff_spl_axis,    sty_err_w{:});
plot(ax2, z_plot_axis, diff_spl_axis_wo, sty_err_wo{:});

set(ax2, 'XScale', 'log');
xlim(ax2, [min(z_plot_axis), max(z_plot_axis)]);

local_set_clean_axis_style(ax2);
ax2.XAxis.FontSize = fs_x_ax;
ax2.YAxis.FontSize = fs_y_ax;
local_set_log_xticks_line(ax2, z_plot_axis);
local_set_y_margin_from_zero(ax2, [diff_spl_axis; diff_spl_axis_wo]);

xlabel(ax2, line_xlab_text, ...
    'FontSize', fs_x_lab, ...
    'Interpreter', line_xlab_interpreter);

ylabel(ax2, diff_ylab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_y_lab);

ax2.ActivePositionProperty = 'position';
ax2.Position = ax_pos;

lgd2 = legend(ax2, {'(W/ local)', '(W/o local)'}, 'Location', 'best');
lgd2.Box = 'on';
lgd2.FontSize = fs_leg;
lgd2.Interpreter = 'tex';
lgd2.AutoUpdate = 'off';

fig3 = local_create_figure(show_figures, 'Color','w', 'Position', fig_pos + [40 40 0 0]);
ax3 = axes('Parent', fig3);
hold(ax3, 'on');

q1 = plot(ax3, rho_plot_radial, spl_d_radial_wo, sty_dim_wo{:});
q2 = plot(ax3, rho_plot_radial, spl_d_radial,    sty_dim_w{:});
q3 = plot(ax3, rho_plot_radial, spl_t_radial_wo, sty_king_wo{:});
q4 = plot(ax3, rho_plot_radial, spl_t_radial,    sty_king_w{:});

set(ax3, 'XScale', 'log');
xlim(ax3, [rho_show_min, rho_show_max]);

local_set_clean_axis_style(ax3);
ax3.XAxis.FontSize = fs_x_ax;
ax3.YAxis.FontSize = fs_y_ax;
local_set_log_xticks_line(ax3, rho_plot_radial(rho_plot_radial > 0));

local_set_y_margin_general(ax3, ...
    [spl_d_radial_wo(idx_radial_show); ...
     spl_d_radial(idx_radial_show); ...
     spl_t_radial_wo(idx_radial_show); ...
     spl_t_radial(idx_radial_show)]);

xlabel(ax3, radial_xlab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_x_lab);

ylabel(ax3, spl_ylab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_y_lab);

ax3.ActivePositionProperty = 'position';
ax3.Position = ax_pos;

lgd3 = legend(ax3, [q1 q2 q3 q4], ...
    {'DIM (W/o local)', 'DIM (W/ local)', 'King (W/o local)', 'King (W/ local)'}, ...
    'NumColumns', 2, ...
    'Location', 'northoutside');
lgd3.Box = 'on';
lgd3.FontSize = fs_leg;
lgd3.Interpreter = 'tex';
lgd3.AutoUpdate = 'off';

ax3.Position = ax_pos;

fig4 = local_create_figure(show_figures, 'Color','w', 'Position', fig_pos + [60 60 0 0]);
ax4 = axes('Parent', fig4);
hold(ax4, 'on');

plot(ax4, rho_plot_radial, diff_spl_radial,    sty_err_w{:});
plot(ax4, rho_plot_radial, diff_spl_radial_wo, sty_err_wo{:});

set(ax4, 'XScale', 'log');
xlim(ax4, [rho_show_min, rho_show_max]);

local_set_clean_axis_style(ax4);
ax4.XAxis.FontSize = fs_x_ax;
ax4.YAxis.FontSize = fs_y_ax;
local_set_log_xticks_line(ax4, rho_plot_radial(rho_plot_radial > 0));

local_set_y_margin_from_zero(ax4, ...
    [diff_spl_radial(idx_radial_show); ...
     diff_spl_radial_wo(idx_radial_show)]);

xlabel(ax4, radial_xlab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_x_lab);

ylabel(ax4, diff_ylab_text, ...
    'Interpreter', 'latex', ...
    'FontSize', fs_y_lab);

ax4.ActivePositionProperty = 'position';
ax4.Position = ax_pos;

lgd4 = legend(ax4, {'(W/ local)', '(W/o local)'}, 'Location', 'best');
lgd4.Box = 'on';
lgd4.FontSize = fs_leg;
lgd4.Interpreter = 'tex';
lgd4.AutoUpdate = 'off';

%% ============================================================
% Save figures and results
%% ============================================================
if save_figures
    save_one(fig1, save_root, [line_file_prefix 'SPL_replot']);
    save_one(fig2, save_root, [line_file_prefix 'SPLDifference_replot']);
    save_one(fig3, save_root, 'RadialSPL_replot');
    save_one(fig4, save_root, 'RadialSPLDifference_replot');
end

if ~show_figures
    close(fig1);
    close(fig2);
    close(fig3);
    close(fig4);
end

if save_calc_results
    save(fullfile(save_root, 'fixed_case_compare_results_localCorr.mat'), ...
        'out_fixed', 'axis_res', 'radial_res', ...
        'N_FHT', 'delta', 'rho_max', 'zu_max', 'za_max', ...
        'rho_axis_req', 'z_axis_req', 'rho_radial_req', 'z_radial_req', ...
        'time_direct_rz_once', 'time_fixed_case', ...
        'loaded_from_cache', 'loaded_from_split_cache', ...
        'pu_cache_dir', '-v7.3');
end

fprintf('\nAll done.\n');
fprintf('Direct ultrasound precompute/load time = %.2f s\n', time_direct_rz_once);
fprintf('Fixed case elapsed time = %.2f s\n', time_fixed_case);
fprintf('Direct ultrasound cache dir = %s\n', pu_cache_dir);
fprintf('Results folder: %s\n', save_root);

if use_oblique_line
    fprintf('Line profile display: oblique, theta = %.1f deg, rho = z tan(theta)\n', theta_line_deg);
else
    fprintf('Line profile display: axial\n');
end

fprintf('Radial display range: rho = %.6g to %.6g m\n', rho_show_min, rho_show_max);

end

%% ============================================================
% Fixed single case: compute transform / direct and apply local correction
%% ============================================================
function out = evaluate_fixed_case_transform_and_direct_with_local_effect( ...
    N_FHT, delta, za_max, ...
    source_cfg, medium_cfg, ...
    c, rho0, beta, pref, f1, f2, fa, ...
    rho_max, zu_max, ...
    green_R_min, ma, ...
    rho_rz, z_rz, q_direct_rz, direct2d, ...
    audio_q_block_size, ...
    audio_point_batch, ...
    use_parallel, use_z_mirror, ...
    rho_axis_req, z_axis_req, rho_radial_req, z_radial_req)

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

%% ---- audio transform (without local effect) ----
Ga_raw = build_green_space_g_transform_raw( ...
    xh, z, ka, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);
Ga = (-4*pi*1j) * Ga_raw;

[pa_transform_full_wo, ~] = compute_paW_from_Gr00( ...
    q_transform, z, z_audio, Ga, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa);

%% ---- extract wo-local-effect line points ----
axis_res = extract_line_from_transform_and_direct( ...
    xh, z_audio, pa_transform_full_wo, ...
    rho_axis_req, z_axis_req, ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d, ma, ka, rho0, wa, pref, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror);

radial_res = extract_line_from_transform_and_direct( ...
    xh, z_audio, pa_transform_full_wo, ...
    rho_radial_req, z_radial_req, ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d, ma, ka, rho0, wa, pref, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror);

%% ---- local-effect correction on target points only ----
calc_cfg_ultra = build_calc_cfg_for_ultra_both( ...
    N_FHT, rho_max, zu_max, za_max, delta, use_parallel, medium_cfg, source_cfg);

[~, fht_tmp] = make_source_velocity(source_cfg, medium_cfg, calc_cfg_ultra);
rho_king_grid = fht_tmp.xh(:);
z_king_grid   = fht_tmp.z_ultra(:);

if numel(unique(rho_axis_req)) == 1
    axis_mode = 'axial';
else
    axis_mode = 'oblique';
end

axis_res = apply_local_effect_correction_to_line_bothstyle( ...
    axis_res, axis_mode, ...
    rho_axis_req, z_axis_req, ...
    rho_king_grid, z_king_grid, ...
    source_cfg, medium_cfg, calc_cfg_ultra);

radial_res = apply_local_effect_correction_to_line_bothstyle( ...
    radial_res, 'radial', ...
    rho_radial_req, z_radial_req, ...
    rho_king_grid, z_king_grid, ...
    source_cfg, medium_cfg, calc_cfg_ultra);

%% ---- 2D transform SPL norm (without local effect) ----
p_abs_2d = abs(pa_transform_full_wo);
p_ref_2d = max(p_abs_2d(:));
spl_transform_2d_norm_wo = 20*log10(p_abs_2d / max(p_ref_2d, eps) + eps);

out = struct();
out.axis_res = axis_res;
out.radial_res = radial_res;

out.rho_audio_grid = xh(:);
out.z_audio_grid   = z_audio(:);
out.pa_transform_full_wo = pa_transform_full_wo;
out.spl_transform_2d_norm_wo = spl_transform_2d_norm_wo;
end

%% ============================================================
% local-effect correction
%
% Revised version:
%   make_source_velocity is not required to support m1/m2 directly.
%   The f1 field is computed separately with source.m = m1, then f1 is extracted.
%   The f2 field is computed separately with source.m = m2, then f2 is extracted.
%% ============================================================
function line_res = apply_local_effect_correction_to_line_bothstyle( ...
    line_res, mode, ...
    rho_pts_req, z_pts_req, ...
    rho_king_grid, z_king_grid, ...
    source_cfg, medium_cfg, calc_cfg_ultra)

rho_pts_req = rho_pts_req(:);
z_pts_req   = z_pts_req(:);

c    = medium_cfg.c0;
rho0 = medium_cfg.rho0;
f1   = source_cfg.f1;
f2   = source_cfg.f2;

%% ---- Build obs_grid ----
switch lower(mode)
    case 'axial'
        rho_fixed_req = rho_pts_req(1);
        [~, ir0] = min(abs(rho_king_grid - rho_fixed_req));
        rho_fixed_use = rho_king_grid(ir0);

        z_use = zeros(size(z_pts_req));
        for ii = 1:numel(z_pts_req)
            [~, iz] = min(abs(z_king_grid - z_pts_req(ii)));
            z_use(ii) = z_king_grid(iz);
        end

        obs_grid = struct();
        obs_grid.dim.x = rho_fixed_use;
        obs_grid.dim.y = 0;
        obs_grid.dim.z = z_use(:).';
        obs_grid.dim.block_size = 200000;

        compare_info = struct();
        compare_info.mode = 'axial';
        compare_info.rho_use = rho_fixed_use;
        compare_info.z_use = z_use(:);

    case 'radial'
        z_fixed_req = z_pts_req(1);
        [~, iz0] = min(abs(z_king_grid - z_fixed_req));
        z_fixed_use = z_king_grid(iz0);

        rho_use = zeros(size(rho_pts_req));
        for ii = 1:numel(rho_pts_req)
            [~, ir] = min(abs(rho_king_grid - rho_pts_req(ii)));
            rho_use(ii) = rho_king_grid(ir);
        end

        obs_grid = struct();
        obs_grid.dim.x = rho_use(:).';
        obs_grid.dim.y = 0;
        obs_grid.dim.z = z_fixed_use;
        obs_grid.dim.block_size = 200000;

        compare_info = struct();
        compare_info.mode = 'radial';
        compare_info.rho_use = rho_use(:);
        compare_info.z_use = z_fixed_use;

    case 'oblique'
        rho_use = zeros(size(rho_pts_req));
        z_use   = zeros(size(z_pts_req));

        for ii = 1:numel(rho_pts_req)
            [~, ir] = min(abs(rho_king_grid - rho_pts_req(ii)));
            [~, iz] = min(abs(z_king_grid   - z_pts_req(ii)));
            rho_use(ii) = rho_king_grid(ir);
            z_use(ii)   = z_king_grid(iz);
        end

        obs_grid = struct();
        obs_grid.dim.x = rho_use(:).';
        obs_grid.dim.y = 0;
        obs_grid.dim.z = z_use(:).';
        obs_grid.dim.block_size = 200000;

        compare_info = struct();
        compare_info.mode = 'oblique';
        compare_info.rho_use = rho_use(:);
        compare_info.z_use   = z_use(:);

    otherwise
        error('Unknown mode in apply_local_effect_correction_to_line_bothstyle.');
end

%% ---- Compute ultrasound p/v for f1/m1 and f2/m2 separately ----
source_cfg_f1 = source_cfg;
source_cfg_f2 = source_cfg;

source_cfg_f1.m  = source_cfg.m1;
source_cfg_f1.m1 = source_cfg.m1;
source_cfg_f1.m2 = source_cfg.m1;

source_cfg_f2.m  = source_cfg.m2;
source_cfg_f2.m1 = source_cfg.m2;
source_cfg_f2.m2 = source_cfg.m2;

fprintf('    [local correction] computing p/v for f1 with m = %d ...\n', source_cfg.m1);
res_p_f1 = calc_ultrasound_field(source_cfg_f1, medium_cfg, calc_cfg_ultra, obs_grid, 'both');
res_v_f1 = calc_ultrasound_velocity_field(source_cfg_f1, medium_cfg, calc_cfg_ultra, obs_grid, 'both');

fprintf('    [local correction] computing p/v for f2 with m = %d ...\n', source_cfg.m2);
res_p_f2 = calc_ultrasound_field(source_cfg_f2, medium_cfg, calc_cfg_ultra, obs_grid, 'both');
res_v_f2 = calc_ultrasound_velocity_field(source_cfg_f2, medium_cfg, calc_cfg_ultra, obs_grid, 'both');

%% ---- Extract f1 and f2 fields separately ----
switch lower(compare_info.mode)
    case 'axial'
        [p1_k, v1r_k, v1phi_k, v1z_k, ...
         p1_d, v1r_d, v1phi_d, v1z_d, ...
         rho_actual_localcorr, z_actual_localcorr] = extract_ultra_line_fields_onefreq_axial_points( ...
            res_p_f1, res_v_f1, 'f1', compare_info.rho_use, compare_info.z_use);

        [p2_k, v2r_k, v2phi_k, v2z_k, ...
         p2_d, v2r_d, v2phi_d, v2z_d, ...
         ~, ~] = extract_ultra_line_fields_onefreq_axial_points( ...
            res_p_f2, res_v_f2, 'f2', compare_info.rho_use, compare_info.z_use);

    case 'radial'
        [p1_k, v1r_k, v1phi_k, v1z_k, ...
         p1_d, v1r_d, v1phi_d, v1z_d, ...
         rho_actual_localcorr, z_actual_localcorr] = extract_ultra_line_fields_onefreq_radial_points( ...
            res_p_f1, res_v_f1, 'f1', compare_info.rho_use, compare_info.z_use);

        [p2_k, v2r_k, v2phi_k, v2z_k, ...
         p2_d, v2r_d, v2phi_d, v2z_d, ...
         ~, ~] = extract_ultra_line_fields_onefreq_radial_points( ...
            res_p_f2, res_v_f2, 'f2', compare_info.rho_use, compare_info.z_use);

    case 'oblique'
        [p1_k, v1r_k, v1phi_k, v1z_k, ...
         p1_d, v1r_d, v1phi_d, v1z_d, ...
         rho_actual_localcorr, z_actual_localcorr] = extract_ultra_line_fields_onefreq_oblique_points( ...
            res_p_f1, res_v_f1, 'f1', compare_info.rho_use, compare_info.z_use);

        [p2_k, v2r_k, v2phi_k, v2z_k, ...
         p2_d, v2r_d, v2phi_d, v2z_d, ...
         ~, ~] = extract_ultra_line_fields_onefreq_oblique_points( ...
            res_p_f2, res_v_f2, 'f2', compare_info.rho_use, compare_info.z_use);

    otherwise
        error('Unknown compare mode.');
end

%% ---- local-effect correction ----
p_loc_transform = compute_local_effect_term( ...
    p1_k, p2_k, v1r_k, v1phi_k, v1z_k, v2r_k, v2phi_k, v2z_k, ...
    rho0, c, f1, f2);

p_loc_direct = compute_local_effect_term( ...
    p1_d, p2_d, v1r_d, v1phi_d, v1z_d, v2r_d, v2phi_d, v2z_d, ...
    rho0, c, f1, f2);

line_res.rho_actual_localcorr = rho_actual_localcorr(:);
line_res.z_actual_localcorr   = z_actual_localcorr(:);

line_res.pa_transform_wo = line_res.pa_transform;
line_res.pa_direct_wo    = line_res.pa_direct;

line_res.spl_transform_wo = line_res.spl_transform;
line_res.spl_direct_wo    = line_res.spl_direct;

line_res.phase_transform_deg_wo = line_res.phase_transform_deg;
line_res.phase_direct_deg_wo    = line_res.phase_direct_deg;

line_res.p_loc_transform = p_loc_transform(:);
line_res.p_loc_direct    = p_loc_direct(:);

line_res.pa_transform = line_res.pa_transform_wo - line_res.p_loc_transform;
line_res.pa_direct    = line_res.pa_direct_wo    - line_res.p_loc_direct;

for ii = 1:numel(line_res.pa_transform)
    line_res.spl_transform(ii,1) = local_pressure_to_spl(abs(line_res.pa_transform(ii)), medium_pref_from_cfg(medium_cfg));
    line_res.spl_direct(ii,1)    = local_pressure_to_spl(abs(line_res.pa_direct(ii)),    medium_pref_from_cfg(medium_cfg));
end

line_res.phase_transform_deg = rad2deg(unwrap(angle(line_res.pa_transform)));
line_res.phase_direct_deg    = rad2deg(unwrap(angle(line_res.pa_direct)));

line_res.rel_err_complex = abs(line_res.pa_transform - line_res.pa_direct) ./ max(abs(line_res.pa_direct), 1e-16);
line_res.rel_err_amp     = abs(abs(line_res.pa_transform) - abs(line_res.pa_direct)) ./ max(abs(line_res.pa_direct), 1e-16);

line_res.ultra_transform = struct( ...
    'p1', p1_k, 'p2', p2_k, ...
    'vr1', v1r_k, 'vphi1', v1phi_k, 'vz1', v1z_k, ...
    'vr2', v2r_k, 'vphi2', v2phi_k, 'vz2', v2z_k);

line_res.ultra_direct = struct( ...
    'p1', p1_d, 'p2', p2_d, ...
    'vr1', v1r_d, 'vphi1', v1phi_d, 'vz1', v1z_d, ...
    'vr2', v2r_d, 'vphi2', v2phi_d, 'vz2', v2z_d);
end

%% ============================================================
% Extract one line from the full transform field and direct point values, without local effect
%% ============================================================
function line_res = extract_line_from_transform_and_direct( ...
    rho_grid_audio, z_grid_audio, pa_transform_full, ...
    rho_req, z_req, ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d, ma, ka, rho0, wa, pref, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror)

nPt = numel(rho_req);

rho_actual = zeros(nPt,1);
z_actual   = zeros(nPt,1);
pa_transform = complex(zeros(nPt,1));

for ii = 1:nPt
    [~, id_rho] = min(abs(rho_grid_audio - rho_req(ii)));
    [~, id_z]   = min(abs(z_grid_audio   - z_req(ii)));

    rho_actual(ii) = rho_grid_audio(id_rho);
    z_actual(ii)   = z_grid_audio(id_z);
    pa_transform(ii) = pa_transform_full(id_rho, id_z);
end

pa_direct = compute_audio_direct_points_from_qrzphi_parallel( ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d.dr, direct2d.dz, direct2d.Nphi, ma, ...
    ka, rho0, wa, ...
    rho_actual, z_actual, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror);

spl_transform = zeros(nPt,1);
spl_direct    = zeros(nPt,1);

for ii = 1:nPt
    spl_transform(ii) = local_pressure_to_spl(abs(pa_transform(ii)), pref);
    spl_direct(ii)    = local_pressure_to_spl(abs(pa_direct(ii)), pref);
end

phase_transform_deg = rad2deg(unwrap(angle(pa_transform)));
phase_direct_deg    = rad2deg(unwrap(angle(pa_direct)));

rel_err_complex = abs(pa_transform - pa_direct) ./ max(abs(pa_direct), 1e-16);
rel_err_amp     = abs(abs(pa_transform) - abs(pa_direct)) ./ max(abs(pa_direct), 1e-16);

line_res = struct();
line_res.rho_req = rho_req(:);
line_res.z_req   = z_req(:);

line_res.rho_actual = rho_actual;
line_res.z_actual   = z_actual;

line_res.pa_transform = pa_transform;
line_res.pa_direct    = pa_direct;

line_res.spl_transform = spl_transform;
line_res.spl_direct    = spl_direct;

line_res.phase_transform_deg = phase_transform_deg;
line_res.phase_direct_deg    = phase_direct_deg;

line_res.rel_err_complex = rel_err_complex;
line_res.rel_err_amp     = rel_err_amp;
end

%% ============================================================
% Extract single-frequency ultrasound p/v at axial target points
%% ============================================================
function [p_k, vr_k, vphi_k, vz_k, ...
          p_d, vr_d, vphi_d, vz_d, ...
          rho_actual, z_actual] = extract_ultra_line_fields_onefreq_axial_points( ...
          res_p, res_v, freq_tag, rho_ax, z_use)

z_use = z_use(:);
nPt = numel(z_use);

zK = res_p.king.z(:);
rhoK = res_p.king.rho(:);
[~, ixK] = min(abs(rhoK - rho_ax));

rho_actual = rhoK(ixK) * ones(nPt,1);
z_actual   = zeros(nPt,1);

[p_names, vr_names, vphi_names, vz_names] = local_get_freq_field_names(freq_tag);

pK_full    = get_field_with_aliases(res_p.king, p_names);
vrK_full   = get_field_with_aliases(res_v.king, vr_names);
vphiK_full = get_field_with_aliases(res_v.king, vphi_names);
vzK_full   = get_field_with_aliases(res_v.king, vz_names);

p_k    = complex(zeros(nPt,1));
vr_k   = complex(zeros(nPt,1));
vphi_k = complex(zeros(nPt,1));
vz_k   = complex(zeros(nPt,1));

for ii = 1:nPt
    [~, izK] = min(abs(zK - z_use(ii)));
    z_actual(ii) = zK(izK);

    p_k(ii)    = pK_full(ixK, izK);
    vr_k(ii)   = vrK_full(ixK, izK);
    vphi_k(ii) = vphiK_full(ixK, izK);
    vz_k(ii)   = vzK_full(ixK, izK);
end

if strcmpi(res_p.calc.dim.method, 'rayleigh')
    zD = res_p.dim.z(:);

    pD_full    = get_field_with_aliases(res_p.dim, p_names);
    vrD_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiD_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzD_full   = get_field_with_aliases(res_v.dim, vz_names);

    p_d    = complex(zeros(nPt,1));
    vr_d   = complex(zeros(nPt,1));
    vphi_d = complex(zeros(nPt,1));
    vz_d   = complex(zeros(nPt,1));

    for ii = 1:nPt
        [~, izD] = min(abs(zD - z_actual(ii)));

        p_d(ii)    = pD_full(1,1,izD);
        vr_d(ii)   = vrD_full(1,1,izD);
        vphi_d(ii) = vphiD_full(1,1,izD);
        vz_d(ii)   = vzD_full(1,1,izD);
    end
else
    xA = res_p.dim.x(:);
    yA = res_p.dim.y(:);
    zA = res_p.dim.z(:);

    [~, ix0] = min(abs(xA - rho_actual(1)));
    [~, iy0] = min(abs(yA - 0));

    pA_full    = get_field_with_aliases(res_p.dim, p_names);
    vrA_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiA_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzA_full   = get_field_with_aliases(res_v.dim, vz_names);

    pA    = squeeze(pA_full(iy0,ix0,:));
    vrA   = squeeze(vrA_full(iy0,ix0,:));
    vphiA = squeeze(vphiA_full(iy0,ix0,:));
    vzA   = squeeze(vzA_full(iy0,ix0,:));

    p_d    = interp1(zA, pA(:),    z_actual, 'linear', 0); p_d = p_d(:);
    vr_d   = interp1(zA, vrA(:),   z_actual, 'linear', 0); vr_d = vr_d(:);
    vphi_d = interp1(zA, vphiA(:), z_actual, 'linear', 0); vphi_d = vphi_d(:);
    vz_d   = interp1(zA, vzA(:),   z_actual, 'linear', 0); vz_d = vz_d(:);
end
end

%% ============================================================
% Extract single-frequency ultrasound p/v at radial target points
%% ============================================================
function [p_k, vr_k, vphi_k, vz_k, ...
          p_d, vr_d, vphi_d, vz_d, ...
          rho_actual, z_actual] = extract_ultra_line_fields_onefreq_radial_points( ...
          res_p, res_v, freq_tag, rho_use, z_use)

rho_use = rho_use(:);
nPt = numel(rho_use);

zK = res_p.king.z(:);
rhoK = res_p.king.rho(:);
[~, izK] = min(abs(zK - z_use));

rho_actual = zeros(nPt,1);
z_actual   = zK(izK) * ones(nPt,1);

[p_names, vr_names, vphi_names, vz_names] = local_get_freq_field_names(freq_tag);

pK_full    = get_field_with_aliases(res_p.king, p_names);
vrK_full   = get_field_with_aliases(res_v.king, vr_names);
vphiK_full = get_field_with_aliases(res_v.king, vphi_names);
vzK_full   = get_field_with_aliases(res_v.king, vz_names);

p_k    = complex(zeros(nPt,1));
vr_k   = complex(zeros(nPt,1));
vphi_k = complex(zeros(nPt,1));
vz_k   = complex(zeros(nPt,1));

for ii = 1:nPt
    [~, irK] = min(abs(rhoK - rho_use(ii)));
    rho_actual(ii) = rhoK(irK);

    p_k(ii)    = pK_full(irK, izK);
    vr_k(ii)   = vrK_full(irK, izK);
    vphi_k(ii) = vphiK_full(irK, izK);
    vz_k(ii)   = vzK_full(irK, izK);
end

if strcmpi(res_p.calc.dim.method, 'rayleigh')
    xD = res_p.dim.x(:);

    pD_full    = get_field_with_aliases(res_p.dim, p_names);
    vrD_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiD_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzD_full   = get_field_with_aliases(res_v.dim, vz_names);

    p_d    = complex(zeros(nPt,1));
    vr_d   = complex(zeros(nPt,1));
    vphi_d = complex(zeros(nPt,1));
    vz_d   = complex(zeros(nPt,1));

    for ii = 1:nPt
        [~, ixD] = min(abs(xD - rho_actual(ii)));

        p_d(ii)    = pD_full(1,ixD,1);
        vr_d(ii)   = vrD_full(1,ixD,1);
        vphi_d(ii) = vphiD_full(1,ixD,1);
        vz_d(ii)   = vzD_full(1,ixD,1);
    end
else
    xA = res_p.dim.x(:);
    yA = res_p.dim.y(:);
    zA = res_p.dim.z(:);

    [~, iy0] = min(abs(yA - 0));
    [~, iz0] = min(abs(zA - z_actual(1)));

    pA_full    = get_field_with_aliases(res_p.dim, p_names);
    vrA_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiA_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzA_full   = get_field_with_aliases(res_v.dim, vz_names);

    pA    = squeeze(pA_full(iy0,:,iz0)).';
    vrA   = squeeze(vrA_full(iy0,:,iz0)).';
    vphiA = squeeze(vphiA_full(iy0,:,iz0)).';
    vzA   = squeeze(vzA_full(iy0,:,iz0)).';

    p_d    = interp1(xA, pA,    rho_actual, 'linear', 0); p_d = p_d(:);
    vr_d   = interp1(xA, vrA,   rho_actual, 'linear', 0); vr_d = vr_d(:);
    vphi_d = interp1(xA, vphiA, rho_actual, 'linear', 0); vphi_d = vphi_d(:);
    vz_d   = interp1(xA, vzA,   rho_actual, 'linear', 0); vz_d = vz_d(:);
end
end

%% ============================================================
% Extract single-frequency ultrasound p/v at oblique target points
%% ============================================================
function [p_k, vr_k, vphi_k, vz_k, ...
          p_d, vr_d, vphi_d, vz_d, ...
          rho_actual, z_actual] = extract_ultra_line_fields_onefreq_oblique_points( ...
          res_p, res_v, freq_tag, rho_use, z_use)

rho_use = rho_use(:);
z_use   = z_use(:);
nPt = numel(rho_use);

rhoK = res_p.king.rho(:);
zK   = res_p.king.z(:);

[p_names, vr_names, vphi_names, vz_names] = local_get_freq_field_names(freq_tag);

pK_full    = get_field_with_aliases(res_p.king, p_names);
vrK_full   = get_field_with_aliases(res_v.king, vr_names);
vphiK_full = get_field_with_aliases(res_v.king, vphi_names);
vzK_full   = get_field_with_aliases(res_v.king, vz_names);

p_k    = complex(zeros(nPt,1));
vr_k   = complex(zeros(nPt,1));
vphi_k = complex(zeros(nPt,1));
vz_k   = complex(zeros(nPt,1));

rho_actual = zeros(nPt,1);
z_actual   = zeros(nPt,1);

for ii = 1:nPt
    [~, irK] = min(abs(rhoK - rho_use(ii)));
    [~, izK] = min(abs(zK   - z_use(ii)));

    rho_actual(ii) = rhoK(irK);
    z_actual(ii)   = zK(izK);

    p_k(ii)    = pK_full(irK, izK);
    vr_k(ii)   = vrK_full(irK, izK);
    vphi_k(ii) = vphiK_full(irK, izK);
    vz_k(ii)   = vzK_full(irK, izK);
end

if strcmpi(res_p.calc.dim.method, 'rayleigh')
    xD = res_p.dim.x(:);
    zD = res_p.dim.z(:);

    pD_full    = get_field_with_aliases(res_p.dim, p_names);
    vrD_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiD_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzD_full   = get_field_with_aliases(res_v.dim, vz_names);

    p_d    = complex(zeros(nPt,1));
    vr_d   = complex(zeros(nPt,1));
    vphi_d = complex(zeros(nPt,1));
    vz_d   = complex(zeros(nPt,1));

    for ii = 1:nPt
        [~, ixD] = min(abs(xD - rho_actual(ii)));
        [~, izD] = min(abs(zD - z_actual(ii)));

        p_d(ii)    = pD_full(1,ixD,izD);
        vr_d(ii)   = vrD_full(1,ixD,izD);
        vphi_d(ii) = vphiD_full(1,ixD,izD);
        vz_d(ii)   = vzD_full(1,ixD,izD);
    end
else
    xA = res_p.dim.x(:);
    yA = res_p.dim.y(:);
    zA = res_p.dim.z(:);

    [~, iy0] = min(abs(yA - 0));

    pA_full    = get_field_with_aliases(res_p.dim, p_names);
    vrA_full   = get_field_with_aliases(res_v.dim, vr_names);
    vphiA_full = get_field_with_aliases(res_v.dim, vphi_names);
    vzA_full   = get_field_with_aliases(res_v.dim, vz_names);

    p_d    = complex(zeros(nPt,1));
    vr_d   = complex(zeros(nPt,1));
    vphi_d = complex(zeros(nPt,1));
    vz_d   = complex(zeros(nPt,1));

    for ii = 1:nPt
        [~, ixA] = min(abs(xA - rho_actual(ii)));
        [~, izA] = min(abs(zA - z_actual(ii)));

        p_d(ii)    = pA_full(iy0,ixA,izA);
        vr_d(ii)   = vrA_full(iy0,ixA,izA);
        vphi_d(ii) = vphiA_full(iy0,ixA,izA);
        vz_d(ii)   = vzA_full(iy0,ixA,izA);
    end
end
end

%% ============================================================
% Return field aliases for f1/f2
%% ============================================================
function [p_names, vr_names, vphi_names, vz_names] = local_get_freq_field_names(freq_tag)

switch lower(freq_tag)
    case 'f1'
        p_names    = {'p_f1','p1','pressure_f1'};
        vr_names   = {'v_rho_f1','vr_f1','v_r_f1'};
        vphi_names = {'v_phi_f1','vphi_f1'};
        vz_names   = {'v_z_f1','vz_f1','uz_f1'};

    case 'f2'
        p_names    = {'p_f2','p2','pressure_f2'};
        vr_names   = {'v_rho_f2','vr_f2','v_r_f2'};
        vphi_names = {'v_phi_f2','vphi_f2'};
        vz_names   = {'v_z_f2','vz_f2','uz_f2'};

    otherwise
        error('Unknown freq_tag: %s', freq_tag);
end
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
% Zeroth-order Hankel transform of spatial-domain g
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
    q_block_size, point_batch_size, ...
    use_parallel, use_z_mirror)

if nargin < 16
    point_batch_size = [];
end


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
% Direct audio calculation at one point, blockwise over z / phi
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

calc_cfg.king = struct();
calc_cfg.king.gspec_method = 'transform';
calc_cfg.king.eps_kzz = 1e-3;
calc_cfg.king.eps_phase = calc_cfg.king.eps_kzz;
calc_cfg.king.kz_min = 1e-12;
calc_cfg.king.band_refine.enable = false;

calc_cfg.asm = struct();
calc_cfg.asm.pad_factor = 16;
calc_cfg.asm.kzz_eps = 1e-12;

calc_cfg.internal = struct();
end

%% ============================================================
% Build calc cfg for ultrasound pressure/velocity calls with both methods
%% ============================================================
function calc_cfg = build_calc_cfg_for_ultra_both( ...
    N_FHT, rho_max, zu_max, za_max, delta, use_parallel, medium_cfg, source_cfg)

if nargin < 7
    medium_cfg = [];
end
if nargin < 8
    source_cfg = [];
end


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
calc_cfg.dim.use_parallel = use_parallel;
calc_cfg.dim.num_workers = 20;
calc_cfg.dim.use_freq = 'f2';
calc_cfg.dim.dis_coe = 32;
calc_cfg.dim.margin = 1;
calc_cfg.dim.src_discretization = 'polar';

calc_cfg.dim.grid_mode = 'uniform';
calc_cfg.dim.ds_rho = 32;
calc_cfg.dim.ds_rho_src = 16;
calc_cfg.dim.ds_rho_obs = 32;
calc_cfg.dim.uniform_dx = delta / 2;
calc_cfg.dim.uniform_dz = delta / 1;

calc_cfg.king = struct();
calc_cfg.king.gspec_method = 'transform';
calc_cfg.king.eps_kzz = 1e-3;
calc_cfg.king.eps_phase = calc_cfg.king.eps_kzz;
calc_cfg.king.kz_min = 1e-12;
calc_cfg.king.band_refine.enable = false;

calc_cfg.asm = struct();
calc_cfg.asm.pad_factor = 16;
calc_cfg.asm.kzz_eps = 1e-12;

calc_cfg.internal = struct();
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
% Return pref
%% ============================================================
function pref = medium_pref_from_cfg(medium_cfg)
pref = medium_cfg.pref;
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
% Build cache metadata for direct ultrasound caches
%% ============================================================
function cache_meta = local_build_direct_cache_meta( ...
    source_cfg, medium_cfg, ...
    fu, fa, f1, f2, ...
    m1, m2, ma, ...
    rho_max, zu_max, ...
    N_FHT, delta, dis_coe, Nphi)

cache_meta = struct();
cache_meta.cache_version = 2;

cache_meta.source_profile = char(string(source_cfg.profile));
cache_meta.a = source_cfg.a;
cache_meta.v0 = source_cfg.v0;
cache_meta.v_ratio = source_cfg.v_ratio;
cache_meta.F = source_cfg.F;

cache_meta.c = medium_cfg.c0;
cache_meta.rho0 = medium_cfg.rho0;
cache_meta.beta = medium_cfg.beta;
cache_meta.pref = medium_cfg.pref;

cache_meta.fu = fu;
cache_meta.fa = fa;
cache_meta.f1 = f1;
cache_meta.f2 = f2;

cache_meta.rho_max = rho_max;
cache_meta.zu_max = zu_max;
cache_meta.N_FHT = N_FHT;
cache_meta.delta = delta;
cache_meta.dis_coe = dis_coe;
cache_meta.Nphi = Nphi;

cache_meta.m1 = m1;
cache_meta.m2 = m2;
cache_meta.ma = ma;

cache_meta.time_tag = datestr(now, 'yyyy-mm-dd HH:MM:SS');
end

%% ============================================================
% Validate cache metadata and return the v0 scaling factor
%% ============================================================
function alpha_v0 = local_validate_cache_meta_and_get_v0_scale( ...
    cache_meta, ...
    a_now, v0_now, c_now, rho0_now, beta_now, f1_now, f2_now, fa_now, ...
    profile_now, ...
    m1_now, m2_now, rho_max_now, zu_max_now, ...
    N_FHT_now, delta_now, dis_coe_now, Nphi_now)

if nargin < 1 || isempty(cache_meta) || ~isstruct(cache_meta)
    error(['Cache metadata is missing or invalid. ' ...
           'Old caches without cache_meta cannot be safely reused. ' ...
           'Please delete the cache and regenerate it.']);
end

required_fields = { ...
    'v0', 'a', 'c', 'rho0', 'beta', ...
    'f1', 'f2', 'fa', ...
    'm1', 'm2', ...
    'rho_max', 'zu_max', ...
    'N_FHT', 'delta', 'dis_coe', 'Nphi'};

for ii = 1:numel(required_fields)
    fn = required_fields{ii};
    if ~isfield(cache_meta, fn)
        error(['Cache metadata field "%s" is missing. ' ...
               'This is likely an old cache. Please delete the cache and regenerate it.'], fn);
    end
end

local_assert_cache_value_close(cache_meta.a,       a_now,       'a');
local_assert_cache_value_close(cache_meta.c,       c_now,       'c');
local_assert_cache_value_close(cache_meta.rho0,    rho0_now,    'rho0');
local_assert_cache_value_close(cache_meta.beta,    beta_now,    'beta');
local_assert_cache_value_close(cache_meta.f1,      f1_now,      'f1');
local_assert_cache_value_close(cache_meta.f2,      f2_now,      'f2');
local_assert_cache_value_close(cache_meta.fa,      fa_now,      'fa');
local_assert_cache_value_close(cache_meta.rho_max, rho_max_now, 'rho_max');
local_assert_cache_value_close(cache_meta.zu_max,  zu_max_now,  'zu_max');
local_assert_cache_value_close(cache_meta.N_FHT,   N_FHT_now,   'N_FHT');
local_assert_cache_value_close(cache_meta.delta,   delta_now,   'delta');
local_assert_cache_value_close(cache_meta.dis_coe, dis_coe_now, 'dis_coe');
local_assert_cache_value_close(cache_meta.Nphi,    Nphi_now,    'Nphi');

if isfield(cache_meta, 'source_profile')
    profile_cache = lower(strtrim(char(string(cache_meta.source_profile))));
    profile_now   = lower(strtrim(char(string(profile_now))));
    if ~strcmp(profile_cache, profile_now)
        error('Cache source.profile mismatch: cache = %s, current = %s.', ...
            profile_cache, profile_now);
    end
else
    fprintf('Cache metadata has no source_profile field. Current profile is assumed for this legacy cache.\n');
end

if cache_meta.m1 ~= m1_now || cache_meta.m2 ~= m2_now
    error('Cache mode mismatch: cache m1=%d, m2=%d; current m1=%d, m2=%d.', ...
        cache_meta.m1, cache_meta.m2, m1_now, m2_now);
end

if cache_meta.v0 == 0
    error('Cache v0 is zero, so v0 scaling cannot be applied.');
end

alpha_v0 = v0_now / cache_meta.v0;
end

function local_assert_cache_value_close(x_cache, x_now, name)
tol_abs = 1e-12;
tol_rel = 1e-10;

x_cache = double(x_cache);
x_now   = double(x_now);

err_abs = abs(x_cache - x_now);
err_rel = err_abs / max(abs(x_now), tol_abs);

if err_abs > tol_abs && err_rel > tol_rel
    error('Cache parameter mismatch for %s: cache = %.16g, current = %.16g.', ...
        name, x_cache, x_now);
end
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
    for ii = 1:numel(fn)
        local_dump_any(fid, sprintf('%s.%s', name, fn{ii}), val.(fn{ii}), indent + 2);
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
            fprintf(fid, '%s%s = %.16g %+ .16gi\n', sp, name, real(val), imag(val));
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
            for kk = 1:nshow
                fprintf(fid, ' %.8g%+.8gi', real(vec(kk)), imag(vec(kk)));
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
    for ii = 1:nshow
        local_dump_any(fid, sprintf('%s{%d}', name, ii), val{ii}, indent + 2);
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
% Get a field from a structure using aliases
%% ============================================================
function A = get_field_with_aliases(S, names)
for kk = 1:numel(names)
    if isfield(S, names{kk})
        A = S.(names{kk});
        return;
    end
end
error('Cannot find any of the fields: %s', strjoin(names, ', '));
end

%% ============================================================
% Compute local-effect term
%% ============================================================
function p_loc = compute_local_effect_term( ...
    p1, p2, v1r, v1phi, v1z, v2r, v2phi, v2z, ...
    rho0, c, f1, f2)

w1 = 2*pi*f1;
w2 = 2*pi*f2;

coef = (w1/w2 + w2/w1 - 1);

vv = conj(v1r).*v2r + conj(v1phi).*v2phi + conj(v1z).*v2z;
pp = conj(p1).*p2;

p_loc = rho0/2 * vv - coef * pp ./ (2*rho0*c^2);
end

%% ============================================================
% Find a reusable cache file for p1 or p2
%
% Revised version:
%   Only pure-mode caches are allowed:
%       m1_here == target_m && m2_here == target_m
%
% Thus the mixed case:
%       m1=0,m2=3
% enforces:
%       p1 <- cache with m1 = m2 = 0
%       p2 <- cache with m1 = m2 = 3
%% ============================================================
function cache_info = find_reusable_single_field_cache( ...
    cache_root, which_field, target_m, ...
    rho_max, zu_max, N_FHT, delta, dis_coe, Nphi)

if nargin < 2
    which_field = '';
end


cache_info = struct();
cache_info.found = false;
cache_info.cache_dir = '';
cache_info.cache_mat = '';
cache_info.m1 = [];
cache_info.m2 = [];

if ~exist(cache_root, 'dir')
    return;
end

d = dir(cache_root);
d = d([d.isdir]);

s_rho   = local_num_to_tag(rho_max);
s_zu    = local_num_to_tag(zu_max);
s_delta = local_num_to_tag(delta);

must_have_1 = sprintf('rm%s_zu%s', s_rho, s_zu);
must_have_2 = sprintf('nf%d', N_FHT);
must_have_3 = sprintf('dt%s', s_delta);
must_have_4 = sprintf('dc%d_np%d', dis_coe, Nphi);

for kk = 1:numel(d)
    name = d(kk).name;
    if strcmp(name,'.') || strcmp(name,'..')
        continue;
    end

    if ~contains(name, must_have_1), continue; end
    if ~contains(name, must_have_2), continue; end
    if ~contains(name, must_have_3), continue; end
    if ~contains(name, must_have_4), continue; end

    tok = regexp(name, '^m1(-?\d+)m2(-?\d+)', 'tokens', 'once');
    if isempty(tok)
        continue;
    end

    m1_here = str2double(tok{1});
    m2_here = str2double(tok{2});

    if ~(m1_here == target_m && m2_here == target_m)
        continue;
    end

    cache_mat = fullfile(cache_root, name, 'pu_direct_rz_cache.mat');
    if ~exist(cache_mat, 'file')
        continue;
    end

    cache_info.found = true;
    cache_info.cache_dir = fullfile(cache_root, name);
    cache_info.cache_mat = cache_mat;
    cache_info.m1 = m1_here;
    cache_info.m2 = m2_here;
    return;
end
end


%% ============================================================
% Save one figure as PNG, FIG, and PDF
%% ============================================================
function save_one(fig, save_dir, base_name)

png_path = fullfile(save_dir, [base_name '.png']);
fig_path = fullfile(save_dir, [base_name '.fig']);
pdf_path = fullfile(save_dir, [base_name '.pdf']);

set(fig, 'Color', 'w');
set(fig, 'InvertHardcopy', 'off');

drawnow;

try
    savefig(fig, fig_path);
catch ME
    warning('Failed to save FIG: %s\n%s', fig_path, ME.message);
end

try
    exportgraphics(fig, png_path, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');
catch
    try
        print(fig, png_path, '-dpng', '-r300');
    catch ME
        warning('Failed to save PNG: %s\n%s', png_path, ME.message);
    end
end

try
    exportgraphics(fig, pdf_path, ...
        'ContentType', 'vector', ...
        'BackgroundColor', 'white');
catch
    try
        set(fig, 'Renderer', 'painters');
        print(fig, pdf_path, '-dpdf', '-painters', '-r300');
    catch ME
        warning('Failed to save PDF: %s\n%s', pdf_path, ME.message);
    end
end

end

%% ============================================================
% Set clean axis style with box and no grid
%% ============================================================
function local_set_clean_axis_style(ax)

ax.Box = 'on';
ax.Layer = 'top';

ax.FontSize = 18;
ax.LineWidth = 1.2;

ax.TickDir = 'in';
ax.TickLength = [0.015 0.015];
ax.XMinorTick = 'on';
ax.YMinorTick = 'off';

ax.XColor = [0 0 0];
ax.YColor = [0 0 0];

grid(ax, 'off');

end

%% ============================================================
% Set sparse major ticks for logarithmic x-axis
%% ============================================================
function local_set_log_xticks_line(ax, x)

x = x(:);
x = x(isfinite(x) & x > 0);

if isempty(x)
    return;
end

xlim_now = xlim(ax);

pmin = floor(log10(xlim_now(1)));
pmax = ceil(log10(xlim_now(2)));

major_ticks = 10.^(pmin:pmax);
major_ticks = major_ticks(major_ticks >= xlim_now(1) & major_ticks <= xlim_now(2));

if isempty(major_ticks)
    major_ticks = unique([xlim_now(1); xlim_now(2)]);
elseif numel(major_ticks) == 1
    major_ticks = unique([xlim_now(1); major_ticks(:); xlim_now(2)]);
end

ax.XTick = major_ticks(:).';

end

%% ============================================================
% Set y-axis margin for SPL curves
%% ============================================================
function local_set_y_margin_general(ax, y)

y = y(:);
y = y(isfinite(y));

if isempty(y)
    return;
end

ymin = min(y);
ymax = max(y);

if ymin == ymax
    ymin = ymin - 1;
    ymax = ymax + 1;
end

yrange = ymax - ymin;
pad = 0.08 * yrange;

ylim(ax, [ymin - pad, ymax + pad]);

end

%% ============================================================
% Set y-axis margin from zero for SPL-difference curves
%% ============================================================
function local_set_y_margin_from_zero(ax, y)

y = y(:);
y = y(isfinite(y));

if isempty(y)
    return;
end

ymax = max(y);

if ymax <= 0
    ymax = 1;
end

pad = 0.08 * ymax;

ylim(ax, [0, ymax + pad]);

end
