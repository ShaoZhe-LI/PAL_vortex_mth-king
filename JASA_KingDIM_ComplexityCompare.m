%% ============================================================
% Benchmark of vortex PAL audio-field calculation
%   Proposed transform/King method vs DIM
%   Accuracy check + article-ready timing / memory / sample-count output
%
% Main changes relative to the original fixed-case comparison script:
%   1) DIM ultrasound field is recomputed on the full direct rz grid.
%      No existing data_pu cache is used for the reported DIM ultrasound time.
%   2) DIM audio direct calculation uses the newly recomputed q(rho,z)
%      and is evaluated at 10 benchmark observation points.
%   3) The script prints and saves article-ready data:
%        - source / virtual-source / field-point sample numbers
%        - wall-clock time
%        - memory estimates
%      for the ultrasound and audible steps of both methods.
%   4) Dense DIM mode forms full Green's-function matrices without blockwise summation.
%   5) Local-effect correction is still evaluated on the same 10 audio points.
%
% External dependencies:
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%   - make_source_velocity.m
%   - calc_ultrasound_field.m
%   - calc_ultrasound_velocity_field.m
%% ============================================================

clear; clc; close all;

%% ===================== Global switches =====================
save_figures      = true;
save_calc_results = true;
save_params_txt   = true;
show_figures      = true;

% IMPORTANT:
% This benchmark always recomputes the DIM ultrasound rz field used by
% the subsequent DIM audible integration. Existing data_pu caches are not
% used for the reported DIM ultrasound timing.
force_recompute_direct_ultrasound_cache = true;

% ===================== Dense DIM settings =====================
% This version removes blockwise DIM summation and forms the full dense
% Green's-function matrices in one shot. With the default grid below, this
% can require tens of TB for the ultrasonic DIM step and hundreds of GB to
% TB for the audible DIM step. Therefore, the script performs a memory
% pre-check before allocation.
%
% Set allow_extreme_allocation = true only if the estimated memory is
% available on your machine. Otherwise, reduce rho_max, zu_max, Nphi, or
% the number of observation points.
dense_dim_cfg = struct();
% Try one-shot dense DIM first. If the estimated peak memory exceeds
% max_estimated_memory_gb, automatically fall back to blockwise DIM.
dense_dim_cfg.allow_extreme_allocation = false;
dense_dim_cfg.max_estimated_memory_gb  = 512;   % GB, threshold for one-shot allocation

% Fallback settings. The block size is chosen automatically so that the
% estimated temporary Green-function arrays per active block are close to
% this value. If parallel workers are used, the target is divided by the
% number of workers to keep the total active block memory near this value.
dense_dim_cfg.auto_fallback_to_blocks = true;
dense_dim_cfg.fallback_target_block_memory_gb = 512;   % GB
%dense_dim_cfg.fallback_target_block_memory_gb = 256;   % safer alternative

dense_dim_cfg.share_fallback_memory_across_workers = true;
dense_dim_cfg.block_memory_safety_factor = 0.85;
dense_dim_cfg.stop_if_exceeds_limit = false;

%% ===================== Benchmark points =====================
% DIM audible benchmark: 10 selected observation points.
% The DIM ultrasound step is recomputed on the full direct rz grid used
% to build the virtual-source distribution.

rho_bench_req = [0.0005; 0.005; 0.015; 0.030; 0.060; 0.100; 0.160; 0.240; 0.350; 0.480];
z_bench_req   = [0.0200; 0.050; 0.100; 0.200; 0.400; 0.700; 1.000; 1.800; 3.000; 4.500];
nBenchPts = numel(rho_bench_req);

all_case_summary = struct([]);

for i = 1:1

clearvars -except i save_figures save_calc_results save_params_txt show_figures ...
    force_recompute_direct_ultrasound_cache dense_dim_cfg ...
    rho_bench_req z_bench_req nBenchPts ...
    all_case_summary
close all;

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

% Benchmark case used for the article-ready cost table.
% Change these two values if another modal pair is required.
m1 = 0;
m2 = 3;
ma = m2 - m1;


%% ===================== Fixed FHT parameters =====================
N_FHT = 16384;
delta = 0.001;

%% ===================== Other calculation parameters =====================
rho_max = 0.5;
zu_max  = 15.0;
za_max  = 5 + delta;
green_R_min = 1e-12;

%% ===================== Parallel settings =====================
use_parallel = true;
num_workers  = 20;

%% ===================== Direct-reference settings =====================
direct_ref_cfg = struct();
direct_ref_cfg.src_block_size_ultra = 120000;
direct_ref_cfg.obs_block_size_rz    = 2048;
direct_ref_cfg.audio_q_block_size   = 200000;
direct_ref_cfg.audio_point_batch    = 16;

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

fprintf('\n============================================================\n');
fprintf('Case %d: m1 = %d, m2 = %d, ma = %d\n', i, m1, m2, ma);
fprintf('DIM ultrasound is recomputed on the full direct rz grid.\n');
fprintf('DIM audible benchmark points    = %d\n', nBenchPts);
fprintf('============================================================\n');

%% ===================== Save path =====================
time_tag  = datestr(now, 'mmdd_HHMMSS');
save_root_parent = 'result_benchmark_timing_memory';
save_root = fullfile(save_root_parent, sprintf('Benchmark10_%s_m1_%d_m2_%d', time_tag, m1, m2));

if (save_figures || save_calc_results || save_params_txt) && ~exist(save_root, 'dir')
    mkdir(save_root);
end

%% ===================== Parallel pool =====================
if use_parallel
    p = gcp('nocreate');
    if isempty(p)
        parpool('local', num_workers);
    elseif p.NumWorkers ~= num_workers
        delete(p);
        parpool('local', num_workers);
    end
end

%% ===================== Build source / medium / calc cfg =====================
source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2);
medium_cfg = build_medium_cfg(c, rho0, beta, pref);

%% ===================== Save parameters to txt =====================
if save_params_txt
    meta = struct();

    meta.script_name = mfilename;
    meta.case_index = i;
    meta.time_tag = time_tag;
    meta.save_root = save_root;

    meta.save_figures = save_figures;
    meta.save_calc_results = save_calc_results;
    meta.save_params_txt = save_params_txt;
    meta.show_figures = show_figures;
    meta.force_recompute_direct_ultrasound_cache = force_recompute_direct_ultrasound_cache;

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

    meta.use_parallel = use_parallel;
    meta.num_workers = num_workers;

    meta.direct_ref_cfg = direct_ref_cfg;
    meta.direct2d = direct2d;
    meta.audio_direct_use_z_mirror = audio_direct_use_z_mirror;
    meta.dense_dim_cfg = dense_dim_cfg;

    meta.rho_bench_req = rho_bench_req;
    meta.z_bench_req = z_bench_req;
    meta.nBenchPts = nBenchPts;

    meta.source_cfg = source_cfg;
    meta.medium_cfg = medium_cfg;

    meta.cache_policy = ['This benchmark recomputes the DIM ultrasound rz field. ' ...
        'Existing data_pu caches are not used for the reported DIM ultrasound timing.'];
    meta.local_effect_formula = ['p_corr = p_wo - [rho0/2 * conj(v1) dot v2 ' ...
        '- (w1/w2 + w2/w1 - 1) * conj(p1)p2 / (2*rho0*c^2)]'];

    local_write_all_params_txt(fullfile(save_root, 'all_parameters.txt'), meta);
end

%% ===================== Direct-reference source discretization =====================
fprintf('Preparing direct-reference source discretization for DIM...\n');
calc_cfg_ref = build_calc_cfg_for_direct_reference( ...
    direct_ref_cfg.N_FHT, rho_max, zu_max, za_max, direct_ref_cfg.delta, ...
    direct_ref_cfg.dis_coe, direct_ref_cfg.num_workers);

% For mixed modes, make_source_velocity is called separately for the two
% primary fields so that the two source topological charges are assigned
% correctly in the direct Rayleigh calculation.
[source_prep_ref, Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref] = ...
    prepare_direct_reference_sources_two_modes(source_cfg, medium_cfg, calc_cfg_ref);

N_src_eff = numel(Xs_ref);  % actual discrete source samples used by the Rayleigh integral

%% ===================== Direct-rz virtual-source grid =====================
rho_rz = build_uniform_grid_0_to_max(rho_max, direct2d.dr);
z_rz   = build_uniform_grid_0_to_max(zu_max,  direct2d.dz);

N_rho_v = numel(rho_rz);
N_z_v_pos = numel(z_rz);
if audio_direct_use_z_mirror
    N_z_v_audio = 2*N_z_v_pos - 1;
else
    N_z_v_audio = N_z_v_pos;
end
N_phi_v = direct2d.Nphi;
N_o_bench = nBenchPts;

%% ============================================================
% DIM ultrasound full-rz recomputation
% This q_direct_rz is then used by the DIM audible direct integration.
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Recomputing DIM ultrasound field on uniform (rho,z) grid ...\n');
fprintf('This step does NOT load existing data_pu caches.\n');
fprintf('Nsrc = %d, Nrho_v = %d, Nz_v = %d, field points = %d\n', ...
    N_src_eff, N_rho_v, N_z_v_pos, N_rho_v*N_z_v_pos);
fprintf('============================================================\n');

loaded_from_cache = false;
loaded_from_split_cache = false;
pu_cache_dir = '';
pu_cache_mat = '';

t_direct_rz = tic;
mem_before_dim_ultra = local_get_matlab_memory_used_gb();

[p1_direct_rz, p2_direct_rz, dim_ultra_dense_info] = compute_ultrasound_direct_rz_grid_parallel_from_ref_sources( ...
    rho_rz, z_rz, ...
    Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref, ...
    source_prep_ref.k1, source_prep_ref.k2, ...
    source_prep_ref.w1, source_prep_ref.w2, ...
    source_prep_ref.medium.rho0, ...
    direct_ref_cfg.src_block_size_ultra, ...
    direct_ref_cfg.obs_block_size_rz, ...
    use_parallel, dense_dim_cfg);

q_direct_rz = conj(p1_direct_rz) .* p2_direct_rz ...
    * beta*(2*pi*fa)/(1j*rho0^2*c^4);

time_direct_rz_once = toc(t_direct_rz);
mem_after_dim_ultra = local_get_matlab_memory_used_gb();

S_dim_ultra_full = whos('Xs_ref','Ys_ref','Zs_ref','q1_ref','q2_ref', ...
    'p1_direct_rz','p2_direct_rz','q_direct_rz','rho_rz','z_rz');
mem_dim_ultra_arrays_gb = local_whos_bytes_gb(S_dim_ultra_full);
mem_dim_ultra_delta_gb  = local_memory_delta_gb(mem_before_dim_ultra, mem_after_dim_ultra);
if exist('dim_ultra_dense_info','var') && isfield(dim_ultra_dense_info,'reported_peak_gb')
    mem_dim_ultra_dim_est_gb = dim_ultra_dense_info.reported_peak_gb;
elseif exist('dim_ultra_dense_info','var') && isfield(dim_ultra_dense_info,'estimated_peak_gb')
    mem_dim_ultra_dim_est_gb = dim_ultra_dense_info.estimated_peak_gb;
else
    mem_dim_ultra_dim_est_gb = NaN;
end
mem_dim_ultra_report_gb = max([mem_dim_ultra_arrays_gb, mem_dim_ultra_delta_gb, mem_dim_ultra_dim_est_gb], [], 'omitnan');

fprintf('DIM ultrasound full-rz recomputation finished. Time = %.3f s, memory = %.4f GB\n', ...
    time_direct_rz_once, mem_dim_ultra_report_gb);

if save_calc_results
    save(fullfile(save_root, 'dim_recomputed_q_direct_rz.mat'), ...
        'p1_direct_rz','p2_direct_rz','q_direct_rz','rho_rz','z_rz', ...
        'time_direct_rz_once','mem_dim_ultra_report_gb','-v7.3');
end

%% ============================================================
% Benchmark calculation: 10 points + timing/memory/article data
%% ============================================================
fprintf('\n============================================================\n');
fprintf('Benchmark calculation WITH local-effect correction: N_FHT = %d, delta = %.6f m\n', N_FHT, delta);
fprintf('DIM audio direct points = %d\n', nBenchPts);
fprintf('============================================================\n');

t_case = tic;

out_fixed = evaluate_benchmark_transform_and_direct_with_timing( ...
    N_FHT, delta, za_max, ...
    source_cfg, medium_cfg, ...
    c, rho0, beta, pref, f1, f2, fa, ...
    rho_max, zu_max, ...
    green_R_min, ma, ...
    rho_rz, z_rz, q_direct_rz, direct2d, ...
    direct_ref_cfg.audio_q_block_size, ...
    direct_ref_cfg.audio_point_batch, ...
    use_parallel, audio_direct_use_z_mirror, dense_dim_cfg, ...
    rho_bench_req, z_bench_req);

time_fixed_case = toc(t_case);
out_fixed.timing.total_case = time_fixed_case;
out_fixed.timing.dim_ultrasound_rz_cache_or_precompute = time_direct_rz_once;
out_fixed.timing.loaded_from_cache = loaded_from_cache;
out_fixed.timing.loaded_from_split_cache = loaded_from_split_cache;

bench_res = out_fixed.bench_res;
timing_stats = out_fixed.timing;
memory_stats = out_fixed.memory;

% Use the full direct-rz recomputation as the DIM ultrasound timing.
timing_stats.dim_ultrasound_full_rz = time_direct_rz_once;
memory_stats.dim_ultrasound_full_rz_gb = mem_dim_ultra_report_gb;

%% ===================== Article-ready sample/time/memory table =====================
article_table = build_article_cost_table( ...
    N_src_eff, N_rho_v, N_phi_v, N_z_v_pos, N_z_v_audio, N_o_bench, ...
    N_FHT, out_fixed.sample_numbers.N_z_v_prop, out_fixed.sample_numbers.N_z_o_prop, out_fixed.sample_numbers.N_z_c_prop, ...
    timing_stats, memory_stats);

fprintf('\n==================== Article-ready computational cost data ====================\n');
print_article_cost_table(article_table);
fprintf('=============================================================================\n');

if save_calc_results
    write_article_cost_table_csv(fullfile(save_root, 'article_cost_table.csv'), article_table);
    write_article_cost_table_txt(fullfile(save_root, 'article_cost_table.txt'), article_table);
end

%% ===================== Simple accuracy plot for 10 points =====================
fig1 = figure('Color','w', 'Position', [80 80 780 560]);
idx = (1:numel(bench_res.rho_actual)).';

subplot(2,1,1);
plot(idx, bench_res.spl_direct_wo, 'r-o', 'LineWidth', 1.1); hold on;
plot(idx, bench_res.spl_transform_wo, 'b--s', 'LineWidth', 1.1);
plot(idx, bench_res.spl_direct, 'm-o', 'LineWidth', 1.3);
plot(idx, bench_res.spl_transform, 'k--s', 'LineWidth', 1.3);
xlabel('Benchmark point index');
ylabel('SPL (dB)');
legend('DIM (wo local)', 'Proposed (wo local)', ...
       'DIM (w local)',  'Proposed (w local)', ...
       'Location', 'best');
title(sprintf('Benchmark points, m_1=%d, m_2=%d', m1, m2));
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'Box', 'on');

subplot(2,1,2);
err_wo = abs(bench_res.pa_transform_wo - bench_res.pa_direct_wo) ...
       ./ max(abs(bench_res.pa_direct_wo), 1e-16);
err_w  = bench_res.rel_err_complex;
plot(idx, log10(max(err_wo, 1e-16)), 'Color', [0 0.45 0.74], 'LineWidth', 1.1); hold on;
plot(idx, log10(max(err_w, 1e-16)),  'Color', [0.85 0.33 0.10], 'LineWidth', 1.1);
xlabel('Benchmark point index');
ylabel('log_{10} E_{rel}');
legend('wo local', 'w local', 'Location', 'best');
set(gca, 'FontSize', 11, 'LineWidth', 0.8, 'Box', 'on');

%% ===================== Save results =====================
if save_figures
    saveas(fig1, fullfile(save_root, 'Benchmark10_wo_and_w_LocalCorr.png'));
    savefig(fig1, fullfile(save_root, 'Benchmark10_wo_and_w_LocalCorr.fig'));
end

if ~show_figures
    close(fig1);
end

if save_calc_results
    save(fullfile(save_root, 'benchmark10_compare_results_localCorr_timing_memory.mat'), ...
        'out_fixed', 'bench_res', 'timing_stats', 'memory_stats', 'article_table', ...
        'N_FHT', 'delta', 'rho_max', 'zu_max', 'za_max', ...
        'rho_bench_req', 'z_bench_req', ...
        'time_direct_rz_once', 'time_fixed_case', ...
        'loaded_from_cache', 'loaded_from_split_cache', ...
        'pu_cache_dir', '-v7.3');
end

fprintf('\nAll done for case m1=%d, m2=%d.\n', m1, m2);
fprintf('Results folder: %s\n', save_root);

all_case_summary(i).m1 = m1;
all_case_summary(i).m2 = m2;
all_case_summary(i).ma = ma;
all_case_summary(i).save_root = save_root;
all_case_summary(i).article_table = article_table;
all_case_summary(i).time_dim_ultra_full_rz = time_direct_rz_once;
all_case_summary(i).mem_dim_ultra_full_rz_gb = mem_dim_ultra_report_gb;

end

%% Save global summary
summary_root = 'result_benchmark_timing_memory';
if ~exist(summary_root, 'dir')
    mkdir(summary_root);
end
save(fullfile(summary_root, ['all_case_summary_' datestr(now,'mmdd_HHMMSS') '.mat']), 'all_case_summary', '-v7.3');

%% ============================================================
% Benchmark function
%% ============================================================
function out = evaluate_benchmark_transform_and_direct_with_timing( ...
    N_FHT, delta, za_max, ...
    source_cfg, medium_cfg, ...
    c, rho0, beta, pref, f1, f2, fa, ...
    rho_max, zu_max, ...
    green_R_min, ma, ...
    rho_rz, z_rz, q_direct_rz, direct2d, ...
    audio_q_block_size, ...
    audio_point_batch, ...
    use_parallel, use_z_mirror, dense_dim_cfg, ...
    rho_bench_req, z_bench_req)

timing = struct();
memory = struct();
sample_numbers = struct();

mem_before = local_get_matlab_memory_used_gb();

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
Nza = numel(z_audio);

sample_numbers.N_z_v_prop = Nz;
sample_numbers.N_z_o_prop = Nza;
% Match the convolution length in compute_paW_from_Gr00.
Nz1 = numel([-fliplr(z(2:end)) z]);
sample_numbers.N_z_c_prop = Nz1 + Nza - 1;

%% ============================================================
% Proposed ultrasound step
%% ============================================================
t_prop_ultra = tic;

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

timing.proposed_ultrasound = toc(t_prop_ultra);
S_ultra = whos('G1_raw','G2_raw','G1','G2','p1_transform','p2_transform','q_transform','Vs1','Vs2');
memory.proposed_ultrasound_arrays_gb = local_whos_bytes_gb(S_ultra);

%% ============================================================
% Proposed audible step
%% ============================================================
t_prop_audio = tic;

Ga_raw = build_green_space_g_transform_raw( ...
    xh, z, ka, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);
Ga = (-4*pi*1j) * Ga_raw;

[pa_transform_full_wo, ~] = compute_paW_from_Gr00( ...
    q_transform, z, z_audio, Ga, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0, ma, ...
    delta, rho0, wa);

timing.proposed_audible = toc(t_prop_audio);
S_audio = whos('Ga_raw','Ga','pa_transform_full_wo');
memory.proposed_audible_arrays_gb = local_whos_bytes_gb(S_audio);
S_qrz = whos('q_direct_rz','rho_rz','z_rz');
memory.direct_qrz_arrays_gb = local_whos_bytes_gb(S_qrz);

%% ============================================================
% DIM audible step: only 10 benchmark points
%% ============================================================
t_dim_audio = tic;
bench_res = extract_line_from_transform_and_direct( ...
    xh, z_audio, pa_transform_full_wo, ...
    rho_bench_req, z_bench_req, ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d, ma, ka, rho0, wa, pref, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror, dense_dim_cfg);
timing.dim_audible_10pts = toc(t_dim_audio);
if isfield(bench_res, 'dim_audio_info') && isfield(bench_res.dim_audio_info, 'reported_peak_gb')
    memory.direct_qrz_arrays_gb = max(memory.direct_qrz_arrays_gb, bench_res.dim_audio_info.reported_peak_gb);
elseif isfield(bench_res, 'dim_audio_info') && isfield(bench_res.dim_audio_info, 'estimated_peak_gb')
    memory.direct_qrz_arrays_gb = max(memory.direct_qrz_arrays_gb, bench_res.dim_audio_info.estimated_peak_gb);
end

%% ============================================================
% Local-effect correction on the same 10 points
%% ============================================================
t_local = tic;
calc_cfg_ultra = build_calc_cfg_for_ultra_both( ...
    N_FHT, rho_max, zu_max, za_max, delta, use_parallel, medium_cfg, source_cfg);

[~, fht_tmp] = make_source_velocity(source_cfg, medium_cfg, calc_cfg_ultra);
rho_king_grid = fht_tmp.xh(:);
z_king_grid   = fht_tmp.z_ultra(:);

bench_res = apply_local_effect_correction_to_line_bothstyle( ...
    bench_res, 'oblique', ...
    rho_bench_req, z_bench_req, ...
    rho_king_grid, z_king_grid, ...
    source_cfg, medium_cfg, calc_cfg_ultra);
timing.local_correction_total = toc(t_local);

memory.matlab_mem_before_gb = mem_before;
memory.matlab_mem_after_gb  = local_get_matlab_memory_used_gb();

%% ============================================================
% Output
%% ============================================================
p_abs_2d = abs(pa_transform_full_wo);
p_ref_2d = max(p_abs_2d(:));
spl_transform_2d_norm_wo = 20*log10(p_abs_2d / max(p_ref_2d, eps) + eps);

out = struct();
out.bench_res = bench_res;
out.rho_audio_grid = xh(:);
out.z_audio_grid   = z_audio(:);
out.pa_transform_full_wo = pa_transform_full_wo;
out.spl_transform_2d_norm_wo = spl_transform_2d_norm_wo;
out.timing = timing;
out.memory = memory;
out.sample_numbers = sample_numbers;
end

%% ============================================================
% Local-effect correction
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

switch lower(mode)
    case {'oblique','points'}
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
        compare_info = struct('mode','oblique','rho_use',rho_use(:),'z_use',z_use(:));
    otherwise
        error('This benchmark version uses paired/oblique points only.');
end

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

[p1_k, v1r_k, v1phi_k, v1z_k, ...
 p1_d, v1r_d, v1phi_d, v1z_d, ...
 rho_actual_localcorr, z_actual_localcorr] = extract_ultra_line_fields_onefreq_oblique_points( ...
    res_p_f1, res_v_f1, 'f1', compare_info.rho_use, compare_info.z_use);

[p2_k, v2r_k, v2phi_k, v2z_k, ...
 p2_d, v2r_d, v2phi_d, v2z_d, ...
 ~, ~] = extract_ultra_line_fields_onefreq_oblique_points( ...
    res_p_f2, res_v_f2, 'f2', compare_info.rho_use, compare_info.z_use);

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
% Extract transform and DIM audio values at benchmark points, wo local
%% ============================================================
function line_res = extract_line_from_transform_and_direct( ...
    rho_grid_audio, z_grid_audio, pa_transform_full, ...
    rho_req, z_req, ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d, ma, ka, rho0, wa, pref, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror, dense_dim_cfg)

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

[pa_direct, dim_audio_info] = compute_audio_direct_points_from_qrzphi_parallel( ...
    rho_rz, z_rz, q_direct_rz, ...
    direct2d.dr, direct2d.dz, direct2d.Nphi, ma, ...
    ka, rho0, wa, ...
    rho_actual, z_actual, ...
    audio_q_block_size, audio_point_batch, ...
    use_parallel, use_z_mirror, dense_dim_cfg);

spl_transform = zeros(nPt,1);
spl_direct    = zeros(nPt,1);
for ii = 1:nPt
    spl_transform(ii) = local_pressure_to_spl(abs(pa_transform(ii)), pref);
    spl_direct(ii)    = local_pressure_to_spl(abs(pa_direct(ii)), pref);
end

line_res = struct();
line_res.rho_req = rho_req(:);
line_res.z_req   = z_req(:);
line_res.rho_actual = rho_actual;
line_res.z_actual   = z_actual;
line_res.pa_transform = pa_transform;
line_res.pa_direct    = pa_direct;
line_res.spl_transform = spl_transform;
line_res.spl_direct    = spl_direct;
line_res.phase_transform_deg = rad2deg(unwrap(angle(pa_transform)));
line_res.phase_direct_deg    = rad2deg(unwrap(angle(pa_direct)));
line_res.rel_err_complex = abs(pa_transform - pa_direct) ./ max(abs(pa_direct), 1e-16);
line_res.rel_err_amp     = abs(abs(pa_transform) - abs(pa_direct)) ./ max(abs(pa_direct), 1e-16);
line_res.dim_audio_info = dim_audio_info;
end

%% ============================================================
% Extract paired oblique target points from ultrasound p/v results
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
% Field aliases
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
% Direct ultrasound on rz grid, auto dense/blockwise version
%
% The function first estimates the memory required by a one-shot dense
% source--field Green-function matrix. If the estimate is below the limit,
% it uses the dense calculation. Otherwise, it automatically switches to
% blockwise summation and chooses the observation-block size so that the
% estimated temporary arrays per active block are close to the configured
% fallback memory target.
%% ============================================================
function [p1_rz, p2_rz, dense_info] = compute_ultrasound_direct_rz_grid_parallel_from_ref_sources( ...
    rho_grid, z_grid, ...
    Xs, Ys, Zs, q1, q2, ...
    k1, k2, w1, w2, rho0, ...
    src_block_size, obs_block_size, use_parallel, dense_dim_cfg)

if nargin < 18 || isempty(dense_dim_cfg)
    dense_dim_cfg = local_default_dense_dim_cfg();
end

[RHO, Z] = ndgrid(rho_grid, z_grid);
xo = RHO(:);
yo = zeros(size(xo));
zo = Z(:);

Nobs = numel(xo);
Nsrc = numel(Xs);

dense_info = local_estimate_dense_ultra_memory_gb(Nobs, Nsrc);
dense_info.mode_requested = 'dense_one_shot_first';

fprintf('    [ultra-direct-rz-auto] Nsrc = %d, Nobs = %d\n', Nsrc, Nobs);
fprintf('    [ultra-direct-rz-auto] dense one-shot estimated peak memory = %.3g GB\n', dense_info.estimated_peak_gb);

if local_should_use_dense_one_shot(dense_info.estimated_peak_gb, dense_dim_cfg)
    fprintf('    [ultra-direct-rz-auto] using dense one-shot Rayleigh integral.\n');
    dense_info.mode_used = 'dense_one_shot';
    dense_info.reported_peak_gb = dense_info.estimated_peak_gb;

    dx = xo - Xs(:).';
    dy = yo - Ys(:).';
    dz = zo - Zs(:).';
    R = sqrt(dx.^2 + dy.^2 + dz.^2);
    R(R < 1e-12) = 1e-12;

    G1 = exp(1j * k1 .* R) ./ (4*pi*R);
    G2 = exp(1j * k2 .* R) ./ (4*pi*R);

    p1_vec = (-2j * rho0 * w1) * (G1 * q1(:));
    p2_vec = (-2j * rho0 * w2) * (G2 * q2(:));

else
    fprintf('    [ultra-direct-rz-auto] dense estimate exceeds %.3g GB; switching to blockwise summation.\n', ...
        dense_dim_cfg.max_estimated_memory_gb);

    target_block_gb_total = local_get_fallback_target_gb(dense_dim_cfg);
    n_active_blocks = local_get_active_block_count(use_parallel, dense_dim_cfg);
    target_block_gb_each = target_block_gb_total / max(n_active_blocks, 1);
    target_block_gb_each = target_block_gb_each * local_get_block_safety_factor(dense_dim_cfg);

    obs_block_size_auto = local_choose_ultra_obs_block_size(Nsrc, target_block_gb_each, Nobs);
    block_info = local_estimate_dense_ultra_memory_gb(obs_block_size_auto, Nsrc);

    dense_info.mode_used = 'blockwise_fallback';
    dense_info.fallback_target_total_gb = target_block_gb_total;
    dense_info.fallback_target_each_active_block_gb = target_block_gb_each;
    dense_info.fallback_active_blocks = n_active_blocks;
    dense_info.obs_block_size_auto = obs_block_size_auto;
    dense_info.block_estimated_peak_gb = block_info.estimated_peak_gb;
    dense_info.reported_peak_gb = block_info.estimated_peak_gb * max(n_active_blocks, 1);

    fprintf('    [ultra-direct-rz-auto] fallback obs_block_size = %d, nBlocks = %d\n', ...
        obs_block_size_auto, ceil(Nobs/obs_block_size_auto));
    fprintf('    [ultra-direct-rz-auto] estimated peak per active block = %.3g GB; total active estimate = %.3g GB\n', ...
        dense_info.block_estimated_peak_gb, dense_info.reported_peak_gb);

    [p1_vec, p2_vec] = local_ultra_direct_blockwise_vectors( ...
        xo, yo, zo, Xs, Ys, Zs, q1, q2, k1, k2, w1, w2, rho0, ...
        min(src_block_size, Nsrc), obs_block_size_auto, use_parallel);
end

p1_rz = reshape(p1_vec, size(RHO));
p2_rz = reshape(p2_vec, size(RHO));
end

function [p1_vec, p2_vec] = local_ultra_direct_blockwise_vectors( ...
    xo, yo, zo, Xs, Ys, Zs, q1, q2, k1, k2, w1, w2, rho0, ...
    src_block_size, obs_block_size, use_parallel)

Nobs = numel(xo);
nBlocks = ceil(Nobs / obs_block_size);
p1_cell = cell(nBlocks,1);
p2_cell = cell(nBlocks,1);

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
        fprintf('        ultra block progress: %5d / %5d (%.1f%%), elapsed = %.1f s\n', ...
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

    function local_update_progress(~)
        progress_count = progress_count + 1;
        fprintf('        ultra block progress: %5d / %5d (%.1f%%), elapsed = %.1f s\n', ...
            progress_count, nBlocks, 100*progress_count/nBlocks, toc(t_progress));
    end
end

%% ============================================================
% Direct ultrasound, single observation block
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
% Proposed ultrasound pressure from spectral Green's function
%% ============================================================
function p_out = compute_ultrasound_pressure_from_G( ...
    G_spec, Vs, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use, rho0, c0, k_use)
F = G_spec .* Vs;
phi = -1j * m_FHT(F, N_FHT, Nz, NH, Nh, a_solve, x0, x1, k0, m_use);
p_out = 1j * rho0 * c0 * real(k_use) .* phi;
end

%% ============================================================
% Space-domain Green's function followed by zeroth-order Hankel transform
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
% Proposed audible field from Gr00
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
% DIM audio direct, auto dense/blockwise version
%
% The function first estimates the memory required by one-shot dense
% virtual-source--observation integration. If it exceeds the configured
% memory limit, it falls back to z/phi blockwise summation with automatically
% selected block sizes.
%% ============================================================
function [pa_pts, dense_info] = compute_audio_direct_points_from_qrzphi_parallel( ...
    rho_grid, z_grid_pos, q_rz_pos, ...
    dr, dz, Nphi, ma, ...
    ka, rho0, wa, ...
    rho_eval, z_eval, ...
    q_block_size, point_batch_size, ...
    use_parallel, use_z_mirror, dense_dim_cfg)

if nargin < 18 || isempty(dense_dim_cfg)
    dense_dim_cfg = local_default_dense_dim_cfg();
end

if use_z_mirror
    [z_grid_use, q_rz_use] = build_mirrored_qrz_about_z0(z_grid_pos, q_rz_pos);
    fprintf('    [audio-direct-rzphi-auto] z-mirror enabled.\n');
else
    z_grid_use = z_grid_pos(:);
    q_rz_use = q_rz_pos;
    fprintf('    [audio-direct-rzphi-auto] z-mirror disabled.\n');
end

rho_grid   = rho_grid(:);
z_grid_use = z_grid_use(:);
phi_all = (0:Nphi-1).' * (2*pi/Nphi);
dphi = 2*pi / Nphi;

Nr = numel(rho_grid);
Nz = numel(z_grid_use);
Npt = numel(rho_eval);
Nvirt = Nr * Nz * Nphi;

dense_info = local_estimate_dense_audio_memory_gb(Npt, Nvirt);
dense_info.mode_requested = 'dense_one_shot_first';

fprintf('    [audio-direct-rzphi-auto] Nr = %d, Nz = %d, Nphi = %d, Nvirt = %.6g, eval pts = %d\n', ...
    Nr, Nz, Nphi, Nvirt, Npt);
fprintf('    [audio-direct-rzphi-auto] dense one-shot estimated peak memory = %.3g GB\n', dense_info.estimated_peak_gb);

if local_should_use_dense_one_shot(dense_info.estimated_peak_gb, dense_dim_cfg)
    fprintf('    [audio-direct-rzphi-auto] using dense one-shot volume integral.\n');
    dense_info.mode_used = 'dense_one_shot';
    dense_info.reported_peak_gb = dense_info.estimated_peak_gb;

    [RHOv, Zv, PHIv] = ndgrid(rho_grid, z_grid_use, phi_all);
    Xv = RHOv .* cos(PHIv);
    Yv = RHOv .* sin(PHIv);

    Qbase = q_rz_use .* (rho_grid * ones(1, Nz)) * (dr * dz);
    Qv = reshape(Qbase, [Nr, Nz, 1]) .* reshape(exp(1j * ma * phi_all) * dphi, [1, 1, Nphi]);

    x_obs = rho_eval(:);
    y_obs = zeros(size(x_obs));
    z_obs = z_eval(:);

    Xv_vec = Xv(:).';
    Yv_vec = Yv(:).';
    Zv_vec = Zv(:).';
    Qv_vec = Qv(:);

    R = sqrt((x_obs - Xv_vec).^2 + (y_obs - Yv_vec).^2 + (z_obs - Zv_vec).^2);
    R(R < 1e-12) = 1e-12;

    G = exp(1j * ka .* R) ./ (4*pi*R);
    pa_pts = -1j * rho0 * wa * (G * Qv_vec);

else
    fprintf('    [audio-direct-rzphi-auto] dense estimate exceeds %.3g GB; switching to blockwise summation.\n', ...
        dense_dim_cfg.max_estimated_memory_gb);

    target_block_gb = local_get_fallback_target_gb(dense_dim_cfg) * local_get_block_safety_factor(dense_dim_cfg);
    [z_block_size, phi_block_size, block_peak_gb] = local_choose_audio_z_phi_block_size( ...
        Nr, Nz, Nphi, target_block_gb);

    dense_info.mode_used = 'blockwise_fallback';
    dense_info.fallback_target_block_memory_gb = target_block_gb;
    dense_info.z_block_size = z_block_size;
    dense_info.phi_block_size = phi_block_size;
    dense_info.block_estimated_peak_gb = block_peak_gb;
    dense_info.reported_peak_gb = max(block_peak_gb, local_whos_bytes_gb(whos('q_rz_use','rho_grid','z_grid_use')));

    fprintf('    [audio-direct-rzphi-auto] fallback z_block_size = %d, phi_block_size = %d\n', ...
        z_block_size, phi_block_size);
    fprintf('    [audio-direct-rzphi-auto] estimated peak per block = %.3g GB\n', block_peak_gb);

    pa_pts = complex(zeros(Npt,1));
    for ip = 1:Npt
        pa_pts(ip) = calc_audio_direct_single_point_from_qrz_blockwise( ...
            rho_grid, z_grid_use, q_rz_use, ...
            dr, dz, Nphi, ma, ...
            ka, rho0, wa, ...
            rho_eval(ip), z_eval(ip), ...
            z_block_size, phi_block_size);
        fprintf('        audio point progress: %4d / %4d (%.1f%%)\n', ip, Npt, 100*ip/Npt);
    end
end
end

function [z_out, q_out] = build_mirrored_qrz_about_z0(z_in, q_in)
z_in = z_in(:);
mask_pos = z_in > 1e-12;
z_neg = -flipud(z_in(mask_pos));
q_neg = fliplr(q_in(:, mask_pos));
z_out = [z_neg; z_in];
q_out = [q_neg, q_in];
end

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
% Config builders
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

function calc_cfg = build_calc_cfg_for_ultra_both( ...
    N_FHT, rho_max, zu_max, za_max, delta, use_parallel, medium_cfg, source_cfg)
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
% General helpers
%% ============================================================
function g = build_uniform_grid_0_to_max(maxv, d)
n = floor(maxv / d);
g = (0:n) * d;
if isempty(g) || abs(g(end) - maxv) > 1e-12
    g = [g, maxv];
end
g = unique(g, 'stable');
end

function spl = local_pressure_to_spl(p_amp, pref)
spl = 20 * log10(p_amp / max(pref, eps) / sqrt(2) + eps);
end

function pref = medium_pref_from_cfg(medium_cfg)
pref = medium_cfg.pref;
end

function folder_name = build_pu_cache_folder_name(m1, m2, rho_max, zu_max, N_FHT, delta, dis_coe, Nphi)
s_rho   = local_num_to_tag(rho_max);
s_zu    = local_num_to_tag(zu_max);
s_delta = local_num_to_tag(delta);
folder_name = sprintf('m1%dm2%drm%s_zu%s_nf%d_dt%s_dc%d_np%d', ...
    m1, m2, s_rho, s_zu, N_FHT, s_delta, dis_coe, Nphi);
end

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

function local_write_all_params_txt(txt_path, S)
fid = fopen(txt_path, 'w');
if fid < 0
    error('Cannot open txt file for writing: %s', txt_path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '============================================================\n');
fprintf(fid, 'ALL PARAMETERS\n');
fprintf(fid, 'Generated time: %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, '============================================================\n\n');
local_dump_any(fid, 'params', S, 0);
end

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
fprintf(fid, '%s%s = <unprintable type: %s>\n', sp, name, class(val));
end

function A = get_field_with_aliases(S, names)
for kk = 1:numel(names)
    if isfield(S, names{kk})
        A = S.(names{kk});
        return;
    end
end
error('Cannot find any of the fields: %s', strjoin(names, ', '));
end

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

function cache_info = find_reusable_single_field_cache( ...
    cache_root, which_field, target_m, ...
    rho_max, zu_max, N_FHT, delta, dis_coe, Nphi)
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
    if isempty(tok), continue; end
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
% Prepare direct-reference source samples for two possibly different modes
%
% IMPORTANT:
%   make_source_velocity uses source.m as a single modal order for both
%   f1 and f2. For mixed cases such as m1 ~= m2, calling it separately
%   may generate different polar source discretizations. Therefore, here
%   we generate ONE common source geometry and then manually assign the
%   two vortex phase factors on the same source points.
%% ============================================================
function [source_prep_ref, Xs_ref, Ys_ref, Zs_ref, q1_ref, q2_ref] = ...
    prepare_direct_reference_sources_two_modes(source_cfg, medium_cfg, calc_cfg_ref)

% Choose a common polar source discretization that is compatible with both
% modal orders. In make_source_velocity, the polar discretization enforces
% the number of angular sectors to be a multiple of source.m_used. Using a
% common nonzero modal order avoids inconsistent source grids in mixed cases.
m_abs = unique(abs([source_cfg.m1, source_cfg.m2]));
m_abs = m_abs(m_abs > 0);
if isempty(m_abs)
    m_common = 0;
else
    m_common = m_abs(1);
    for kk = 2:numel(m_abs)
        m_common = lcm(m_common, m_abs(kk));
    end
end

source_cfg_common = source_cfg;
source_cfg_common.m  = m_common;
source_cfg_common.m1 = m_common;
source_cfg_common.m2 = m_common;

[source_prep_ref, ~, ~] = make_source_velocity(source_cfg_common, medium_cfg, calc_cfg_ref);

Xs_ref = source_prep_ref.dim.Xs(:);
Ys_ref = source_prep_ref.dim.Ys(:);
Zs_ref = source_prep_ref.dim.Zs(:);
Nsrc = numel(Xs_ref);

if isfield(source_prep_ref.dim, 'dA_pts') && ~isempty(source_prep_ref.dim.dA_pts)
    dA_ref = source_prep_ref.dim.dA_pts(:);
else
    dA_ref = source_prep_ref.dim.dA;
end

if isscalar(dA_ref)
    dA_vec = dA_ref * ones(Nsrc, 1);
else
    dA_vec = dA_ref(:);
    if numel(dA_vec) ~= Nsrc
        error('Length of dA_pts does not match the number of source samples.');
    end
end

rho_s = hypot(Xs_ref, Ys_ref);
phi_s = atan2(Ys_ref, Xs_ref);
aperture = double(rho_s <= source_cfg.a + 1e-12);

% Build radial velocity amplitudes consistently with make_source_velocity.
% See make_source_velocity: Vortex-m uses a uniform radial core and the
% azimuthal factor exp(1i*m*phi); Focus uses exp(-1i*Re(k)*(sqrt(rho^2+F^2)-F)).
profile = lower(strtrim(char(source_cfg.profile)));
if ~isfield(source_cfg, 'v_ratio') || isempty(source_cfg.v_ratio)
    source_cfg.v_ratio = 1;
end

switch profile
    case 'uniform'
        amp1 = source_cfg.v0 .* aperture;
        amp2 = source_cfg.v_ratio * source_cfg.v0 .* aperture;

    case 'vortex-m'
        amp1 = source_cfg.v0 .* aperture;
        amp2 = source_cfg.v_ratio * source_cfg.v0 .* aperture;

    case 'focus'
        F = source_cfg.F;
        amp1 = source_cfg.v0 .* exp(-1i*real(source_prep_ref.k1).*(sqrt(rho_s.^2 + F.^2) - F)) .* aperture;
        amp2 = source_cfg.v_ratio * source_cfg.v0 .* exp(-1i*real(source_prep_ref.k2).*(sqrt(rho_s.^2 + F.^2) - F)) .* aperture;

    case 'poly'
        if isfield(source_cfg,'poly_n') && ~isempty(source_cfg.poly_n)
            n = source_cfg.poly_n;
        else
            n = 2;
        end
        amp1 = source_cfg.v0 .* (max(rho_s,0)./source_cfg.a).^n .* aperture;
        amp2 = source_cfg.v_ratio * source_cfg.v0 .* (max(rho_s,0)./source_cfg.a).^n .* aperture;

    case 'custom'
        % For arbitrary custom source functions, use make_source_velocity's
        % common-geometry handles if available. This branch assumes that any
        % required modal phase is already included in the custom handles.
        if isfield(source_prep_ref, 'dim') && isfield(source_prep_ref.dim, 'Vn_pts_f1')
            vn1_ref = source_prep_ref.dim.Vn_pts_f1(:);
            vn2_ref = source_prep_ref.dim.Vn_pts_f2(:);
            q1_ref = vn1_ref .* dA_vec;
            q2_ref = vn2_ref .* dA_vec;
            source_prep_ref.dim.Xs = Xs_ref;
            source_prep_ref.dim.Ys = Ys_ref;
            source_prep_ref.dim.Zs = Zs_ref;
            source_prep_ref.dim.Vn_pts_f1 = vn1_ref;
            source_prep_ref.dim.Vn_pts_f2 = vn2_ref;
            source_prep_ref.dim.dA_pts_common = dA_vec;
            fprintf('    [source prep] custom source uses make_source_velocity output on common geometry.\n');
            fprintf('    [source prep] Nsrc = %d, m_common = %d\n', Nsrc, m_common);
            return;
        else
            error('Custom source profile is not available in source_prep_ref.dim.Vn_pts_f1/2.');
        end

    otherwise
        error('Unsupported source profile in prepare_direct_reference_sources_two_modes: %s', source_cfg.profile);
end

% Assign the two requested vortex orders on the same source geometry.
vn1_ref = amp1 .* exp(1i * source_cfg.m1 * phi_s);
vn2_ref = amp2 .* exp(1i * source_cfg.m2 * phi_s);

q1_ref = vn1_ref .* dA_vec;
q2_ref = vn2_ref .* dA_vec;

% Optional self-check: for the common order, manual construction should
% agree with make_source_velocity. This validates the phase convention.
try
    if source_cfg.m1 == m_common && isfield(source_prep_ref.dim, 'Vn_pts_f1')
        vn0 = source_prep_ref.dim.Vn_pts_f1(:);
        rel_err = norm(vn0 - vn1_ref) / max(norm(vn0), eps);
        fprintf('    [source prep check] rel_err manual vs make_source f1 = %.3e\n', rel_err);
    end
catch
end

source_prep_ref.dim.Xs = Xs_ref;
source_prep_ref.dim.Ys = Ys_ref;
source_prep_ref.dim.Zs = Zs_ref;
source_prep_ref.dim.Vn_pts_f1 = vn1_ref;
source_prep_ref.dim.Vn_pts_f2 = vn2_ref;
source_prep_ref.dim.dA_pts_common = dA_vec;
source_prep_ref.dim.m_common_for_source_grid = m_common;

fprintf('    [source prep] common source geometry used for both modes.\n');
fprintf('    [source prep] Nsrc = %d, m1 = %d, m2 = %d, m_common = %d\n', ...
    Nsrc, source_cfg.m1, source_cfg.m2, m_common);
end

%% ============================================================
% Direct ultrasound at selected points
%% ============================================================
function [p1_pts, p2_pts] = compute_ultrasound_direct_points_from_ref_sources( ...
    rho_pts, z_pts, ...
    Xs, Ys, Zs, q1, q2, ...
    k1, k2, w1, w2, rho0, src_block_size)

rho_pts = rho_pts(:);
z_pts   = z_pts(:);

if numel(rho_pts) ~= numel(z_pts)
    error('rho_pts and z_pts must have the same number of elements.');
end

xo = rho_pts;
yo = zeros(size(xo));
zo = z_pts;

[p1_pts, p2_pts] = calc_ultra_direct_obs_block_from_xyz( ...
    1, numel(xo), xo, yo, zo, ...
    Xs, Ys, Zs, q1, q2, ...
    k1, k2, w1, w2, rho0, src_block_size);
end

%% ============================================================
% Dense-DIM memory / fallback helpers
%% ============================================================
function cfg = local_default_dense_dim_cfg()
cfg = struct();
cfg.allow_extreme_allocation = false;
cfg.max_estimated_memory_gb  = 512;
cfg.auto_fallback_to_blocks = true;
cfg.fallback_target_block_memory_gb = 512;
cfg.share_fallback_memory_across_workers = true;
cfg.block_memory_safety_factor = 0.85;
cfg.stop_if_exceeds_limit = false;
end

function tf = local_should_use_dense_one_shot(estimated_gb, cfg)
cfg = local_complete_dense_cfg(cfg);
if cfg.allow_extreme_allocation
    tf = true;
else
    tf = estimated_gb <= cfg.max_estimated_memory_gb;
end
if ~tf && ~cfg.auto_fallback_to_blocks && cfg.stop_if_exceeds_limit
    error(['Dense one-shot DIM requires %.3g GB, exceeding the configured limit %.3g GB, ', ...
           'and auto_fallback_to_blocks is disabled.'], estimated_gb, cfg.max_estimated_memory_gb);
end
end

function cfg = local_complete_dense_cfg(cfg)
if nargin < 1 || isempty(cfg)
    cfg = local_default_dense_dim_cfg();
    return;
end
if ~isfield(cfg,'allow_extreme_allocation') || isempty(cfg.allow_extreme_allocation)
    cfg.allow_extreme_allocation = false;
end
if ~isfield(cfg,'max_estimated_memory_gb') || isempty(cfg.max_estimated_memory_gb)
    cfg.max_estimated_memory_gb = 512;
end
if ~isfield(cfg,'auto_fallback_to_blocks') || isempty(cfg.auto_fallback_to_blocks)
    cfg.auto_fallback_to_blocks = true;
end
if ~isfield(cfg,'fallback_target_block_memory_gb') || isempty(cfg.fallback_target_block_memory_gb)
    cfg.fallback_target_block_memory_gb = cfg.max_estimated_memory_gb;
end
if ~isfield(cfg,'share_fallback_memory_across_workers') || isempty(cfg.share_fallback_memory_across_workers)
    cfg.share_fallback_memory_across_workers = true;
end
if ~isfield(cfg,'block_memory_safety_factor') || isempty(cfg.block_memory_safety_factor)
    cfg.block_memory_safety_factor = 0.85;
end
if ~isfield(cfg,'stop_if_exceeds_limit') || isempty(cfg.stop_if_exceeds_limit)
    cfg.stop_if_exceeds_limit = false;
end
end

function target_gb = local_get_fallback_target_gb(cfg)
cfg = local_complete_dense_cfg(cfg);
target_gb = cfg.fallback_target_block_memory_gb;
end

function safety = local_get_block_safety_factor(cfg)
cfg = local_complete_dense_cfg(cfg);
safety = cfg.block_memory_safety_factor;
end

function n = local_get_active_block_count(use_parallel, cfg)
cfg = local_complete_dense_cfg(cfg);
if use_parallel && cfg.share_fallback_memory_across_workers
    p = gcp('nocreate');
    if isempty(p)
        n = 1;
    else
        n = max(1, p.NumWorkers);
    end
else
    n = 1;
end
end

function obs_block_size = local_choose_ultra_obs_block_size(Nsrc, target_gb_each, Nobs)
bytes_per_pair = 4*8 + 2*16;  % dx,dy,dz,R + G1,G2
max_pairs = floor(target_gb_each * 1024^3 / bytes_per_pair);
obs_block_size = floor(max_pairs / max(double(Nsrc),1));
obs_block_size = max(1, min(Nobs, obs_block_size));
end

function [z_block_size, phi_block_size, peak_gb] = local_choose_audio_z_phi_block_size(Nr, Nz, Nphi, target_gb)
% Estimate block arrays in calc_audio_direct_single_point_from_qrz_blockwise.
% Use all azimuth samples when possible and reduce z_block first.
bytes_per_triplet = 8 + 8 + 16 + 16; % R2/R + G + Qfull, conservative lower estimate
phi_block_size = Nphi;
max_triplets = floor(target_gb * 1024^3 / bytes_per_triplet);
z_block_size = floor(max_triplets / max(double(Nr)*double(phi_block_size),1));
if z_block_size < 1
    z_block_size = 1;
    phi_block_size = floor(max_triplets / max(double(Nr),1));
    phi_block_size = max(1, min(Nphi, phi_block_size));
else
    z_block_size = min(Nz, z_block_size);
end
peak_gb = double(Nr) * double(z_block_size) * double(phi_block_size) * bytes_per_triplet / 1024^3;
end

function info = local_estimate_dense_ultra_memory_gb(Nobs, Nsrc)
bytes_real = 8;
bytes_complex = 16;
M = double(Nobs) * double(Nsrc);
bytes = M * (4*bytes_real + 2*bytes_complex);
info = struct();
info.Nobs = Nobs;
info.Nsrc = Nsrc;
info.matrix_entries = M;
info.estimated_peak_gb = bytes / 1024^3;
info.reported_peak_gb = info.estimated_peak_gb;
end

function info = local_estimate_dense_audio_memory_gb(Npt, Nvirt)
bytes_real = 8;
bytes_complex = 16;
Mobs = double(Npt) * double(Nvirt);
Mvirt = double(Nvirt);
bytes = Mobs * (bytes_real + bytes_complex) + Mvirt * (5*bytes_real + 2*bytes_complex);
info = struct();
info.Npt = Npt;
info.Nvirt = Nvirt;
info.matrix_entries = Mobs;
info.estimated_peak_gb = bytes / 1024^3;
info.reported_peak_gb = info.estimated_peak_gb;
end

%% ============================================================
% Memory and article-table helpers
%% ============================================================
function gb = local_whos_bytes_gb(S)
if isempty(S)
    gb = 0;
else
    gb = sum([S.bytes]) / 1024^3;
end
end

function gb = local_get_matlab_memory_used_gb()
gb = NaN;
try
    if ispc
        u = memory;
        gb = u.MemUsedMATLAB / 1024^3;
    end
catch
    gb = NaN;
end
end

function dgb = local_memory_delta_gb(gb_before, gb_after)
if isfinite(gb_before) && isfinite(gb_after)
    dgb = max(gb_after - gb_before, 0);
else
    dgb = 0;
end
end

function article_table = build_article_cost_table( ...
    N_src_eff, N_rho_v, N_phi_v, N_z_v_pos, N_z_v_audio, N_o_bench, ...
    N_FHT, N_z_v_prop, N_z_o_prop, N_z_c_prop, ...
    timing, memory)

article_table = struct([]);

article_table(1).Step = 'Ultrasound';
article_table(1).Method = 'DIM';
article_table(1).Complexity = 'O(Nrho_s Nphi_s Nrho_v Nphi_v Nz_v)';
article_table(1).Samples = sprintf('Nsrc=%d; virtual rz=%d x %d; field=%d', ...
    N_src_eff, N_rho_v, N_z_v_pos, N_rho_v*N_z_v_pos);
article_table(1).Time_s = timing.dim_ultrasound_full_rz;
article_table(1).Memory_GB = memory.dim_ultrasound_full_rz_gb;
article_table(1).Note = 'full DIM ultrasound recomputation on rz grid; no cache used';

article_table(2).Step = 'Ultrasound';
article_table(2).Method = 'Proposed';
article_table(2).Complexity = 'O(Nz_v N_FHT log N_FHT)';
article_table(2).Samples = sprintf('N_FHT=%d; Nz_v=%d', N_FHT, N_z_v_prop);
article_table(2).Time_s = timing.proposed_ultrasound;
article_table(2).Memory_GB = memory.proposed_ultrasound_arrays_gb;
article_table(2).Note = 'transform ultrasound + q construction';

article_table(3).Step = 'Audible';
article_table(3).Method = 'DIM';
article_table(3).Complexity = 'O(Nrho_o Nphi_o Nz_o Nrho_v Nphi_v Nz_v)';
article_table(3).Samples = sprintf('virtual=%d x %d x %d; obs=%d points', ...
    N_rho_v, N_phi_v, N_z_v_audio, N_o_bench);
article_table(3).Time_s = timing.dim_audible_10pts;
article_table(3).Memory_GB = memory.direct_qrz_arrays_gb;
article_table(3).Note = 'DIM audio direct integration at 10 points';

article_table(4).Step = 'Audible';
article_table(4).Method = 'Proposed';
article_table(4).Complexity = 'O((Nz_v+Nz_o)N_FHT log N_FHT + N_FHT Nz_c log Nz_c)';
article_table(4).Samples = sprintf('N_FHT=%d; Nz_v=%d; Nz_o=%d; Nz_c=%d', ...
    N_FHT, N_z_v_prop, N_z_o_prop, N_z_c_prop);
article_table(4).Time_s = timing.proposed_audible;
article_table(4).Memory_GB = memory.proposed_audible_arrays_gb;
article_table(4).Note = 'FHT + axial FFT convolution';
end

function print_article_cost_table(T)
fprintf('%-12s %-10s %-62s %-42s %12s %12s  %s\n', ...
    'Step','Method','Complexity','Samples','Time(s)','Memory(GB)','Note');
fprintf('%s\n', repmat('-',1,180));
for k = 1:numel(T)
    fprintf('%-12s %-10s %-62s %-42s %12.3f %12.4f  %s\n', ...
        T(k).Step, T(k).Method, T(k).Complexity, T(k).Samples, ...
        T(k).Time_s, T(k).Memory_GB, T(k).Note);
end
end

function write_article_cost_table_csv(path, T)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot open CSV file for writing: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, 'Step,Method,Complexity,Samples,Time_s,Memory_GB,Note\n');
for k = 1:numel(T)
    fprintf(fid, '"%s","%s","%s","%s",%.9g,%.9g,"%s"\n', ...
        T(k).Step, T(k).Method, T(k).Complexity, T(k).Samples, ...
        T(k).Time_s, T(k).Memory_GB, T(k).Note);
end
end

function write_article_cost_table_txt(path, T)
fid = fopen(path, 'w');
if fid < 0
    error('Cannot open TXT file for writing: %s', path);
end
cleanupObj = onCleanup(@() fclose(fid));
fprintf(fid, '%-12s %-10s %-62s %-42s %12s %12s  %s\n', ...
    'Step','Method','Complexity','Samples','Time(s)','Memory(GB)','Note');
fprintf(fid, '%s\n', repmat('-',1,180));
for k = 1:numel(T)
    fprintf(fid, '%-12s %-10s %-62s %-42s %12.3f %12.4f  %s\n', ...
        T(k).Step, T(k).Method, T(k).Complexity, T(k).Samples, ...
        T(k).Time_s, T(k).Memory_GB, T(k).Note);
end
end
