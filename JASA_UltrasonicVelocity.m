%% ============================================================
% MAIN: King vs DIM velocity comparison (two profiles automatically)
%
% OUTPUT RULE:
%   if m == 0:
%       1) radial profile @ z ~= z_target
%       2) axial profile   @ rho = first positive rho on King grid
%   else:
%       1) radial profile @ z ~= z_target
%       2) oblique profile along rho(z)=z*tan(theta_line_deg)
%
% FEATURES
% - output all three velocity components: vrho / vphi / vz
% - output 1D particle-velocity-level figures and PVL-difference figures
% - x-axis of 1D figures is logarithmic
% - all comparison points are snapped to King grid first
% - z grid is fully determined by make_source_velocity() via fht.z_ultra
% - radial profile display is limited to rho_max_radial_use
% - for m = 3, oblique profile can use an expanded rho range to reach z = 10 m
% - figures are saved as PNG / PDF / FIG with unified size
% ============================================================

clear; clc; close all;

%% -------------------- save figures (GLOBAL control) --------------------
global SAVE_FIGS SAVE_DIR
SAVE_FIGS = false;   % <<< true to save, false to disable saving

%% -------------------- comparison settings --------------------
z_target = 1;              % radial profile taken at z ~= z_target
ds = 16 * 2;               % radial downsample factor on King rho-grid

z_show_max = 10;           % x-axis upper limit for axial / oblique profiles
rho_show_min = 1e-2;       % x-axis lower limit for radial profiles

rho_max_radial_use = 0.5;  % radial profile sampling / display upper limit

%% -------------------- all velocity components --------------------
vel_comp_list = {'vrho','vphi','vz'};

%% -------------------- medium --------------------
medium.c0 = 343;
medium.rho0 = 1.21;
medium.beta = 1.2;
medium.pref = 2e-5;
medium.use_absorp = true;
medium.atten_handle = @(f) AbsorpAttenCoef(f);

%% -------------------- source --------------------
source.profile = 'Vortex-m'; % 'Uniform' | 'Focus' | 'Vortex-m' | 'Poly' | 'Custom'
source.a = 0.05;
source.v0 = 0.108;
source.v_ratio = 1;
source.m = 0;       % 0 or 3
source.F = 0.2;
% source.poly_n = 0;

source.f1 = 40e3;
source.fa = 0.5e3;

theta_line_deg = 0;   % default for m = 0
if source.m == 1
    theta_line_deg = 3;
end
if source.m == 2
    theta_line_deg = 6;
end
if source.m == 3
    theta_line_deg = 9;
end

%% -------------------- calc parameters --------------------
calc = struct();

% --- FHT / King ---
calc.fht.N_FHT = 16384 * 1;
calc.fht.rho_max = rho_max_radial_use;
calc.fht.Nh_scale = 2;

if source.m == 0
    calc.fht.NH_scale = 1.2;
end
if source.m == 1
    calc.fht.NH_scale = 4;
end
if source.m == 3
    calc.fht.NH_scale = 4;
end

calc.fht.Nh_v_scale = 1.1;
calc.fht.zu_max = 10;
calc.fht.za_max = 0.5;
calc.fht.delta = medium.c0 / source.f1 / 0.5;
calc.fht.delta = 0.01;

% For m = 3, only expand the computation range for the oblique profile.
% Radial profile remains limited to rho_max_radial_use.
if source.m == 3
    rho_needed_for_oblique = z_show_max * tand(theta_line_deg);
    calc.fht.rho_max = max(calc.fht.rho_max, 1.05 * rho_needed_for_oblique);
end

% --- z grid mode determined inside make_source_velocity ---
calc.fht.z_sampling = 'log';
calc.fht.Nz_ultra = 400;
calc.fht.z_log_min = calc.fht.delta;

% --- DIM source discretization ---
calc.dim.use_freq = 'f2';
calc.dim.dis_coe = 16 * 1;
calc.dim.margin = 1;
calc.dim.src_discretization = 'polar';   % 'cart' | 'polar'

% --- King analytic spectrum stability ---
calc.king.gspec_method = 'transform';
calc.king.eps_kzz = 1e-3;
calc.king.eps_phase = calc.king.eps_kzz;
calc.king.kz_min = 1e-12;
calc.king.band_refine.enable = false;

% --- ASM settings ---
calc.asm.pad_factor = 16;
calc.asm.kzz_eps = 1e-12;

%% -------------------- save setup --------------------
if SAVE_FIGS
    tstr = datestr(datetime('now'), 'mmdd_HHMM');
    a_str = sprintf('%.2fm', source.a);
    m_str = sprintf('m=%d', source.m);
    SAVE_DIR = fullfile('JASA_ultra_velocity_compare', ...
        sprintf('%s_%s__%s', a_str, m_str, tstr));

    if ~exist(SAVE_DIR, 'dir')
        mkdir(SAVE_DIR);
    end

    local_write_runinfo_txt(SAVE_DIR, medium, source, calc, ...
        z_target, ds, theta_line_deg, vel_comp_list, ...
        z_show_max, rho_show_min, rho_max_radial_use);
else
    SAVE_DIR = '';
end

%% ==================== PROFILER helpers ====================
get_mem_mb = @() local_get_mem_mb();
fmt_mem = @(x) local_fmt_mem(x);

%% ============================================================
% PREPARE KING GRID (z grid determined by make_source_velocity)
%% ============================================================
[~, fht_tmp] = make_source_velocity(source, medium, calc);

rho_full = fht_tmp.xh(:);
zK_full  = fht_tmp.z_ultra(:);

rho_max_oblique_use = calc.fht.rho_max;   % expanded only when needed, e.g., m = 3

% ----- radial rho points: clipped to radial display / sampling range -----
rho_rad = rho_full(1:ds:end);
rho_rad = rho_rad(:);
rho_rad = rho_rad(rho_rad > 0 & rho_rad <= rho_max_radial_use);

% ----- radial z plane -----
[~, iz_tar] = min(abs(zK_full - z_target));
z_use = zK_full(iz_tar);

% ----- positive z points already determined by make_source_velocity -----
z_pos = zK_full(zK_full > 0);

%% ============================================================
% PROFILE 1: radial @ z = z_use
%% ============================================================
fprintf('\n==================== PROFILE 1 / RADIAL ====================\n');
fprintf('z_target = %.3f m, snapped z_use = %.6f m\n', z_target, z_use);
fprintf('Radial rho range for sampling/display: %.6g to %.6g m\n', ...
    rho_show_min, rho_max_radial_use);

obs_grid_rad.dim.x = rho_rad(:).';
obs_grid_rad.dim.y = 0;
obs_grid_rad.dim.z = z_use;
obs_grid_rad.dim.block_size = 200000;

mem0 = get_mem_mb(); t0 = tic;
res_rad = calc_ultrasound_velocity_field(source, medium, calc, obs_grid_rad, 'both');
t_rad = toc(t0); mem1 = get_mem_mb();

fprintf('Time (radial): %.3f s\n', t_rad);
fprintf('Memory: start %s, end %s, delta %s\n', ...
    fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

zK = res_rad.king.z(:);
rhoK = res_rad.king.rho(:);

[~, izK_rad] = min(abs(zK - z_use));
idx_rho_rad = local_snap_to_grid_indices(rho_rad, rhoK);

%% ----- plot all 3 velocity components on radial profile -----
for ic = 1:numel(vel_comp_list)
    vel_comp = vel_comp_list{ic};
    fprintf('  -> radial component: %s\n', vel_comp);

    [vK_all_rad, vD_all_rad, ref_name] = local_get_velocity_pair(res_rad, vel_comp);

    vK_rad = vK_all_rad(idx_rho_rad, izK_rad);
    vK_rad = vK_rad(:);

    vD_rad = local_extract_dim_radial(vD_all_rad, res_rad, rho_rad, z_use);
    vD_rad = vD_rad(:);

    spec_rad = struct();
    spec_rad.coord_vec   = rho_rad(:);
    spec_rad.coord_label = '\rho (m)';
    spec_rad.coord_name  = 'rho';
    spec_rad.coord_title = sprintf('radial profile @ z = %.6f m', z_use);
    spec_rad.file_tag    = sprintf('%s_radial', vel_comp);
    spec_rad.use_logx    = true;
    spec_rad.is_oblique  = false;
    spec_rad.x_show      = [rho_show_min, rho_max_radial_use];

    local_plot_compare_1d_velocity(spec_rad, vK_rad, vD_rad, ref_name, vel_comp, medium);
end

%% ============================================================
% PROFILE 2:
%   m==0  -> axial
%   m~=0  -> oblique
%% ============================================================
if source.m == 0
    fprintf('\n==================== PROFILE 2 / AXIAL ====================\n');

    rho_ax = rho_rad(1);   % first positive snapped rho

    obs_grid_line.dim.x = rho_ax;
    obs_grid_line.dim.y = 0;
    obs_grid_line.dim.z = z_pos(:).';
    obs_grid_line.dim.block_size = 200000;

    mem0 = get_mem_mb(); t0 = tic;
    res_line = calc_ultrasound_velocity_field(source, medium, calc, obs_grid_line, 'both');
    t_line = toc(t0); mem1 = get_mem_mb();

    fprintf('Axial rho = %.6g m\n', rho_ax);
    fprintf('Time (axial): %.3f s\n', t_line);
    fprintf('Memory: start %s, end %s, delta %s\n', ...
        fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

    zK2 = res_line.king.z(:);
    rhoK2 = res_line.king.rho(:);

    [~, ixK_ax] = min(abs(rhoK2 - rho_ax));
    idx_z_ax = local_snap_to_grid_indices(z_pos, zK2);

    for ic = 1:numel(vel_comp_list)
        vel_comp = vel_comp_list{ic};
        fprintf('  -> axial component: %s\n', vel_comp);

        [vK_all_line, vD_all_line, ref_name] = local_get_velocity_pair(res_line, vel_comp);

        vK_line = vK_all_line(ixK_ax, idx_z_ax).';
        vK_line = vK_line(:);

        vD_line = local_extract_dim_axial(vD_all_line, res_line, rho_ax, z_pos);
        vD_line = vD_line(:);

        spec_line = struct();
        spec_line.coord_vec   = z_pos(:);
        spec_line.coord_label = 'z (m)';
        spec_line.coord_name  = 'z';
        spec_line.coord_title = sprintf('axial profile @ \\rho = %.6g m', rho_ax);
        spec_line.file_tag    = sprintf('%s_axial', vel_comp);
        spec_line.use_logx    = true;
        spec_line.is_oblique  = false;
        spec_line.x_show      = [1e-2, z_show_max];

        local_plot_compare_1d_velocity(spec_line, vK_line, vD_line, ref_name, vel_comp, medium);
    end

else
    fprintf('\n==================== PROFILE 2 / OBLIQUE ====================\n');

    z_line = z_pos(:);
    rho_tar = z_line * tand(theta_line_deg);

    keep = z_line > 0 & z_line <= z_show_max & ...
           rho_tar > 0 & rho_tar <= rho_max_oblique_use;

    z_line = z_line(keep);
    rho_tar = rho_tar(keep);

    idx_rho_line = local_snap_to_grid_indices(rho_tar, rhoK);
    rho_snap = rhoK(idx_rho_line);

    % build DIM observation grid from snapped King points
    x_line_unique = unique(rho_snap(:).', 'stable');

    obs_grid_line.dim.x = x_line_unique;
    obs_grid_line.dim.y = 0;
    obs_grid_line.dim.z = z_line(:).';
    obs_grid_line.dim.block_size = 200000;

    mem0 = get_mem_mb(); t0 = tic;
    res_line = calc_ultrasound_velocity_field(source, medium, calc, obs_grid_line, 'both');
    t_line = toc(t0); mem1 = get_mem_mb();

    fprintf('Oblique angle to z-axis = %.2f deg\n', theta_line_deg);
    fprintf('Oblique target: z <= %.6f m\n', z_show_max);
    fprintf('Oblique rho at z = %.2f m: rho = %.6f m\n', ...
        z_show_max, z_show_max * tand(theta_line_deg));
    fprintf('rho_max_oblique_use = %.6f m\n', rho_max_oblique_use);
    fprintf('rho_max_radial_use  = %.6f m\n', rho_max_radial_use);
    fprintf('Number of snapped line points = %d\n', numel(z_line));
    fprintf('Time (oblique): %.3f s\n', t_line);
    fprintf('Memory: start %s, end %s, delta %s\n', ...
        fmt_mem(mem0), fmt_mem(mem1), fmt_mem(mem1-mem0));

    zK2 = res_line.king.z(:);

    idx_z_line = local_snap_to_grid_indices(z_line, zK2);

    for ic = 1:numel(vel_comp_list)
        vel_comp = vel_comp_list{ic};
        fprintf('  -> oblique component: %s\n', vel_comp);

        [vK_all_line, vD_all_line, ref_name] = local_get_velocity_pair(res_line, vel_comp);

        vK_line = zeros(numel(z_line),1);
        for k = 1:numel(z_line)
            vK_line(k) = vK_all_line(idx_rho_line(k), idx_z_line(k));
        end

        vD_line = local_extract_dim_oblique(vD_all_line, res_line, rho_snap, z_line);
        vD_line = vD_line(:);

        spec_line = struct();
        spec_line.coord_vec   = z_line(:);
        spec_line.coord_label = 'z (m)';
        spec_line.coord_name  = 'z';
        spec_line.coord_title = sprintf(['oblique profile @ \\theta = %.2f^\\circ to z-axis, ', ...
                                         '\\rho(z)=z\\tan\\theta (snapped to King grid)'], ...
                                         theta_line_deg);
        spec_line.file_tag    = sprintf('%s_oblique_%gdeg', vel_comp, theta_line_deg);
        spec_line.use_logx    = true;
        spec_line.is_oblique  = true;
        spec_line.theta_line_deg = theta_line_deg;
        spec_line.x_show      = [1e-2, z_show_max];

        local_plot_compare_1d_velocity(spec_line, vK_line, vD_line, ref_name, vel_comp, medium);
    end
end

fprintf('\n==================== SUMMARY ====================\n');
fprintf('velocity components = %s, %s, %s\n', vel_comp_list{:});
fprintf('Profile-1: radial @ z = %.6f m\n', z_use);
fprintf('Radial display range: rho = %.6g to %.6g m\n', ...
    rho_show_min, rho_max_radial_use);

if source.m == 0
    fprintf('Profile-2: axial @ rho = first positive snapped King rho\n');
else
    fprintf('Profile-2: oblique @ theta = %.2f deg to z-axis\n', theta_line_deg);
    fprintf('Oblique display target: z <= %.6f m\n', z_show_max);
    fprintf('Oblique computation rho limit = %.6g m\n', rho_max_oblique_use);
end

fprintf('The second profile is chosen according to m, exactly as in the pressure code.\n');
fprintf('z-grid is taken directly from make_source_velocity() -> fht.z_ultra\n');
fprintf('All comparison points are snapped to King grid before DIM comparison.\n');

%% ==================== local functions ====================

function [vK_all, vD_all, ref_name] = local_get_velocity_pair(res, vel_comp)
switch lower(strtrim(string(vel_comp)))
    case "vrho"
        vK_all = res.king.v_rho_f1;
        vD_all = res.dim.v_rho_f1;
    case "vphi"
        vK_all = res.king.v_phi_f1;
        vD_all = res.dim.v_phi_f1;
    case "vz"
        vK_all = res.king.v_z_f1;
        vD_all = res.dim.v_z_f1;
    otherwise
        error('vel_comp must be ''vrho'', ''vphi'', or ''vz''.');
end

ref_name = 'DIM';
if isfield(res, 'dim') && isfield(res.dim, 'method') && ~isempty(res.dim.method)
    method_str = lower(char(string(res.dim.method)));
    if contains(method_str, 'rayleigh')
        ref_name = 'Rayleigh';
    elseif contains(method_str, 'asm')
        ref_name = 'ASM';
    else
        ref_name = 'DIM';
    end
end
end

function idx = local_snap_to_grid_indices(vals, grid)
vals = vals(:);
grid = grid(:);
idx = zeros(size(vals));
for k = 1:numel(vals)
    [~, idx(k)] = min(abs(grid - vals(k)));
end
end

function vD = local_extract_dim_radial(vD_all, res, rho_vec, z_use)
if local_is_dim_rayleigh(res)
    vD = squeeze(vD_all(1,:,1)).';
else
    xA = res.dim.x(:);
    zA = res.dim.z(:);
    [~, iy0] = min(abs(res.dim.y(:) - 0));
    [~, iz0] = min(abs(zA - z_use));
    vA_line = squeeze(vD_all(iy0,:,iz0)).';
    vD = interp1(xA, vA_line, rho_vec(:), 'linear', 0);
end
vD = vD(:);
end

function vD = local_extract_dim_axial(vD_all, res, rho_ax, z_vec)
if local_is_dim_rayleigh(res)
    vD = squeeze(vD_all(1,1,:));
    vD = vD(:);

    zA = res.dim.z(:);
    if numel(zA) ~= numel(z_vec) || max(abs(zA(:) - z_vec(:))) > 1e-12
        vD = interp1(zA, vD, z_vec(:), 'linear', 0);
    end
else
    xA = res.dim.x(:);
    yA = res.dim.y(:);
    zA = res.dim.z(:);

    [~, ix0] = min(abs(xA - rho_ax));
    [~, iy0] = min(abs(yA - 0));

    vA_line = squeeze(vD_all(iy0,ix0,:));
    vA_line = vA_line(:);

    if numel(zA) == numel(z_vec) && max(abs(zA(:) - z_vec(:))) < 1e-12
        vD = vA_line;
    else
        vD = interp1(zA, vA_line, z_vec(:), 'linear', 0);
    end
end
vD = vD(:);
end

function vD = local_extract_dim_oblique(vD_all, res, rho_snap, z_line)
xA = res.dim.x(:);
zA = res.dim.z(:);

if local_is_dim_rayleigh(res)
    P = squeeze(vD_all(1,:,:));   % Nx x Nz
else
    yA = res.dim.y(:);
    [~, iy0] = min(abs(yA - 0));
    P = squeeze(vD_all(iy0,:,:)); % Nx x Nz
end

if size(P,1) ~= numel(xA)
    P = P.';
end

vD = zeros(numel(z_line),1);
for k = 1:numel(z_line)
    [~, ix] = min(abs(xA - rho_snap(k)));
    [~, iz] = min(abs(zA - z_line(k)));
    vD(k) = P(ix, iz);
end
end

function local_plot_compare_1d_velocity(spec, vK, vD, ref_name, vel_comp, medium)
coord_vec = spec.coord_vec(:);
vK = vK(:);
vD = vD(:);

mask = isfinite(coord_vec) & coord_vec > 0 & isfinite(vK) & isfinite(vD);
coord_vec = coord_vec(mask);
vK = vK(mask);
vD = vD(mask);

% ============================================================
% Particle velocity level definition
%   PVL = 20 log10((|v|/sqrt(2)) / v_ref)
%   v_ref = p_ref / (rho0 c0)
% ============================================================
v_ref   = medium.pref / (medium.rho0 * medium.c0);
v_floor = 1e-16;

PVLK = 20*log10(max(abs(vK), v_floor) / sqrt(2) / v_ref);
PVLD = 20*log10(max(abs(vD), v_floor) / sqrt(2) / v_ref);

diff_PVL = abs(PVLD - PVLK);

switch lower(strtrim(string(vel_comp)))
    case "vrho"
        ylab_level = '$\mathrm{PVL}_{\rho}\ \mathrm{(dB)}$';
        ylab_diff  = '$|\Delta \mathrm{PVL}_{\rho}|\ \mathrm{(dB)}$';
    case "vphi"
        ylab_level = '$\mathrm{PVL}_{\varphi}\ \mathrm{(dB)}$';
        ylab_diff  = '$|\Delta \mathrm{PVL}_{\varphi}|\ \mathrm{(dB)}$';
    case "vz"
        ylab_level = '$\mathrm{PVL}_{z}\ \mathrm{(dB)}$';
        ylab_diff  = '$|\Delta \mathrm{PVL}_{z}|\ \mathrm{(dB)}$';
    otherwise
        ylab_level = '$\mathrm{PVL}\ \mathrm{(dB)}$';
        ylab_diff  = '$|\Delta \mathrm{PVL}|\ \mathrm{(dB)}$';
end

% ---------- display range ----------
switch lower(spec.coord_name)
    case 'rho'
        if isfield(spec, 'x_show') && ~isempty(spec.x_show)
            x_show = spec.x_show;
        else
            x_show = [1e-2, max(coord_vec)];
        end

        xlab_text = '$\mathrm{Radial\ distance,}\ \rho\ \mathrm{(m)}$';
        xlab_interpreter = 'latex';

    case 'z'
        if isfield(spec, 'x_show') && ~isempty(spec.x_show)
            x_show = spec.x_show;
        else
            x_show = [1e-2, 10];
        end

        if isfield(spec, 'is_oblique') && spec.is_oblique
            if isfield(spec, 'theta_line_deg') && ~isempty(spec.theta_line_deg)
                xlab_text = sprintf(['$\\mathrm{Axial\\ coordinate,}\\ z\\ \\mathrm{(m)},\\ ', ...
                                     '\\rho=z\\tan %.0f^\\circ$'], ...
                                     spec.theta_line_deg);
            else
                xlab_text = '$\mathrm{Axial\ coordinate,}\ z\ \mathrm{(m)},\ \rho=z\tan\theta$';
            end
        else
            xlab_text = '$\mathrm{Axial\ distance,}\ z\ \mathrm{(m)}$';
        end

        xlab_interpreter = 'latex';

    otherwise
        x_show = [min(coord_vec), max(coord_vec)];
        xlab_text = spec.coord_label;
        xlab_interpreter = 'latex';
end

% ============================================================
% y-limits are determined only by points inside the displayed x-range.
% ============================================================
show_mask = coord_vec >= x_show(1) & coord_vec <= x_show(2) & isfinite(coord_vec);

PVL_for_ylim = [PVLD(show_mask); PVLK(show_mask)];
diff_for_ylim = diff_PVL(show_mask);

% ---------- unified style ----------
fig_pos = [100 100 800 320];
lw_main = 2.0;
lw_ref  = 2.0;
fs_ax   = 18;
fs_lab  = 24;
fs_leg  = 19;

fs_x_ax  = fs_ax  * 1.25;
fs_y_ax  = fs_ax  * 1.0;
fs_x_lab = fs_lab * 1.25;
fs_y_lab = fs_lab * 1.15;

% ============================================================
% 1) PARTICLE VELOCITY LEVEL FIGURE
% ============================================================
f1 = figure('Position', fig_pos, 'Color', 'w');
ax1 = axes('Parent', f1, 'Units', 'normalized', 'Position', [0.13 0.24 0.83 0.60]);
hold(ax1, 'on');

plot(coord_vec, PVLD, '-',  'LineWidth', lw_ref,  'Color', [0 0.2 1.0], ...
    'DisplayName', ref_name);
plot(coord_vec, PVLK, '--', 'LineWidth', lw_main, 'Color', [1 0 0], ...
    'DisplayName', 'King');

set(ax1, 'XScale', 'log');
xlim(ax1, x_show);
local_set_y_margin_general(ax1, PVL_for_ylim);
box(ax1, 'on');

xlabel(ax1, xlab_text, ...
    'FontSize', fs_x_lab, ...
    'Interpreter', xlab_interpreter);

ylabel(ax1, ylab_level, ...
    'FontSize', fs_y_lab, ...
    'Interpreter', 'latex');

set(ax1, 'LineWidth', 1.2, ...
    'FontSize', fs_y_ax, ...
    'TickDir', 'in', ...
    'TickLength', [0.015 0.015], ...
    'Layer', 'top');

ax1.XAxis.FontSize = fs_x_ax;
ax1.YAxis.FontSize = fs_y_ax;

leg1 = legend(ax1, 'Orientation', 'horizontal');
set(leg1, 'Units', 'normalized');
set(leg1, 'Position', [0.36 0.88 0.28 0.06]);
set(leg1, 'Box', 'on', 'FontSize', fs_leg, 'Interpreter', 'tex');
set(leg1, 'AutoUpdate', 'off');

local_save_fig_all(f1, sprintf('King_vs_%s_%s_PVL', ref_name, spec.file_tag));

% ============================================================
% 2) PARTICLE VELOCITY LEVEL DIFFERENCE FIGURE
% ============================================================
f2 = figure('Position', fig_pos, 'Color', 'w');
ax2 = axes('Parent', f2, 'Units', 'normalized', 'Position', [0.13 0.24 0.83 0.60]);
hold(ax2, 'on');

plot(coord_vec, diff_PVL, '-', 'LineWidth', 2.0, 'Color', [0 0.2 1.0]);

set(ax2, 'XScale', 'log');
xlim(ax2, x_show);

if isfield(spec, 'diff_ylim') && ~isempty(spec.diff_ylim)
    ylim(ax2, spec.diff_ylim);
else
    local_set_y_margin_from_zero(ax2, diff_for_ylim);
end

box(ax2, 'on');

xlabel(ax2, xlab_text, ...
    'FontSize', fs_x_lab, ...
    'Interpreter', xlab_interpreter);

ylabel(ax2, ylab_diff, ...
    'FontSize', fs_y_lab, ...
    'Interpreter', 'latex');

set(ax2, 'LineWidth', 1.2, ...
    'FontSize', fs_y_ax, ...
    'TickDir', 'in', ...
    'TickLength', [0.015 0.015], ...
    'Layer', 'top');

ax2.XAxis.FontSize = fs_x_ax;
ax2.YAxis.FontSize = fs_y_ax;

local_save_fig_all(f2, sprintf('King_vs_%s_%s_PVLdiff', ref_name, spec.file_tag));
end

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

function mem_mb = local_get_mem_mb()
mem_mb = NaN;

if ispc
    try
        m = memory();
        mem_mb = double(m.MemUsedMATLAB) / (1024^2);
        return;
    catch
    end
end

try
    rt = java.lang.Runtime.getRuntime();
    used = double(rt.totalMemory() - rt.freeMemory());
    mem_mb = used / (1024^2);
catch
    mem_mb = NaN;
end
end

function s = local_fmt_mem(x)
if isnan(x)
    s = 'NaN';
else
    s = sprintf('%.1f MB', x);
end
end

function local_save_fig_all(fig_handle, base_name)
global SAVE_FIGS SAVE_DIR
if isempty(SAVE_FIGS) || ~SAVE_FIGS
    return;
end
if isempty(SAVE_DIR) || ~exist(SAVE_DIR,'dir')
    return;
end

fname = local_sanitize_filename(base_name);

fp_png = fullfile(SAVE_DIR, [fname, '.png']);
fp_pdf = fullfile(SAVE_DIR, [fname, '.pdf']);
fp_fig = fullfile(SAVE_DIR, [fname, '.fig']);

set(fig_handle, 'Color', 'w');
set(fig_handle, 'PaperPositionMode', 'auto');
drawnow;

try
    exportgraphics(fig_handle, fp_png, 'Resolution', 300);
catch
    print(fig_handle, fp_png, '-dpng', '-r300');
end

try
    exportgraphics(fig_handle, fp_pdf, 'ContentType', 'vector');
catch
    print(fig_handle, fp_pdf, '-dpdf', '-painters');
end

try
    savefig(fig_handle, fp_fig);
catch
    warning('Failed to save FIG file: %s', fp_fig);
end
end

function s = local_sanitize_filename(s)
s = char(string(s));
s = strrep(s, ' ', '_');
bad = '<>:"/\|?*';
for k = 1:numel(bad)
    s = strrep(s, bad(k), '_');
end
end

function local_write_runinfo_txt(save_dir, medium, source, calc, ...
    z_target, ds, theta_line_deg, vel_comp_list, ...
    z_show_max, rho_show_min, rho_max_radial_use)

fp = fullfile(save_dir, 'run_info.txt');
fid = fopen(fp, 'w');
if fid < 0
    warning('Cannot create run_info.txt in %s', save_dir);
    return;
end

fprintf(fid, '===== RUN INFO =====\n');
fprintf(fid, 'Time: %s\n\n', datestr(datetime('now'), 'yyyy-mm-dd HH:MM:SS'));
fprintf(fid, 'z_target = %.15g\n', z_target);
fprintf(fid, 'ds = %.15g\n', ds);
fprintf(fid, 'theta_line_deg = %.15g\n', theta_line_deg);
fprintf(fid, 'z_show_max = %.15g\n', z_show_max);
fprintf(fid, 'rho_show_min = %.15g\n', rho_show_min);
fprintf(fid, 'rho_max_radial_use = %.15g\n', rho_max_radial_use);
fprintf(fid, 'vel_comp_list = %s, %s, %s\n', vel_comp_list{:});
fprintf(fid, '\n');

fprintf(fid, 'medium.c0 = %.15g\n', medium.c0);
fprintf(fid, 'medium.rho0 = %.15g\n', medium.rho0);
fprintf(fid, 'medium.pref = %.15g\n', medium.pref);
fprintf(fid, 'v_ref = pref / (rho0*c0) = %.15g\n', medium.pref/(medium.rho0*medium.c0));
fprintf(fid, '\n');

fprintf(fid, 'source.m = %d\n', source.m);
fprintf(fid, 'calc.fht.rho_max = %.15g\n', calc.fht.rho_max);
fprintf(fid, 'calc.fht.zu_max = %.15g\n', calc.fht.zu_max);
fprintf(fid, 'calc.fht.delta = %.15g\n', calc.fht.delta);

if isfield(calc.fht,'z_sampling')
    fprintf(fid, 'calc.fht.z_sampling = %s\n', char(string(calc.fht.z_sampling)));
end
if isfield(calc.fht,'Nz_ultra') && ~isempty(calc.fht.Nz_ultra)
    fprintf(fid, 'calc.fht.Nz_ultra = %.15g\n', calc.fht.Nz_ultra);
end
if isfield(calc.fht,'z_log_min') && ~isempty(calc.fht.z_log_min)
    fprintf(fid, 'calc.fht.z_log_min = %.15g\n', calc.fht.z_log_min);
end

fclose(fid);
end

function tf = local_is_dim_rayleigh(res)
tf = false;
if isfield(res, 'dim') && isfield(res.dim, 'method') && ~isempty(res.dim.method)
    method_str = lower(char(string(res.dim.method)));
    tf = contains(method_str, 'rayleigh');
end
end