%% ============================================================
% King-only PAL audio field with King local-effect correction
%
% Sweep test:
%   m1, m2 = -2:2, total 5 x 5 = 25 cases
%
% Calculation:
%   1) Compute the audio field pa_wo_local(rho,z) using the king method.
%   2) Call calc_ultrasound_field(...,'king') and
%      calc_ultrasound_velocity_field(...,'king')
%      calc_ultrasound_velocity_field(...,king) to compute ultrasonic p/v on the King/FHT rho-z grid.
%   3) Compute the local-effect correction:
%
%        p_loc = rho0/2 * conj(v1)·v2
%              - (w1/w2 + w2/w1 - 1) * conj(p1)p2 / (2 rho0 c^2)
%
%      pa_w_local = pa_wo_local - p_loc
%
%   4) Save compact replot-style PNG images on the z = 1 m plane and the xOz plane.
%
% Saved figures:
%   W_local_z1_SPL.png
%   W_local_z1_phs.png
%   W_local_xOz_SPL.png
%   W_o_local_xOz_SPL.png
%
% Notes:
%   - W_o_local = without local effect
%   - W_local   = with local effect
%   - The local effect uses only King/FHT results and does not use DIM.
%   - Phase figures use colormap('hsv').
%   - Figures have no title inside the axes.
%   - plot_data.mat saves only the data needed for later replotting.
%
% External dependencies:
%   - AbsorpAttenCoef.m
%   - solve_kappa0.m
%   - m_FHT.m
%   - MyColor.m
%   - calc_ultrasound_field.m
%   - calc_ultrasound_velocity_field.m
%   - make_source_velocity.m
%% ============================================================

clear; clc; close all;

%% ===================== Modal sweep settings =====================
m_list = -2:2;     % [-2, -1, 0, 1, 2]

%% ===================== Save and display control =====================
save_figures = true;     % save compact PNG figures
save_results = true;     % save plot_data.mat and summary files
show_figures = false;    % false: close figure windows after each case is saved

save_root = 'AudioTransform_KingLocalEffect_mSweep';
out_root  = fullfile(save_root, 'Replot_original_style_exportgraphics_png');

if save_figures || save_results
    if ~exist(save_root, 'dir')
        mkdir(save_root);
    end
    if ~exist(out_root, 'dir')
        mkdir(out_root);
    end
end

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

%% ===================== FHT / calculation parameters =====================
N_FHT = 16384 / 1;
delta = 0.002 * 1;

rho_max = 0.5;
zu_max  = 12.0;
za_max  = 3.0 + delta;

Nh_scale   = 2.0;
NH_scale   = 4.0;
Nh_v_scale = 1.1;

green_R_min = 1e-12;

%% ===================== Plotting parameters =====================
plot_cfg = struct();

% z = 1 m plane
plot_cfg.z_plane = 1.0;
plot_cfg.xy_max  = 0.5;
plot_cfg.n_xy    = 401;
plot_cfg.mask_radius_z1 = 0.5;

% xOz plane
plot_cfg.xoz_x_max = 0.5;
plot_cfg.xoz_z_max = 3.0;
plot_cfg.max_nx_xoz = 401;
plot_cfg.max_nz_xoz = 401;

% Original SPL limits kept for summary records
plot_cfg.dynamic_range_dB = 30;
plot_cfg.clim_step_dB = 5;

% Compact replot-style display limits
plot_cfg.spl_clim = [-30, 0];
plot_cfg.phase_clim = [-pi, pi];

% Figure size and axes position for compact images
plot_cfg.fig_pos_z1  = [100 80 800 760];
plot_cfg.fig_pos_xoz = [80 80 900 700];
plot_cfg.ax_pos_z1   = [0.08 0.08 0.84 0.84];
plot_cfg.ax_pos_xoz  = [0.08 0.18 0.84 0.64];

% Original full-style parameters kept for compatibility
plot_cfg.xy_ticks = [-0.5, -0.25, 0, 0.25, 0.5];
plot_cfg.xoz_x_ticks = [-0.5, -0.25, 0, 0.25, 0.5];
plot_cfg.font_size_axis  = 20;
plot_cfg.font_size_label = 22;
plot_cfg.font_size_cb    = 20;
plot_cfg.line_width_axis = 2;

%% ===================== Source velocity profile =====================
source_profile = 'Vortex-m';    % 'Uniform' | 'Focus' | 'Vortex-m'
F_focus = 0.2;

%% ===================== Overall summary =====================
summary = struct();
summary.case_id = [];
summary.m1 = [];
summary.m2 = [];
summary.ma = [];
summary.time_all = [];
summary.z_plane_actual = [];
summary.clim_wo_local_z1_low = [];
summary.clim_wo_local_z1_high = [];
summary.clim_w_local_z1_low = [];
summary.clim_w_local_z1_high = [];
summary.clim_wo_local_xOz_low = [];
summary.clim_wo_local_xOz_high = [];
summary.clim_w_local_xOz_low = [];
summary.clim_w_local_xOz_high = [];
summary.save_dir = {};

case_id = 0;
total_cases = numel(m_list)^2;

%% ============================================================
% Sweep loop
%% ============================================================
for i_m1 = 1:numel(m_list)
    for i_m2 = 1:numel(m_list)

        m1 = m_list(i_m1);
        m2 = m_list(i_m2);
        ma = m2 - m1;

        case_id = case_id + 1;

        %% ===================== case save dir =====================
        case_dir_name = sprintf('m1_%+d_m2_%+d', m1, m2);
        case_dir_name = strrep(case_dir_name, '+', 'p');
        case_dir_name = strrep(case_dir_name, '-', 'm');

        save_dir = fullfile(save_root, case_dir_name);
        fig_out_dir = fullfile(out_root, case_dir_name);

        if save_figures || save_results
            if ~exist(save_dir, 'dir')
                mkdir(save_dir);
            end
            if ~exist(fig_out_dir, 'dir')
                mkdir(fig_out_dir);
            end
        end

        fprintf('\n============================================================\n');
        fprintf('Case %d / %d\n', case_id, total_cases);
        fprintf('m1 = %d, m2 = %d, ma = %d\n', m1, m2, ma);
        fprintf('Save dir: %s\n', save_dir);
        fprintf('Figure dir: %s\n', fig_out_dir);
        fprintf('============================================================\n');

        %% ===================== King local-effect calculation parameters =====================
        calc_king = struct();

        calc_king.fht = struct();
        calc_king.fht.N_FHT = N_FHT;
        calc_king.fht.rho_max = rho_max;
        calc_king.fht.Nh_scale = Nh_scale;
        calc_king.fht.NH_scale = NH_scale;
        calc_king.fht.Nh_v_scale = Nh_v_scale;
        calc_king.fht.zu_max = zu_max;
        calc_king.fht.za_max = za_max;
        calc_king.fht.delta = delta;

        % The dim field is kept only for compatibility with existing function structures; this script does not call DIM.
        calc_king.dim = struct();
        calc_king.dim.method = 'rayleigh';
        calc_king.dim.use_parallel = false;
        calc_king.dim.num_workers = 1;
        calc_king.dim.use_freq = 'f2';
        calc_king.dim.dis_coe = 32;
        calc_king.dim.margin = 1;
        calc_king.dim.src_discretization = 'polar';

        calc_king.dim.grid_mode = 'uniform';
        calc_king.dim.ds_rho = 32;
        calc_king.dim.ds_rho_src = 16;
        calc_king.dim.ds_rho_obs = 32;
        calc_king.dim.uniform_dx = delta / 2;
        calc_king.dim.uniform_dz = delta / 1;

        calc_king.king = struct();
        calc_king.king.gspec_method = 'transform';
        calc_king.king.eps_kzz = 1e-3;
        calc_king.king.eps_phase = calc_king.king.eps_kzz;
        calc_king.king.kz_min = 1e-12;
        calc_king.king.band_refine.enable = false;

        calc_king.asm = struct();
        calc_king.asm.pad_factor = 16;
        calc_king.asm.kzz_eps = 1e-12;

        calc_king.internal = struct();

        source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2, source_profile, F_focus);
        medium_cfg = build_medium_cfg(c, rho0, beta, pref);

        %% ============================================================
        % 1) Compute the wo-local audio field using the transform method
        %% ============================================================
        t_all = tic;

        out_audio = calc_audio_transform_wo_local( ...
            a, v0, c, rho0, beta, ...
            f1, f2, fa, ...
            m1, m2, ...
            source_profile, F_focus, ...
            N_FHT, delta, ...
            rho_max, zu_max, za_max, ...
            Nh_scale, NH_scale, Nh_v_scale, ...
            green_R_min);

        rho_audio = out_audio.rho_audio(:);
        z_audio   = out_audio.z_audio(:).';
        pa_wo_local_rz_coeff = out_audio.pa_audio;

        fprintf('\nTransform wo-local calculation finished.\n');
        fprintf('Audio coefficient field size: [%d, %d]\n', ...
            size(pa_wo_local_rz_coeff,1), size(pa_wo_local_rz_coeff,2));

        %% ============================================================
        % 2) Compute the local-effect term p_loc(rho,z) using King/FHT
        %% ============================================================
        fprintf('\n============================================================\n');
        fprintf('Computing King local-effect term on rho-z grid...\n');
        fprintf('============================================================\n');

        [p_loc_king, rho_king, z_king, local_debug] = compute_local_effect_king_rz( ...
            source_cfg, medium_cfg, calc_king, ...
            f1, f2, m1, m2, ...
            rho0, c);

        fprintf('King local-effect field size: [%d, %d]\n', size(p_loc_king,1), size(p_loc_king,2));

        %% ============================================================
        % 3) Align p_loc to the audio rho-z grid
        %% ============================================================
        if isequal(size(p_loc_king), size(pa_wo_local_rz_coeff)) && ...
           numel(rho_king) == numel(rho_audio) && ...
           numel(z_king) == numel(z_audio) && ...
           max(abs(rho_king(:) - rho_audio(:))) < 1e-12 && ...
           max(abs(z_king(:)   - z_audio(:)))   < 1e-12

            fprintf('King local grid matches audio grid exactly.\n');
            p_loc_rz_coeff = p_loc_king;

        else
            fprintf('Interpolating King local-effect term onto audio grid...\n');

            [RHO_A, Z_A] = ndgrid(rho_audio, z_audio);

            p_loc_rz_coeff = interp_complex_rz( ...
                rho_king, z_king, p_loc_king, RHO_A, Z_A);

            p_loc_rz_coeff(~isfinite(real(p_loc_rz_coeff))) = 0;
            p_loc_rz_coeff(~isfinite(imag(p_loc_rz_coeff))) = 0;
        end

        pa_w_local_rz_coeff = pa_wo_local_rz_coeff - p_loc_rz_coeff;

        %% ============================================================
        % 4) z = 1 m plane: wo-local and w-local
        %% ============================================================
        [~, iz_plane] = min(abs(z_audio - plot_cfg.z_plane));
        z_plane_actual = z_audio(iz_plane);

        x_plane = linspace(-plot_cfg.xy_max, plot_cfg.xy_max, plot_cfg.n_xy);
        y_plane = linspace(-plot_cfg.xy_max, plot_cfg.xy_max, plot_cfg.n_xy);

        [Xp, Yp] = meshgrid(x_plane, y_plane);
        RHO_p = sqrt(Xp.^2 + Yp.^2);
        PHI_p = atan2(Yp, Xp);

        pa_wo_local_z1_coeff = pa_wo_local_rz_coeff(:, iz_plane);
        pa_w_local_z1_coeff  = pa_w_local_rz_coeff(:, iz_plane);

        pa_wo_local_z1_base = interp1( ...
            rho_audio, pa_wo_local_z1_coeff, RHO_p, 'linear', NaN);

        pa_w_local_z1_base = interp1( ...
            rho_audio, pa_w_local_z1_coeff, RHO_p, 'linear', NaN);

        pa_wo_local_z1_base(RHO_p > max(rho_audio)) = NaN;
        pa_w_local_z1_base(RHO_p > max(rho_audio))  = NaN;

        % Full azimuthal phase exp(j ma phi)
        pa_wo_local_z1 = pa_wo_local_z1_base .* exp(1j * ma * PHI_p);
        pa_w_local_z1  = pa_w_local_z1_base  .* exp(1j * ma * PHI_p);

        %% ============================================================
        % 5) xOz plane: wo-local and w-local
        %% ============================================================
        rho_xoz_max = min(plot_cfg.xoz_x_max, max(rho_audio));
        z_xoz_max   = min(plot_cfg.xoz_z_max, max(z_audio));

        idx_r_max = find(rho_audio <= rho_xoz_max, 1, 'last');
        idx_z_max = find(z_audio   <= z_xoz_max,   1, 'last');

        rho_use = rho_audio(1:idx_r_max);
        z_use   = z_audio(1:idx_z_max);

        pa_wo_xoz_coeff_use = pa_wo_local_rz_coeff(1:idx_r_max, 1:idx_z_max);
        pa_w_xoz_coeff_use  = pa_w_local_rz_coeff(1:idx_r_max,  1:idx_z_max);

        [x_full, pa_wo_xoz_full] = build_xoz_complex_map_from_rz( ...
            rho_use, pa_wo_xoz_coeff_use, ma);

        [~, pa_w_xoz_full] = build_xoz_complex_map_from_rz( ...
            rho_use, pa_w_xoz_coeff_use, ma);

        idx_x_ds = make_ds_index(numel(x_full), plot_cfg.max_nx_xoz);
        idx_z_ds = make_ds_index(numel(z_use),  plot_cfg.max_nz_xoz);

        x_xoz = x_full(idx_x_ds);
        z_xoz = z_use(idx_z_ds);

        pa_wo_local_xOz = pa_wo_xoz_full(idx_x_ds, idx_z_ds);
        pa_w_local_xOz  = pa_w_xoz_full(idx_x_ds,  idx_z_ds);

        time_all = toc(t_all);

        fprintf('\nTotal calculation finished. elapsed = %.2f s\n', time_all);

        %% ============================================================
        % 6) SPL / phase
        %% ============================================================
        spl_wo_local_z1 = pressure_to_spl(abs(pa_wo_local_z1), pref);
        spl_w_local_z1  = pressure_to_spl(abs(pa_w_local_z1),  pref);

        phs_wo_local_z1 = angle(pa_wo_local_z1);
        phs_w_local_z1  = angle(pa_w_local_z1);

        spl_wo_local_xOz = pressure_to_spl(abs(pa_wo_local_xOz), pref);
        spl_w_local_xOz  = pressure_to_spl(abs(pa_w_local_xOz),  pref);

        phs_wo_local_xOz = angle(pa_wo_local_xOz);
        phs_w_local_xOz  = angle(pa_w_local_xOz);

        clim_wo_local_z1 = auto_clim_ceil_step( ...
            spl_wo_local_z1, plot_cfg.dynamic_range_dB, plot_cfg.clim_step_dB);

        clim_w_local_z1 = auto_clim_ceil_step( ...
            spl_w_local_z1, plot_cfg.dynamic_range_dB, plot_cfg.clim_step_dB);

        clim_wo_local_xOz = auto_clim_ceil_step( ...
            spl_wo_local_xOz, plot_cfg.dynamic_range_dB, plot_cfg.clim_step_dB);

        clim_w_local_xOz = auto_clim_ceil_step( ...
            spl_w_local_xOz, plot_cfg.dynamic_range_dB, plot_cfg.clim_step_dB);

        fprintf('\nColor limits:\n');
        fprintf('  W/o local z1  SPL clim = [%.1f, %.1f] dB\n', clim_wo_local_z1(1), clim_wo_local_z1(2));
        fprintf('  W   local z1  SPL clim = [%.1f, %.1f] dB\n', clim_w_local_z1(1),  clim_w_local_z1(2));
        fprintf('  W/o local xOz SPL clim = [%.1f, %.1f] dB\n', clim_wo_local_xOz(1), clim_wo_local_xOz(2));
        fprintf('  W   local xOz SPL clim = [%.1f, %.1f] dB\n', clim_w_local_xOz(1),  clim_w_local_xOz(2));

        %% ============================================================
        % 7) Build compact replot-style figures after saving plot data
        %% ============================================================

        %% ============================================================
        % 9) Save lightweight plotting data
        %% ============================================================
        if save_results
            plot_data = struct();

            % ===== case info =====
            plot_data.m1 = m1;
            plot_data.m2 = m2;
            plot_data.ma = ma;

            % ===== coordinates for z = 1 m plane =====
            plot_data.x_plane = x_plane;
            plot_data.y_plane = y_plane;
            plot_data.z_plane_actual = z_plane_actual;

            % ===== coordinates for xOz plane =====
            plot_data.x_xoz = x_xoz;
            plot_data.z_xoz = z_xoz;

            % ===== SPL data for plotting =====
            plot_data.spl_wo_local_z1  = spl_wo_local_z1;
            plot_data.spl_w_local_z1   = spl_w_local_z1;
            plot_data.spl_wo_local_xOz = spl_wo_local_xOz;
            plot_data.spl_w_local_xOz  = spl_w_local_xOz;

            % ===== phase data for plotting =====
            plot_data.phs_wo_local_z1  = phs_wo_local_z1;
            plot_data.phs_w_local_z1   = phs_w_local_z1;
            plot_data.phs_wo_local_xOz = phs_wo_local_xOz;
            plot_data.phs_w_local_xOz  = phs_w_local_xOz;

            % ===== color limits =====
            plot_data.clim_wo_local_z1  = clim_wo_local_z1;
            plot_data.clim_w_local_z1   = clim_w_local_z1;
            plot_data.clim_wo_local_xOz = clim_wo_local_xOz;
            plot_data.clim_w_local_xOz  = clim_w_local_xOz;

            % ===== plotting config =====
            plot_data.plot_cfg = plot_cfg;

            % ===== physical and numerical parameters, for record only =====
            plot_data.param = struct();
            plot_data.param.a = a;
            plot_data.param.v0 = v0;
            plot_data.param.c = c;
            plot_data.param.rho0 = rho0;
            plot_data.param.beta = beta;
            plot_data.param.pref = pref;
            plot_data.param.fu = fu;
            plot_data.param.fa = fa;
            plot_data.param.f1 = f1;
            plot_data.param.f2 = f2;

            plot_data.param.N_FHT = N_FHT;
            plot_data.param.delta = delta;
            plot_data.param.rho_max = rho_max;
            plot_data.param.zu_max = zu_max;
            plot_data.param.za_max = za_max;
            plot_data.param.Nh_scale = Nh_scale;
            plot_data.param.NH_scale = NH_scale;
            plot_data.param.Nh_v_scale = Nh_v_scale;

            plot_data.param.source_profile = source_profile;
            plot_data.param.F_focus = F_focus;
            plot_data.param.time_all = time_all;

            save(fullfile(save_dir, 'plot_data.mat'), 'plot_data');
        end

        %% ============================================================
        % 8) Save compact replot-style PNG figures
        %% ============================================================
        [Xp_plot, Yp_plot] = meshgrid(x_plane, y_plane);
        mask_z1 = sqrt(Xp_plot.^2 + Yp_plot.^2) <= plot_cfg.mask_radius_z1;

        spl_w_local_z1_plot = spl_w_local_z1;
        spl_w_local_z1_plot(~mask_z1) = NaN;
        spl_w_local_z1_plot = normalize_spl_to_0dB(spl_w_local_z1_plot);

        [fig1, ax1] = plot_z1_spl_no_axis( ...
            x_plane, y_plane, spl_w_local_z1_plot, plot_cfg, ...
            'W local z1 SPL', show_figures);
        save_png_if_needed(ax1, save_figures, fullfile(fig_out_dir, 'W_local_z1_SPL.png'));

        phs_w_local_z1_plot = phs_w_local_z1;
        phs_w_local_z1_plot(~mask_z1) = NaN;

        [fig2, ax2] = plot_z1_phase_no_axis( ...
            x_plane, y_plane, phs_w_local_z1_plot, plot_cfg, ...
            'W local z1 phase', show_figures);
        save_png_if_needed(ax2, save_figures, fullfile(fig_out_dir, 'W_local_z1_phs.png'));

        spl_w_local_xOz_plot = normalize_spl_to_0dB(spl_w_local_xOz);

        [fig3, ax3] = plot_xOz_spl_no_axis( ...
            z_xoz, x_xoz, spl_w_local_xOz_plot, plot_cfg, ...
            'W local xOz SPL', show_figures);
        save_png_if_needed(ax3, save_figures, fullfile(fig_out_dir, 'W_local_xOz_SPL.png'));

        spl_wo_local_xOz_plot = normalize_spl_to_0dB(spl_wo_local_xOz);

        [fig4, ax4] = plot_xOz_spl_no_axis( ...
            z_xoz, x_xoz, spl_wo_local_xOz_plot, plot_cfg, ...
            'W/o local xOz SPL', show_figures);
        save_png_if_needed(ax4, save_figures, fullfile(fig_out_dir, 'W_o_local_xOz_SPL.png'));

        %% ============================================================
        % summary
        %% ============================================================
        summary.case_id(end+1,1) = case_id;
        summary.m1(end+1,1) = m1;
        summary.m2(end+1,1) = m2;
        summary.ma(end+1,1) = ma;
        summary.time_all(end+1,1) = time_all;
        summary.z_plane_actual(end+1,1) = z_plane_actual;

        summary.clim_wo_local_z1_low(end+1,1) = clim_wo_local_z1(1);
        summary.clim_wo_local_z1_high(end+1,1) = clim_wo_local_z1(2);
        summary.clim_w_local_z1_low(end+1,1) = clim_w_local_z1(1);
        summary.clim_w_local_z1_high(end+1,1) = clim_w_local_z1(2);

        summary.clim_wo_local_xOz_low(end+1,1) = clim_wo_local_xOz(1);
        summary.clim_wo_local_xOz_high(end+1,1) = clim_wo_local_xOz(2);
        summary.clim_w_local_xOz_low(end+1,1) = clim_w_local_xOz(1);
        summary.clim_w_local_xOz_high(end+1,1) = clim_w_local_xOz(2);

        summary.save_dir{end+1,1} = save_dir;

        if ~show_figures
            close(fig1); close(fig2); close(fig3); close(fig4);
        end

        %% ===================== clear large variables =====================
        clear out_audio local_debug
        clear p_loc_king p_loc_rz_coeff pa_w_local_rz_coeff pa_wo_local_rz_coeff
        clear pa_wo_local_z1 pa_w_local_z1 pa_wo_local_xOz pa_w_local_xOz
        clear spl_wo_local_z1 spl_w_local_z1 phs_wo_local_z1 phs_w_local_z1
        clear spl_wo_local_xOz spl_w_local_xOz phs_wo_local_xOz phs_w_local_xOz
        clear pa_wo_local_z1_base pa_w_local_z1_base
        clear pa_wo_xoz_full pa_w_xoz_full
        clear p_loc_king rho_king z_king
        clear Xp Yp RHO_p PHI_p
        clear Xp_plot Yp_plot mask_z1
        clear spl_w_local_z1_plot phs_w_local_z1_plot
        clear spl_w_local_xOz_plot spl_wo_local_xOz_plot

        fprintf('\nCase %d / %d done.\n', case_id, total_cases);
    end
end

%% ============================================================
% Save overall summary
%% ============================================================
if save_results
    save(fullfile(save_root, 'summary_all_cases.mat'), 'summary');

    T_summary = table( ...
        summary.case_id, ...
        summary.m1, ...
        summary.m2, ...
        summary.ma, ...
        summary.time_all, ...
        summary.z_plane_actual, ...
        summary.clim_wo_local_z1_low, ...
        summary.clim_wo_local_z1_high, ...
        summary.clim_w_local_z1_low, ...
        summary.clim_w_local_z1_high, ...
        summary.clim_wo_local_xOz_low, ...
        summary.clim_wo_local_xOz_high, ...
        summary.clim_w_local_xOz_low, ...
        summary.clim_w_local_xOz_high, ...
        summary.save_dir, ...
        'VariableNames', { ...
        'case_id', 'm1', 'm2', 'ma', 'time_all', 'z_plane_actual', ...
        'clim_wo_local_z1_low', 'clim_wo_local_z1_high', ...
        'clim_w_local_z1_low', 'clim_w_local_z1_high', ...
        'clim_wo_local_xOz_low', 'clim_wo_local_xOz_high', ...
        'clim_w_local_xOz_low', 'clim_w_local_xOz_high', ...
        'save_dir'});

    writetable(T_summary, fullfile(save_root, 'summary_all_cases.csv'));
end

fprintf('\n============================================================\n');
fprintf('All cases done.\n');
fprintf('Total cases = %d\n', total_cases);
fprintf('Root save folder = %s\n', save_root);
fprintf('============================================================\n');

%% ============================================================
% Local functions
%% ============================================================

function source_cfg = build_source_cfg(a, v0, m1, m2, f1, fa, f2, source_profile, F_focus)

source_cfg = struct();

source_cfg.profile = source_profile;
source_cfg.a = a;
source_cfg.v0 = v0;
source_cfg.v_ratio = 1;

source_cfg.m1 = m1;
source_cfg.m2 = m2;
source_cfg.m  = m1;

source_cfg.F = F_focus;

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

function out = calc_audio_transform_wo_local( ...
    a, v0, c, rho0, beta, ...
    f1, f2, fa, ...
    m1, m2, ...
    source_profile, F_focus, ...
    N_FHT, delta, ...
    rho_max, zu_max, za_max, ...
    Nh_scale, NH_scale, Nh_v_scale, ...
    green_R_min)

w1 = 2*pi*f1;
w2 = 2*pi*f2;
wa = 2*pi*fa;

k1 = w1/c + 1j*AbsorpAttenCoef(f1);
k2 = w2/c + 1j*AbsorpAttenCoef(f2);
ka = wa/c + 1j*AbsorpAttenCoef(fa);

ma = m2 - m1;

Nh = Nh_scale * rho_max;
NH = NH_scale * w2 / c;

n_FHT = 0:N_FHT-1;
[a_solve, k0, x1, x0] = solve_kappa0(N_FHT, n_FHT);

rho_audio = (x1 * Nh).';
z_ultra   = 0:delta:zu_max;
z_audio   = 0:delta:za_max;

Nz_ultra = numel(z_ultra);

fprintf('\n[FHT grid]\n');
fprintf('N_FHT = %d\n', N_FHT);
fprintf('delta = %.6g m\n', delta);
fprintf('rho points = %d\n', numel(rho_audio));
fprintf('z_ultra points = %d\n', numel(z_ultra));
fprintf('z_audio points = %d\n', numel(z_audio));
fprintf('Nh = %.6g, NH = %.6g\n', Nh, NH);

%% -------------------- source velocity spectra --------------------
fprintf('\nBuilding source velocity spectra...\n');

Nh_v = Nh_v_scale * a;
NH_v = NH;

rho_v_grid = (x1 * Nh_v).';

switch source_profile
    case 'Uniform'
        vs1 = v0 * double(rho_v_grid <= a);
        vs2 = vs1;

    case 'Vortex-m'
        vs1 = v0 * double(rho_v_grid <= a);
        vs2 = vs1;

    case 'Focus'
        vs1 = v0 * exp(-1j*real(k1)*sqrt(rho_v_grid.^2 + F_focus^2)) ...
            .* double(rho_v_grid <= a);

        vs2 = v0 * exp(-1j*real(k2)*sqrt(rho_v_grid.^2 + F_focus^2)) ...
            .* double(rho_v_grid <= a);

    otherwise
        error('Unknown source_profile: %s', source_profile);
end

Vs1 = m_FHT( ...
    vs1, ...
    N_FHT, 1, ...
    Nh_v, NH_v, ...
    a_solve, x0, x1, k0, ...
    m1);

Vs2 = m_FHT( ...
    vs2, ...
    N_FHT, 1, ...
    Nh_v, NH_v, ...
    a_solve, x0, x1, k0, ...
    m2);

%% -------------------- ultrasonic Green transforms --------------------
fprintf('\nComputing ultrasonic Green transforms...\n');

G1_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, k1, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G2_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, k2, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

G1 = (-4*pi*1j) * G1_raw;
G2 = (-4*pi*1j) * G2_raw;

%% -------------------- ultrasonic pressure fields --------------------
fprintf('\nComputing ultrasonic pressure fields...\n');

p1_transform = compute_ultrasound_pressure_from_G( ...
    G1, Vs1, ...
    N_FHT, Nz_ultra, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    m1, rho0, c, k1);

p2_transform = compute_ultrasound_pressure_from_G( ...
    G2, Vs2, ...
    N_FHT, Nz_ultra, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    m2, rho0, c, k2);

%% -------------------- nonlinear source term --------------------
fprintf('\nComputing nonlinear source term q...\n');

q_transform = conj(p1_transform) .* p2_transform ...
    * beta*wa/(1j*rho0^2*c^4);

mask_rho = double(rho_audio(:) <= rho_max);
q_transform = q_transform .* repmat(mask_rho, 1, size(q_transform, 2));

%% -------------------- audio Green transform --------------------
fprintf('\nComputing audio Green transform...\n');

Ga_raw = build_green_space_g_transform_raw( ...
    rho_audio, z_ultra, ka, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0);

Ga = (-4*pi*1j) * Ga_raw;

%% -------------------- audio pressure field without local effect --------------------
fprintf('\nComputing audio pressure field without local effect...\n');

[pa_audio, phia_audio] = compute_paW_from_Gr00( ...
    q_transform, z_ultra, z_audio, Ga, ...
    N_FHT, Nh, NH, ...
    a_solve, x0, x1, k0, ...
    ma, delta, rho0, wa);

out = struct();

out.rho_audio = rho_audio(:);
out.z_ultra   = z_ultra(:);
out.z_audio   = z_audio(:);

out.p1_transform = p1_transform;
out.p2_transform = p2_transform;
out.q_transform  = q_transform;

out.pa_audio   = pa_audio;
out.phia_audio = phia_audio;

out.k1 = k1;
out.k2 = k2;
out.ka = ka;

out.w1 = w1;
out.w2 = w2;
out.wa = wa;

out.m1 = m1;
out.m2 = m2;
out.ma = ma;
end

function [p_loc_king, rho_king, z_king, dbg] = compute_local_effect_king_rz( ...
    source_cfg, medium_cfg, calc_cfg_king, ...
    f1, f2, m1, m2, ...
    rho0, c)

w1 = 2*pi*f1;
w2 = 2*pi*f2;

obs_grid = struct();

%% -------------------- f1 / m1, King only --------------------
source_cfg_f1 = source_cfg;
source_cfg_f1.m  = m1;
source_cfg_f1.m1 = m1;
source_cfg_f1.m2 = m1;

calc_cfg_f1 = calc_cfg_king;
calc_cfg_f1.dim.use_freq = 'f1';

fprintf('  Computing KING ultrasound p/v for f1 with m = %d ...\n', m1);

res_p_f1 = calc_ultrasound_field( ...
    source_cfg_f1, medium_cfg, calc_cfg_f1, obs_grid, 'king');

res_v_f1 = calc_ultrasound_velocity_field( ...
    source_cfg_f1, medium_cfg, calc_cfg_f1, obs_grid, 'king');

%% -------------------- f2 / m2, King only --------------------
source_cfg_f2 = source_cfg;
source_cfg_f2.m  = m2;
source_cfg_f2.m1 = m2;
source_cfg_f2.m2 = m2;

calc_cfg_f2 = calc_cfg_king;
calc_cfg_f2.dim.use_freq = 'f2';

fprintf('  Computing KING ultrasound p/v for f2 with m = %d ...\n', m2);

res_p_f2 = calc_ultrasound_field( ...
    source_cfg_f2, medium_cfg, calc_cfg_f2, obs_grid, 'king');

res_v_f2 = calc_ultrasound_velocity_field( ...
    source_cfg_f2, medium_cfg, calc_cfg_f2, obs_grid, 'king');

%% -------------------- extract King fields --------------------
rho_king = res_p_f1.king.rho(:);
z_king   = res_p_f1.king.z(:).';

p1 = extract_king_field_as_rz(res_p_f1.king, {'p_f1','p1','pressure_f1'}, rho_king, z_king);
p2 = extract_king_field_as_rz(res_p_f2.king, {'p_f2','p2','pressure_f2'}, rho_king, z_king);

v1r = extract_king_field_as_rz(res_v_f1.king, {'v_rho_f1','vr_f1','v_r_f1'}, rho_king, z_king);
v1p = extract_king_field_as_rz(res_v_f1.king, {'v_phi_f1','vphi_f1'}, rho_king, z_king);
v1z = extract_king_field_as_rz(res_v_f1.king, {'v_z_f1','vz_f1','uz_f1'}, rho_king, z_king);

v2r = extract_king_field_as_rz(res_v_f2.king, {'v_rho_f2','vr_f2','v_r_f2'}, rho_king, z_king);
v2p = extract_king_field_as_rz(res_v_f2.king, {'v_phi_f2','vphi_f2'}, rho_king, z_king);
v2z = extract_king_field_as_rz(res_v_f2.king, {'v_z_f2','vz_f2','uz_f2'}, rho_king, z_king);

%% -------------------- local-effect term on King rho-z grid --------------------
coef = (w1/w2 + w2/w1 - 1);

vv = conj(v1r).*v2r + conj(v1p).*v2p + conj(v1z).*v2z;
pp = conj(p1).*p2;

p_loc_king = rho0/2 * vv - coef * pp ./ (2*rho0*c^2);

dbg = struct();
dbg.p1 = p1;
dbg.p2 = p2;
dbg.v1r = v1r;
dbg.v1p = v1p;
dbg.v1z = v1z;
dbg.v2r = v2r;
dbg.v2p = v2p;
dbg.v2z = v2z;
end

function A = extract_king_field_as_rz(king_struct, names, rho_vec, z_vec)

A = get_field_with_aliases(king_struct, names);
A = squeeze(A);

Nr = numel(rho_vec);
Nz = numel(z_vec);

if isequal(size(A), [Nr, Nz])
    return;
end

if isequal(size(A), [Nz, Nr])
    A = A.';
    return;
end

if isvector(A)
    if numel(A) == Nr && Nz == 1
        A = A(:);
        return;
    elseif numel(A) == Nz && Nr == 1
        A = A(:).';
        return;
    end
end

error('King field size mismatch. Field size = [%s], expected [%d, %d].', ...
    num2str(size(A)), Nr, Nz);
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

function Vq = interp_complex_rz(rho_vec, z_vec, V, RHO_q, Z_q)

rho_vec = rho_vec(:);
z_vec = z_vec(:).';

Fr = griddedInterpolant({rho_vec, z_vec}, real(V), 'linear', 'none');
Fi = griddedInterpolant({rho_vec, z_vec}, imag(V), 'linear', 'none');

Vq = Fr(RHO_q, Z_q) + 1j * Fi(RHO_q, Z_q);
end

function p_out = compute_ultrasound_pressure_from_G( ...
    G_spec, Vs, ...
    N_FHT, Nz, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    m_use, rho0, c0, k_use)

F = G_spec .* Vs;

phi = -1j * m_FHT( ...
    F, ...
    N_FHT, Nz, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    m_use);

p_out = 1j * rho0 * c0 * real(k_use) .* phi;
end

function G_raw = build_green_space_g_transform_raw( ...
    rho_vec, z_vec, k_use, green_R_min, ...
    N_FHT, Nh, NH, a_solve, x0, x1, k0)

[RHO, Z] = ndgrid(rho_vec, z_vec);

RR = sqrt(RHO.^2 + Z.^2);
RR_use = max(RR, green_R_min);

g_space = exp(1j * k_use * RR_use) ./ (4*pi * RR_use);

G_raw = m_FHT( ...
    g_space, ...
    N_FHT, numel(z_vec), ...
    Nh, NH, ...
    a_solve, x0, x1, k0, ...
    0);
end

function [pa_W, phia_W] = compute_paW_from_Gr00( ...
    q_full, z, z_audio, Gr00, ...
    N_FHT, Nh, NH, ...
    a_solve, x0, x1, k0, ...
    ma, delta, rho0, wa)

Nz  = numel(z);
Nza = numel(z_audio);

absz = [-fliplr(z(2:end)) z];

Nz1 = numel(absz);
N_conv = Nz1 + Nza - 1;

Qr00 = m_FHT( ...
    q_full, ...
    N_FHT, Nz, ...
    Nh, NH, ...
    a_solve, x0, x1, k0, ...
    ma);

Qr0 = [fliplr(Qr00(:,2:end)) Qr00];

[M_0, N_0] = size(Qr0);

Qr = [Qr0, zeros(M_0, N_conv-N_0)];
Q  = (fft(Qr.')).';

Gr0 = [fliplr(Gr00(:,2:end)) Gr00];
Gr  = [Gr0, zeros(M_0, N_conv-N_0)];
G   = (fft(Gr.')).';

Pa = Q .* G;

par0 = (ifft(Pa.')).';
par  = par0(:, N_conv-Nza+1:N_conv);

phia0 = m_FHT( ...
    par, ...
    N_FHT, Nza, ...
    NH, Nh, ...
    a_solve, x0, x1, k0, ...
    ma);

phia_W = -phia0 * delta * 1j / 2;
pa_W   = 1j * rho0 * wa * phia_W;
end

function [x_full, p_xoz_full] = build_xoz_complex_map_from_rz(rho_vec, p_rz_coeff, ma)

rho_vec = rho_vec(:);

if abs(rho_vec(1)) < 1e-15
    x_full = [flipud(-rho_vec(2:end)); rho_vec];

    p_neg = flipud(p_rz_coeff(2:end,:)) .* exp(1j * ma * pi);
    p_pos = p_rz_coeff;

    p_xoz_full = [p_neg; p_pos];
else
    x_full = [flipud(-rho_vec); rho_vec];

    p_neg = flipud(p_rz_coeff) .* exp(1j * ma * pi);
    p_pos = p_rz_coeff;

    p_xoz_full = [p_neg; p_pos];
end
end

function spl = pressure_to_spl(p_amp, pref)

spl = 20 * log10(p_amp ./ max(pref, eps) ./ sqrt(2) + eps);
end

function clim_use = auto_clim_ceil_step(SPL_data, dynamic_range_dB, step_dB)

valid_data = SPL_data(isfinite(SPL_data));

if isempty(valid_data)
    clim_use = [0 dynamic_range_dB];
    return;
end

max_val = max(valid_data);

clim_high = ceil(max_val / step_dB) * step_dB;
clim_low  = clim_high - dynamic_range_dB;

clim_use = [clim_low, clim_high];
end

function idx = make_ds_index(n, n_keep)

if n <= n_keep
    idx = 1:n;
else
    idx = unique(round(linspace(1, n, n_keep)));
end

idx = idx(:);
end

function fig = plot_z1_spl_figure(x, y, spl_data, clim_use, plot_cfg, fig_name)

fig = figure( ...
    'Name', fig_name, ...
    'Position', [100 80 800 760], ...
    'Color', 'w');

pcolor(x, y, spl_data);
shading flat;
colormap(MyColor('vik'));

cb = colorbar;
cb.Title.Interpreter = 'none';
cb.Title.String = 'SPL (dB)';
set(cb, 'FontSize', plot_cfg.font_size_cb);

caxis(clim_use);

xlim([-plot_cfg.xy_max plot_cfg.xy_max]);
ylim([-plot_cfg.xy_max plot_cfg.xy_max]);

xticks(plot_cfg.xy_ticks);
yticks(plot_cfg.xy_ticks);

xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);

daspect([1 1 1]);

set(gca, ...
    'LineWidth', plot_cfg.line_width_axis, ...
    'FontSize', plot_cfg.font_size_axis, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');
end

function fig = plot_z1_phase_figure(x, y, phase_data, plot_cfg, fig_name)

fig = figure( ...
    'Name', fig_name, ...
    'Position', [140 100 800 760], ...
    'Color', 'w');

pcolor(x, y, phase_data);
shading flat;
colormap('hsv');

cb = colorbar;
cb.Title.Interpreter = 'latex';
cb.Title.String = 'Phase (rad)';
set(cb, 'FontSize', plot_cfg.font_size_cb);

caxis(plot_cfg.phase_clim);

cb.Ticks = [-pi, -pi/2, 0, pi/2, pi];
cb.TickLabels = {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'};
cb.TickLabelInterpreter = 'latex';

xlim([-plot_cfg.xy_max plot_cfg.xy_max]);
ylim([-plot_cfg.xy_max plot_cfg.xy_max]);

xticks(plot_cfg.xy_ticks);
yticks(plot_cfg.xy_ticks);

xlabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);
ylabel('$y$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);

daspect([1 1 1]);

set(gca, ...
    'LineWidth', plot_cfg.line_width_axis, ...
    'FontSize', plot_cfg.font_size_axis, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');
end

function fig = plot_xOz_spl_figure(z, x, spl_data, clim_use, plot_cfg, fig_name)

fig = figure( ...
    'Name', fig_name, ...
    'Position', [80 80 900 700], ...
    'Color', 'w');

pcolor(z, x, spl_data);
shading flat;
colormap(MyColor('vik'));

cb = colorbar;
cb.Title.Interpreter = 'none';
cb.Title.String = 'SPL (dB)';
set(cb, 'FontSize', plot_cfg.font_size_cb);

caxis(clim_use);

xlim([0 plot_cfg.xoz_z_max]);
ylim([-plot_cfg.xoz_x_max plot_cfg.xoz_x_max]);

xticks(0:0.5:plot_cfg.xoz_z_max);
yticks(plot_cfg.xoz_x_ticks);

xlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);
ylabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);

daspect([1 1 1]);

set(gca, ...
    'LineWidth', plot_cfg.line_width_axis, ...
    'FontSize', plot_cfg.font_size_axis, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');
end

function fig = plot_xOz_phase_figure(z, x, phase_data, plot_cfg, fig_name)

fig = figure( ...
    'Name', fig_name, ...
    'Position', [120 100 900 700], ...
    'Color', 'w');

pcolor(z, x, phase_data);
shading flat;
colormap('hsv');

cb = colorbar;
cb.Title.Interpreter = 'latex';
cb.Title.String = 'Phase (rad)';
set(cb, 'FontSize', plot_cfg.font_size_cb);

caxis(plot_cfg.phase_clim);

cb.Ticks = [-pi, -pi/2, 0, pi/2, pi];
cb.TickLabels = {'$-\pi$', '$-\pi/2$', '$0$', '$\pi/2$', '$\pi$'};
cb.TickLabelInterpreter = 'latex';

xlim([0 plot_cfg.xoz_z_max]);
ylim([-plot_cfg.xoz_x_max plot_cfg.xoz_x_max]);

xticks(0:0.5:plot_cfg.xoz_z_max);
yticks(plot_cfg.xoz_x_ticks);

xlabel('$z$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);
ylabel('$x$ (m)', 'Interpreter', 'latex', 'FontSize', plot_cfg.font_size_label);

daspect([1 1 1]);

set(gca, ...
    'LineWidth', plot_cfg.line_width_axis, ...
    'FontSize', plot_cfg.font_size_axis, ...
    'TickLabelInterpreter', 'latex', ...
    'Box', 'on');
end


function spl_norm = normalize_spl_to_0dB(spl_data)

valid_data = spl_data(isfinite(spl_data));

if isempty(valid_data)
    spl_norm = spl_data;
    return;
end

max_val = max(valid_data(:));
spl_norm = spl_data - max_val;
end

function [fig, ax] = plot_z1_spl_no_axis(x, y, spl_data, plot_cfg, fig_name, show_figures)

fig = create_figure(show_figures, ...
    'Name', fig_name, ...
    'Position', plot_cfg.fig_pos_z1, ...
    'Color', 'w');

ax = axes(fig);

pcolor(ax, x, y, spl_data);
shading(ax, 'flat');
colormap(ax, MyColor('vik'));
caxis(ax, plot_cfg.spl_clim);

xlim(ax, [-plot_cfg.xy_max, plot_cfg.xy_max]);
ylim(ax, [-plot_cfg.xy_max, plot_cfg.xy_max]);

daspect(ax, [1 1 1]);
axis(ax, 'off');

set(ax, ...
    'Position', plot_cfg.ax_pos_z1, ...
    'Box', 'off');
end

function [fig, ax] = plot_z1_phase_no_axis(x, y, phase_data, plot_cfg, fig_name, show_figures)

fig = create_figure(show_figures, ...
    'Name', fig_name, ...
    'Position', plot_cfg.fig_pos_z1, ...
    'Color', 'w');

ax = axes(fig);

pcolor(ax, x, y, phase_data);
shading(ax, 'flat');
colormap(ax, 'hsv');
caxis(ax, plot_cfg.phase_clim);

xlim(ax, [-plot_cfg.xy_max, plot_cfg.xy_max]);
ylim(ax, [-plot_cfg.xy_max, plot_cfg.xy_max]);

daspect(ax, [1 1 1]);
axis(ax, 'off');

set(ax, ...
    'Position', plot_cfg.ax_pos_z1, ...
    'Box', 'off');
end

function [fig, ax] = plot_xOz_spl_no_axis(z, x, spl_data, plot_cfg, fig_name, show_figures)

fig = create_figure(show_figures, ...
    'Name', fig_name, ...
    'Position', plot_cfg.fig_pos_xoz, ...
    'Color', 'w');

ax = axes(fig);

pcolor(ax, z, x, spl_data);
shading(ax, 'flat');
colormap(ax, MyColor('vik'));
caxis(ax, plot_cfg.spl_clim);

xlim(ax, [0, plot_cfg.xoz_z_max]);
ylim(ax, [-plot_cfg.xoz_x_max, plot_cfg.xoz_x_max]);

daspect(ax, [1 1 1]);
axis(ax, 'off');

set(ax, ...
    'Position', plot_cfg.ax_pos_xoz, ...
    'Box', 'off');
end

function fig = create_figure(show_figures, varargin)

if show_figures
    fig = figure(varargin{:});
else
    fig = figure(varargin{:}, 'Visible', 'off');
end
end

function save_png_if_needed(ax_handle, save_png, png_file)

if ~save_png
    return;
end

try
    exportgraphics(ax_handle, png_file, ...
        'Resolution', 300, ...
        'BackgroundColor', 'white');
catch
    fig_handle = ancestor(ax_handle, 'figure');
    print(fig_handle, png_file, '-dpng', '-r300');
end
end

function save_if_needed(fig_handle, save_figures, save_fig_file, save_dir, base_name)

if save_figures
    save_fig_all_formats(fig_handle, fullfile(save_dir, base_name), save_fig_file);
end
end

function save_fig_all_formats(fig_handle, file_base, save_fig_file)

set(fig_handle, 'Color', 'w');
set(fig_handle, 'PaperPositionMode', 'auto');

fp_fig = [file_base, '.fig'];
fp_png = [file_base, '.png'];
fp_pdf = [file_base, '.pdf'];

if save_fig_file
    savefig(fig_handle, fp_fig);
end

try
    exportgraphics(fig_handle, fp_png, 'Resolution', 300);
catch
    print(fig_handle, fp_png, '-dpng', '-r300');
end

try
    exportgraphics(fig_handle, fp_pdf, ...
        'ContentType', 'image', ...
        'Resolution', 300);
catch
    print(fig_handle, fp_pdf, '-dpdf', '-opengl', '-r300');
end
end