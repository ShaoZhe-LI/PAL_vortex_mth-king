function result = calc_ultrasound_field(source, medium, calc, grid, method)
% =========================================================================
% CALC_ULTRASOUND_FIELD
% - ULTRASOUND fields at f1/f2 using:
%   (A) King/FHT axisymmetric
%       Green spectrum G(k_rho,z): 'analytic' (default) | 'transform'
%   (B) DIM:
%       'rayleigh' (default, memory-lean block integration)
%       'asm'
%
% UNIFIED GRID POLICY
%   - KING/FHT always uses its own internally generated FHT grid.
%   - DIM uses:
%       * grid.dim.x / y / z, if grid is provided
%       * otherwise, auto-snaps to FHT grid by using:
%           x = fht.xh, y = 0, z = fht.z_ultra
%
% INPUT
%   source, medium, calc : passed to make_source_velocity()
%   grid : optional
%   method : 'king' | 'dim' | 'both' (default 'both')
%
% OUTPUT
%   result.king.* and/or result.dim.*
%
% NOTE ON DIM-RAYLEIGH NEAR-FIELD REGULARIZATION
%   To improve robustness when z is very small, the Rayleigh point-kernel
%   1/R is regularized by an equivalent panel radius:
%       R_reg = sqrt(R^2 + a_eq^2),   a_eq = sqrt(dA/pi)
%   This mainly affects the very near source plane and leaves far-field
%   behavior almost unchanged.
% =========================================================================

if nargin < 5 || isempty(method)
    method = 'both';
end
if nargin < 4
    grid = [];
end
method = lower(strtrim(string(method)));

if nargin < 3 || isempty(calc) || ~isstruct(calc)
    calc = struct();
end

if ~isfield(calc,'dim') || isempty(calc.dim); calc.dim = struct(); end
if ~isfield(calc.dim,'method') || isempty(calc.dim.method)
    calc.dim.method = 'rayleigh';
end
dim_method = lower(strtrim(string(calc.dim.method)));

if ~isfield(calc,'king') || isempty(calc.king); calc.king = struct(); end
if ~isfield(calc.king,'gspec_method') || isempty(calc.king.gspec_method)
    calc.king.gspec_method = 'transform';
end
king_gspec_method = lower(strtrim(string(calc.king.gspec_method)));

do_king = any(method == ["king","both"]);
do_dim  = any(method == ["dim","both"]);

[source, fht, dim_src] = make_source_velocity(source, medium, calc);

result = struct();
result.source  = source;
result.medium  = source.medium;
result.calc    = calc;
result.fht     = fht;
result.dim_src = dim_src;

if ~isfield(calc,'king') || isempty(calc.king); calc.king = struct(); end
if ~isfield(calc.king,'eps_phase') || isempty(calc.king.eps_phase)
    calc.king.eps_phase = 1e-3;
end
if ~isfield(calc.king,'kz_min') || isempty(calc.king.kz_min)
    calc.king.kz_min = 1e-12;
end

if ~isfield(calc.king,'band_refine') || isempty(calc.king.band_refine)
    calc.king.band_refine = struct();
end
if ~isfield(calc.king.band_refine,'enable') || isempty(calc.king.band_refine.enable)
    calc.king.band_refine.enable = false;
end
if ~isfield(calc.king.band_refine,'N_mult') || isempty(calc.king.band_refine.N_mult)
    calc.king.band_refine.N_mult = 2;
end
if ~isfield(calc.king.band_refine,'delta_k_factor') || isempty(calc.king.band_refine.delta_k_factor)
    calc.king.band_refine.delta_k_factor = 8;
end
if ~isfield(calc.king.band_refine,'taper_k_factor') || isempty(calc.king.band_refine.taper_k_factor)
    calc.king.band_refine.taper_k_factor = 2;
end

if ~isfield(calc,'dim') || isempty(calc.dim); calc.dim = struct(); end
if ~isfield(calc.dim,'use_parallel') || isempty(calc.dim.use_parallel)
    calc.dim.use_parallel = true;
end
if ~isfield(calc.dim,'num_workers') || isempty(calc.dim.num_workers)
    calc.dim.num_workers = 4;
end

% ===== 新增：DIM-Rayleigh 近源正则化默认参数 =====
if ~isfield(calc.dim,'nearfield_regularization') || isempty(calc.dim.nearfield_regularization)
    calc.dim.nearfield_regularization = struct();
end
if ~isfield(calc.dim.nearfield_regularization,'enable') || isempty(calc.dim.nearfield_regularization.enable)
    calc.dim.nearfield_regularization.enable = true;
end
if ~isfield(calc.dim.nearfield_regularization,'mode') || isempty(calc.dim.nearfield_regularization.mode)
    calc.dim.nearfield_regularization.mode = 'panel_radius'; % only supported mode now
end
if ~isfield(calc.dim.nearfield_regularization,'min_radius_factor') || isempty(calc.dim.nearfield_regularization.min_radius_factor)
    calc.dim.nearfield_regularization.min_radius_factor = 0.5;
end

if ~isfield(calc,'asm') || isempty(calc.asm); calc.asm = struct(); end
if ~isfield(calc.asm,'pad_factor') || isempty(calc.asm.pad_factor)
    calc.asm.pad_factor = 8;
end
if ~isfield(calc.asm,'kzz_eps') || isempty(calc.asm.kzz_eps)
    calc.asm.kzz_eps = 1e-12;
end

% =========================================================================
% (A) KING / FHT
% =========================================================================
if do_king
    rho_k = fht.xh;
    z_k   = fht.z_ultra;
    Nz    = numel(z_k);

    Vs1 = m_FHT(source.fht.vs1_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    Vs2 = m_FHT(source.fht.vs2_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    kr   = fht.xH;
    zrow = reshape(z_k, 1, []);

    switch king_gspec_method
        case "analytic"
            G1 = local_green_spec_analytic(source.k1, kr, zrow, calc.king.eps_phase, calc.king.kz_min);
            G2 = local_green_spec_analytic(source.k2, kr, zrow, calc.king.eps_phase, calc.king.kz_min);

        case "transform"
            Nr = numel(rho_k);

            G1 = complex(zeros(Nr, Nz));
            G2 = complex(zeros(Nr, Nz));

            for iz = 1:Nz
                z0 = z_k(iz);
                r = hypot(rho_k, z0);
                r(r < 1e-9) = 1e-9;

                g1 = exp(1j * source.k1 .* r) ./ (4*pi*r);
                g2 = exp(1j * source.k2 .* r) ./ (4*pi*r);

                G1(:,iz) = m_FHT(g1, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);

                G2(:,iz) = m_FHT(g2, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);
            end

        otherwise
            error('calc_ultrasound_field:BadKingGspecMethod', ...
                'calc.king.gspec_method must be ''analytic'' or ''transform''.');
    end

    F1_c = G1 .* Vs1;
    F2_c = G2 .* Vs2;

    phi1_c = -4*pi * m_FHT(F1_c, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    phi2_c = -4*pi * m_FHT(F2_c, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

    phi1 = phi1_c;
    phi2 = phi2_c;

    if isfield(calc,'king') && isfield(calc.king,'band_refine') ...
            && calc.king.band_refine.enable

        kr_c = kr(:);
        zrow = reshape(z_k, 1, []);

        k1r = real(source.k1); a1 = abs(imag(source.k1));
        k2r = real(source.k2); a2 = abs(imag(source.k2));

        dk1 = max(calc.king.band_refine.delta_k_factor * a1, 1e-12);
        dk2 = max(calc.king.band_refine.delta_k_factor * a2, 1e-12);
        tw1 = calc.king.band_refine.taper_k_factor * dk1;
        tw2 = calc.king.band_refine.taper_k_factor * dk2;

        W1_c = local_band_window(kr_c, k1r, dk1, tw1);
        W2_c = local_band_window(kr_c, k2r, dk2, tw2);

        F1_cb = F1_c .* W1_c;
        F2_cb = F2_c .* W2_c;

        phi1_cb = -4*pi * m_FHT(F1_cb, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
            fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

        phi2_cb = -4*pi * m_FHT(F2_cb, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
            fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);

        calc_f = calc;
        calc_f.fht.N_FHT = calc.fht.N_FHT * calc.king.band_refine.N_mult;

        [source_f, fht_f] = make_source_velocity(source, medium, calc_f);

        kr_f = fht_f.xH(:);

        Nh_rho = fht.Nh;
        rho_c  = fht.xh(:);
        rho_f  = fht_f.x1(:) * Nh_rho;

        Vs1_f = m_FHT(source_f.fht.vs1_on_xh_v, ...
            fht_f.N_FHT, 1, fht_f.Nh_v, fht_f.NH_v, ...
            fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, source_f.m_used);

        Vs2_f = m_FHT(source_f.fht.vs2_on_xh_v, ...
            fht_f.N_FHT, 1, fht_f.Nh_v, fht_f.NH_v, ...
            fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, source_f.m_used);

        switch king_gspec_method
            case "analytic"
                G1_f = local_green_spec_analytic(source.k1, kr_f, zrow, ...
                    calc.king.eps_phase, calc.king.kz_min);
                G2_f = local_green_spec_analytic(source.k2, kr_f, zrow, ...
                    calc.king.eps_phase, calc.king.kz_min);

            case "transform"
                rho_k_f = fht_f.xh(:);
                Nr_f = numel(rho_k_f);

                G1_f = complex(zeros(Nr_f, Nz));
                G2_f = complex(zeros(Nr_f, Nz));

                for iz = 1:Nz
                    z0 = z_k(iz);
                    r = hypot(rho_k_f, z0);
                    r(r < 1e-9) = 1e-9;

                    g1 = exp(1j * source.k1 .* r) ./ (4*pi*r);
                    g2 = exp(1j * source.k2 .* r) ./ (4*pi*r);

                    G1_f(:,iz) = m_FHT(g1, fht_f.N_FHT, 1, fht_f.Nh, fht_f.NH, ...
                        fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, 0);

                    G2_f(:,iz) = m_FHT(g2, fht_f.N_FHT, 1, fht_f.Nh, fht_f.NH, ...
                        fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, 0);
                end

            otherwise
                error('calc_ultrasound_field:BadKingGspecMethod', ...
                    'calc.king.gspec_method must be ''analytic'' or ''transform''.');
        end

        F1_f = G1_f .* Vs1_f;
        F2_f = G2_f .* Vs2_f;

        W1_f = local_band_window(kr_f, k1r, dk1, tw1);
        W2_f = local_band_window(kr_f, k2r, dk2, tw2);

        F1_fb = F1_f .* W1_f;
        F2_fb = F2_f .* W2_f;

        phi1_fb = -4*pi * m_FHT(F1_fb, fht_f.N_FHT, Nz, fht_f.NH, Nh_rho, ...
            fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, source_f.m_used);

        phi2_fb = -4*pi * m_FHT(F2_fb, fht_f.N_FHT, Nz, fht_f.NH, Nh_rho, ...
            fht_f.a_solve, fht_f.x0, fht_f.x1, fht_f.k0, source_f.m_used);

        phi1_fb_on_c = interp1(log(rho_f), phi1_fb, log(rho_c), 'linear', 0);
        phi2_fb_on_c = interp1(log(rho_f), phi2_fb, log(rho_c), 'linear', 0);

        phi1 = phi1_c - phi1_cb + phi1_fb_on_c;
        phi2 = phi2_c - phi2_cb + phi2_fb_on_c;
    end

    p1 = 1j * source.medium.rho0 * source.medium.c0 * real(source.k1) .* phi1;
    p2 = 1j * source.medium.rho0 * source.medium.c0 * real(source.k2) .* phi2;

    result.king = struct();
    result.king.method = "King-FHT";
    result.king.green_spec_method = king_gspec_method;
    result.king.eps_phase = calc.king.eps_phase;
    result.king.kz_min = calc.king.kz_min;
    result.king.rho = rho_k;
    result.king.z = z_k;
    result.king.p_f1 = p1;
    result.king.p_f2 = p2;
    result.king.f1 = source.f1;
    result.king.f2 = source.f2;
    result.king.G1 = G1;
    result.king.Vs1 = Vs1;
end

% =========================================================================
% (B) DIM
% =========================================================================
if do_dim
    gd = local_get_dim_grid_or_snap(grid, fht);

    [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd);
    Nx = numel(x_obs); Ny = numel(y_obs); Nz = numel(z_obs); %#ok<NASGU>

    if strcmpi(dim_method, "asm")
        if isfield(source,'dim') && isfield(source.dim,'src_discretization')
            if strcmpi(string(source.dim.src_discretization), "polar")
                error('calc_ultrasound_field:ASMNeedsCartesian', ...
                    'DIM-ASM requires Cartesian source sampling. Set calc.dim.src_discretization = ''cart''.');
            end
        end
    end

    switch dim_method
        case "rayleigh"
            blk = 200000;
            if isfield(gd,'block_size') && ~isempty(gd.block_size)
                blk = gd.block_size;
            end

            src_blk = 5000;
            if isfield(gd,'src_block_size') && ~isempty(gd.src_block_size)
                src_blk = gd.src_block_size;
            end

            Xs = source.dim.Xs(:);
            Ys = source.dim.Ys(:);
            Zs = source.dim.Zs(:);

            vn1 = source.dim.Vn_pts_f1(:);
            vn2 = source.dim.Vn_pts_f2(:);

            if isfield(source.dim,'dA_pts') && ~isempty(source.dim.dA_pts)
                dA_w = source.dim.dA_pts;
            else
                dA_w = source.dim.dA;
            end
            dA_w = dA_w(:);

            if isscalar(dA_w)
                q1 = vn1 * dA_w;
                q2 = vn2 * dA_w;
                dA_vec = repmat(dA_w, numel(vn1), 1);
            else
                if numel(dA_w) ~= numel(vn1)
                    error('calc_ultrasound_field:BadDA', 'dA_pts size mismatch with Vn_pts.');
                end
                q1 = vn1 .* dA_w;
                q2 = vn2 .* dA_w;
                dA_vec = dA_w;
            end

            Nsrc = numel(Xs);

            % ===== 新增：每个源点对应的等效面元半径 =====
            reg = calc.dim.nearfield_regularization;
            aeq = local_get_equiv_panel_radius(dA_vec, reg.min_radius_factor);

            P1 = complex(zeros(Ny, Nx, Nz));
            P2 = complex(zeros(Ny, Nx, Nz));

            [X2, Y2] = meshgrid(x_obs, y_obs);
            xo = X2(:); yo = Y2(:);
            nP = numel(xo);

            use_par = logical(calc.dim.use_parallel);
            nW_req  = calc.dim.num_workers;

            if Nz < nW_req
                use_par = false;
            end

            p = gcp('nocreate');
            if use_par
                if ~isempty(p)
                    nW_eff = p.NumWorkers;
                    if Nz < nW_eff
                        use_par = false;
                    end
                else
                    try
                        parpool('local', nW_req);
                    catch ME
                        warning('calc_ultrasound_field:ParpoolFailed', ...
                            'Failed to start parpool (%s). Falling back to serial.', ME.message);
                        use_par = false;
                    end
                end
            end

            if use_par
                parfor iz = 1:Nz
                    z0 = z_obs(iz);
                    p1_vec = complex(zeros(nP,1));
                    p2_vec = complex(zeros(nP,1));

                    for s = 1:blk:nP
                        id = s:min(s+blk-1, nP);

                        phi1 = complex(zeros(numel(id),1));
                        phi2 = complex(zeros(numel(id),1));

                        for js = 1:src_blk:Nsrc
                            jd = js:min(js+src_blk-1, Nsrc);

                            dXj = xo(id) - Xs(jd).';
                            dYj = yo(id) - Ys(jd).';
                            dZj = z0     - Zs(jd).';
                            Rj  = sqrt(dXj.^2 + dYj.^2 + dZj.^2);

                            if reg.enable
                                Rj_use = local_regularize_distance(Rj, aeq(jd).', reg.mode);
                            else
                                Rj_use = Rj;
                                Rj_use(Rj_use < 1e-9) = 1e-9;
                            end

                            Gj = exp(1j*source.k1.*Rj_use) ./ (4*pi*Rj_use);
                            phi1 = phi1 + Gj * q1(jd);

                            Gj = exp(1j*source.k2.*Rj_use) ./ (4*pi*Rj_use);
                            phi2 = phi2 + Gj * q2(jd);
                        end

                        p1_vec(id) = -2j * source.medium.rho0 * source.w1 .* phi1;
                        p2_vec(id) = -2j * source.medium.rho0 * source.w2 .* phi2;
                    end

                    P1(:,:,iz) = reshape(p1_vec, Ny, Nx);
                    P2(:,:,iz) = reshape(p2_vec, Ny, Nx);
                end
            else
                for iz = 1:Nz
                    z0 = z_obs(iz);
                    p1_vec = complex(zeros(nP,1));
                    p2_vec = complex(zeros(nP,1));

                    for s = 1:blk:nP
                        id = s:min(s+blk-1, nP);

                        phi1 = complex(zeros(numel(id),1));
                        phi2 = complex(zeros(numel(id),1));

                        for js = 1:src_blk:Nsrc
                            jd = js:min(js+src_blk-1, Nsrc);

                            dXj = xo(id) - Xs(jd).';
                            dYj = yo(id) - Ys(jd).';
                            dZj = z0     - Zs(jd).';
                            Rj  = sqrt(dXj.^2 + dYj.^2 + dZj.^2);

                            if reg.enable
                                Rj_use = local_regularize_distance(Rj, aeq(jd).', reg.mode);
                            else
                                Rj_use = Rj;
                                Rj_use(Rj_use < 1e-9) = 1e-9;
                            end

                            Gj = exp(1j*source.k1.*Rj_use) ./ (4*pi*Rj_use);
                            phi1 = phi1 + Gj * q1(jd);

                            Gj = exp(1j*source.k2.*Rj_use) ./ (4*pi*Rj_use);
                            phi2 = phi2 + Gj * q2(jd);
                        end

                        p1_vec(id) = -2j * source.medium.rho0 * source.w1 .* phi1;
                        p2_vec(id) = -2j * source.medium.rho0 * source.w2 .* phi2;
                    end

                    P1(:,:,iz) = reshape(p1_vec, Ny, Nx);
                    P2(:,:,iz) = reshape(p2_vec, Ny, Nx);
                end
            end

            result.dim = struct();
            result.dim.method = "DIM-Rayleigh";
            result.dim.dim_method = "rayleigh";
            result.dim.x = x_obs;
            result.dim.y = y_obs;
            result.dim.z = z_obs;
            result.dim.p_f1 = P1;
            result.dim.p_f2 = P2;
            result.dim.f1 = source.f1;
            result.dim.f2 = source.f2;
            result.dim.block_size = blk;
            result.dim.src_block_size = src_blk;
            result.dim.use_parallel = use_par;
            result.dim.nearfield_regularization = reg;
            if use_par
                result.dim.num_workers = calc.dim.num_workers;
            end

        case "asm"
            dx = source.dim.dx; dy = source.dim.dy;

            xg = source.dim.x_grid(:).';
            yg = source.dim.y_grid(:).';
            if numel(xg) < 2 || numel(yg) < 2
                error('ASM requires x_grid/y_grid length >= 2.');
            end
            if max(abs(diff(xg) - dx)) > 1e-12 || max(abs(diff(yg) - dy)) > 1e-12
                error('ASM requires uniform x_grid/y_grid consistent with dx/dy.');
            end

            if ~isfield(source.dim,'dx') || ~isfield(source.dim,'dy')
                error('calc_ultrasound_field:ASMNeedDxDy', 'ASM needs source.dim.dx and source.dim.dy.');
            end
            if ~isfield(source.dim,'x_grid') || ~isfield(source.dim,'y_grid')
                error('calc_ultrasound_field:ASMNeedXYGrid', 'ASM needs source.dim.x_grid and source.dim.y_grid.');
            end

            Vn0_f1 = source.dim.Vn_grid_f1;
            Vn0_f2 = source.dim.Vn_grid_f2;
            x0 = source.dim.x_grid(:).';
            y0 = source.dim.y_grid(:).';
            [Ny0, Nx0] = size(Vn0_f1);

            pf = calc.asm.pad_factor;
            NyP = pf * Ny0;
            NxP = pf * Nx0;

            VnP_f1 = zeros(NyP, NxP);
            VnP_f2 = zeros(NyP, NxP);

            iy0 = floor((NyP - Ny0)/2) + (1:Ny0);
            ix0 = floor((NxP - Nx0)/2) + (1:Nx0);

            VnP_f1(iy0, ix0) = Vn0_f1;
            VnP_f2(iy0, ix0) = Vn0_f2;

            x0c = x0( floor((Nx0+1)/2) );
            y0c = y0( floor((Ny0+1)/2) );

            x_pad = ((0:NxP-1) - floor((NxP+1)/2)) * dx + x0c;
            y_pad = ((0:NyP-1) - floor((NyP+1)/2)) * dy + y0c;

            P1_full = local_asm_pressure_from_vn(VnP_f1, dx, dy, source.k1, source.medium.rho0, source.w1, z_obs, ...
                calc.asm.kzz_eps);
            P2_full = local_asm_pressure_from_vn(VnP_f2, dx, dy, source.k2, source.medium.rho0, source.w2, z_obs, ...
                calc.asm.kzz_eps);

            [x_out, y_out] = deal(x_pad, y_pad);
            [Xo, Yo] = meshgrid(x_out, y_out);

            result.dim = struct();
            result.dim.method = "DIM-ASM";
            result.dim.dim_method = "asm";
            result.dim.x = x_out;
            result.dim.y = y_out;
            result.dim.z = z_obs;
            result.dim.X = Xo;
            result.dim.Y = Yo;
            result.dim.p_f1 = P1_full;
            result.dim.p_f2 = P2_full;
            result.dim.f1 = source.f1;
            result.dim.f2 = source.f2;
            result.dim.dx = dx;
            result.dim.dy = dy;
            result.dim.pad_factor = pf;

        otherwise
            error('calc_ultrasound_field:BadDimMethod', ...
                'calc.dim.method must be ''rayleigh'' or ''asm''.');
    end
end

end

% ======================================================================
% helpers
% ======================================================================
function G = local_green_spec_analytic(k, kr_col, z_row, eps_phase, kz_min)
kz = sqrt(k.^2 - kr_col.^2);
idx = imag(kz) < 0;
kz(idx) = -kz(idx);

Nr = numel(kz);
Nz = numel(z_row);

KZ = kz * ones(1, Nz);
Z  = ones(Nr, 1) * reshape(z_row, 1, []);

KZs = KZ;
mask_kz0 = abs(KZs) < kz_min;
if any(mask_kz0(:))
    KZs(mask_kz0) = kz_min .* exp(1j * angle(KZs(mask_kz0)));
end

G = (1j/(4*pi)) * exp(1j * KZ .* Z) ./ KZs;

phase  = abs(KZ .* Z);
rel_im = abs(imag(KZ)) ./ max(abs(KZ), kz_min);

mask = (phase < eps_phase) & (rel_im < 1e-6);
if any(mask(:))
    Gt = (1./KZs) + 1j*Z - (KZs .* (Z.^2))/2;
    G(mask) = (1j/(4*pi)) * Gt(mask);
end
end

function [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd)
if isfield(gd,'x') && ~isempty(gd.x)
    x_obs = gd.x(:).';
else
    error('calc_ultrasound_field:BadGridX', 'grid.dim.x is required (vector).');
end

if isfield(gd,'y') && ~isempty(gd.y)
    y_obs = gd.y(:).';
else
    error('calc_ultrasound_field:BadGridY', 'grid.dim.y is required (vector).');
end

if isfield(gd,'z') && ~isempty(gd.z)
    z_obs = gd.z(:).';
else
    error('calc_ultrasound_field:BadGridZ', 'grid.dim.z is required (vector).');
end
end

function gd = local_get_dim_grid_or_snap(grid, fht)
if nargin >= 1 && ~isempty(grid) && isstruct(grid) && isfield(grid,'dim') && ~isempty(grid.dim)
    gd = grid.dim;
else
    gd = struct();
    gd.x = fht.xh(:).';
    gd.y = 0;
    gd.z = fht.z_ultra(:).';
end
end

function P = local_asm_pressure_from_vn(Vn0, dx, dy, k, rho0, w, z_list, kzz_eps)
[Ny, Nx] = size(Vn0);

kx = (2*pi) * [0:floor(Nx/2), -ceil(Nx/2)+1:-1] / (Nx*dx);
ky = (2*pi) * [0:floor(Ny/2), -ceil(Ny/2)+1:-1] / (Ny*dy);
[KX, KY] = meshgrid(kx, ky);

kz = sqrt(k.^2 - KX.^2 - KY.^2);
kz(imag(kz)<0) = -kz(imag(kz)<0);

Vhat = fft2(Vn0);

coef = rho0 * w;

Nz = numel(z_list);
P = complex(zeros(Ny, Nx, Nz));

for iz = 1:Nz
    z = z_list(iz);
    H = exp(1j*kz*z) ./ (kz + kzz_eps);
    P(:,:,iz) = ifft2( coef * (Vhat .* H) );
end
end

function W = local_band_window(kr, k0, dk, tw)
d = abs(kr - k0);
W = zeros(size(kr));

id1 = d <= dk;
W(id1) = 1;

id2 = (d > dk) & (d < dk + tw);
xi = (d(id2) - dk) / tw;
W(id2) = 0.5 * (1 + cos(pi * xi));
end

% ===== 新增 helper 1：等效面元半径 =====
function aeq = local_get_equiv_panel_radius(dA_vec, min_radius_factor)
dA_vec = abs(dA_vec(:));
aeq = sqrt(max(dA_vec, 0) / pi);

% 避免出现极小/零面积导致正则化半径失效
a_floor = min_radius_factor * sqrt(max(median(dA_vec(dA_vec > 0)), eps) / pi);
if isempty(a_floor) || ~isfinite(a_floor) || a_floor <= 0
    a_floor = 1e-12;
end
aeq = max(aeq, a_floor);
end

% ===== 新增 helper 2：近源正则化距离 =====
function R_use = local_regularize_distance(R, aeq_row, mode)
switch lower(string(mode))
    case "panel_radius"
        R_use = sqrt(R.^2 + aeq_row.^2);
    otherwise
        error('calc_ultrasound_field:BadNearfieldRegularizationMode', ...
            'Unsupported nearfield regularization mode.');
end

R_use(R_use < 1e-12) = 1e-12;
end