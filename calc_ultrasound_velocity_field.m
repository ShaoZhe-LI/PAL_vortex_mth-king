function result = calc_ultrasound_velocity_field(source, medium, calc, grid, method)
% =========================================================================
% CALC_ULTRASOUND_VELOCITY_FIELD
% - Compute ultrasonic particle velocity field.
%
% UNIFIED GRID POLICY
%   - KING/FHT always uses its own internally generated FHT grid.
%   - DIM uses:
%       * grid.dim.x / y / z, if grid is provided
%       * otherwise, auto-snaps to FHT grid by using:
%           x = fht.xh, y = 0, z = fht.z_ultra
%
% METHODS
%   method = 'king'  : King/FHT velocity on (rho,z)
%   method = 'dim'   : DIM velocity on Cartesian grid.dim
%   method = 'both'  : both
%
% OUTPUT
%   result.king : axisymmetric velocity on (rho,z)
%   result.dim  : Cartesian + cylindrical components on (x,y,z)
%
% KING/FHT GREEN-SPECTRUM OPTIONS
%   calc.king.gspec_method = 'analytic'  : use analytic spectrum
%   calc.king.gspec_method = 'transform' : use 0th-order FHT of Green func
% =========================================================================

if nargin < 5 || isempty(method), method = 'king'; end
if nargin < 4, grid = []; end
method = lower(strtrim(string(method)));

do_king = any(method == ["king","both"]);
do_dim  = any(method == ["dim","both"]);

if ~do_king && ~do_dim
    error('calc_ultrasound_velocity_field:BadMethod', ...
        'method must be ''king'', ''dim'', or ''both''.');
end

if nargin < 3 || isempty(calc) || ~isstruct(calc), calc = struct(); end
if ~isfield(calc,'king') || isempty(calc.king), calc.king = struct(); end

if ~isfield(calc.king,'eps_phase') || isempty(calc.king.eps_phase)
    calc.king.eps_phase = 1e-3;
end
if ~isfield(calc.king,'kz_min') || isempty(calc.king.kz_min)
    calc.king.kz_min = 1e-12;
end
if ~isfield(calc.king,'gspec_method') || isempty(calc.king.gspec_method)
    calc.king.gspec_method = 'transform';
end
king_gspec_method = lower(strtrim(string(calc.king.gspec_method)));

if ~isfield(calc,'dim') || isempty(calc.dim), calc.dim = struct(); end
if ~isfield(calc.dim,'src_discretization') || isempty(calc.dim.src_discretization)
    calc.dim.src_discretization = 'cart';
end
calc.dim.src_discretization = lower(strtrim(string(calc.dim.src_discretization)));
if ~any(calc.dim.src_discretization == ["cart","polar"])
    error('calc_ultrasound_velocity_field:BadSrcDiscretization', ...
        'calc.dim.src_discretization must be ''cart'' or ''polar''.');
end
if ~isfield(calc.dim,'Nphi') || isempty(calc.dim.Nphi)
    calc.dim.Nphi = 256;
end
if ~isfield(calc.dim,'block_size') || isempty(calc.dim.block_size)
    calc.dim.block_size = 20000;
end
if ~isfield(calc.dim,'src_block_size') || isempty(calc.dim.src_block_size)
    calc.dim.src_block_size = 5000;
end

[source, fht, ~] = make_source_velocity(source, medium, calc);

result = struct();
result.source = source;
result.medium = source.medium;
result.calc   = calc;
result.fht    = fht;

rho = fht.xh(:);
z   = fht.z_ultra(:).';
Nrho = numel(rho);
Nz   = numel(z);
m = source.m_used;

% =========================================================================
% KING / FHT
% =========================================================================
if do_king
    kr = fht.xH(:);
    Nkr = numel(kr);

    Vs1 = m_FHT(source.fht.vs1_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vs2 = m_FHT(source.fht.vs2_on_xh_v, ...
        fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vs1 = Vs1(:);
    Vs2 = Vs2(:);

    Zabs = abs(z);
    tol_z0 = 1e-15;
    sgnZ = sign(z);
    sgnZ(abs(z) < tol_z0) = 0;

    KZ1 = local_kz_branch(source.k1, kr, Nz);
    KZ2 = local_kz_branch(source.k2, kr, Nz);

    zrow = reshape(Zabs, 1, []);

    switch king_gspec_method
        case "analytic"
            G1 = local_green_spec_analytic(source.k1, kr, zrow, ...
                calc.king.eps_phase, calc.king.kz_min);
            G2 = local_green_spec_analytic(source.k2, kr, zrow, ...
                calc.king.eps_phase, calc.king.kz_min);

        case "transform"
            G1 = complex(zeros(Nkr, Nz));
            G2 = complex(zeros(Nkr, Nz));

            for iz = 1:Nz
                z0 = Zabs(iz);

                r = hypot(rho, z0);
                r(r < 1e-9) = 1e-9;

                g1 = exp(1j * source.k1 .* r) ./ (4*pi*r);
                g2 = exp(1j * source.k2 .* r) ./ (4*pi*r);

                G1(:,iz) = m_FHT(g1, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);

                G2(:,iz) = m_FHT(g2, fht.N_FHT, 1, fht.Nh, fht.NH, ...
                    fht.a_solve, fht.x0, fht.x1, fht.k0, 0);
            end

        otherwise
            error('calc_ultrasound_velocity_field:BadKingGspecMethod', ...
                'calc.king.gspec_method must be ''analytic'' or ''transform''.');
    end

    % exp(i*kz*z)/kz = -4*pi*i * G
    Hcom1 = -4*pi*1i * G1;
    Hcom2 = -4*pi*1i * G2;

    % exp(i*kz*z)
    Hz1 = KZ1 .* Hcom1;
    Hz2 = KZ2 .* Hcom2;

    % kr * exp(i*kz*z)/kz
    Hr1 = (kr * ones(1,Nz)) .* Hcom1;
    Hr2 = (kr * ones(1,Nz)) .* Hcom2;

    Scom1 = (Vs1 * ones(1,Nz)) .* Hcom1;
    Scom2 = (Vs2 * ones(1,Nz)) .* Hcom2;

    Sz1   = (Vs1 * ones(1,Nz)) .* Hz1;
    Sz2   = (Vs2 * ones(1,Nz)) .* Hz2;

    Sr1   = (Vs1 * ones(1,Nz)) .* Hr1;
    Sr2   = (Vs2 * ones(1,Nz)) .* Hr2;

    Vz1 = m_FHT(Sz1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);
    Vz2 = m_FHT(Sz2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    Vz1 = Vz1 .* (ones(Nrho,1) * sgnZ);
    Vz2 = Vz2 .* (ones(Nrho,1) * sgnZ);

    Vphi_core1 = m_FHT(Scom1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);
    Vphi_core2 = m_FHT(Scom2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m);

    if m == 0
        Vphi1 = complex(zeros(Nrho, Nz));
        Vphi2 = complex(zeros(Nrho, Nz));
    else
        rho_safe = rho;
        pos = rho_safe(rho_safe>0);
        if isempty(pos)
            rho_safe(:) = 1e-12;
        else
            rho_safe(rho_safe==0) = min(pos) * 1e-6 + 1e-12;
        end
        Vphi1 = (m ./ rho_safe) * ones(1,Nz) .* Vphi_core1;
        Vphi2 = (m ./ rho_safe) * ones(1,Nz) .* Vphi_core2;
    end

    Inv_m1_f1 = m_FHT(Sr1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m-1);
    Inv_p1_f1 = m_FHT(Sr1, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m+1);

    Inv_m1_f2 = m_FHT(Sr2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m-1);
    Inv_p1_f2 = m_FHT(Sr2, fht.N_FHT, Nz, fht.NH, fht.Nh, ...
        fht.a_solve, fht.x0, fht.x1, fht.k0, m+1);

    Vr1 = -(1i/2) * (Inv_m1_f1 - Inv_p1_f1);
    Vr2 = -(1i/2) * (Inv_m1_f2 - Inv_p1_f2);

    result.king = struct();
    result.king.method            = "King-FHT velocity (Eq.13)";
    result.king.green_spec_method = king_gspec_method;
    result.king.m                 = m;
    result.king.rho               = rho;
    result.king.z                 = z;

    result.king.v_rho_f1 = Vr1;
    result.king.v_phi_f1 = Vphi1;
    result.king.v_z_f1   = Vz1;

    result.king.v_rho_f2 = Vr2;
    result.king.v_phi_f2 = Vphi2;
    result.king.v_z_f2   = Vz2;

    result.king.f1 = source.f1;
    result.king.f2 = source.f2;

    result.king.G1 = G1;
    result.king.G2 = G2;
    result.king.Vs1 = Vs1;
    result.king.Vs2 = Vs2;
end

% =========================================================================
% DIM
% =========================================================================
if do_dim
    gd = local_get_dim_grid_or_snap(grid, fht);
    [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd);

    blk     = max(1000, round(calc.dim.block_size));
    src_blk = max(500,  round(calc.dim.src_block_size));

    if ~isfield(source,'dim') || ~isfield(source.dim,'Xs') || isempty(source.dim.Xs)
        error('calc_ultrasound_velocity_field:DIMNeedSourceDim', ...
            'source.dim.* not found. Ensure make_source_velocity builds DIM source points.');
    end

    Xs = source.dim.Xs(:);
    Ys = source.dim.Ys(:);
    Zs = source.dim.Zs(:);

    vn1 = source.dim.Vn_pts_f1(:);
    vn2 = source.dim.Vn_pts_f2(:);

    src_mode = string(calc.dim.src_discretization);

    if src_mode == "polar"
        if ~isfield(source.dim,'dA_pts') || isempty(source.dim.dA_pts)
            error('calc_ultrasound_velocity_field:PolarNeedsDApts', ...
                'Polar source discretization requires source.dim.dA_pts.');
        end
        dA_w = source.dim.dA_pts(:);
        if numel(dA_w) ~= numel(vn1)
            error('calc_ultrasound_velocity_field:BadDA', ...
                'Polar: numel(dA_pts) must equal numel(Vn_pts).');
        end
        q1 = vn1 .* dA_w;
        q2 = vn2 .* dA_w;
    else
        if isfield(source.dim,'dA_pts') && ~isempty(source.dim.dA_pts)
            dA_w = source.dim.dA_pts(:);
            if numel(dA_w) ~= numel(vn1)
                error('calc_ultrasound_velocity_field:BadDA', ...
                    'Cart: numel(dA_pts) must equal numel(Vn_pts).');
            end
            q1 = vn1 .* dA_w;
            q2 = vn2 .* dA_w;
        else
            if ~isfield(source.dim,'dA') || isempty(source.dim.dA)
                error('calc_ultrasound_velocity_field:BadDA', ...
                    'Cart: need source.dim.dA or dA_pts.');
            end
            q1 = vn1 * source.dim.dA;
            q2 = vn2 * source.dim.dA;
        end
    end

    Nsrc = numel(Xs);
    [X2, Y2] = meshgrid(x_obs, y_obs);
    Ny = size(X2,1);
    Nx = size(X2,2);
    xo = X2(:);
    yo = Y2(:);
    Nobs = numel(xo);
    Nz_dim = numel(z_obs);

    VX1 = complex(zeros(Ny, Nx, Nz_dim));
    VY1 = complex(zeros(Ny, Nx, Nz_dim));
    VZ1 = complex(zeros(Ny, Nx, Nz_dim));
    VR1 = complex(zeros(Ny, Nx, Nz_dim));
    VPHI1 = complex(zeros(Ny, Nx, Nz_dim));

    VX2 = complex(zeros(Ny, Nx, Nz_dim));
    VY2 = complex(zeros(Ny, Nx, Nz_dim));
    VZ2 = complex(zeros(Ny, Nx, Nz_dim));
    VR2 = complex(zeros(Ny, Nx, Nz_dim));
    VPHI2 = complex(zeros(Ny, Nx, Nz_dim));

    k1 = source.k1;
    k2 = source.k2;
    C_DIM = -2;

    PHI = atan2(Y2, X2);
    cP = cos(PHI);
    sP = sin(PHI);

    for iz = 1:Nz_dim
        zo = z_obs(iz);

        vx1_vec = complex(zeros(Nobs,1));
        vy1_vec = complex(zeros(Nobs,1));
        vz1_vec = complex(zeros(Nobs,1));

        vx2_vec = complex(zeros(Nobs,1));
        vy2_vec = complex(zeros(Nobs,1));
        vz2_vec = complex(zeros(Nobs,1));

        for s = 1:blk:Nobs
            id = s:min(s+blk-1, Nobs);
            nb = numel(id);

            vx1_blk = complex(zeros(nb,1));
            vy1_blk = complex(zeros(nb,1));
            vz1_blk = complex(zeros(nb,1));

            vx2_blk = complex(zeros(nb,1));
            vy2_blk = complex(zeros(nb,1));
            vz2_blk = complex(zeros(nb,1));

            for js = 1:src_blk:Nsrc
                jd = js:min(js+src_blk-1, Nsrc);

                dX = xo(id) - Xs(jd).';
                dY = yo(id) - Ys(jd).';
                dZ = zo     - Zs(jd).';

                R  = sqrt(dX.^2 + dY.^2 + dZ.^2);
                R(R < 1e-9) = 1e-9;

                common1 = exp(1j*k1.*R) ./ (4*pi) .* (1j*k1./(R.^2) - 1./(R.^3));
                common2 = exp(1j*k2.*R) ./ (4*pi) .* (1j*k2./(R.^2) - 1./(R.^3));

                Gx1 = dX .* common1;  Gy1 = dY .* common1;  Gz1 = dZ .* common1;
                Gx2 = dX .* common2;  Gy2 = dY .* common2;  Gz2 = dZ .* common2;

                vx1_blk = vx1_blk + Gx1 * q1(jd);
                vy1_blk = vy1_blk + Gy1 * q1(jd);
                vz1_blk = vz1_blk + Gz1 * q1(jd);

                vx2_blk = vx2_blk + Gx2 * q2(jd);
                vy2_blk = vy2_blk + Gy2 * q2(jd);
                vz2_blk = vz2_blk + Gz2 * q2(jd);
            end

            vx1_vec(id) = C_DIM * vx1_blk;
            vy1_vec(id) = C_DIM * vy1_blk;
            vz1_vec(id) = C_DIM * vz1_blk;

            vx2_vec(id) = C_DIM * vx2_blk;
            vy2_vec(id) = C_DIM * vy2_blk;
            vz2_vec(id) = C_DIM * vz2_blk;
        end

        vx1p = reshape(vx1_vec, Ny, Nx);
        vy1p = reshape(vy1_vec, Ny, Nx);
        vz1p = reshape(vz1_vec, Ny, Nx);

        vx2p = reshape(vx2_vec, Ny, Nx);
        vy2p = reshape(vy2_vec, Ny, Nx);
        vz2p = reshape(vz2_vec, Ny, Nx);

        VX1(:,:,iz) = vx1p;
        VY1(:,:,iz) = vy1p;
        VZ1(:,:,iz) = vz1p;
        VR1(:,:,iz) = vx1p .* cP + vy1p .* sP;
        VPHI1(:,:,iz) = -vx1p .* sP + vy1p .* cP;

        VX2(:,:,iz) = vx2p;
        VY2(:,:,iz) = vy2p;
        VZ2(:,:,iz) = vz2p;
        VR2(:,:,iz) = vx2p .* cP + vy2p .* sP;
        VPHI2(:,:,iz) = -vx2p .* sP + vy2p .* cP;
    end

    result.dim = struct();
    result.dim.method         = "DIM-Rayleigh velocity";
    result.dim.x              = x_obs;
    result.dim.y              = y_obs;
    result.dim.z              = z_obs;
    result.dim.X              = X2;
    result.dim.Y              = Y2;

    result.dim.v_x_f1         = VX1;
    result.dim.v_y_f1         = VY1;
    result.dim.v_z_f1         = VZ1;
    result.dim.v_rho_f1       = VR1;
    result.dim.v_phi_f1       = VPHI1;

    result.dim.v_x_f2         = VX2;
    result.dim.v_y_f2         = VY2;
    result.dim.v_z_f2         = VZ2;
    result.dim.v_rho_f2       = VR2;
    result.dim.v_phi_f2       = VPHI2;

    result.dim.block_size     = blk;
    result.dim.src_block_size = src_blk;
    result.dim.f1             = source.f1;
    result.dim.f2             = source.f2;
end

end

% ======================================================================
% helpers
% ======================================================================
function KZ = local_kz_branch(k, kr_col, Nz)
kz = sqrt(k.^2 - kr_col.^2);
idx = imag(kz) < 0;
kz(idx) = -kz(idx);
KZ = kz * ones(1, Nz);
end

function G = local_green_spec_analytic(k, kr_col, z_row, eps_phase, kz_min)
kz = sqrt(k.^2 - kr_col.^2);
idx = imag(kz) < 0;
kz(idx) = -kz(idx);

Nr = numel(kz);
Nz = numel(z_row);

KZ = kz * ones(1, Nz);
Z  = ones(Nr,1) * reshape(z_row, 1, []);

KZs = KZ;
mask_kz0 = abs(KZs) < kz_min;
if any(mask_kz0(:))
    KZs(mask_kz0) = kz_min .* exp(1j * angle(KZs(mask_kz0)));
end

G = (1j/(4*pi)) * exp(1j * KZ .* Z) ./ KZs;

phase = abs(KZ .* Z);
rel_im = abs(imag(KZ)) ./ max(abs(KZ), kz_min);

mask = (phase < eps_phase) & (rel_im < 1e-6);
if any(mask(:))
    Gt = (1./KZs) + 1j*Z - (KZs .* (Z.^2))/2;
    G(mask) = (1j/(4*pi)) * Gt(mask);
end
end

function w = local_trapz_weights_1d(x)
x = x(:);
n = numel(x);
if n < 2
    w = 0;
    return;
end
dx = diff(x);
w = zeros(n,1);
w(1)   = dx(1)/2;
w(end) = dx(end)/2;
if n > 2
    w(2:end-1) = (dx(1:end-1) + dx(2:end))/2;
end
end

function [x_obs, y_obs, z_obs] = local_parse_dim_grid(gd)
if isfield(gd,'x') && ~isempty(gd.x)
    x_obs = gd.x(:).';
else
    error('calc_ultrasound_velocity_field:BadGridX', 'grid.dim.x is required.');
end

if isfield(gd,'y') && ~isempty(gd.y)
    y_obs = gd.y(:).';
else
    error('calc_ultrasound_velocity_field:BadGridY', 'grid.dim.y is required.');
end

if isfield(gd,'z') && ~isempty(gd.z)
    z_obs = gd.z(:).';
else
    error('calc_ultrasound_velocity_field:BadGridZ', 'grid.dim.z is required.');
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