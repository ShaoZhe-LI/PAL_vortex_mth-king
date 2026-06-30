function [source, fht, dim] = make_source_velocity(source, medium, calc)
% =========================================================================
% MAKE_SOURCE_VELOCITY
% - Build consistent aperture normal-velocity models for:
%       (1) FHT (radial samples on xh_v)   and
%       (2) DIM/ASM (2D samples on (x,y))
% - Meanwhile build parameter structs for FHT/DIM that are directly usable
%   in computation.
%
% IMPORTANT UPDATE (multi-frequency, PAL-style):
%   source contains THREE frequencies:
%       source.f1 : ultrasonic carrier (Hz)
%       source.f2 : ultrasonic sideband (Hz) (typically f1 + fa)
%       source.fa : difference/audio (Hz)
%
%   We separate two "concepts":
%   (A) FHT SAMPLING/GRIDS: use f2 as the ultrasonic reference frequency
%       - NH (k_rho truncation) uses w2/c0
%       - DIM discretization (by default) uses lambda(f2)
%   (B) VELOCITY PROFILES: MUST be consistent with BOTH f1 and f2
%       - Build two velocity handles: vs_rho_1 / vs_rho_2  (radial, m/s)
%       - Build two 2D handles:       vn_xy_1  / vn_xy_2   (x,y, m/s)
%       - Sample vectors for FHT:     vs1_on_xh_v, vs2_on_xh_v
%       - Transform for FHT:          Vs1, Vs2 (you can compute outside)
%
% This ensures you can compute:
%   - pressure field at f1 using (vs_rho_1, k1)
%   - pressure field at f2 using (vs_rho_2, k2)
% while still using a single, unified FHT grid (based on f2) for efficiency.
% -------------------------------------------------------------------------
% INPUT
% source (struct):
%   source.profile : 'Uniform' | 'Focus' | 'Vortex-m' | 'Poly' | 'Custom'
%   source.a       : aperture radius (m)
%   source.v0      : velocity amplitude (m/s)
%   source.m       : topological charge (integer), used by 'Vortex-m'
%   source.F       : focus distance (m), used by 'Focus' (optional)
%   source.poly_n  : power, used by 'Poly' (optional)
%   source.custom_vrho_handle_1 : @(rho) radial velocity at f1 (m/s), for 'Custom' (optional)
%   source.custom_vrho_handle_2 : @(rho) radial velocity at f2 (m/s), for 'Custom' (optional)
%
%   source.f1, source.f2, source.fa:
%       - if only source.f1 and source.fa are provided, f2 = f1 + fa
%       - if source.f2 is missing, it will be generated
%
% medium (struct, optional):
%   medium.c0    : sound speed (m/s), default 343
%   medium.rho0  : density (kg/m^3), default 1.21
%   medium.beta  : nonlinearity coefficient, default 1.2
%   medium.pref  : reference pressure (Pa), default 2e-5
%   medium.temp/humi/p0 : optional
%
% calc (struct):
%   calc.fht (struct): parameters for FHT sampling / grids
%       .N_FHT       : number of FHT points (default 32768)
%       .rho_max     : rho truncation max for field computation (m) (default 2)
%       .Nh_scale    : Nh = Nh_scale * rho_max (default 1.2)
%       .NH_scale    : NH = NH_scale * (w2/c0) (default 4)  <-- based on f2
%       .Nh_v_scale  : Nh_v = Nh_v_scale * a (default 1.1)
%       .zu_max      : max z for ultrasonic (m) (default 15)
%       .za_max      : max z for audio (m) (default 4)
%       .solve_handle: function handle to solve_kappa0 (default @solve_kappa0)
%       .z_sampling  : 'uniform' | 'log' (default 'uniform')
%       .Nz_ultra    : total number of z_ultra points when using 'log'
%                      default = numel(0:delta:zu_max)
%       .z_log_min   : minimum positive z used in log sampling
%                      default = delta
%
%       delta rule:
%       - ultrasonic field delta is set by audio wavelength: delta = (c0/fa)/8
%         (your requirement)
%
%   calc.dim (struct): parameters for DIM source discretization
%       .use_freq    : 'f1' | 'f2' | 'fa' (default 'f2')
%       .dis_coe     : points-per-wavelength (default 16) => dx=lambda/dis_coe
%       .dx/.dy      : override dx/dy directly (m) (optional)
%       .margin      : extra margin added to bounding box (default 1 grid step)
%
% OUTPUT
% source (struct) adds:
%   source.medium
%   source.w1,w2,wa
%   source.k1,k2,ka
%   source.v_ratio               : ratio of v0 between at f2 and f1 (default 1)
%   source.m_used                : actually-used m (0 except 'Vortex-m')
%   source.vs_rho_1 / vs_rho_2   : @(rho) radial normal velocity at f1/f2 [m/s], disk-masked
%   source.vn_xy_1  / vn_xy_2    : @(X,Y) 2D normal velocity at f1/f2 [m/s]
%   source.meta                  : record of resolved parameters
%
% fht (struct) adds:
%   fht.N_FHT, fht.delta, fht.Nh, fht.NH, fht.Nh_v, fht.NH_v
%   fht.n_FHT, fht.a_solve, fht.k0, fht.x1, fht.x0
%   fht.xh, fht.xH, fht.xh_v, fht.xH_v
%   fht.z_ultra, fht.z_audio
%   fht.freq_ultra = f2, fht.freq_audio = fa
%
% dim (struct) adds:
%   dim.use_freq, dim.lambda, dim.dx, dim.dy, dim.dA
%   dim.make_grid : helper to create rectangular grid covering the disk
%
% NOTES
% - Velocity outputs (vs_rho_1/2, vn_xy_1/2) are in m/s and already include v0.
% - In DIM Rayleigh integral, do NOT multiply v0 again.
% - This function only prepares inputs; it does not compute fields.
% =========================================================================

% -------------------- validate inputs --------------------
if nargin < 1 || ~isstruct(source)
    error('make_source_velocity:BadInput', 'Input must be a struct "source".');
end
if ~isfield(source,'profile') || isempty(source.profile)
    error('make_source_velocity:MissingProfile', 'source.profile is required.');
end

if nargin < 2 || isempty(medium); medium = struct(); end
if ~isstruct(medium); error('make_source_velocity:BadMedium', 'medium must be a struct.'); end

if nargin < 3 || isempty(calc); calc = struct(); end
if ~isstruct(calc); error('make_source_velocity:BadCalc', 'calc must be a struct.'); end

% -------------------- for bookkeeping --------------------
meta = struct();
meta.source_in = source;

% -------------------- medium defaults --------------------
if ~isfield(medium,'c0')           || isempty(medium.c0);           medium.c0 = 343; end
if ~isfield(medium,'rho0')         || isempty(medium.rho0);         medium.rho0 = 1.21; end
if ~isfield(medium,'beta')         || isempty(medium.beta);         medium.beta = 1.2; end
if ~isfield(medium,'pref')         || isempty(medium.pref);         medium.pref = 2e-5; end
if ~isfield(medium,'use_absorp')   || isempty(medium.use_absorp);   medium.use_absorp = false; end
if ~isfield(medium,'atten_handle') || isempty(medium.atten_handle); medium.atten_handle = []; end
source.medium = medium;

% -------------------- source defaults (geometry / mode) --------------------
if ~isfield(source,'a')        || isempty(source.a);        source.a = 0.1; end
if ~isfield(source,'v0')       || isempty(source.v0);       source.v0 = 0.172; end
if ~isfield(source,'v_ratio')  || isempty(source.v_ratio);  source.v_ratio = 1; end
if ~isfield(source,'m')        || isempty(source.m);        source.m = 3; end
if ~isfield(source,'F')        || isempty(source.F);        source.F = 1; end
if ~isfield(source,'poly_n')   || isempty(source.poly_n);   source.poly_n = 2; end

% -------------------- frequencies: f1, f2, fa --------------------
if ~isfield(source,'fa') || isempty(source.fa)
    error('make_source_velocity:MissingFa', 'source.fa (difference/audio frequency) is required.');
end

if ~isfield(source,'f1') || isempty(source.f1)
    if isfield(source,'fu') && ~isempty(source.fu)
        source.f1 = source.fu;
    else
        error('make_source_velocity:MissingF1', 'Provide source.f1 (or source.fu).');
    end
end

if ~isfield(source,'f2') || isempty(source.f2)
    source.f2 = source.f1 + source.fa;
end

% Angular frequencies
source.w1 = 2*pi*source.f1;
source.w2 = 2*pi*source.f2;
source.wa = 2*pi*source.fa;

% --- absorption / attenuation (optional) ---
alpha1 = 0; alpha2 = 0; alphaa = 0;

if isfield(medium,'use_absorp') && medium.use_absorp
    if isfield(medium,'atten_handle') && isa(medium.atten_handle,'function_handle')
        alpha1 = medium.atten_handle(source.f1);
        alpha2 = medium.atten_handle(source.f2);
        alphaa = medium.atten_handle(source.fa);
    else
        error('make_source_velocity:MissingAttenHandle', ...
            'medium.use_absorp=true but medium.atten_handle is missing/not a function handle.');
    end
end

source.k1 = source.w1/medium.c0 + 1j*alpha1;
source.k2 = source.w2/medium.c0 + 1j*alpha2;
source.ka = source.wa/medium.c0 + 1j*alphaa;


% -------------------- normalize profile label --------------------
% Make profile robust to string type / spaces / capitalization
profile = lower(strtrim(char(source.profile)));

a  = source.a;
v0 = source.v0;
m0 = source.m;     % original user-provided m
F  = source.F;

% =========================================================================
% (1) Build TWO radial velocity cores: one for f1, one for f2
%     - For uniform/vortex/poly: same at both frequencies
%     - For focus: phase uses Re(k1) for f1 and Re(k2) for f2
%     - For custom: allow separate handles for f1/f2; if only one provided, reuse it
% =========================================================================

switch profile
    case 'uniform'
        % v_n(rho) = v0 (both f1 and f2)
        vrho_core_1 = @(rho) v0 + 0*rho;
        vrho_core_2 = @(rho) source.v_ratio * v0 + 0*rho;
        m_use = 0;

    case 'vortex-m'
        % radial part uniform; azimuth exp(i*m*phi) applied later
        vrho_core_1 = @(rho) v0 + 0*rho;
        vrho_core_2 = @(rho) source.v_ratio * v0 + 0*rho;
        m_use = source.m;

    case 'focus'
        % v_n(rho;fi) = v0 * exp(-i*Re(ki)*sqrt(rho^2+F^2))
        kr1 = real(source.k1);
        kr2 = real(source.k2);
        vrho_core_1 = @(rho) v0 .* exp(-1i * kr1 .* (sqrt(rho.^2 + F.^2) - F));
        vrho_core_2 = @(rho) source.v_ratio * v0 .* exp(-1i * kr2 .* (sqrt(rho.^2 + F.^2) - F));
        m_use = 0;

    case 'poly'
        % v_n(rho) = v0 * (rho/a)^n (both f1 and f2)
        n = source.poly_n;
        vrho_core_1 = @(rho) v0 .* (max(rho,0)./a).^n;
        vrho_core_2 = @(rho) source.v_ratio * v0 .* (max(rho,0)./a).^n;
        m_use = 0;

        case 'custom'
        % =============================================================
        % Custom supports TWO kinds of user inputs:
        %   (A) radial handles: custom_vrho_handle_1/2 (or custom_vrho_handle)
        %   (B) 2D handles:     custom_vn_xy_handle_1/2
        %
        % Priority:
        %   - If vn_xy handles exist -> use them for DIM (and set radial to 0 placeholder)
        %   - Else -> require vrho handle(s) as before
        % =============================================================

        % ---- check 2D handles ----
        has_vn1 = isfield(source,'custom_vn_xy_handle_1') && isa(source.custom_vn_xy_handle_1,'function_handle');
        has_vn2 = isfield(source,'custom_vn_xy_handle_2') && isa(source.custom_vn_xy_handle_2,'function_handle');

        % ---- check radial handles ----
        has1 = isfield(source,'custom_vrho_handle_1') && isa(source.custom_vrho_handle_1,'function_handle');
        has2 = isfield(source,'custom_vrho_handle_2') && isa(source.custom_vrho_handle_2,'function_handle');
        has_single = isfield(source,'custom_vrho_handle') && isa(source.custom_vrho_handle,'function_handle');

        if (has_vn1 || has_vn2)
            % ========== Mode B: arbitrary vn(x,y) ==========
            % Provide radial "placeholder" (zeros) so FHT-related code won't crash.
            vrho_core_1 = @(rho) 0*rho;
            vrho_core_2 = @(rho) 0*rho;

            % Choose m
            if isfield(source,'m_custom') && ~isempty(source.m_custom)
                m_use = source.m_custom;
            else
                m_use = 0;
            end

            % Build vn_xy directly from user handles, masked by disk
            if has_vn1
                vn_xy_1 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_1);
            else
                % if only vn2 provided, reuse it
                vn_xy_1 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_2);
            end

            if has_vn2
                vn_xy_2 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_2);
            else
                % if only vn1 provided, reuse it
                vn_xy_2 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_1);
            end

        else
            % ========== Mode A: radial vrho(r) (legacy behavior) ==========
            if ~(has1 || has2 || has_single)
                error('make_source_velocity:MissingCustomHandle', ...
                    ['For profile="Custom", provide one of:\n', ...
                     '  (A) source.custom_vrho_handle (single), or\n', ...
                     '      source.custom_vrho_handle_1 and/or source.custom_vrho_handle_2, or\n', ...
                     '  (B) source.custom_vn_xy_handle_1 and/or source.custom_vn_xy_handle_2.\n']);
            end

            if has_single
                vrho_core_1 = source.custom_vrho_handle;
                vrho_core_2 = source.custom_vrho_handle;
            else
                if has1; vrho_core_1 = source.custom_vrho_handle_1; end
                if has2; vrho_core_2 = source.custom_vrho_handle_2; end
                if ~has1 && has2; vrho_core_1 = vrho_core_2; end
                if has1 && ~has2; vrho_core_2 = vrho_core_1; end
            end

            if isfield(source,'m_custom') && ~isempty(source.m_custom)
                m_use = source.m_custom;
            else
                m_use = 0;
            end

            % vn_xy will be constructed later by the common code:
            % vn_xy = vs_rho(rho) .* exp(1i*m*phi)
        end

    otherwise
        error('make_source_velocity:UnknownProfile', ...
            'Unknown source.profile "%s". Use Uniform|Focus|Vortex-m|Poly|Custom.', profile);
end

% Store actually-used m (for clarity)
source.m_used = m_use;

% =========================================================================
% (2) Build LOCAL velocity handles (f1/f2) for later FHT/DIM sampling
%   - Default: vn(x,y)=vs_rho(rho)*exp(j*m*phi)  (axisymmetric / vortex)
%   - If provided: use arbitrary 2D handle custom_vn_xy_handle_1/2
% =========================================================================
vs_rho_1 = @(rho) local_vrho_masked(rho, a, vrho_core_1);
vs_rho_2 = @(rho) local_vrho_masked(rho, a, vrho_core_2);

% If vn_xy_1/2 already defined in Custom-mode-B, keep them.
if ~exist('vn_xy_1','var') || ~exist('vn_xy_2','var')
    has_vn1 = isfield(source,'custom_vn_xy_handle_1') && isa(source.custom_vn_xy_handle_1,'function_handle');
    has_vn2 = isfield(source,'custom_vn_xy_handle_2') && isa(source.custom_vn_xy_handle_2,'function_handle');

    if has_vn1
        vn_xy_1 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_1);
    else
        vn_xy_1 = @(X,Y) local_vn_xy_eval(X, Y, vs_rho_1, source.m_used);
    end

    if has_vn2
        vn_xy_2 = @(X,Y) local_vn_xy_masked(X, Y, a, source.custom_vn_xy_handle_2);
    else
        vn_xy_2 = @(X,Y) local_vn_xy_eval(X, Y, vs_rho_2, source.m_used);
    end
end

% =========================================================================
% (3) Build FHT parameter struct
%     - Sampling/Truncation follows f2 (your requirement)
%     - delta follows fa wavelength / 8 (your requirement)
% =========================================================================
if ~isfield(calc,'fht') || isempty(calc.fht); calc.fht = struct(); end
cf = calc.fht;

if ~isfield(cf,'N_FHT')      || isempty(cf.N_FHT);      cf.N_FHT = 32768; end
if ~isfield(cf,'rho_max')    || isempty(cf.rho_max);    cf.rho_max = 2; end
if ~isfield(cf,'Nh_scale')   || isempty(cf.Nh_scale);   cf.Nh_scale = 1.2; end
if ~isfield(cf,'NH_scale')   || isempty(cf.NH_scale);   cf.NH_scale = 1.2; end
if ~isfield(cf,'Nh_v_scale') || isempty(cf.Nh_v_scale); cf.Nh_v_scale = cf.Nh_scale; end
if ~isfield(cf,'zu_max')     || isempty(cf.zu_max);     cf.zu_max = 15; end
if ~isfield(cf,'za_max')     || isempty(cf.za_max);     cf.za_max = 4; end
if ~isfield(cf,'solve_handle') || isempty(cf.solve_handle)
    cf.solve_handle = @solve_kappa0;
end

fht = struct();
fht.N_FHT     = cf.N_FHT;
fht.Nh_scale  = cf.Nh_scale;
fht.NH_scale  = cf.NH_scale;
fht.Nh_v_scale= cf.Nh_v_scale;
fht.rho_max   = cf.rho_max;
fht.zu_max    = cf.zu_max;
fht.za_max    = cf.za_max;

% delta: priority = external calc.fht.delta > default lambda_a/8
lambda_2 = medium.c0 / source.f2;
if isfield(cf,'delta') && ~isempty(cf.delta)
    fht.delta = cf.delta;
    fht.delta_source = 'external calc.fht.delta';
else
    fht.delta = lambda_2 / 4;
    fht.delta_source = 'default lambda_2/4';
end

% Truncation lengths:
%   Nh by rho_max; NH by ultrasonic reference f2 (w2/c0)
fht.Nh = cf.Nh_scale * cf.rho_max;
fht.NH = cf.NH_scale * (source.w2/medium.c0);

% Indices
fht.n_FHT = (0:fht.N_FHT-1);

% Parameters from solve_kappa0
[fht.a_solve, fht.k0, fht.x1, fht.x0] = cf.solve_handle(fht.N_FHT, fht.n_FHT);

% Sampling grids in rho and k_rho
fht.xh = (fht.x1 * fht.Nh).';
fht.xH = (fht.x1 * fht.NH).';

% Velocity transform grids (tied to aperture)
fht.Nh_v = cf.Nh_v_scale * source.a;
fht.NH_v = fht.NH;  % keep same k_rho truncation as field (based on f2)
fht.xh_v = (fht.x1 * fht.Nh_v).';
fht.xH_v = (fht.x1 * fht.NH_v).';

% -------------------- z sampling mode --------------------
if ~isfield(cf,'z_sampling') || isempty(cf.z_sampling)
    cf.z_sampling = 'uniform';   % 'uniform' | 'log'
end
fht.z_sampling = lower(strtrim(char(cf.z_sampling)));

if ~isfield(cf,'Nz_ultra') || isempty(cf.Nz_ultra)
    % default total point count consistent with old uniform grid
    cf.Nz_ultra = numel(0:fht.delta:fht.zu_max);
end
fht.Nz_ultra_target = cf.Nz_ultra;

if ~isfield(cf,'z_log_min') || isempty(cf.z_log_min)
    cf.z_log_min = fht.delta;   % default positive start for log grid
end
fht.z_log_min = cf.z_log_min;

switch fht.z_sampling
    case 'uniform'
        % keep exactly the old behavior
        fht.z_ultra = 0:fht.delta:fht.zu_max;

    case 'log'
        if fht.zu_max <= 0
            error('make_source_velocity:BadZuMax', ...
                'For log z sampling, calc.fht.zu_max must be > 0.');
        end

        if fht.Nz_ultra_target < 2
            error('make_source_velocity:BadNzUltra', ...
                'For log z sampling, calc.fht.Nz_ultra must be >= 2.');
        end

        zmin = max(fht.z_log_min, eps);
        if zmin >= fht.zu_max
            error('make_source_velocity:BadZLogMin', ...
                'calc.fht.z_log_min must be smaller than calc.fht.zu_max.');
        end

        zpos = logspace(log10(zmin), log10(fht.zu_max), fht.Nz_ultra_target - 1);
        fht.z_ultra = unique([0, zpos]);

    otherwise
        error('make_source_velocity:BadZSampling', ...
            'calc.fht.z_sampling must be ''uniform'' or ''log''.');
end

% audio z grid remains unchanged
fht.z_audio = 0:fht.delta:fht.za_max;

fht.Nz_ultra = numel(fht.z_ultra);
fht.Nz_audio = numel(fht.z_audio);

% record which frequencies control the grids
fht.freq_ultra = source.f2; % sampling based on f2
fht.freq_audio = source.fa; % delta based on fa

% =========================================================================
% (4) Build DIM discretization struct  (and PRE-SAMPLED source plane points)
% - Mimic the FHT part: provide not only step sizes, but also sampled (X,Y)
%   coordinates and corresponding velocity samples at f1/f2.
% =========================================================================
if ~isfield(calc,'dim') || isempty(calc.dim); calc.dim = struct(); end
cd = calc.dim;

% ---- init DIM struct ONCE (do not re-init later) ----
dim = struct();

% --- source discretization type for DIM-Rayleigh ---
if ~isfield(cd,'src_discretization') || isempty(cd.src_discretization)
    cd.src_discretization = 'cart';   % 'cart' (default) | 'polar'
end
dim.src_discretization = lower(strtrim(char(cd.src_discretization)));

% --- wavelengths for reference ---
dim.lambda_f1 = medium.c0 / source.f1;
dim.lambda_f2 = medium.c0 / source.f2;
dim.lambda_fa = medium.c0 / source.fa;

if ~isfield(calc.dim,'dis_coe') || isempty(calc.dim.dis_coe)
    calc.dim.dis_coe = 16;   % default
end
dim.dis_coe = calc.dim.dis_coe;

% --- choose which frequency controls discretization step (default: f2) ---
if ~isfield(cd,'use_freq') || isempty(cd.use_freq)
    cd.use_freq = 'f2';   % 'f1' | 'f2' | 'fa'
end
dim.use_freq = cd.use_freq;

switch lower(strtrim(cd.use_freq))
    case 'f1'
        dim.lambda = dim.lambda_f1;
    case 'fa'
        dim.lambda = dim.lambda_fa;
    otherwise
        dim.lambda = dim.lambda_f2;  % f2
end

% --- step size selection (dx,dy) ---
if isfield(cd,'dx') && ~isempty(cd.dx)
    dim.dx = cd.dx;
else
    if ~isfield(cd,'dis_coe') || isempty(cd.dis_coe); cd.dis_coe = 16; end
    dim.dx = dim.lambda / cd.dis_coe;   % points-per-wavelength
end
if isfield(cd,'dy') && ~isempty(cd.dy)
    dim.dy = cd.dy;
else
    dim.dy = dim.dx;
end
dim.dA = dim.dx * dim.dy;

% --- source plane center and z (default consistent with your script) ---
if ~isfield(cd,'center_x') || isempty(cd.center_x); cd.center_x = 0; end
if ~isfield(cd,'center_y') || isempty(cd.center_y); cd.center_y = 0; end
if ~isfield(cd,'z0')       || isempty(cd.z0);       cd.z0 = 0; end
dim.center_x = cd.center_x;
dim.center_y = cd.center_y;
dim.z0       = cd.z0;

% --- margin used when creating rectangular grid covering the disk ---
if ~isfield(cd,'margin') || isempty(cd.margin); cd.margin = 1; end
dim.margin = cd.margin;

% --- 1D grids covering the disk (bbox) ---
[x_grid, y_grid] = local_make_disk_bbox_grid(dim.center_x, dim.center_y, ...
    source.a, dim.dx, dim.dy, dim.margin);
dim.x_grid = x_grid;
dim.y_grid = y_grid;

% --- full meshgrid (for matrix-style ASM/DIM) ---
[dim.X, dim.Y] = meshgrid(dim.x_grid, dim.y_grid);      % size Ny x Nx
dim.RHO = hypot(dim.X - dim.center_x, dim.Y - dim.center_y);
dim.PHI = atan2(dim.Y - dim.center_y, dim.X - dim.center_x);

% --- aperture mask on the grid ---
dim.mask = (dim.RHO <= source.a);

% --- velocity samples on the GRID (matrix form, Ny x Nx) ---
dim.Vn_grid_f1 = vn_xy_1(dim.X - dim.center_x, dim.Y - dim.center_y);
dim.Vn_grid_f2 = vn_xy_2(dim.X - dim.center_x, dim.Y - dim.center_y);

% enforce exact zero outside disk (numerical cleanliness)
dim.Vn_grid_f1(~dim.mask) = 0;
dim.Vn_grid_f2(~dim.mask) = 0;

% --- sampled source points inside the disk (vector form) ---
switch dim.src_discretization
    case 'cart'
        dim.Xs = dim.X(dim.mask);
        dim.Ys = dim.Y(dim.mask);
        dim.Zs = dim.z0 * ones(size(dim.Xs));
        dim.num_source_points = numel(dim.Xs);

        dim.Vn_pts_f1 = dim.Vn_grid_f1(dim.mask);
        dim.Vn_pts_f2 = dim.Vn_grid_f2(dim.mask);

        dim.dA_pts = dim.dA * ones(dim.num_source_points, 1);

    case 'polar'
        dr = dim.dx;  % dr = lambda/dis_coe, e.g., lambda/16
        [Xs, Ys, dA_pts] = local_make_disk_polar_cells(dim.center_x, dim.center_y, source.a, dr, source.m_used);

        dim.Xs = Xs;
        dim.Ys = Ys;
        dim.Zs = dim.z0 * ones(size(dim.Xs));
        dim.num_source_points = numel(dim.Xs);

        % velocity evaluated directly on polar cell centers (NO interpolation)
        dim.Vn_pts_f1 = vn_xy_1(dim.Xs - dim.center_x, dim.Ys - dim.center_y);
        dim.Vn_pts_f2 = vn_xy_2(dim.Xs - dim.center_x, dim.Ys - dim.center_y);

        dim.dA_pts = dA_pts;

    otherwise
        error('make_source_velocity:BadSrcDiscretization', ...
            'calc.dim.src_discretization must be ''cart'' or ''polar''.');
end

% --- keep helper for custom centers if you need later ---
dim.make_grid = @(cx,cy) local_make_disk_bbox_grid(cx, cy, source.a, dim.dx, dim.dy, dim.margin);

% =========================================================================
% (5) FULL ISOLATION PACK
% - FHT namespace: ONLY rho-handle + rho sampling vectors (no 2D)
% - DIM namespace: ONLY 2D handle + 2D sampled grids/velocities (no rho)
% NOTE:
%   - "fht" (output argument) stores ALL numerical grids/parameters
%   - "source.fht" stores ONLY source-related quantities (velocity, sampling)
% =========================================================================

% -------------------- pack FHT: rho only --------------------
% ensure NH_v exists (robustness for future refactor)
if ~isfield(fht,'NH_v') || isempty(fht.NH_v)
    fht.NH_v = fht.NH;
end
source.fht = struct();

% (a) rho velocity handles (f1/f2)
source.fht.vs_rho_1 = vs_rho_1;     % @(rho) masked radial v_n at f1
source.fht.vs_rho_2 = vs_rho_2;     % @(rho) masked radial v_n at f2

% (b) rho sampling grid used by FHT velocity transform
source.fht.xh_v = fht.xh_v;

% (c) sampled radial velocity vectors on xh_v (ready for m_FHT)
source.fht.vs1_on_xh_v = source.fht.vs_rho_1(source.fht.xh_v);
source.fht.vs2_on_xh_v = source.fht.vs_rho_2(source.fht.xh_v);

% (optional bookkeeping)
source.fht.Nh_v = fht.Nh_v;
source.fht.NH_v = fht.NH_v;

% -------------------- pack DIM: 2D only --------------------
source.dim = struct();

source.dim.src_discretization = dim.src_discretization;

% (a) 2D velocity handles (f1/f2)
source.dim.vn_xy_1 = vn_xy_1;       % @(X,Y) v_n(x,y) for f1
source.dim.vn_xy_2 = vn_xy_2;       % @(X,Y) v_n(x,y) for f2

% (b) sampled 2D grids (matrix form)
source.dim.X    = dim.X;
source.dim.Y    = dim.Y;
source.dim.mask = dim.mask;
source.dim.x_grid = dim.x_grid;
source.dim.y_grid = dim.y_grid;

% (c) sampled 2D velocities (matrix form)
source.dim.Vn_grid_f1 = dim.Vn_grid_f1;   % 2D Cartesian grid velocity (for ASM / FFT)
source.dim.Vn_grid_f2 = dim.Vn_grid_f2;   % same velocity as Vn_pts_f2, but stored on grid (FFT-friendly)

% sampled disk points (vector form, for direct integration)
source.dim.Xs = dim.Xs;                   % x-coordinates of disk points
source.dim.Ys = dim.Ys;                   % y-coordinates of disk points
source.dim.Zs = dim.Zs;                   % z-coordinates of disk points (source plane)
source.dim.Vn_pts_f1 = dim.Vn_pts_f1;     % velocity at disk points (for Rayleigh integral)
source.dim.Vn_pts_f2 = dim.Vn_pts_f2;     % velocity at disk points (for Rayleigh integral)
source.dim.dis_coe = dim.dis_coe;
source.dim.num_source_points = dim.num_source_points; % number of disk points

% (e) discretization parameters
source.dim.dx = dim.dx;
source.dim.dy = dim.dy;
source.dim.dA = dim.dA;
source.dim.dA_pts = dim.dA_pts;   % Nsrc x 1 (for Rayleigh, supports polar)
source.dim.lambda = dim.lambda;
source.dim.use_freq = dim.use_freq;
source.dim.make_grid = dim.make_grid;

% =========================================================================
% Remove any legacy/top-level fields to enforce isolation
% (Your current code no longer creates top-level vs_rho/vn_xy/vs*_on_xh_v,
%  but keep this for safety if you iterate.)
% =========================================================================
to_remove = {};
% legacy cleanup (kept for safety during refactoring; normally empty)
legacy = {'vs_rho_1','vs_rho_2','vn_xy_1','vn_xy_2','vs1_on_xh_v','vs2_on_xh_v', ...
    'dim_X','dim_Y','dim_mask','dim_Vn_grid_f1','dim_Vn_grid_f2', ...
    'dim_Xs','dim_Ys','dim_Zs','dim_Vn_pts_f1','dim_Vn_pts_f2'};
for i = 1:numel(legacy)
    if isfield(source, legacy{i}); to_remove{end+1} = legacy{i}; end %#ok<AGROW>
end
if ~isempty(to_remove)
    source = rmfield(source, unique(to_remove));
end

% (You can immediately do:)
%   Vs1 = m_FHT(source.fht.vs1_on_xh_v, ...
%               fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
%               fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);
%   Vs2 = m_FHT(source.fht.vs2_on_xh_v, ...
%               fht.N_FHT, 1, fht.Nh_v, fht.NH_v, ...
%               fht.a_solve, fht.x0, fht.x1, fht.k0, source.m_used);


% =========================================================================
% meta (for reproducibility & method-level decisions)
% =========================================================================

% ---- profile & geometry ----
meta.profile_in = char(source.profile);   % user input (raw)
meta.profile    = profile;                % normalized label
meta.a          = source.a;               % aperture radius
meta.v0         = source.v0;              % velocity amplitude
meta.m_original = m0;                     % user-provided m
meta.m_used     = source.m_used;           % actually used m

% ---- frequencies (PAL) ----
meta.f1 = source.f1;                      % carrier
meta.f2 = source.f2;                      % sideband
meta.fa = source.fa;                      % difference frequency
meta.w1 = source.w1;
meta.w2 = source.w2;
meta.wa = source.wa;
meta.k1 = source.k1;
meta.k2 = source.k2;
meta.ka = source.ka;

% ---- medium ----
meta.medium = medium;

% ---- numerical design decisions ----
meta.fht_grid_based_on   = 'f2';           % FHT sampling/truncation uses f2
meta.delta_based_on      = 'fa/8';         % z-step based on audio wavelength
meta.dim_discretize_freq = dim.use_freq;   % which freq controls DIM dx,dy
meta.use_absorp = medium.use_absorp;       % whether absorption is used
meta.fht_delta        = fht.delta;
meta.fht_delta_source = fht.delta_source;

% ---- method separation (important) ----
meta.method_separation = struct( ...
    'FHT', 'radial-only (rho)', ...
    'DIM', '2D-only (x,y)' );

source.meta = meta;


end

% ======================================================================
% helpers (local functions)
% ======================================================================
function v = local_vrho_masked(rho, a, vrho_core)
rho = abs(rho);
v = zeros(size(rho));
idx = rho <= a;
if any(idx(:))
    v(idx) = vrho_core(rho(idx));
end
end

function vxy = local_vn_xy_eval(X, Y, vs_rho, m)
rho = hypot(X, Y);
phi = atan2(Y, X);
vxy = vs_rho(rho) .* exp(1i * m * phi);
end

function [x_grid, y_grid] = local_make_disk_bbox_grid(cx, cy, a, dx, dy, margin_steps)
mx = margin_steps * dx;
my = margin_steps * dy;
x_min = cx - a - mx;
x_max = cx + a + mx;
y_min = cy - a - my;
y_max = cy + a + my;
x_grid = x_min:dx:x_max;
y_grid = y_min:dy:y_max;
end

function [Xs, Ys, dA_pts] = local_make_disk_polar_cells(cx, cy, a, dr, m)
% ring-by-ring polar discretization:
% r: fixed step dr
% theta: adaptive sectors by outer radius (arc length <= dr),
%        AND enforced to be a multiple of m (for exact angular cancellation)
%
% outputs:
%   Xs, Ys   : cell-center coordinates (Ns x 1)
%   dA_pts   : corresponding cell areas (Ns x 1)

Nr = ceil(a/dr);
r_edges = (0:Nr) * dr;
r_edges(end) = a;  % clamp to aperture radius

Xs = [];
Ys = [];
dA_pts = [];

m_use = abs(m);

for i = 1:Nr
    r_in  = r_edges(i);
    r_out = r_edges(i+1);
    if r_out <= r_in
        continue;
    end

    % ring center
    r_c = 0.5 * (r_in + r_out);

    % ----------------------------------------------------------
    % 1) geometric requirement: arc length at outer edge <= dr
    % 2) minimum angular resolution
    % 3) enforce Ntheta to be a multiple of m (if m ~= 0)
    % ----------------------------------------------------------
    Nth_geo = ceil(2*pi*r_out/dr);
    Nth0    = max(32, Nth_geo);   % base requirement

    if m_use == 0
        Nth = Nth0;        
    else
        Nth = ceil(Nth0 / abs(m_use)) * m_use;        
    end

    dth = 2*pi / Nth;

    % cell-center angles
    th_c = ((0:Nth-1) + 0.5) * dth;

    % Cartesian coordinates of cell centers
    xs = cx + r_c * cos(th_c(:));
    ys = cy + r_c * sin(th_c(:));

    % exact sector area
    dA = 0.5 * (r_out^2 - r_in^2) * dth;
    dA = dA * ones(size(xs));

    % accumulate
    Xs = [Xs; xs];
    Ys = [Ys; ys];
    dA_pts = [dA_pts; dA];
end
end

function vxy = local_vn_xy_masked(X, Y, a, vn_handle)
rho = hypot(X, Y);
vxy = zeros(size(rho));
mask = (rho <= a);
if any(mask(:))
    vxy(mask) = vn_handle(X(mask), Y(mask));
end
end
