function H = m_FHT(h, N, Nz, Nh, NH, alpha, x0, x1, k0, miu)
% =========================================================================
% m_FHT
% - miu-order Hankel transform (fast Hankel transform) on columns of h
% =========================================================================

% ---- handle negative integer order ----
if miu < 0
    miu_abs = abs(miu);

    H = m_FHT(h, N, Nz, Nh, NH, alpha, x0, x1, k0, miu_abs);

    % J_{-m}(x) = (-1)^m J_m(x)
    H = (-1)^miu_abs * H;

    return;
end

% ---- shape checks ----
x1 = x1(:).';  % ensure row
assert(numel(x1) == N, ...
    'm_FHT: x1 must have length N (%d), but got %d.', N, numel(x1));
assert(size(h,1) == N, ...
    'm_FHT: h must have N rows (%d).', N);

% ---- constants ----
Nf = Nh * NH;

% exp(alpha*(n+1-N)), n=0..2N-1 or 0..N-2
n2 = (0:2*N-1).' + 1 - N;         % (2N x 1)
n1 = (0:N-2).'   + 1 - N;         % (N-1 x 1)

E2 = exp(alpha * n2);             % (2N x 1)
E1 = exp(alpha * n1);             % (N-1 x 1)

% ---- allocate ----
phi = zeros(2*N, Nz);

% ---- build phi (difference form) ----
dh1 = (E1 ./ x1(1:end-1).').^miu;  % (N-1 x 1)
dh2 = (E1 ./ x1(2:end).').^miu;    % (N-1 x 1)

h1 = h(1:end-1,:) .* E1 .* dh1;
h2 = h(2:end,:)   .* E1 .* dh2;

phi(1:N-1,:) = h1 - h2;
phi(1,:)     = k0 * phi(1,:);

phi(N,:)     = h(end,:) ./ (x1(end)^miu);

% ---- FFT pipeline ----
Phi = fft(phi);                    % (2N x Nz)

jnn = Nf * x0 * E2;                 % (2N x 1)
jn0 = besselj(miu + 1, jnn);        % (2N x 1)

J   = ifft(jn0 * ones(1, Nz));      % (2N x Nz)

G   = fft(Phi .* J);               % (2N x Nz)

% ---- output scaling ----
H = G(1:N,:) ./ (Nf * x1.') * Nh^2;

end