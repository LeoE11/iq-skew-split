function [Z, T, truth] = simulate_link(varargin)
%SIMULATE_LINK  Synthetic coherent link with known Tx and Rx I/Q imbalance.
%
%   [Z, T, truth] = simulate_link('nPol', 2, 'nCap', 20)
%
% Builds K captures of a 16-QAM signal (one tributary per polarisation),
% applies a skew and gain imbalance at each transmitter, a random carrier phase
% (single polarisation) or polarisation state (dual) per capture, then a skew,
% gain and quadrature imbalance at each receive pair, and returns symbol-rate
% samples ready for iq_skew_split.
%
% The defects are applied at 2 samples per symbol with exact frequency-domain
% fractional delays, so nothing here shares algebra with the estimator.
%
% OPTIONS (name, value; defaults in brackets)
%   'nPol'    1 or 2                                   [2]
%   'nCap'    number of captures                       [20]
%   'Ns'      symbols per capture                      [32768]
%   'Br'      symbol rate, Hz                          [64e9]
%   'tauRx'   Rx skew per receive pair, s              [3.6e-12 -2.1e-12]
%   'gRx'     Rx gain imbalance per pair, dB           [0.35 -0.12]
%   'qRx'     Rx quadrature error per pair, deg        [0 0]
%   'tauTx'   Tx skew per tributary, s                 [1.4e-12 -0.85e-12]
%   'gTx'     Tx gain imbalance per tributary, dB      [0.22 -0.31]
%   'snr'     SNR at the receiver, dB (Inf: no noise)  [30]
%   'seed'    random seed                              [7]
%
% Convention: skew > 0 delays Q relative to I; gain is 20*log10(|Q|/|I|);
% a quadrature error phi makes the Q axis read Q*cos(phi) + I*sin(phi).

p = struct('nPol', 2, 'nCap', 20, 'Ns', 32768, 'Br', 64e9, ...
           'tauRx', [3.6e-12 -2.1e-12], 'gRx', [0.35 -0.12], 'qRx', [0 0], ...
           'tauTx', [1.4e-12 -0.85e-12], 'gTx', [0.22 -0.31], ...
           'snr', 30, 'seed', 7);
for i = 1:2:numel(varargin)
    if ~isfield(p, varargin{i})
        error('simulate_link:option', 'Unknown option ''%s''.', varargin{i});
    end
    p.(varargin{i}) = varargin{i+1};
end
n = p.nPol;
p.tauRx = p.tauRx(1:n); p.gRx = p.gRx(1:n); p.qRx = p.qRx(1:n);
p.tauTx = p.tauTx(1:n); p.gTx = p.gTx(1:n);

rng(p.seed);
ovr = 2;
pam = @() 2*randi([0 3], p.Ns, n) - 3;
T = pam() + 1i*pam();
T = T / sqrt(mean(abs(T(:)).^2));

Z = cell(1, p.nCap);
for k = 1:p.nCap
    % Transmitter: one I/Q modulator per tributary.
    xt = zeros(p.Ns*ovr, n);
    for r = 1:n
        s = upfir(T(:, r), ovr);
        g = 10^(p.gTx(r)/20);
        xt(:, r) = real(s) + 1i*g*fracdel(imag(s), p.tauTx(r)*p.Br*ovr);
    end

    % Fibre: a carrier phase or a polarisation state that changes per capture.
    % The common phase c is what spreads lambda; the mixing angle a is kept
    % between 23 and 47 deg so that every tributary reaches every pair.
    c = deg2rad(29*k);
    if n == 1
        J = exp(1i*c);
    else
        a = deg2rad(35 + 12*sin(0.7*k));
        b = deg2rad(17*k);
        J = exp(1i*c) * [cos(a), sin(a)*exp(1i*b); -sin(a)*exp(-1i*b), cos(a)];
    end
    y = xt * J.';

    if isfinite(p.snr)
        sig = mean(abs(y(:)).^2);
        y = y + sqrt(sig/2 * 10^(-p.snr/10)) * (randn(size(y)) + 1i*randn(size(y)));
    end

    % Receiver: skew, gain and quadrature error per I/Q pair.
    z = zeros(size(y));
    for q = 1:n
        g = 10^(p.gRx(q)/20);
        ph = deg2rad(p.qRx(q));
        I = real(y(:, q));
        Q = imag(y(:, q))*cos(ph) + I*sin(ph);
        z(:, q) = I + 1i*g*fracdel(Q, p.tauRx(q)*p.Br*ovr);
    end

    Z{k} = z(ovr:ovr:end, :);
end

truth = p;
end

% ------------------------------------------------------------------
function y = upfir(s, ovr)
% Sinc interpolation: zero-insert, then brick-wall at the symbol-rate Nyquist
% frequency, so decimating at ovr:ovr:end returns the symbols exactly.
N = numel(s) * ovr;
u = zeros(N, 1);
u(ovr:ovr:end) = s;
U = fft(u);
f = (0:N-1).' / N; f(f > 0.5) = f(f > 0.5) - 1;
U(abs(f) > 0.5/ovr) = 0;
y = ifft(U) * ovr;
end

function y = fracdel(x, d)
% Exact fractional delay of d samples on a real signal.
N = numel(x);
f = (0:N-1).' / N; f(f > 0.5) = f(f > 0.5) - 1;
y = real(ifft(fft(x) .* exp(-1i*2*pi*f*d)));
end
