function S = iq_skew_split(Z, T, Fs, varargin)
%IQ_SKEW_SPLIT  Separate transmitter and receiver I/Q imbalance, per frequency,
%from ordinary data captures with known symbols.
%
%   S = iq_skew_split(Z, T, Fs)
%   S = iq_skew_split(Z, T, Fs, 'fmax', 28e9, 'plot', true)
%
% INPUTS
%   Z    1 x K cell, one entry per capture. Z{k} is Ns x nPair complex: the
%        received I + iQ of each receive pair (one pair per polarisation), one
%        sample per symbol, BEFORE any I/Q imbalance correction.
%   T    Ns x nTrib complex: the known transmitted symbols of each tributary,
%        the same for every capture. Pass a 1 x K cell if they differ.
%   Fs   symbol rate in Hz (the sample rate of Z).
%
% OPTIONS (name, value)
%   'fmax'         band used for the fit, |f| <= fmax.        default 0.44*Fs
%   'Nfft'         block length of the Welch channel estimate.  default 256
%   'align'        find each capture's integer delay to T by circular
%                  cross-correlation (looping test pattern).   default false
%   'minElemFrac'  drop a Jones element carrying less than this fraction of
%                  its receive pair's power.                   default 0.02
%   'minRows'      fewest captures a tributary needs in a pair to be solved
%                  there.                                       default 3
%   'tribName'     cell of names for the tributaries.
%   'verbose'      print the report.                            default true
%   'plot'         plot B(f) and the line fits.                 default false
%
% OUTPUT
%   S.f            frequencies of the fitted bins (Hz)
%   S.pair{q}      per receive pair q:
%       .rx        receiver side: gain_dB, quad_deg, skew_ps, mag, resid, echo
%       .tx{r}     transmitter side for tributary r (NaN where unidentifiable)
%       .null{r}   the null term, which nothing physical drives: the noise
%                  floor of the split
%       .Rx, .Tx{r}, .Null{r}   the per-bin curves the lines are fitted to
%       .lever     |mean(lambda)|, 0 is ideal, 1 means no phase spread
%       .cond      condition number of the regression
%   S.txSkew_ps    nTrib x nPair: each transmitter's skew seen through each pair
%
% SIGN CONVENTIONS
%   gain_dB  = 20*log10(|Q| / |I|)
%   quad_deg = phi, where the Q branch reads Q*cos(phi) + I*sin(phi)
%   skew_ps  > 0 means Q lags I. Undo it by delaying I by skew_ps.
%
% METHOD
%   1. Per capture, per receive pair q and tributary r, a widely linear channel
%      estimate Z(f) = A(f) T(f) + B(f) conj(T(-f)) gives the signal path A and
%      the image path B. Their ratio
%
%         rho(f) = B(f) / conj(A(-f)) = B_rx(f) + B_tx(f) * lambda
%
%      cancels anything common to I and Q (delay, loss, real filters). The
%      receiver term is fixed to the ADC axes. The transmitter term is seen
%      through the optical phase theta, once directly and once conjugated, so
%      it carries lambda = J / conj(J), which is exp(2i*theta) for one
%      polarisation and the phase of each Jones-matrix element J(q,r) for two.
%   2. theta drifts randomly between captures, so regressing rho against
%      [1, lambda, conj(lambda)] in every frequency bin returns the receiver
%      and transmitter images separately. The conj(lambda) term has no
%      physical source; it is kept as the noise floor of the split.
%   3. Each side's B(f) is fitted with a straight line in f. To first order a
%      gain imbalance is real and flat, a quadrature error imaginary and flat,
%      and a skew imaginary and proportional to f (i*tan(pi*f*tau)), so the
%      real intercept, imaginary intercept and imaginary slope give all three.
%      A band average would read the skew as zero because it is odd in f.
%   4. What the line leaves is ripple; its inverse transform gives the delay
%      of the reflection causing it.
%
% REQUIREMENTS AND LIMITS
%   - theta (or the polarisation state) must change between captures. The
%     printed lever says whether it did; above 0.5 the split is unreliable.
%   - theta must be static within a capture: the same laser for transmitter
%     and LO, or a frequency offset small enough to ignore over one capture.
%   - The link between the two I/Q pairs must have a real impulse response
%     apart from the carrier phase. Uncompensated chromatic dispersion breaks
%     the cancellation of the Tx term.
%   - The Tx quadrature-error reading is biased by the receiver's own phase;
%     the Tx skew and Tx gain readings are not.

opt = struct('fmax', 0.44*Fs, 'Nfft', 256, 'align', false, ...
             'minElemFrac', 0.02, 'minRows', 3, 'tribName', {{}}, ...
             'verbose', true, 'plot', false);
for i = 1:2:numel(varargin)
    if ~isfield(opt, varargin{i})
        error('iq_skew_split:option', 'Unknown option ''%s''.', varargin{i});
    end
    opt.(varargin{i}) = varargin{i+1};
end

if ~iscell(Z), Z = {Z}; end
K = numel(Z);
if ~iscell(T), T = repmat({T}, 1, K); end
if numel(T) ~= K
    error('iq_skew_split:input', 'T must be one matrix or one cell per capture.');
end
nPair = size(Z{1}, 2);
nTrib = size(T{1}, 2);
if isempty(opt.tribName)
    opt.tribName = arrayfun(@(r) sprintf('trib %d', r), 1:nTrib, 'UniformOutput', false);
end

Nfft = opt.Nfft;
kk = (0:Nfft-1).';
fk = kk / Nfft * Fs;
fk(fk >= Fs/2) = fk(fk >= Fs/2) - Fs;
mir = mod(Nfft - kk, Nfft) + 1;              % index of the -f bin
band = abs(fk) <= opt.fmax & fk ~= 0;
fb = fk(band);
nB = sum(band);

%% rho(f) and its lever for every capture, pair and tributary
rho = cell(1, nPair); lam = cell(1, nPair); trb = cell(1, nPair);
nSkip = 0;
for k = 1:K
    t = T{k};
    z = Z{k};
    n = min(size(z, 1), size(t, 1));
    z = z(1:n, :);
    t = t(1:n, :);
    if opt.align
        z = align_to(z, t);
    end
    J = (t \ z).';                           % nPair x nTrib
    rowP = sum(abs(J).^2, 2);
    for q = 1:nPair
        [Af, Bf] = wl_channel(z(:, q), t, Nfft, mir);
        for r = 1:nTrib
            % An element the tributary barely reaches gives rho = noise/noise.
            if abs(J(q, r))^2 < opt.minElemFrac * rowP(q)
                nSkip = nSkip + 1;
                continue;
            end
            rho{q} = [rho{q}, Bf(band, r) ./ conj(Af(mir(band), r))];
            lam{q} = [lam{q}; J(q, r) / conj(J(q, r))];
            trb{q} = [trb{q}; r];
        end
    end
end

%% One regression per frequency bin, per receive pair
% Rx has one column (it belongs to the pair); Tx and Null have one per
% tributary (they belong to the transmitter).
S = struct('K', K, 'Fs', Fs, 'fmax', opt.fmax, 'f', fb, ...
           'nPair', nPair, 'nTrib', nTrib, 'nSkipped', nSkip);
S.pair = cell(1, nPair);
for q = 1:nPair
    L = lam{q}(:);
    t = trb{q}(:);
    nRow = numel(L);
    nRowTrib = arrayfun(@(r) sum(t == r), 1:nTrib);
    okTrib = nRowTrib >= opt.minRows;

    Atx = zeros(nRow, nTrib);
    Anl = zeros(nRow, nTrib);
    for r = 1:nTrib
        s = (t == r);
        Atx(s, r) = L(s);
        Anl(s, r) = conj(L(s));
    end
    use = [true, okTrib, okTrib];
    A = [ones(nRow, 1), Atx, Anl];
    A = A(:, use);

    Cfull = nan(1 + 2*nTrib, nB);
    if nRow > size(A, 2) && rank(A) == size(A, 2)
        C = A \ rho{q}.';
        res = rho{q}.' - A*C;
        Cfull(use, :) = C;
    else
        res = nan(1, nB);
    end

    P = struct();
    P.nRowTrib = nRowTrib;
    P.okTrib = okTrib;
    P.lever = abs(mean(L));
    P.cond = cond(A);
    P.Rx = Cfull(1, :).';
    P.rx = fit_side(fb, P.Rx);
    P.rx.echo = echo_profile(P.rx.resid, fb, band, Nfft, Fs);
    P.Tx = cell(1, nTrib); P.Null = cell(1, nTrib);
    P.tx = cell(1, nTrib); P.null = cell(1, nTrib);
    for r = 1:nTrib
        P.Tx{r} = Cfull(1 + r, :).';
        P.Null{r} = Cfull(1 + nTrib + r, :).';
        P.tx{r} = fit_side(fb, P.Tx{r});
        P.tx{r}.echo = echo_profile(P.tx{r}.resid, fb, band, Nfft, Fs);
        P.null{r} = fit_side(fb, P.Null{r});
    end
    P.residRMS = sqrt(mean(abs(res(:)).^2));
    S.pair{q} = P;
end

S.txSkew_ps = nan(nTrib, nPair);
for q = 1:nPair
    for r = 1:nTrib
        if S.pair{q}.okTrib(r)
            S.txSkew_ps(r, q) = S.pair{q}.tx{r}.skew_ps;
        end
    end
end

if opt.verbose, report(S, opt); end
if opt.plot, plot_split(S, opt); end
end

% ------------------------------------------------------------------
function [Af, Bf] = wl_channel(z, t, Nfft, mir)
% Widely linear channel of one receive pair, block by block:
%   Z(f) = sum_r A_r(f) T_r(f) + B_r(f) conj(T_r(-f)),
% solved by least squares over blocks in each bin, with every tributary and
% its mirror in the same solve. Separate cross-spectra would leave the
% finite-sample correlation between regressors in the result, and since the
% same symbols are sent in every capture that leakage would not average out
% across captures: it shows up as a fixed Tx gain error.
nTrib = size(t, 2);
nBlk = floor(numel(z) / Nfft);
idx = (1:Nfft).' + (0:nBlk-1)*Nfft;
Zb = fft(z(idx));                            % Nfft x nBlk
Tb = zeros(Nfft, nBlk, nTrib);
for r = 1:nTrib
    tr = t(:, r);
    Tb(:, :, r) = fft(tr(idx));
end
X = cat(3, Tb, conj(Tb(mir, :, :)));         % Nfft x nBlk x 2*nTrib
Af = zeros(Nfft, nTrib); Bf = zeros(Nfft, nTrib);
for k = 1:Nfft
    Xk = reshape(X(k, :, :), nBlk, 2*nTrib);
    c = Xk \ Zb(k, :).';
    Af(k, :) = c(1:nTrib);
    Bf(k, :) = c(nTrib+1:end);
end
end

% ------------------------------------------------------------------
function s = fit_side(f, Bf)
% Straight line through B(f): real intercept -> gain imbalance, imaginary
% intercept -> quadrature error, imaginary slope -> skew. The slope is then
% re-solved against the exact i*tan(pi*f*tau), which a tangent line overstates.
% An unidentifiable column comes in as NaN and goes out as NaN.
if isempty(Bf) || any(~isfinite(Bf))
    s = struct('coef', [NaN; NaN], 'mag', NaN, 'gain_dB', NaN, 'quad_deg', NaN, ...
               'skew_ps', NaN, 'skewLin_ps', NaN, 'resid', nan(numel(f), 1), ...
               'residRMS', NaN);
    return
end
D = [ones(numel(f), 1), f];
c = D \ Bf;
r = Bf - D*c;

eps_g = real(c(1));
q0 = imag(c(1));
slope = imag(c(2));

ff = f(:);
lsSlope = @(tau) (ff' * tan(pi*ff*tau)) / (ff' * ff);
tau0 = slope / pi;
tau = tau0;
if abs(tau0) > 1e-15 && abs(pi*max(abs(ff))*tau0) < 1.2
    try
        tau = fzero(@(x) lsSlope(x) - slope, tau0);
    catch
        tau = tau0;
    end
end

s = struct();
s.coef = c;
s.mag = sqrt(mean(abs(Bf).^2));
s.gain_dB = 20*log10((1 - eps_g) / (1 + eps_g));
s.quad_deg = rad2deg(2 * q0);
s.skew_ps = tau * 1e12;
s.skewLin_ps = tau0 * 1e12;
s.resid = r;
s.residRMS = sqrt(mean(abs(r).^2));
end

% ------------------------------------------------------------------
function e = echo_profile(r, f, band, Nfft, Fs)
% A reflection with round-trip delay D makes B(f) ripple with period 1/D, so
% the inverse transform of the line-fit residual peaks at D. Delays inside one
% band-limited main lobe are the line fit's own residual and are excluded.
e = struct('t', [], 'h', [], 'peak', NaN, 'peak_ps', NaN, 'tMin_ps', NaN);
if any(~isfinite(r)), return; end
n = numel(r);
[~, ord] = sort(f);                          % bins are in FFT order
w = zeros(n, 1);
w(ord) = 0.5 - 0.5*cos(2*pi*(0:n-1).' / (n - 1));
full = zeros(Nfft, 1);
full(band) = r(:) .* w;
h = fftshift(ifft(full));
t = (-Nfft/2 : Nfft/2-1).' / Fs;
tMin = 2 / (Fs * sum(band) / Nfft);
[pk, i] = max(abs(h) .* (abs(t) > tMin));
e = struct('t', t, 'h', h, 'peak', pk, 'peak_ps', t(i)*1e12, 'tMin_ps', tMin*1e12);
end

% ------------------------------------------------------------------
function z = align_to(z, t)
% Integer delay of each receive pair to the tributary it correlates best with,
% by circular cross-correlation. Assumes a looping pattern, as from an AWG.
Tf = fft(t);
for q = 1:size(z, 2)
    Zf = fft(z(:, q));
    best = 0; lag = 0;
    for r = 1:size(t, 2)
        x = ifft(Zf .* conj(Tf(:, r)));
        [m, i] = max(abs(x));
        if m > best, best = m; lag = i - 1; end
    end
    z(:, q) = circshift(z(:, q), -lag);
end
end

% ------------------------------------------------------------------
function report(S, opt)
if S.nPair == 2, pn = {'X', 'Y'}; else, pn = arrayfun(@num2str, 1:S.nPair, 'UniformOutput', false); end
fprintf('\n--- iq_skew_split: %d captures, %d bins, |f| <= %.1f GHz ---\n', ...
    S.K, numel(S.f), S.fmax/1e9);
if S.nSkipped > 0
    fprintf('dropped %d weak Jones element(s) below %.1f%% of their pair''s power\n', ...
        S.nSkipped, 100*opt.minElemFrac);
end
for q = 1:S.nPair
    P = S.pair{q};
    fprintf('\nreceive pair %s   lever %.3f, cond %.1f', pn{q}, P.lever, P.cond);
    if P.lever > 0.5 || P.cond > 20
        fprintf('   POOR: too little phase spread, split unreliable\n');
    else
        fprintf('\n');
    end
    fprintf('%-18s %9s %11s %12s %11s\n', 'side', '|B| rms', 'gain dB', 'quad deg', 'skew ps');
    fprintf('%s\n', repmat('-', 1, 65));
    print_side('Rx', P.rx);
    for r = 1:S.nTrib
        if P.okTrib(r)
            print_side(['Tx ' opt.tribName{r}], P.tx{r});
        else
            fprintf('%-18s %9s %11s %12s %11s   (%d rows, not identifiable)\n', ...
                ['Tx ' opt.tribName{r}], '--', '--', '--', '--', P.nRowTrib(r));
        end
    end
    for r = 1:S.nTrib
        if P.okTrib(r), print_side(['null ' opt.tribName{r}], P.null{r}); end
    end
    fprintf('%-18s %9.4f   per-bin scatter across captures\n', 'unexplained', P.residRMS);
    fprintf('%-18s %9.4f   ripple the Rx line does not explain (see rx.echo)\n', ...
        'Rx ripple', P.rx.residRMS);
end
if S.nPair > 1
    fprintf('\nTx skew (ps) of each transmitter, seen through each receive pair:\n');
    for r = 1:S.nTrib
        fprintf('  %-10s', opt.tribName{r});
        fprintf(' %9.3f', S.txSkew_ps(r, :));
        fprintf('\n');
    end
end
fprintf(['\nThe null row has no physical source: read the Tx row only where it\n' ...
         'stands clear of it. Skew > 0 means Q lags I.\n']);
end

function print_side(name, s)
fprintf('%-18s %9.4f %11.3f %12.3f %11.3f\n', name, s.mag, s.gain_dB, s.quad_deg, s.skew_ps);
end

% ------------------------------------------------------------------
function plot_split(S, opt)
figure('Name', 'iq_skew_split: image vs frequency');
nCol = 1 + S.nTrib;
for q = 1:S.nPair
    P = S.pair{q};
    subplot(S.nPair, nCol, (q-1)*nCol + 1);
    plot_side(sprintf('pair %d, Rx', q), S.f, P.Rx, P.rx);
    for r = 1:S.nTrib
        subplot(S.nPair, nCol, (q-1)*nCol + 1 + r);
        if P.okTrib(r)
            plot_side(sprintf('pair %d, Tx %s', q, opt.tribName{r}), S.f, P.Tx{r}, P.tx{r});
        else
            axis off;
            title(sprintf('pair %d, Tx %s: not identifiable', q, opt.tribName{r}));
        end
    end
end
end

function plot_side(name, f, Bf, s)
[f, i] = sort(f);
Bf = Bf(i);
fit = s.coef(1) + s.coef(2)*f;
plot(f/1e9, real(Bf), '.', 'Color', [0 0.35 1], 'MarkerSize', 5); hold on;
plot(f/1e9, imag(Bf), '.', 'Color', [0.96 0.5 0.19], 'MarkerSize', 5);
plot(f/1e9, real(fit), '-', 'Color', [0 0.35 1], 'LineWidth', 1.5);
plot(f/1e9, imag(fit), '-', 'Color', [0.96 0.5 0.19], 'LineWidth', 1.5);
grid on; xlabel('f (GHz)'); ylabel('B(f)');
legend('real (gain)', 'imag (quadrature + skew)', 'Location', 'best');
title(sprintf('%s: %+.2f dB, %+.2f deg, %+.2f ps', name, s.gain_dB, s.quad_deg, s.skew_ps));
end
