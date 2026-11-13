function run_tests()
%RUN_TESTS  Ground-truth checks for iq_skew_split.
%
% Known defects are injected on each side and each must come back on the
% right side at the right size. Two things cannot be checked by inspection:
%   1. Which side is which. A lever with the wrong exponent swaps the Tx term
%      with the null term and still fits the data; only known values catch it.
%   2. That a skew is recovered at all. A band-averaged estimate reads an
%      odd-in-f image as zero and would pass a test that only checked gain.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(root, fullfile(root, 'examples'));
nFail = 0;

    function chk(name, got, want, tol, unit)
        ok = abs(got - want) <= tol;
        if ~ok, nFail = nFail + 1; end
        tag = 'PASS'; if ~ok, tag = 'FAIL'; end
        fprintf('  [%s] %-28s got %+8.3f  want %+8.3f  (tol %.3f %s)\n', ...
            tag, name, got, want, tol, unit);
    end

%% single polarisation
fprintf('\n== single polarisation ==\n');
Fs = 64e9;
[Z, T, t] = simulate_link('nPol', 1, 'Br', Fs, ...
    'tauRx', 3.6e-12, 'gRx', 0.35, 'qRx', 1.0, 'tauTx', 1.4e-12, 'gTx', 0.22);
S = iq_skew_split(Z, T, Fs, 'verbose', false);
P = S.pair{1};
chk('Rx skew', P.rx.skew_ps, t.tauRx*1e12, 0.15, 'ps');
chk('Rx gain', P.rx.gain_dB, t.gRx, 0.06, 'dB');
chk('Rx quad', P.rx.quad_deg, t.qRx, 0.30, 'deg');
chk('Tx skew', P.tx{1}.skew_ps, t.tauTx*1e12, 0.25, 'ps');
chk('Tx gain', P.tx{1}.gain_dB, t.gTx, 0.12, 'dB');
chk('null |B|', P.null{1}.mag, 0, 0.02, '');
chk('lever', P.lever, 0, 0.5, '');

%% dual polarisation, each transmitter and receive pair different
fprintf('\n== dual polarisation ==\n');
[Z, T, t] = simulate_link('nPol', 2, 'Br', Fs);
S = iq_skew_split(Z, T, Fs, 'verbose', false);
for q = 1:2
    P = S.pair{q};
    chk(sprintf('pair %d Rx skew', q), P.rx.skew_ps, t.tauRx(q)*1e12, 0.15, 'ps');
    chk(sprintf('pair %d Rx gain', q), P.rx.gain_dB, t.gRx(q), 0.06, 'dB');
    chk(sprintf('pair %d Rx quad', q), P.rx.quad_deg, 0, 0.30, 'deg');
    for r = 1:2
        if ~P.okTrib(r)
            nFail = nFail + 1;
            fprintf('  [FAIL] pair %d trib %d not identifiable\n', q, r);
            continue
        end
        chk(sprintf('pair %d Tx%d skew', q, r), P.tx{r}.skew_ps, t.tauTx(r)*1e12, 0.25, 'ps');
        chk(sprintf('pair %d Tx%d gain', q, r), P.tx{r}.gain_dB, t.gTx(r), 0.12, 'dB');
        chk(sprintf('pair %d null%d |B|', q, r), P.null{r}.mag, 0, 0.02, '');
    end
end
% The same transmitter seen through two receive pairs must agree.
for r = 1:2
    chk(sprintf('Tx%d pair-to-pair spread', r), abs(diff(S.txSkew_ps(r, :))), 0, 0.30, 'ps');
end

%% no phase spread: the lever must say the split is not trustworthy
fprintf('\n== no phase spread ==\n');
[Z, T] = simulate_link('nPol', 1, 'nCap', 1, 'Br', Fs);
Z = repmat(Z, 1, 6);
S = iq_skew_split(Z, T, Fs, 'verbose', false);
chk('lever with fixed phase', S.pair{1}.lever, 1, 1e-6, '');

%% integer misalignment of a looping pattern
fprintf('\n== alignment ==\n');
[Z, T, t] = simulate_link('nPol', 1, 'nCap', 8, 'Br', Fs, ...
    'tauRx', 2.5e-12, 'gRx', 0, 'tauTx', 0, 'gTx', 0);
for k = 1:numel(Z), Z{k} = circshift(Z{k}, 137*k); end
S = iq_skew_split(Z, T, Fs, 'align', true, 'verbose', false);
chk('Rx skew after align', S.pair{1}.rx.skew_ps, t.tauRx*1e12, 0.15, 'ps');

fprintf('\n');
if nFail > 0
    error('run_tests:fail', '%d check(s) failed', nFail);
end
fprintf('all checks passed\n');
end
