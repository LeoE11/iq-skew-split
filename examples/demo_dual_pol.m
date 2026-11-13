%DEMO_DUAL_POL  Split Tx and Rx I/Q imbalance on a dual-polarisation link.
%
% Two independent I/Q modulators (X and Y tributaries), each with its own skew
% and gain imbalance, and a receiver whose two I/Q pairs also differ. The
% polarisation state changes between captures. Each transmitter is measured
% once through each receive pair, which gives a consistency check on the Tx
% reading that needs no extra equipment.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

Fs = 64e9;
[Z, T, truth] = simulate_link('nPol', 2, 'nCap', 20, 'Br', Fs);

S = iq_skew_split(Z, T, Fs, 'tribName', {'X', 'Y'}, 'plot', true);

fprintf('\n%-16s %10s %10s %10s\n', '', 'injected', 'pair X', 'pair Y');
for q = 1:2
    fprintf('%-16s %10.2f %10.2f\n', sprintf('Rx skew, pair %d', q), ...
        truth.tauRx(q)*1e12, S.pair{q}.rx.skew_ps);
end
for r = 1:2
    fprintf('%-16s %10.2f %10.2f %10.2f\n', sprintf('Tx skew, trib %d', r), ...
        truth.tauTx(r)*1e12, S.txSkew_ps(r, 1), S.txSkew_ps(r, 2));
end
