%DEMO_SINGLE_POL  Split Tx and Rx I/Q imbalance on a single-polarisation link.
%
% 20 captures of 64 GBd 16-QAM. The transmitter has 1.4 ps of skew and 0.22 dB
% of gain imbalance; the receiver has 3.6 ps, 0.35 dB and 1 deg of quadrature
% error. The carrier phase changes between captures, as it does when a fibre
% drifts.

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));

Fs = 64e9;
[Z, T, truth] = simulate_link('nPol', 1, 'nCap', 20, 'Br', Fs, ...
    'tauRx', 3.6e-12, 'gRx', 0.35, 'qRx', 1.0, 'tauTx', 1.4e-12, 'gTx', 0.22);

S = iq_skew_split(Z, T, Fs, 'plot', true);

P = S.pair{1};
fprintf('\n%-10s %10s %10s\n', '', 'injected', 'measured');
fprintf('%-10s %10.2f %10.2f  ps\n', 'Rx skew', truth.tauRx*1e12, P.rx.skew_ps);
fprintf('%-10s %10.2f %10.2f  dB\n', 'Rx gain', truth.gRx, P.rx.gain_dB);
fprintf('%-10s %10.2f %10.2f  deg\n', 'Rx quad', truth.qRx, P.rx.quad_deg);
fprintf('%-10s %10.2f %10.2f  ps\n', 'Tx skew', truth.tauTx*1e12, P.tx{1}.skew_ps);
fprintf('%-10s %10.2f %10.2f  dB\n', 'Tx gain', truth.gTx, P.tx{1}.gain_dB);
