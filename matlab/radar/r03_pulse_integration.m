%% Integrazione di N impulsi: coerente (FFT) contro non coerente (somma dei moduli quadri)
%  coerente:      | sum_n r_n exp(-j 2 pi fd n T) |   -> guadagno N
%  non coerente:  sum_n |r_n|^2                      -> guadagno minore di N
clear; close all; clc; rng(1);
N = 16; trials = 40000; pfa = 1e-3;
kd = 3;                                         % bin Doppler del bersaglio
steer = exp(1j*2*pi*kd*(0:N-1)/N);              % vettore Doppler (riga)
cnoise = @() (randn(trials, N) + 1j*randn(trials, N))/sqrt(2);

names = {'singolo impulso', 'non coerente (N=16)', 'coerente (N=16)'};
stat = {@(r) abs(r(:,1)), @(r) sum(abs(r).^2, 2), @(r) abs(r*steer')};

% soglie fissate sul solo rumore, per avere la stessa Pfa per tutti
thr = zeros(1,3);
for i = 1:3
    s = sort(stat{i}(cnoise()));
    thr(i) = s(round((1-pfa)*trials));
end

snr_db = -12:1:15;
pd = zeros(3, numel(snr_db));
for k = 1:numel(snr_db)
    r = sqrt(10^(snr_db(k)/10))*steer + cnoise();
    for i = 1:3
        pd(i,k) = mean(stat{i}(r) > thr(i));
    end
end

figure; plot(snr_db, pd, 'LineWidth', 1.2); grid on; yline(0.9, 'k:');
xlabel('SNR per impulso [dB]'); ylabel('Pd'); legend(names, 'Location', 'southeast');
title(sprintf('Integrazione di impulsi (Pfa = %.0e)', pfa));

need = zeros(1,3);
for i = 1:3
    [p, iu] = unique(pd(i,:));                  % interp1 richiede ascisse distinte
    need(i) = interp1(p, snr_db(iu), 0.9);
    fprintf('SNR per impulso per Pd = 0.9, %-20s: %5.1f dB\n', names{i}, need(i));
end
fprintf('Guadagno coerente     : %.1f dB (teorico 10log10(N) = %.1f dB)\n', need(1)-need(3), 10*log10(N));
fprintf('Guadagno non coerente : %.1f dB\n', need(1)-need(2));
