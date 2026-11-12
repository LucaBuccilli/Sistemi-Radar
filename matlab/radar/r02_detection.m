%% Rivelazione a singolo impulso: soglia, Pfa, Pd (Monte Carlo) e formula di Albersheim
% H0: z = n       (inviluppo di Rayleigh)
% H1: z = a + n   (inviluppo di Rice)
% con n complesso gaussiano di potenza 1 e SNR = |a|^2.
clear; close all; clc; rng(0);
cnoise = @(n) (randn(n,1) + 1j*randn(n,1))/sqrt(2);

% 1) Verifica di Pfa = exp(-T^2) con una Pfa misurabile
pfa = 1e-3; T = sqrt(-log(pfa));
fprintf('Pfa teorica = %.1e, Pfa Monte Carlo = %.2e\n', pfa, mean(abs(cnoise(2e6)) > T));

% 2) Curve Pd(SNR) per diverse Pfa
snr_db = 0:0.5:18;
figure; hold on; grid on;
for pfa = [1e-2 1e-4 1e-6]
    T = sqrt(-log(pfa));
    pd = zeros(size(snr_db));
    for i = 1:numel(snr_db)
        pd(i) = mean(abs(sqrt(10^(snr_db(i)/10)) + cnoise(2e5)) > T);
    end
    plot(snr_db, pd, 'LineWidth', 1.2, 'DisplayName', sprintf('Pfa = %.0e', pfa));
    plot(albersheim(pfa, 0.9, 1), 0.9, 'ko', 'HandleVisibility', 'off');
end
xlabel('SNR [dB]'); ylabel('Pd'); legend('Location', 'southeast');
title('Pd in funzione dell''SNR (pallini: Albersheim per Pd = 0.9)');

fprintf('Albersheim: SNR per Pfa=1e-6, Pd=0.9, 1 impulso  = %.1f dB (atteso ~13 dB)\n', albersheim(1e-6, 0.9, 1));
fprintf('Albersheim: SNR per Pfa=1e-6, Pd=0.9, 18 impulsi = %.1f dB per impulso\n', albersheim(1e-6, 0.9, 18));

function s = albersheim(pfa, pd, n)
% SNR richiesto [dB], bersaglio non fluttuante, n impulsi integrati non coerentemente
A = log(0.62/pfa); B = log(pd/(1-pd));
s = -5*log10(n) + (6.2 + 4.54/sqrt(n+0.44))*log10(A + 0.12*A*B + 1.7*B);
end
