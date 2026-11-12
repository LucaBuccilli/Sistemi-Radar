%% Canceller MTI binomiali: attenuazione del clutter e Improvement Factor
% Clutter con spettro gaussiano: rho(kT) = exp(-2 pi^2 sigma_f^2 k^2 T^2).
% Canceller di ordine M: H(f) = (1 - exp(-j 2 pi f T))^M, |H|^2 = (4 sin^2(pi f T))^M.
% Scenario: radar di ricerca a corto raggio in banda S, sigma_v = 0.2 m/s, lambda = 0.1 m, PRF = 1250 Hz.
clear; close all; clc; rng(2);
PRF = 1250; lam = 0.1; sigma_v = 0.2; T = 1/PRF;
sigma_f = 2*sigma_v/lam;                         % 4 Hz
rho = @(k, sf) exp(-2*pi^2*sf^2*(k*T).^2);

fprintf('sigma_f = %.1f Hz, rho(T) = %.4f, rho(2T) = %.4f, rho(3T) = %.4f\n', ...
    sigma_f, rho(1,sigma_f), rho(2,sigma_f), rho(3,sigma_f));
fprintf('Canceller | CA [dB] | IF a PRF/2 | IF a PRF/20 | IF medio\n');
fd = linspace(0, PRF, 1000);
figure; hold on; grid on;
for M = 1:3
    Pout = clutter_out(M, sigma_f, rho);
    IF = @(f) (4*sin(pi*f*T).^2).^M / Pout;
    fprintf('  M = %d   | %6.1f  | %8.1f   | %8.1f    | %6.1f\n', M, -10*log10(Pout), ...
        10*log10(IF(PRF/2)), 10*log10(IF(PRF/20)), 10*log10(nchoosek(2*M, M)/Pout));
    plot(fd, 10*log10(max(IF(fd), 1e-30)), 'LineWidth', 1.2, 'DisplayName', sprintf('ordine %d', M));
end
ylim([-20 120]); xlabel('Doppler del bersaglio [Hz]'); ylabel('Improvement Factor [dB]');
title('IF dei canceller binomiali: velocità cieche a 0 e PRF'); legend('Location', 'south');

% Effetto della rotazione dell'antenna: lo spettro del clutter si allarga.
% Nota: arrotondando rho(T) a 0.999 si otterrebbero circa 87.5 dB; con il valore esatto il risultato è più alto.
ToT = 66.7e-3; sigma_rot = hypot(sigma_f, 0.265/ToT);
fprintf('Con rotazione dell''antenna: sigma_f'' = %.2f Hz -> IF triplo canceller a PRF/2 = %.1f dB\n', ...
    sigma_rot, 10*log10(64/clutter_out(3, sigma_rot, rho)));

% Verifica per simulazione: clutter gaussiano correlato + canceller singolo
n = 2^18;
f = (0:n-1)*PRF/n; f(f >= PRF/2) = f(f >= PRF/2) - PRF;       % asse come fftfreq
shape = exp(-f.^2/(4*sigma_f^2));                             % radice della PSD gaussiana
cl = ifft(shape .* (randn(1,n) + 1j*randn(1,n)));
cl = cl / sqrt(mean(abs(cl).^2));
out = filter([1 -1], 1, cl); out = out(2:end);                % y[n] = x[n] - x[n-1]
fprintf('Simulazione canceller singolo: CA = %.1f dB (teoria 1/(2(1-rho)) = %.1f dB)\n', ...
    -10*log10(mean(abs(out).^2)), -10*log10(clutter_out(1, sigma_f, rho)));

function P = clutter_out(M, sf, rho)
% Potenza di clutter in uscita (ingresso unitario): sum_m (-1)^m C(2M, M+m) rho(mT)
P = 0;
for m = -M:M
    P = P + (-1)^m * nchoosek(2*M, M+m) * rho(m, sf);
end
end
