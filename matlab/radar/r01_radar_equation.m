%% Equazione radar: SNR in funzione della distanza, portata massima, orizzonte radar
% Caso di riferimento: radar di sorveglianza in banda S (lambda = 0.2 m) e dimensionamento della potenza.
clear; close all; clc;
k = 1.380649e-23; T0 = 290; c = 3e8;

% Caso di riferimento: Pt=32 kW, G=35 dB, lambda=0.2 m, RCS=10 m^2, L=10 dB, F=5 dB, B=1 MHz
Pt = 32e3; G = 10^3.5; lam = 0.2; sigma = 10; L = 10; F = 10^0.5; B = 1e6;
snr = @(R) Pt*G^2*lam^2*sigma ./ ((4*pi)^3 * R.^4 * L*k*T0*B*F);

for R = [100 1e3 10e3 100e3 400e3]
    fprintf('R = %7.1f km   SNR = %6.1f dB\n', R/1e3, 10*log10(snr(R)));
end

snr_min = 10^(13.3/10);                 % Pd = 0.9, Pfa = 1e-6, bersaglio non fluttuante
Rmax = (Pt*G^2*lam^2*sigma / ((4*pi)^3*L*k*T0*B*F*snr_min))^(1/4);
fprintf('Portata massima con SNR_min = 13.3 dB: %.1f km\n', Rmax/1e3);

R = logspace(2, log10(400e3), 400);
figure; semilogx(R/1e3, 10*log10(snr(R)), 'LineWidth', 1.2); hold on; grid on;
yline(13.3, 'r--'); xline(Rmax/1e3, 'k:');
xlabel('Distanza [km]'); ylabel('SNR [dB]');
title('Equazione radar: l''SNR cala di 40 dB per decade (1/R^4)');
legend('SNR singolo impulso', 'SNR richiesto (13.3 dB)');

% Dimensionamento: potenza di picco per Rmax = 150 km
f = 2e9; Bx = 1e6; eta = 0.6; Ag = 12; Fx = 10^0.4; sig = 5; Rx = 150e3; snr_req = 10^1.2;
lam_x = c/f;
Gx = 4*pi*eta*Ag/lam_x^2;
Pt_req = snr_req*(4*pi)^3*Rx^4*k*T0*Fx*Bx / (Gx^2*lam_x^2*sig);
fprintf('Dimensionamento: G = %.2f dB, Pt richiesta = %.1f kW (calcolo a mano in dB: 87.7 kW)\n', 10*log10(Gx), Pt_req/1e3);

% Orizzonte radar con raggio terrestre equivalente 4/3
Re = 4/3*6371e3; hT = [10e3 900 100 0];
for hA = [100 4]
    Rh = sqrt(2*Re*hA) + sqrt(2*Re*hT);
    fprintf('h_antenna = %5.0f m -> orizzonte [km] per h_T = 10 km, 900 m, 100 m, 0: %s\n', ...
        hA, num2str(round(Rh/1e3)));
end
