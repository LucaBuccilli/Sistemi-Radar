%% Compressione d'impulso (chirp + filtro adattato) e mappa range-Doppler
% Schema: filtro adattato in tempo veloce (distanza), poi FFT in tempo lento (Doppler).
clear; close all; clc; rng(3);
c = 3e8; fc = 1.5e9; B = 3e6; tau = 10e-6; fs = 12e6;     % prodotto tempo-banda B*tau = 30
PRT = 100e-6; N = 32; lam = c/fc;

t = (0:round(tau*fs)-1)/fs;
chirp = exp(1j*pi*(B/tau)*(t - tau/2).^2);                % LFM in banda base
mf = conj(fliplr(chirp));                                 % filtro adattato h(t) = s*(-t)

% 1) Compressione del singolo impulso
comp = abs(conv(chirp, mf)); comp = comp/max(comp);
fprintf('Durata impulso: %.1f us -> dopo la compressione (-3 dB): %.2f us (teoria 0.886/B = %.2f us)\n', ...
    tau*1e6, sum(comp > 1/sqrt(2))/fs*1e6, 0.886e6/B);
fprintf('Risoluzione in distanza: %.0f m senza compressione, %.0f m con compressione\n', c*tau/2, c/(2*B));
figure; plot(((0:numel(comp)-1) - (numel(chirp)-1))/fs*1e6, 20*log10(max(comp, 1e-4))); grid on;
ylim([-50 1]); xlabel('Ritardo [us]'); ylabel('dB');
title('Uscita del filtro adattato al chirp (lobi laterali a circa -13 dB)');

% 2) Treno di N impulsi: due bersagli + rumore
ns = round(PRT*fs);
targets = [6000 150 1.0; 9000 -400 0.5];                  % distanza [m], velocità [m/s], ampiezza
data = (randn(N, ns) + 1j*randn(N, ns))/sqrt(2);
for it = 1:size(targets, 1)
    d = round(2*targets(it,1)/c*fs);
    fd = 2*targets(it,2)/lam;
    for n = 1:N
        idx = d + (1:numel(chirp));
        data(n, idx) = data(n, idx) + targets(it,3)*chirp*exp(1j*2*pi*fd*(n-1)*PRT);
    end
end

rc = zeros(N, ns);
for n = 1:N
    y = conv(data(n,:), mf);
    rc(n,:) = y(numel(mf):end);                           % compressione in distanza
end
win = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));                 % finestra di Hann (senza toolbox)
rd = fftshift(fft(rc .* win, [], 1), 1);                  % FFT Doppler lungo il tempo lento

rng_axis = (0:ns-1)*c/(2*fs);
vel_axis = (-N/2:N/2-1)/(N*PRT)*lam/2;
P = 20*log10(abs(rd));
[~, imax] = max(P(:)); [i, j] = ind2sub(size(P), imax);
fprintf('Picco della mappa: R = %.0f m, v = %.0f m/s (bersaglio vero: 6000 m, 150 m/s; passo Doppler %.0f m/s)\n', ...
    rng_axis(j), vel_axis(i), lam/(2*N*PRT));
fprintf('Velocità non ambigua: +/- %.0f m/s\n', lam/(4*PRT));

figure; imagesc(rng_axis/1e3, vel_axis, P); axis xy; colorbar; caxis(max(P(:)) + [-50 0]);
xlabel('Distanza [km]'); ylabel('Velocità radiale [m/s]');
title('Mappa range-Doppler (filtro adattato + FFT in tempo lento)');
