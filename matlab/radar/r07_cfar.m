%% CA-CFAR (Cell Averaging CFAR) sulla mappa range-Doppler
% Rivelatore quadratico: senza bersaglio |z|^2 è esponenziale con media Pn (ignota e variabile).
% Il CA-CFAR stima Pn mediando N celle di addestramento (escluse le celle di guardia) e usa
%   T = alpha * Pn_stimata,   alpha = N * (Pfa^(-1/N) - 1)
% 1) verifica della Pfa su rumore esponenziale con potenza variabile
% 2) catena completa: compressione -> (MTI) -> FFT Doppler -> CA-CFAR, con clutter a Doppler zero
clear; close all; clc; rng(9);

% 1) verifica della Pfa
pfa = 1e-4; Nd = 256; Nr = 2048;
level = 10.^(linspace(0, 3, Nr)/10);                      % potenza da 0 a 30 dB lungo la distanza
noise = -log(rand(Nd, Nr)) .* level;                       % esponenziale a media unitaria * livello
[det, ~, Ntr] = ca_cfar_2d(noise, [2 2], [4 8], pfa);
alpha = Ntr*(pfa^(-1/Ntr) - 1);
fixed_thr = -log(pfa)*mean(level(1:Nr/10));
fprintf('Celle di addestramento N = %d, alpha = %.2f (%.1f dB)\n', Ntr, alpha, 10*log10(alpha));
fprintf('Pfa progettata = %.1e, misurata CA-CFAR = %.2e, misurata soglia fissa = %.2e\n', ...
    pfa, mean(det(:)), mean(noise(:) > fixed_thr));

% 2) catena completa
c = 3e8; fc = 1.5e9; B = 3e6; tau = 10e-6; fs = 12e6; PRT = 100e-6; Np = 32; lam = c/fc;
t = (0:round(tau*fs)-1)/fs;
chirp = exp(1j*pi*(B/tau)*(t - tau/2).^2);
mf = conj(fliplr(chirp));
ns = round(PRT*fs);
rng_axis = (0:ns-1)*c/(2*fs);
vel_axis = (-Np/2:Np/2-1)/(Np*PRT)*lam/2;

data = (randn(Np, ns) + 1j*randn(Np, ns))/sqrt(2);
amp_cl = 30*min(1000./max(rng_axis, 1), 1).^1.5;           % clutter fermo, forte a breve distanza
cl_tx = amp_cl.*(randn(1, ns) + 1j*randn(1, ns))/sqrt(2);
cl_rx = conv(cl_tx, chirp);
data = data + repmat(cl_rx(1:ns), Np, 1);
targets = [2500 120 0.15; 6000 150 0.5; 9000 -400 0.1];    % R [m], v [m/s], ampiezza
for it = 1:3
    d = round(2*targets(it,1)/c*fs); fd = 2*targets(it,2)/lam;
    idx = d + (1:numel(chirp));
    for n = 1:Np
        data(n, idx) = data(n, idx) + targets(it,3)*chirp*exp(1j*2*pi*fd*(n-1)*PRT);
    end
end
rc = zeros(Np, ns);
for n = 1:Np
    y = conv(data(n,:), mf); rc(n,:) = y(numel(mf):numel(mf)+ns-1);
end

P = doppler_map(rc, Np);                                   % senza MTI
P_mti = doppler_map(rc(2:end,:) - rc(1:end-1,:), Np);      % canceller singolo prima della FFT

pfa2 = 1e-5;
floor_n = median(reshape(P(:, rng_axis > 12e3), [], 1));
det_fixed = P > -log(pfa2)*floor_n;
det_2d  = ca_cfar_2d(P, [2 4], [4 12], pfa2);              % finestra anche in Doppler
det_rng = ca_cfar_2d(P, [0 6], [0 16], pfa2);              % solo lungo la distanza
det_mti = ca_cfar_2d(P_mti, [0 6], [0 16], pfa2);          % MTI + CFAR in distanza

names = {'soglia fissa', 'CA-CFAR 2D', 'CA-CFAR in distanza', 'MTI + CA-CFAR'};
dets = {det_fixed, det_2d, det_rng, det_mti};
ridge = repmat(abs(vel_axis(:)) < 40, 1, ns);
for k = 1:4
    D = dets{k}; hits = false(1,3);
    for it = 1:3
        [~, j] = min(abs(rng_axis - targets(it,1))); [~, i] = min(abs(vel_axis - targets(it,2)));
        sub = D(max(i-1,1):min(i+1,Np), max(j-3,1):min(j+3,ns)); hits(it) = any(sub(:));
    end
    fprintf('%-22s: bersagli %d/3, celle sul clutter = %d, celle totali = %d\n', names{k}, sum(hits), nnz(D & ridge), nnz(D));
end

figure;
subplot(2,2,1); imagesc(rng_axis/1e3, vel_axis, 10*log10(P)); axis xy; title('Mappa range-Doppler [dB]');
subplot(2,2,2); imagesc(rng_axis/1e3, vel_axis, det_fixed); axis xy; title('Soglia fissa');
subplot(2,2,3); imagesc(rng_axis/1e3, vel_axis, det_rng); axis xy; title('CA-CFAR in distanza');
subplot(2,2,4); imagesc(rng_axis/1e3, vel_axis, det_mti); axis xy; title('MTI + CA-CFAR');
colormap(gray);

function P = doppler_map(x, Np)
n = size(x, 1);
w = 0.5 - 0.5*cos(2*pi*(0:n-1)'/(n-1));                   % finestra di Hann
P = abs(fftshift(fft(x.*w, Np, 1), 1)).^2;
end

function [det, thr, Ntr] = ca_cfar_2d(P, guard, train, pfa)
% guard/train: semiampiezze [Doppler distanza]. Bordi trattati in modo circolare.
gd = guard(1); gr = guard(2); td = train(1); tr = train(2);
K = ones(2*(gd+td)+1, 2*(gr+tr)+1);
K(td+1:td+2*gd+1, tr+1:tr+2*gr+1) = 0;                     % guardia + cella sotto test
Ntr = nnz(K);
alpha = Ntr*(pfa^(-1/Ntr) - 1);
pd = size(K,1)-1; pr = size(K,2)-1;                        % estensione circolare per conv2 'valid'
Pe = P([end-pd/2+1:end, 1:end, 1:pd/2], [end-pr/2+1:end, 1:end, 1:pr/2]);
est = conv2(Pe, K, 'valid')/Ntr;
thr = alpha*est;
det = P > thr;
end
