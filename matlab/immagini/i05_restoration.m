%% Image recovery: y = h * x + n. Inverso, pseudo-inverso, minimi quadrati vincolati e Wiener
% Tutti i filtri sono applicati in frequenza (la sfocatura è una convoluzione circolare):
%   inverso          F = 1 / H
%   pseudo-inverso   F = 1 / H dove |H| > soglia, 0 altrove
%   CLS              F = H* / (|H|^2 + lambda |C|^2),  C = Laplaciano (passa-alto)
%   Wiener           F = H* / (|H|^2 + NSR),           NSR = Pww / Pxx
% Con Image Processing Toolbox le stesse cose si fanno con deconvwnr e deconvreg.
clear; close all; clc; rng(6);
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
x = test_image(256); N = size(x, 1);
psnr_db = @(ref, a) 10*log10(1/mean((ref(:) - a(:)).^2));

% PSF: sfocatura da movimento lungo la diagonale (15 pixel), centrata nell'origine
psf = zeros(N); L = 15;
for i = 0:L-1
    k = mod(i - floor(L/2), N) + 1;
    psf(k, k) = 1/L;
end
H = fft2(psf);

sigma_n = 0.01;
y = real(ifft2(fft2(x).*H)) + sigma_n*randn(N);
Y = fft2(y);

lap = zeros(N); lap(1,1) = 4; lap(1,2) = -1; lap(2,1) = -1; lap(1,N) = -1; lap(N,1) = -1;
C = fft2(lap);
Pxx = abs(fft2(x)).^2/N^2;                % spettro di potenza del segnale (qui noto: caso ideale)
Pww = sigma_n^2;

Hs = H; Hs(abs(H) < 1e-12) = 1e-12;       % evita la divisione per zero nel filtro inverso
names = {'inverso', 'pseudo-inverso', 'CLS (lambda = 0.01)', 'Wiener (NSR costante)', 'Wiener (NSR vero)'};
Fpinv = zeros(N); Fpinv(abs(H) > 0.1) = 1./Hs(abs(H) > 0.1);
filt = {1./Hs, Fpinv, conj(H)./(abs(H).^2 + 0.01*abs(C).^2), ...
        conj(H)./(abs(H).^2 + sigma_n^2/var(x(:))), conj(H)./(abs(H).^2 + Pww./max(Pxx, 1e-12))};

fprintf('PSNR immagine osservata      : %6.2f dB\n', psnr_db(x, y));
res = cell(1, 5);
for i = 1:5
    res{i} = real(ifft2(Y.*filt{i}));
    fprintf('PSNR %-24s: %6.2f dB\n', names{i}, psnr_db(x, res{i}));
end

ims = {x, y, res{1}, res{2}, res{3}, res{5}};
ttl = {'originale', 'sfocata + rumore', names{1}, names{2}, names{3}, names{5}};
figure; colormap gray;
for i = 1:6
    subplot(2,3,i); imagesc(ims{i}, [0 1]); axis image off; title(ttl{i});
end
