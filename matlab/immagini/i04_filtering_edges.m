%% Filtraggio spaziale 2D: smoothing, sharpening, gradiente (Sobel) e Laplaciano
% Verifica anche la risposta in frequenza della maschera a croce 3x3:
% [0 1 0; 1 1 1; 0 1 0]  ->  H = 1 + 2 cos(w1) + 2 cos(w2).
clear; close all; clc; rng(5);
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
img = test_image(256);
noisy = img + 0.08*randn(size(img));

mean3 = ones(3)/9;
lap = [0 1 0; 1 -4 1; 0 1 0];
sob = [-1 0 1; -2 0 2; -1 0 1];

smooth = conv2(noisy, mean3, 'same');                      % passa-basso: toglie rumore
laplacian = conv2(img, lap, 'same');                       % derivata seconda: zero crossing sui bordi
sharp = img - 0.7*laplacian;                               % originale + passa-alto
grad = hypot(conv2(img, sob, 'same'), conv2(img, sob', 'same'));   % picco sui bordi

% il passa-basso toglie rumore nelle zone lisce ma cancella i dettagli fini (le righe in basso)
r1 = 5:150; r2 = 170:250; cc = 5:250;                      % si escludono i bordi dell'immagine
mse = @(a, rr) mean(mean((a(rr,cc) - img(rr,cc)).^2));
fprintf('MSE rumorosa : zona liscia = %.5f, zona a righe = %.5f\n', mse(noisy, r1), mse(noisy, r2));
fprintf('MSE media 3x3: zona liscia = %.5f, zona a righe = %.5f\n', mse(smooth, r1), mse(smooth, r2));

% risposta in frequenza della maschera a croce
cross = [0 1 0; 1 1 1; 0 1 0];
w = linspace(-pi, pi, 65); [W1, W2] = meshgrid(w, w);
H = zeros(size(W1));
for m = -1:1
    for n = -1:1
        H = H + cross(m+2, n+2)*exp(-1j*(W1*n + W2*m));
    end
end
err = abs(H - (1 + 2*cos(W1) + 2*cos(W2)));
fprintf('Maschera a croce: errore rispetto a 1 + 2cos(w1) + 2cos(w2) = %.1e; H(0,0) = %.0f (passa-basso)\n', ...
    max(err(:)), real(H(33,33)));

ims = {noisy, smooth, sharp, grad, abs(laplacian), abs(H)};
ttl = {'con rumore', 'media 3x3 (smoothing)', 'sharpening', 'gradiente (Sobel)', '|Laplaciano|', '|H| maschera a croce'};
figure; colormap gray;
for i = 1:6
    subplot(2,3,i); imagesc(ims{i}); axis image off; title(ttl{i});
end
