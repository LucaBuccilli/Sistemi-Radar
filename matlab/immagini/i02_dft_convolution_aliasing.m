%% DFT e convoluzione (circolare contro lineare, zero padding) e aliasing
% 1) Prodotto delle DFT = convoluzione CIRCOLARE. Con zero padding fino a M1+M2-1 si ottiene quella lineare.
% 2) Sottocampionare senza passa-basso produce aliasing; la decimazione (filtro + sottocampionamento) no.
clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
img = test_image(128);
h = ones(9)/81;                                            % filtro di media 9x9

% 1) convoluzione lineare, tramite DFT con e senza zero padding
lin = conv2(img, h);                                       % dimensione (128+9-1)^2
M = size(img,1) + size(h,1) - 1;
dft_zp = real(ifft2(fft2(img, M, M) .* fft2(h, M, M)));
circ = real(ifft2(fft2(img) .* fft2(h, size(img,1), size(img,2))));
fprintf('Dimensione uscita lineare: %d x %d = N1 + N2 - 1\n', size(lin,1), size(lin,2));
fprintf('Errore massimo DFT con zero padding rispetto a conv2: %.2e\n', max(abs(lin(:) - dft_zp(:))));
d = lin(1:size(img,1), 1:size(img,2)) - circ;
fprintf('Errore massimo DFT senza zero padding (circolare)   : %.2f -> NON è la convoluzione lineare\n', max(abs(d(:))));

% 2) aliasing: zone plate (frequenza crescente dal centro)
n = 512; [x, y] = meshgrid((-n/2:n/2-1)/n);
zone = 0.5 + 0.5*cos(2*pi*120*(x.^2 + y.^2));
Mfac = 4;
down = zone(1:Mfac:end, 1:Mfac:end);                       % solo sottocampionamento
g = exp(-(-8:8).^2/(2*(Mfac/2)^2)); g = g/sum(g);          % passa-basso gaussiano separabile
zf = conv2(g, g, zone, 'same');
dec = zf(1:Mfac:end, 1:Mfac:end);                          % passa-basso, poi sottocampionamento

figure; colormap gray;
subplot(1,3,1); imagesc(zone, [0 1]); axis image off; title('originale 512x512');
subplot(1,3,2); imagesc(down, [0 1]); axis image off; title('sottocampionata x4: aliasing');
subplot(1,3,3); imagesc(dec,  [0 1]); axis image off; title('decimata x4: nessun falso anello');

% 3) 1D: un tono oltre Fs/2 si ripiega
Fs = 8000; f0 = 5000; Nn = 800;
xs = cos(2*pi*f0*(0:Nn-1)/Fs);
X = abs(fft(xs)); [~, k] = max(X(1:Nn/2));
fprintf('Tono a %d Hz campionato a %d Hz: picco dello spettro a %.0f Hz (atteso |f0 - Fs| = %d Hz)\n', ...
    f0, Fs, (k-1)*Fs/Nn, abs(f0-Fs));
