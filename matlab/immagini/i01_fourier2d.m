%% Trasformata di Fourier 2D: spettro, passa-basso e passa-alto ideali, rotazione
% Proprietà verificata: ruotare l'immagine ruota lo spettro dello stesso angolo.
clear; close all; clc;
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));
img = test_image(256); N = size(img, 1);
spec = @(a) log(1 + abs(fftshift(fft2(a))));          % modulo dello spettro in scala log

F = fftshift(fft2(img));
u = -N/2:N/2-1; [U, V] = meshgrid(u, u); radius = hypot(U, V);
lp = real(ifft2(ifftshift(F .* (radius < 20))));       % passa-basso ideale
hp = real(ifft2(ifftshift(F .* (radius >= 20))));      % passa-alto ideale

% rotazione di 30 gradi con interp2 (senza Image Processing Toolbox)
th = 30*pi/180; [X, Y] = meshgrid(1:N, 1:N); cx = (N+1)/2;
Xr =  cos(th)*(X-cx) + sin(th)*(Y-cx) + cx;
Yr = -sin(th)*(X-cx) + cos(th)*(Y-cx) + cx;
rot = interp2(X, Y, img, min(max(Xr,1),N), min(max(Yr,1),N), 'linear');

% le righe oblique danno due picchi in (+-24, +-24) cicli per immagine
Fa = abs(F); Fa(radius < 5) = 0;
[~, imax] = max(Fa(:)); [iy, ix] = ind2sub(size(Fa), imax);
fprintf('Picco dello spettro (esclusa la continua): (%d, %d) cicli/immagine; atteso (+-24, +-24)\n', u(ix), u(iy));
fprintf('Energia nelle basse frequenze (raggio < 20): %.1f %%\n', ...
    100*sum(abs(F(radius < 20)).^2)/sum(abs(F(:)).^2));

ims = {img, lp, hp, rot};
ttl = {'originale', 'passa-basso (sfocata, ringing)', 'passa-alto (solo bordi)', 'ruotata di 30 gradi'};
figure; colormap gray;
for i = 1:4
    subplot(2,4,i);   imagesc(ims{i});       axis image off; title(ttl{i});
    subplot(2,4,4+i); imagesc(spec(ims{i})); axis image off; title('spettro (log)');
end
