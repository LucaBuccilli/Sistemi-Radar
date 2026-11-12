%% KLT, DCT e schema JPEG (DCT a blocchi 8x8 + quantizzazione con la matrice di luminanza)
% 1) KLT: autovettori della covarianza, coefficienti scorrelati (esempio 2x2 risolvibile a mano).
% 2) DCT: base ortonormale, approssima la KLT dei processi di Markov molto correlati.
% 3) JPEG: blocchi 8x8 -> DCT -> round(X / Q) -> la maggior parte dei coefficienti va a zero.
clear; close all; clc; rng(7);
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));

% 1) KLT, esempio 2x2: Kx = [2 1; 1 2]
Kx = [2 1; 1 2];
[U, D] = eig(Kx); [lam, order] = sort(diag(D), 'descend'); U = U(:, order);
xs = randn(2e5, 2)*chol(Kx);                      % campioni gaussiani con covarianza Kx
Xk = xs*U;                                        % coefficienti KLT
fprintf('KLT: autovalori = %s (attesi: 3 e 1), primo autovettore = %s\n', mat2str(lam', 4), mat2str(abs(U(:,1))', 4));
disp('     covarianza dei coefficienti (diagonale: scorrelati):'); disp(round(cov(Xk), 3));
fprintf('     energia nel primo coefficiente: %.0f %%\n', 100*lam(1)/sum(lam));

% 2) matrice DCT 8x8: righe = funzioni di base
n = 8; [m, k] = meshgrid(0:n-1, 0:n-1);
Cm = sqrt(2/n)*cos(pi*(2*m+1).*k/(2*n)); Cm(1,:) = 1/sqrt(n);
fprintf('DCT: errore di ortonormalità = %.1e\n', max(max(abs(Cm*Cm' - eye(n)))));

rho = 0.95; Km = rho.^abs((0:7)' - (0:7));        % processo di Markov del primo ordine
[Uk, Dk] = eig(Km); [~, order] = sort(diag(Dk), 'descend'); Uk = Uk(:, order);
fprintf('Somiglianza |<base DCT k, autovettore KLT k>|, k = 0..7: %s\n', mat2str(abs(sum(Cm.*Uk', 2))', 4));

% 3) JPEG semplificato
Q = [16 11 10 16 24 40 51 61; 12 12 14 19 26 58 60 55; 14 13 16 24 40 57 69 56; 14 17 22 29 51 87 80 62;
     18 22 37 56 68 109 103 77; 24 35 55 64 81 104 113 92; 49 64 78 87 103 121 120 101; 72 92 95 98 112 100 103 99];
img = test_image(256)*255; scales = [1 4 12];
figure; colormap gray;
subplot(1,4,1); imagesc(img, [0 255]); axis image off; title('originale');
for s = 1:3
    rec = zeros(size(img)); nz = 0;
    for r = 1:8:size(img,1)
        for c = 1:8:size(img,2)
            X = Cm*(img(r:r+7, c:c+7) - 128)*Cm';          % DCT 2D separabile
            Xq = round(X./(Q*scales(s)));                  % quantizzazione (qui si perde informazione)
            nz = nz + nnz(Xq);
            rec(r:r+7, c:c+7) = Cm'*(Xq.*Q*scales(s))*Cm + 128;
        end
    end
    p = 10*log10(255^2/mean((img(:) - rec(:)).^2));
    fprintf('JPEG-like, Q x %2d: coefficienti non nulli = %5.1f %%, PSNR = %.1f dB\n', scales(s), 100*nz/numel(img), p);
    subplot(1,4,s+1); imagesc(rec, [0 255]); axis image off;
    title(sprintf('Q x %d: %.0f%% coeff., %.0f dB', scales(s), 100*nz/numel(img), p));
end

% le 64 funzioni di base della DCT 8x8
figure; colormap gray;
for k1 = 1:8
    for k2 = 1:8
        subplot(8, 8, (k1-1)*8 + k2); imagesc(Cm(k1,:)'*Cm(k2,:)); axis image off;
    end
end
