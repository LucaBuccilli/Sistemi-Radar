%% Quantizzazione: uniforme (regola dei 6 dB per bit) e ottima di Lloyd-Max
% Uniforme su [-A, A] con b bit: passo D = 2A/2^b, distorsione D^2/12 = sigma_x^2 / 2^(2b).
% Lloyd-Max: soglie a metà tra i livelli, livelli = baricentri degli intervalli.
clear; close all; clc; rng(4);
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils'));

% 1) quantizzatore uniforme
A = 1; x = 2*A*rand(5e5, 1) - A;
fprintf('bit | SNR misurato | teoria 6.02*b\n');
for b = 2:8
    step = 2*A/2^b;
    xq = step*(floor(x/step) + 0.5);                       % centro dell'intervallo
    fprintf(' %d  |   %6.2f dB  |  %6.2f dB\n', b, 10*log10(var(x)/mean((x-xq).^2)), 6.02*b);
end

% 2) Lloyd-Max per una gaussiana N(0,1)
g = randn(1e6, 1);
expected = {'', '+-0.7979', '0, +-1.224', '+-0.4528, +-1.5104'};
for L = 2:4
    [q, theta, dist] = lloyd_max(g, L);
    fprintf('Lloyd-Max %d livelli: q = [%s], soglie = [%s], D = %.4f   (valori tabulati: %s)\n', ...
        L, num2str(q', '%.3f '), num2str(theta', '%.3f '), dist, expected{L});
end
step = 6/4;                                                % confronto: uniforme a 4 livelli su [-3, 3]
gu = min(max(step*(floor(g/step) + 0.5), -2.25), 2.25);
[~, ~, dlm] = lloyd_max(g, 4);
fprintf('Gaussiana, 4 livelli: D uniforme su [-3,3] = %.4f, D Lloyd-Max = %.4f\n', mean((g-gu).^2), dlm);

% 3) luminanza di un'immagine a 8, 4 e 2 bit
img = test_image(256); bits = [8 4 2];
figure; colormap gray;
for i = 1:3
    L = 2^bits(i);
    subplot(1,3,i); imagesc(round(img*(L-1))/(L-1), [0 1]); axis image off;
    title(sprintf('%d bit (%d livelli)', bits(i), L));
end

function [q, theta, dist] = lloyd_max(s, L)
s = sort(s);
q = s(round(((1:L)' - 0.5)/L*numel(s)));                   % inizializzazione con i quantili
for it = 1:200
    theta = (q(1:end-1) + q(2:end))/2;                     % soglie a metà tra i livelli
    edges = [-inf; theta; inf];
    for k = 1:L
        q(k) = mean(s(s >= edges(k) & s < edges(k+1)));    % baricentro dell'intervallo
    end
end
theta = (q(1:end-1) + q(2:end))/2;
edges = [-inf; theta; inf]; sq = zeros(size(s));
for k = 1:L
    sq(s >= edges(k) & s < edges(k+1)) = q(k);
end
dist = mean((s - sq).^2);
end
