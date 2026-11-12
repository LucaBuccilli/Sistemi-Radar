%% Segnali su grafo: Laplaciano L = D - A e trasformata di Fourier su grafo (GFT)
% Base = autovettori di L; gli autovalori fanno la parte delle frequenze.
% Un segnale "liscio" sul grafo ha i coefficienti concentrati sugli autovalori piccoli.
clear; close all; clc; rng(8);

% grafo con 3 gruppi (cluster) di 10 nodi: molti archi dentro i gruppi, pochi tra gruppi
centers = [0 0; 3 0.5; 1.5 2.5];
pos = []; label = [];
for g = 1:3
    pos = [pos; centers(g,:) + 0.5*randn(10, 2)]; %#ok<AGROW>
    label = [label; g*ones(10, 1)];               %#ok<AGROW>
end
n = size(pos, 1); A = zeros(n);
for i = 1:n
    for j = i+1:n
        if label(i) == label(j), p = 0.7; else, p = 0.02; end
        A(i,j) = rand < p; A(j,i) = A(i,j);
    end
end
links = [1 11; 11 21; 21 1];                      % garantisce che il grafo sia connesso
for r = 1:3, A(links(r,1), links(r,2)) = 1; A(links(r,2), links(r,1)) = 1; end

D = diag(sum(A, 2));
Lap = D - A;
[U, E] = eig(Lap); [lam, order] = sort(diag(E)); U = U(:, order);   % "frequenze" crescenti
fprintf('Autovalori più piccoli: %s (il primo è 0: autovettore costante)\n', mat2str(round(abs(lam(1:4))', 3)));

vals = [1; 3; -2];
smooth = vals(label) + 0.1*randn(n, 1);           % quasi costante dentro i gruppi
rough = randn(n, 1);                              % varia a caso tra nodi vicini
sigs = {smooth, rough}; names = {'liscio', 'irregolare'};
for i = 1:2
    s = sigs{i}; S = U'*s;                        % GFT: S_k = <s, u_k>
    fprintf('Segnale %-10s: s''Ls = %7.1f, energia nei primi 3 coefficienti = %5.1f %%\n', ...
        names{i}, s'*Lap*s, 100*sum(S(1:3).^2)/sum(S.^2));
end
fprintf('Ricostruzione s = U S: errore = %.1e\n', max(abs(U*(U'*smooth) - smooth)));

figure;
subplot(1,3,1); hold on;
[ii, jj] = find(triu(A));
for e = 1:numel(ii)
    plot(pos([ii(e) jj(e)], 1), pos([ii(e) jj(e)], 2), 'Color', [0.8 0.8 0.8]);
end
scatter(pos(:,1), pos(:,2), 40, smooth, 'filled'); colormap(gca, jet); colorbar; axis off;
title('Segnale liscio sul grafo');
subplot(1,3,2); stem(lam, abs(U'*smooth)); title('GFT del segnale liscio'); xlabel('autovalore');
subplot(1,3,3); stem(lam, abs(U'*rough));  title('GFT del segnale irregolare'); xlabel('autovalore');
