%% SAR: chirp in azimut e focalizzazione. La risoluzione lungo traccia vale circa d_a/2
% Scenario avionico in banda X: lambda = 3.1 cm, antenna d_a = 1.8 m, R = 20 km.
clear; close all; clc;
lam = 0.031; da = 1.8; Ry = 20e3; v = 100;
Ds = lam*Ry/da;                 % impronta in azimut (apertura sintetica)
Tobs = Ds/v;                    % tempo di osservazione di un punto
Bd = 2*v/da;                    % banda Doppler
PRF = 1.5*Bd;                   % sopra il minimo 2v/d_a
fprintf('Impronta reale: %.0f m, Tobs = %.2f s, banda Doppler = %.1f Hz, PRF = %.1f Hz\n', Ds, Tobs, Bd, PRF);
fprintf('Pendenza del chirp: %.2f Hz/s, prodotto tempo-banda = %.0f\n', 2*v^2/(lam*Ry), Bd*Tobs);

ta = -2*Tobs : 1/PRF : 2*Tobs;  % tempo lento
x = v*ta;                       % posizione della piattaforma
% eco di un punto in x0: fase -4 pi R(t)/lambda, visibile solo dentro il fascio
echo = @(x0) exp(-1j*4*pi/lam*(sqrt(Ry^2 + (x-x0).^2) - Ry)) .* (abs(x-x0) <= Ds/2);

ref = echo(0);
sig = echo(0) + echo(2);        % due punti distanti 2 m lungo traccia

% Filtro adattato al punto in x0 = correlazione con l'eco atteso.
% Il passo tra impulsi è v/PRF = 0.6 m: si valuta l'uscita su una griglia più fitta.
xf = -6:0.02:8;
foc = zeros(size(xf)); single = zeros(size(xf));
for i = 1:numel(xf)
    e = echo(xf(i));
    foc(i) = abs(sum(conj(e).*sig));
    single(i) = abs(sum(conj(e).*ref));
end
foc = foc/max(foc); single = single/max(single);
in3 = xf(single > 1/sqrt(2));
fprintf('Risoluzione reale (senza SAR): %.0f m\n', Ds);
fprintf('Risoluzione focalizzata (-3 dB): %.2f m; teoria 0.886*d_a/2 = %.2f m (ordine d_a/2 = %.2f m)\n', ...
    max(in3) - min(in3), 0.886*da/2, da/2);

figure;
subplot(1,3,1); plot(ta, real(ref)); xlim([-Tobs/2 Tobs/2]);
title('Eco di un punto: chirp in tempo lento'); xlabel('tempo lento [s]');
fi = gradient(unwrap(angle(ref)))*PRF/(2*pi); m = abs(ta) < Tobs/2;
subplot(1,3,2); plot(ta(m), fi(m)); grid on;
title('Frequenza Doppler istantanea'); xlabel('tempo lento [s]'); ylabel('Hz');
subplot(1,3,3); plot(xf, foc); grid on;
title('Dopo la focalizzazione: due punti a 2 m'); xlabel('lungo traccia [m]');
