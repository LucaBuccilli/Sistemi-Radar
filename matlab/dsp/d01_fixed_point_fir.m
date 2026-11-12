%% Filtro FIR in virgola fissa: effetto della lunghezza di parola
% Formato Q con segno: w bit totali, f bit frazionari (passo 2^-f).
%   ingresso/uscita: f = w-2 (range +-2), coefficienti: f = w-1, accumulatore a precisione piena.
% La quantizzazione è simulata con round + saturazione, senza Fixed-Point Designer
% (con il toolbox si userebbe fi(x, 1, w, f)).
clear; close all; clc;
Fs = 8000; fc = 1000; Nt = 64;
n = (0:Nt) - Nt/2;
sincf = @(x) (x == 0) + (x ~= 0).*sin(pi*x)./(pi*x + (x == 0));
h = 2*fc/Fs*sincf(2*fc/Fs*n) .* (0.54 - 0.46*cos(2*pi*(0:Nt)/Nt));   % sinc con finestra di Hamming

t = (0:7999)/Fs;
x = sin(2*pi*300*t) + sin(2*pi*2500*t);
y_ref = filter(h, 1, x);

Nf = 8192; f = (0:Nf/2-1)*Fs/Nf; stop = f > 1500;
Href = abs(fft(h, Nf)); Href = Href(1:Nf/2);
wl = [8 10 12 14 16 20]; snr = zeros(size(wl)); att = zeros(size(wl));
figure; hold on; grid on;
fprintf('bit | SNR uscita [dB] | attenuazione in banda oscura [dB]\n');
for i = 1:numel(wl)
    w = wl(i);
    hq = quantize(h, w, w-1);
    yq = quantize(filter(hq, 1, quantize(x, w, w-2)), w, w-2);
    snr(i) = 10*log10(mean(y_ref.^2)/mean((y_ref - yq).^2));
    H = abs(fft(hq, Nf)); H = H(1:Nf/2);
    att(i) = -20*log10(max(H(stop)));
    fprintf(' %2d |     %6.1f      |     %6.1f\n', w, snr(i), att(i));
    if any(w == [8 12 16]), plot(f, 20*log10(H + 1e-12), 'DisplayName', sprintf('%d bit', w)); end
end
plot(f, 20*log10(Href + 1e-12), 'k--', 'DisplayName', 'double');
ylim([-110 5]); xlabel('Hz'); ylabel('|H| [dB]'); legend; title('Coefficienti quantizzati');
fprintf('Pendenza misurata: %.2f dB per bit (teoria circa 6.02)\n', (snr(wl==16) - snr(wl==8))/8);
fprintf('Attenuazione in banda oscura in doppia precisione: %.1f dB\n', -20*log10(max(Href(stop))));

figure; plot(wl, snr, 'o-'); hold on; plot(wl, snr(1) + 6.02*(wl - 8), 'k--'); grid on;
xlabel('lunghezza di parola [bit]'); ylabel('SNR uscita [dB]'); legend('misurato', '6.02 dB/bit');

function y = quantize(x, w, f)
q = 2^-f; lo = -2^(w-f-1); hi = 2^(w-f-1) - q;
y = min(max(round(x/q)*q, lo), hi);
end
