%% Elaborazione a frame: Buffer -> finestra -> FFT -> spettro mediato
% Equivalente MATLAB della catena Simulink Buffer -> FFT -> |.|^2 (modello sim_fft_spectrum).
clear; close all; clc; rng(10);
Fs = 8000; N = 256; M = 64; df = Fs/N;
t = (0:N*M-1)/Fs;
noise = 1e-4*randn(size(t));                                % rumore molto basso: rende visibile il leakage
rect = ones(1, N); hann = 0.5 - 0.5*cos(2*pi*(0:N-1)/(N-1));
f = (0:N/2)*Fs/N;

f0s = [1000, 1000 + df/2]; labels = {'sul bin', 'tra due bin'};
wins = {rect, hann}; wnames = {'rettangolare', 'Hann'};
figure; hold on; grid on;
for a = 1:2
    x = sin(2*pi*f0s(a)*t) + noise;
    for b = 1:2
        S = averaged_spectrum(x, wins{b}, N, M);
        [~, k] = max(S);
        far = 10*log10(median(S(abs(f - f0s(a)) > 1000))/max(S));
        fprintf('tono %-11s (%7.1f Hz), finestra %-12s: picco al bin %d (%.1f Hz), livello a >1 kHz = %6.1f dB\n', ...
            labels{a}, f0s(a), wnames{b}, k-1, f(k), far);
        if a == 2, plot(f, 10*log10(S/max(S)), 'DisplayName', wnames{b}); end
    end
end
xlim([0 4000]); ylim([-110 5]); legend; xlabel('Hz'); ylabel('dB'); title('Tono tra due bin: leakage');
fprintf('Risoluzione Fs/N = %.2f Hz; bin atteso per 1000 Hz = %d\n', df, 1000*N/Fs);

single = abs(fft(noise(1:N))).^2/N; single = single(1:N/2+1);
avg = averaged_spectrum(noise, rect, N, M);
r1 = std(single(6:end-5))/mean(single(6:end-5)); rM = std(avg(6:end-5))/mean(avg(6:end-5));
fprintf('Rumore: dev. std relativa dello spettro, 1 frame = %.2f, %d frame = %.3f (atteso 1/sqrt(M) = %.3f)\n', ...
    r1, M, rM, 1/sqrt(M));

function S = averaged_spectrum(x, win, N, M)
frames = reshape(x(1:N*M), N, M).';                        % Buffer: M frame da N campioni
X = fft(frames.*win, [], 2);
S = mean(abs(X(:, 1:N/2+1)).^2, 1)/sum(win.^2);
end
