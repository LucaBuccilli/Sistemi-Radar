%% build_models.m - Crea e simula i modelli Simulink della repository
% Lancia questo script da MATLAB (con Simulink installato): crea i file .slx nella
% cartella simulink/, li simula e mostra i risultati.
% Usa solo blocchi della libreria base di Simulink, così funziona senza toolbox aggiuntivi.
%
% Modelli creati:
%   sim_aliasing.slx      Sine Wave -> Zero-Order Hold: campionamento e aliasing
%   sim_quantizer.slx     Sine Wave -> Quantizer: errore di quantizzazione e regola dei 6 dB/bit
%   sim_fir_lowpass.slx   segnale + rumore -> Discrete FIR Filter (passa-basso)
%   sim_mti_canceller.slx clutter + bersaglio -> canceller singolo (Unit Delay + Sum) e triplo (FIR)
clear; close all; clc;
here = fileparts(mfilename('fullpath'));
cd(here);

%% Parametri comuni (letti dai blocchi dal workspace)
Fs   = 8000;          % frequenza di campionamento [Hz]
f0   = 5000;          % tono oltre Fs/2: deve ricomparire a |f0 - Fs| = 3000 Hz
Tsim = 0.05;          % durata della simulazione [s]

%% 1) Aliasing: Sine Wave continuo -> Zero-Order Hold (Ts = 1/Fs) -> To Workspace
mdl = 'sim_aliasing';
new_model(mdl);
add_block('simulink/Sources/Sine Wave', [mdl '/Sine'], 'Frequency', '2*pi*f0', 'Position', [40 40 80 80]);
add_block('simulink/Discrete/Zero-Order Hold', [mdl '/ZOH'], 'SampleTime', '1/Fs', 'Position', [140 40 180 80]);
add_block('simulink/Sinks/To Workspace', [mdl '/out'], 'VariableName', 'x_sampled', ...
    'SaveFormat', 'Array', 'SampleTime', '1/Fs', 'Position', [240 40 300 80]);   % un campione ogni 1/Fs
add_block('simulink/Sinks/Scope', [mdl '/Scope'], 'Position', [240 110 270 140]);
add_line(mdl, 'Sine/1', 'ZOH/1'); add_line(mdl, 'ZOH/1', 'out/1'); add_line(mdl, 'ZOH/1', 'Scope/1');
set_param(mdl, 'StopTime', 'Tsim', 'Solver', 'ode45', 'MaxStep', '1/(20*f0)');
finish_model(mdl);
so = sim(mdl);
xs = so.x_sampled(:);                                        % campioni a passo 1/Fs
N = numel(xs); X = abs(fft(xs - mean(xs)));
[~, k] = max(X(1:floor(N/2)));
fprintf('[aliasing] tono a %d Hz campionato a %d Hz -> picco a %.0f Hz (atteso %d Hz)\n', ...
    f0, Fs, (k-1)*Fs/N, abs(f0 - Fs));

%% 2) Quantizzatore: Sine Wave discreto -> Quantizer -> errore = uscita - ingresso
mdl = 'sim_quantizer';
b = 6; A = 1; Delta = 2*A/2^b;  %#ok<NASGU>  % passo per b bit su [-A, A]
new_model(mdl);
add_block('simulink/Sources/Sine Wave', [mdl '/Sine'], 'Amplitude', '0.999*A', 'Frequency', '2*pi*37', ...
    'SampleTime', '1/Fs', 'Position', [40 40 80 80]);
add_block('simulink/Discontinuities/Quantizer', [mdl '/Q'], 'QuantizationInterval', 'Delta', 'Position', [140 40 180 80]);
add_block('simulink/Math Operations/Sum', [mdl '/Err'], 'Inputs', '+-', 'Position', [240 50 260 70]);
add_block('simulink/Sinks/To Workspace', [mdl '/xin'], 'VariableName', 'xin', 'SaveFormat', 'Array', 'Position', [240 120 300 150]);
add_block('simulink/Sinks/To Workspace', [mdl '/err'], 'VariableName', 'qerr', 'SaveFormat', 'Array', 'Position', [320 45 380 75]);
add_line(mdl, 'Sine/1', 'Q/1'); add_line(mdl, 'Q/1', 'Err/1'); add_line(mdl, 'Sine/1', 'Err/2');
add_line(mdl, 'Err/1', 'err/1'); add_line(mdl, 'Sine/1', 'xin/1');
set_param(mdl, 'StopTime', '1', 'Solver', 'FixedStepDiscrete', 'FixedStep', '1/Fs');
finish_model(mdl);
so = sim(mdl);
fprintf('[quantizer] b = %d bit: SNR = %.1f dB (teoria per sinusoide 6.02 b + 1.76 = %.1f dB)\n', ...
    b, 10*log10(var(so.xin)/var(so.qerr)), 6.02*b + 1.76);

%% 3) FIR passa-basso: tono a 300 Hz + tono a 2500 Hz -> Discrete FIR Filter (fc = 1000 Hz)
mdl = 'sim_fir_lowpass';
Nt = 64; n = (0:Nt) - Nt/2; fc = 1000;
h = 2*fc/Fs*sinc(2*fc/Fs*n) .* (0.54 - 0.46*cos(2*pi*(0:Nt)/Nt));  %#ok<NASGU> % sinc finestrata (Hamming)
new_model(mdl);
add_block('simulink/Sources/Sine Wave', [mdl '/Low'],  'Frequency', '2*pi*300',  'SampleTime', '1/Fs', 'Position', [40 30 80 70]);
add_block('simulink/Sources/Sine Wave', [mdl '/High'], 'Frequency', '2*pi*2500', 'SampleTime', '1/Fs', 'Position', [40 110 80 150]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum'], 'Inputs', '++', 'Position', [130 80 150 100]);
add_block('simulink/Discrete/Discrete FIR Filter', [mdl '/FIR'], 'Coefficients', 'h', 'SampleTime', '1/Fs', 'Position', [200 70 260 110]);
add_block('simulink/Sinks/To Workspace', [mdl '/in'],  'VariableName', 'fir_in',  'SaveFormat', 'Array', 'Position', [200 150 260 180]);
add_block('simulink/Sinks/To Workspace', [mdl '/out'], 'VariableName', 'fir_out', 'SaveFormat', 'Array', 'Position', [310 75 370 105]);
add_block('simulink/Sinks/Scope', [mdl '/Scope'], 'NumInputPorts', '2', 'Position', [310 130 340 170]);
add_line(mdl, 'Low/1', 'Sum/1'); add_line(mdl, 'High/1', 'Sum/2'); add_line(mdl, 'Sum/1', 'FIR/1');
add_line(mdl, 'Sum/1', 'in/1'); add_line(mdl, 'FIR/1', 'out/1');
add_line(mdl, 'Sum/1', 'Scope/1'); add_line(mdl, 'FIR/1', 'Scope/2');
set_param(mdl, 'StopTime', '0.1', 'Solver', 'FixedStepDiscrete', 'FixedStep', '1/Fs');
finish_model(mdl);
so = sim(mdl);
N = numel(so.fir_out); f = (0:N-1)*Fs/N;
Xin = abs(fft(so.fir_in)); Xout = abs(fft(so.fir_out));
[~, i300] = min(abs(f - 300)); [~, i2500] = min(abs(f - 2500));
fprintf('[FIR] attenuazione a 300 Hz = %.1f dB, a 2500 Hz = %.1f dB\n', ...
    20*log10(Xout(i300)/Xin(i300)), 20*log10(Xout(i2500)/Xin(i2500)));

%% 4) Canceller MTI: clutter correlato (lento) + bersaglio a PRF/4, campione per impulso
mdl = 'sim_mti_canceller';
PRF = 1250; sigma_f = 4; Np = 4096;
rng(2);
fr = (0:Np-1)*PRF/Np; fr(fr >= PRF/2) = fr(fr >= PRF/2) - PRF;
cl = real(ifft(exp(-fr.^2/(4*sigma_f^2)) .* (randn(1, Np) + 1j*randn(1, Np))));
cl = cl/std(cl) * 100;                                   % clutter 40 dB sopra il bersaglio
tg = cos(2*pi*(PRF/4)*(0:Np-1)/PRF);                     % bersaglio con Doppler PRF/4
t = (0:Np-1)'/PRF;
clutter_ts = [t cl(:)]; target_ts = [t tg(:)];           %#ok<NASGU>
new_model(mdl);
add_block('simulink/Sources/From Workspace', [mdl '/Clutter'], 'VariableName', 'clutter_ts', 'SampleTime', '1/PRF', 'Position', [30 30 90 60]);
add_block('simulink/Sources/From Workspace', [mdl '/Target'],  'VariableName', 'target_ts',  'SampleTime', '1/PRF', 'Position', [30 100 90 130]);
add_block('simulink/Math Operations/Sum', [mdl '/Rx'], 'Inputs', '++', 'Position', [130 70 150 90]);
% canceller singolo costruito "a mano": y[n] = x[n] - x[n-1]
add_block('simulink/Discrete/Unit Delay', [mdl '/z^-1'], 'SampleTime', '1/PRF', 'Position', [190 120 220 150]);
add_block('simulink/Math Operations/Sum', [mdl '/SC'], 'Inputs', '+-', 'Position', [260 70 280 90]);
% canceller triplo come filtro FIR binomiale [1 -3 3 -1]
add_block('simulink/Discrete/Discrete FIR Filter', [mdl '/TC'], 'Coefficients', '[1 -3 3 -1]', 'SampleTime', '1/PRF', 'Position', [190 190 250 220]);
add_block('simulink/Sinks/To Workspace', [mdl '/in'],  'VariableName', 'mti_in',  'SaveFormat', 'Array', 'Position', [330 10 390 40]);
add_block('simulink/Sinks/To Workspace', [mdl '/sc'],  'VariableName', 'mti_sc',  'SaveFormat', 'Array', 'Position', [330 65 390 95]);
add_block('simulink/Sinks/To Workspace', [mdl '/tc'],  'VariableName', 'mti_tc',  'SaveFormat', 'Array', 'Position', [330 190 390 220]);
add_line(mdl, 'Clutter/1', 'Rx/1'); add_line(mdl, 'Target/1', 'Rx/2');
add_line(mdl, 'Rx/1', 'SC/1'); add_line(mdl, 'Rx/1', 'z^-1/1'); add_line(mdl, 'z^-1/1', 'SC/2');
add_line(mdl, 'Rx/1', 'TC/1'); add_line(mdl, 'Rx/1', 'in/1');
add_line(mdl, 'SC/1', 'sc/1'); add_line(mdl, 'TC/1', 'tc/1');
set_param(mdl, 'StopTime', num2str((Np-1)/PRF), 'Solver', 'FixedStepDiscrete', 'FixedStep', '1/PRF');
finish_model(mdl);
so = sim(mdl);
skip = 10;                                               % transitorio iniziale dei filtri
fprintf('[MTI] potenza in ingresso: clutter/bersaglio = %.1f dB\n', 10*log10(var(cl)/0.5));
fprintf('[MTI] potenza in uscita: singolo = %.1f dB, triplo = %.1f dB rispetto all''ingresso\n', ...
    10*log10(var(so.mti_sc(skip:end))/var(so.mti_in)), 10*log10(var(so.mti_tc(skip:end))/var(so.mti_in)));
figure;
subplot(3,1,1); plot(so.mti_in); title('Ingresso: clutter + bersaglio (il bersaglio non si vede)');
subplot(3,1,2); plot(so.mti_sc); title('Dopo il canceller singolo');
subplot(3,1,3); plot(so.mti_tc); title('Dopo il canceller triplo: resta il bersaglio');

disp('Modelli creati: sim_aliasing.slx, sim_quantizer.slx, sim_fir_lowpass.slx, sim_mti_canceller.slx');

%% funzioni di supporto
function new_model(name)
if bdIsLoaded(name), close_system(name, 0); end
if exist([name '.slx'], 'file'), delete([name '.slx']); end
new_system(name);
open_system(name);
end

function finish_model(name)
save_system(name);
end
