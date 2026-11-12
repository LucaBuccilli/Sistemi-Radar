%% build_models_dsp.m - Modelli Simulink con DSP System Toolbox e versione in virgola fissa
% Requisiti:
%   sim_fft_spectrum    -> Simulink + DSP System Toolbox (Buffer, FFT, Spectrum Analyzer)
%   sim_fir_fixed_point -> Simulink + Fixed-Point Designer (simulazione con tipi fixdt)
% Lo script crea i modelli, li simula e confronta i risultati con la teoria.
% Gli stessi esperimenti in MATLAB/Python puri sono in matlab/dsp e python/dsp.
clear; close all; clc;
here = fileparts(mfilename('fullpath'));
cd(here);

%% Parametri comuni
Fs = 8000;                 % frequenza di campionamento [Hz]
N  = 256;                  % lunghezza del frame per la FFT
f0 = 1000;                 % tono sul bin k = f0*N/Fs = 32
hannwin = 0.5 - 0.5*cos(2*pi*(0:N-1)'/(N-1));   %#ok<NASGU> finestra di Hann (vettore colonna)

%% 1) Spettro a frame: Sine + rumore -> Buffer -> finestra -> FFT -> |.| ; Spectrum Analyzer in parallelo
mdl = 'sim_fft_spectrum';
new_model(mdl);
add_block('simulink/Sources/Sine Wave', [mdl '/Sine'], 'Frequency', '2*pi*f0', 'SampleTime', '1/Fs', 'Position', [30 30 70 70]);
add_block('simulink/Sources/Random Number', [mdl '/Noise'], 'Variance', '1e-4', 'SampleTime', '1/Fs', 'Position', [30 100 70 140]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum'], 'Inputs', '++', 'Position', [110 70 130 90]);
add_block(lib_block('Buffer'), [mdl '/Buffer'], 'Position', [170 60 230 100]);
set_param([mdl '/Buffer'], 'N', 'N');                                      % frame da N campioni, senza sovrapposizione
add_block('simulink/Sources/Constant', [mdl '/Hann'], 'Value', 'hannwin', 'Position', [170 130 230 160]);
add_block('simulink/Math Operations/Product', [mdl '/Window'], 'Multiplication', 'Element-wise(.*)', 'Position', [270 70 300 110]);
add_block(lib_block('FFT'), [mdl '/FFT'], 'Position', [340 70 400 110]);
add_block('simulink/Math Operations/Abs', [mdl '/Abs'], 'Position', [440 75 480 105]);
add_block('simulink/Sinks/To Workspace', [mdl '/out'], 'VariableName', 'fft_mag', 'SaveFormat', 'Array', 'Position', [520 75 580 105]);
add_block(lib_block('Spectrum Analyzer'), [mdl '/Spectrum Analyzer'], 'Position', [170 190 230 230]);
add_line(mdl, 'Sine/1', 'Sum/1'); add_line(mdl, 'Noise/1', 'Sum/2');
add_line(mdl, 'Sum/1', 'Buffer/1'); add_line(mdl, 'Sum/1', 'Spectrum Analyzer/1');
add_line(mdl, 'Buffer/1', 'Window/1'); add_line(mdl, 'Hann/1', 'Window/2');
add_line(mdl, 'Window/1', 'FFT/1'); add_line(mdl, 'FFT/1', 'Abs/1'); add_line(mdl, 'Abs/1', 'out/1');
set_param(mdl, 'StopTime', num2str(64*N/Fs), 'Solver', 'FixedStepDiscrete');
save_system(mdl);
so = sim(mdl);
X = squeeze(so.fft_mag);                  % N x numero di frame (o viceversa, a seconda della release)
if size(X, 1) ~= N, X = X.'; end
S = mean(X(1:N/2+1, 2:end).^2, 2);        % spettro mediato (si scarta il primo frame, parzialmente vuoto)
[~, k] = max(S);
fprintf('[FFT] picco al bin %d (atteso %d = f0*N/Fs), risoluzione Fs/N = %.2f Hz\n', k-1, f0*N/Fs, Fs/N);
fprintf('[FFT] livello del rumore rispetto al picco: %.1f dB\n', 10*log10(median(S)/max(S)));
figure; plot((0:N/2)*Fs/N, 10*log10(S/max(S))); grid on; xlabel('Hz'); ylabel('dB');
title('Spettro mediato da Simulink (Buffer -> Hann -> FFT)');

%% 2) FIR in virgola fissa: confronto con la versione double al variare della lunghezza di parola
mdl = 'sim_fir_fixed_point';
Nt = 64; n = (0:Nt) - Nt/2; fc = 1000;
sincf = @(x) (x == 0) + (x ~= 0).*sin(pi*x)./(pi*x + (x == 0));
h = 2*fc/Fs*sincf(2*fc/Fs*n) .* (0.54 - 0.46*cos(2*pi*(0:Nt)/Nt));   %#ok<NASGU> passa-basso a 1 kHz
WL = 16;                                                               % lunghezza di parola (modificata nel ciclo)
new_model(mdl);
add_block('simulink/Sources/Sine Wave', [mdl '/Low'],  'Frequency', '2*pi*300',  'SampleTime', '1/Fs', 'Position', [30 30 70 70]);
add_block('simulink/Sources/Sine Wave', [mdl '/High'], 'Frequency', '2*pi*2500', 'SampleTime', '1/Fs', 'Position', [30 100 70 140]);
add_block('simulink/Math Operations/Sum', [mdl '/Sum'], 'Inputs', '++', 'Position', [110 70 130 90]);
% ramo di riferimento in double
add_block('simulink/Discrete/Discrete FIR Filter', [mdl '/FIR double'], 'Coefficients', 'h', 'SampleTime', '1/Fs', 'Position', [200 20 270 60]);
% ramo in virgola fissa: ingresso Q(WL, WL-2), coefficienti Q(WL, WL-1), uscita Q(WL, WL-2)
add_block('simulink/Signal Attributes/Data Type Conversion', [mdl '/To fixed'], 'OutDataTypeStr', 'fixdt(1,WL,WL-2)', ...
    'RndMeth', 'Nearest', 'SaturateOnIntegerOverflow', 'on', 'Position', [170 110 230 140]);
add_block('simulink/Discrete/Discrete FIR Filter', [mdl '/FIR fixed'], 'Coefficients', 'h', 'SampleTime', '1/Fs', 'Position', [270 105 340 145]);
fixed_params = {'CoefDataTypeStr', 'fixdt(1,WL,WL-1)'; 'ProductDataTypeStr', 'fixdt(1,2*WL+8,2*WL-3)'; ...
                'AccumDataTypeStr', 'fixdt(1,2*WL+8,2*WL-3)'; 'OutDataTypeStr', 'fixdt(1,WL,WL-2)'; ...
                'RndMeth', 'Nearest'; 'SaturateOnIntegerOverflow', 'on'};
for i = 1:size(fixed_params, 1)
    try
        set_param([mdl '/FIR fixed'], fixed_params{i,1}, fixed_params{i,2});
    catch err
        warning('Parametro %s non impostato (%s): controllarlo a mano nel blocco.', fixed_params{i,1}, err.message);
    end
end
add_block('simulink/Signal Attributes/Data Type Conversion', [mdl '/To double'], 'OutDataTypeStr', 'double', 'Position', [380 110 440 140]);
add_block('simulink/Math Operations/Sum', [mdl '/Err'], 'Inputs', '+-', 'Position', [490 70 510 90]);
add_block('simulink/Sinks/To Workspace', [mdl '/ref'], 'VariableName', 'y_ref', 'SaveFormat', 'Array', 'Position', [550 15 610 45]);
add_block('simulink/Sinks/To Workspace', [mdl '/err'], 'VariableName', 'y_err', 'SaveFormat', 'Array', 'Position', [550 65 610 95]);
add_line(mdl, 'Low/1', 'Sum/1'); add_line(mdl, 'High/1', 'Sum/2');
add_line(mdl, 'Sum/1', 'FIR double/1'); add_line(mdl, 'Sum/1', 'To fixed/1');
add_line(mdl, 'To fixed/1', 'FIR fixed/1'); add_line(mdl, 'FIR fixed/1', 'To double/1');
add_line(mdl, 'FIR double/1', 'Err/1'); add_line(mdl, 'To double/1', 'Err/2');
add_line(mdl, 'FIR double/1', 'ref/1'); add_line(mdl, 'Err/1', 'err/1');
set_param(mdl, 'StopTime', '1', 'Solver', 'FixedStepDiscrete', 'FixedStep', '1/Fs');
save_system(mdl);

wls = [8 12 16]; snr = zeros(size(wls));
for i = 1:numel(wls)
    WL = wls(i);                                                      %#ok<NASGU> letto dai blocchi
    so = sim(mdl);
    snr(i) = 10*log10(mean(so.y_ref.^2)/mean(so.y_err.^2));
    fprintf('[virgola fissa] WL = %2d bit: SNR uscita = %.1f dB\n', wls(i), snr(i));
end
fprintf('[virgola fissa] pendenza: %.2f dB per bit (teoria circa 6.02)\n', (snr(end) - snr(1))/(wls(end) - wls(1)));
disp('Modelli creati: sim_fft_spectrum.slx, sim_fir_fixed_point.slx');

%% funzioni di supporto
function new_model(name)
if bdIsLoaded(name), close_system(name, 0); end
if exist([name '.slx'], 'file'), delete([name '.slx']); end
new_system(name);
open_system(name);
end

function p = lib_block(name)
% Percorso del blocco nella libreria DSP System Toolbox (i nomi delle librerie interne
% sono stabili da molte release; se cambiano, trascinare il blocco a mano dal Library Browser).
candidates = struct('Buffer', {{'dspbuff3/Buffer'}}, ...
                    'FFT', {{'dspxfrm3/FFT'}}, ...
                    'SpectrumAnalyzer', {{'dspsnks4/Spectrum Analyzer'}});
key = strrep(name, ' ', '');
for c = candidates.(key)
    lib = strtok(c{1}, '/');
    try
        load_system(lib);
        if getSimulinkBlockHandle(c{1}) ~= -1
            p = c{1};
            return
        end
    catch
    end
end
error('Blocco "%s" non trovato: serve DSP System Toolbox.', name);
end
