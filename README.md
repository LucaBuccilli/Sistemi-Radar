# Radar & Image Signal Processing – MATLAB · Simulink · Python

Raccolta di simulazioni che ho sviluppato per consolidare, con il codice, gli argomenti di **elaborazione dei segnali radar** e di **elaborazione delle immagini** studiati durante la laurea magistrale in Ingegneria delle Telecomunicazioni (Sapienza Università di Roma).

Per ogni argomento ho scritto lo **stesso esperimento in Python e in MATLAB**, con parametri identici, e ho confrontato i risultati numerici con i valori attesi dalla teoria. Per gli argomenti più vicini all'implementazione in tempo reale (campionamento, quantizzazione, filtraggio FIR, cancellazione del clutter) ho costruito anche dei **modelli Simulink**, generati da script.

## Obiettivi
- Trasformare le formule in **esperimenti numerici verificabili**: ogni script stampa il risultato misurato accanto a quello teorico.
- Confrontare due ecosistemi (MATLAB/Simulink e Python/NumPy/SciPy) sullo stesso problema.
- Tenere il codice **senza dipendenze pesanti**: niente toolbox MATLAB, solo NumPy, SciPy e Matplotlib in Python, immagini di test sintetiche.

## Struttura
```
├── python/
│   ├── radar/       r01 … r06
│   ├── immagini/    i01 … i07
│   ├── figures/     figure generate
│   └── utils.py     funzioni comuni (salvataggio figure, immagine di test, PSNR)
├── matlab/
│   ├── radar/       r01 … r06
│   ├── immagini/    i01 … i07
│   └── utils/       test_image.m
├── simulink/
│   └── build_models.m   costruisce e simula i modelli .slx con add_block/add_line
├── run_all.py       lancia tutte le simulazioni Python
└── requirements.txt
```

## Avvio rapido
```bash
pip install -r requirements.txt
python run_all.py            # 13 simulazioni, circa 1 minuto
python run_all.py r05 i05    # solo alcune
```
In MATLAB: aprire uno script in `matlab/` e premere **Run**. Per i modelli Simulink: `cd simulink` e poi `build_models`.

---

## Parte 1 – Elaborazione dei segnali radar

### r01 · Equazione radar e portata
SNR in funzione della distanza per un radar in banda S, portata massima per un SNR richiesto, dimensionamento della potenza di picco e orizzonte radar con raggio terrestre equivalente 4/3.

```math
\mathrm{SNR}(R)=\frac{P_t\,G^2\,\lambda^2\,\sigma}{(4\pi)^3\,R^4\,L\,k\,T_0\,B\,F}
\qquad
R_{max}=\left[\frac{P_t\,G^2\,\lambda^2\,\sigma}{(4\pi)^3\,L\,k\,T_0\,B\,F\,\mathrm{SNR}_{min}}\right]^{1/4}
```
```math
R_{oriz}\approx\sqrt{2R_e h_A}+\sqrt{2R_e h_T},\qquad R_e=\tfrac43\cdot 6371\ \mathrm{km}
```
- SNR da **+127 dB a 100 m** a **−17 dB a 400 km**: 40 dB in meno per ogni decade di distanza.
- Potenza di picco richiesta per 150 km su un bersaglio di 5 m²: **88 kW**.

![SNR vs distanza](python/figures/r01_snr_vs_range.png)

### r02 · Rivelazione a singolo impulso
Soglia di Neyman–Pearson sull'inviluppo, Pfa e Pd stimate con Monte Carlo (inviluppo di Rayleigh e di Rice), confronto con l'approssimazione di Albersheim.

Con rumore complesso gaussiano di potenza $\sigma_d^2$ e $\mathrm{SNR}=|a|^2/\sigma_d^2$:
```math
P_{fa}=e^{-T^2/\sigma_d^2}\;\Rightarrow\;T=\sigma_d\sqrt{\ln(1/P_{fa})}
\qquad
P_d=\int_{T/\sigma_d}^{\infty}2t\,e^{-(t^2+\mathrm{SNR})}\,I_0\!\left(2t\sqrt{\mathrm{SNR}}\right)dt
```
Approssimazione di Albersheim, con $A=\ln(0.62/P_{fa})$ e $B=\ln\big(P_d/(1-P_d)\big)$:
```math
\mathrm{SNR}_{dB}\approx-5\log_{10}N+\left(6.2+\frac{4.54}{\sqrt{N+0.44}}\right)\log_{10}\left(A+0.12AB+1.7B\right)
```
- Per **Pd = 0.9, Pfa = 10⁻⁶** servono circa **13 dB** di SNR.

![Pd vs SNR](python/figures/r02_pd_vs_snr.png)

### r03 · Integrazione di impulsi
Rivelatore coerente (somma con compensazione della fase Doppler, cioè un filtro della FFT) contro rivelatore non coerente (energia), a parità di Pfa.

```math
\text{coerente:}\ \left|\sum_{n=0}^{N-1} r_n\,e^{-j2\pi f_d nT}\right|\gtrless\lambda
\qquad
\text{non coerente:}\ \sum_{n=0}^{N-1}|r_n|^2\gtrless\lambda
\qquad
G_{coer}=10\log_{10}N\ \mathrm{dB}
```
- Con N = 16: guadagno coerente **12.0 dB = 10·log₁₀N**, non coerente 8.8 dB.

![Integrazione](python/figures/r03_pulse_integration.png)

### r04 · Cancellazione del clutter (MTI)
Canceller binomiali di ordine 1–3 contro clutter con spettro gaussiano: attenuazione del clutter e Improvement Factor, con verifica tramite simulazione del clutter correlato.

```math
H_M(f)=\left(1-e^{-j2\pi fT}\right)^M,\qquad |H_M(f)|^2=\left[4\sin^2(\pi fT)\right]^M,\qquad \rho(kT)=e^{-2\pi^2\sigma_f^2k^2T^2},\quad \sigma_f=\frac{2\sigma_v}{\lambda}
```
```math
CA_1=\frac{1}{2\,(1-\rho(T))}
\qquad
IF(f_d)=\frac{|H_M(f_d)|^2}{\sum_{m=-M}^{M}(-1)^m\binom{2M}{M+m}\rho(mT)}
```
- IF a PRF/2: **40 / 75 / 108 dB**; la simulazione del canceller singolo coincide con la teoria entro 0.1 dB.

![MTI](python/figures/r04_mti_improvement_factor.png)

### r05 · Compressione d'impulso e mappa range-Doppler
Chirp con prodotto tempo-banda 30, filtro adattato, treno di 32 impulsi con due bersagli. Compressione in tempo veloce, poi FFT in tempo lento.

```math
s(t)=e^{\,j\pi\frac{B}{\tau}t^2},\qquad h(t)=s^*(-t),\qquad \Delta R=\frac{c\,\tau}{2}\;\longrightarrow\;\Delta R=\frac{c}{2B}
```
```math
f_d=\frac{2v}{\lambda},\qquad \Delta v=\frac{\lambda}{2NT},\qquad v_{na}=\pm\frac{\lambda}{4T},\qquad R_{na}=\frac{c\,T}{2}
```
- Risoluzione in distanza da **1500 m a 50 m**; entrambi i bersagli localizzati correttamente in distanza e velocità.

![Range-Doppler](python/figures/r05_range_doppler.png)

### r06 · SAR: focalizzazione in azimut
Eco di un punto visto da una piattaforma in moto (chirp in tempo lento) e compressione con il filtro adattato.

```math
R(t_a)=\sqrt{R_y^2+v^2t_a^2}\approx R_y+\frac{v^2t_a^2}{2R_y}
\;\Rightarrow\;
s(t_a)\approx e^{-j\pi\frac{2v^2}{\lambda R_y}t_a^2}
```
```math
D_s=\frac{\lambda R_y}{d_a},\qquad T_{obs}=\frac{D_s}{v},\qquad B_d=\frac{2v}{d_a},\qquad \delta x=\frac{v}{B_d}=\frac{d_a}{2}
```
- Risoluzione lungo traccia da **344 m** (apertura reale) a **0.76 m**, coerente con il limite teorico ≈ d_a/2.

![SAR](python/figures/r06_sar_azimuth.png)

---

## Parte 2 – Elaborazione delle immagini

### i01 · Trasformata di Fourier 2D
Spettro di un'immagine, filtri ideali passa-basso e passa-alto, proprietà di rotazione.

```math
X(f_1,f_2)=\iint x(t_1,t_2)\,e^{-j2\pi(f_1t_1+f_2t_2)}\,dt_1\,dt_2
\qquad
x(\mathbf{A}\mathbf{t})\;\longleftrightarrow\;\frac{1}{|\det\mathbf{A}|}\,X\!\left(\mathbf{A}^{-T}\mathbf{f}\right)
```
Per una rotazione $\mathbf{A}=\mathbf{Q}(\alpha)$ si ha $\mathbf{A}^{-T}=\mathbf{Q}(\alpha)$: lo spettro ruota dello stesso angolo.

![Fourier 2D](python/figures/i01_fourier2d.png)

### i02 · Convoluzione tramite DFT e aliasing
Convoluzione circolare contro lineare (zero padding a N₁+N₂−1, errore ~10⁻¹⁶), aliasing nel sottocampionamento e decimazione con filtro anti-aliasing.

```math
X[k]=\sum_{n=0}^{N-1}x[n]\,e^{-j\frac{2\pi}{N}kn}=X\!\left(e^{j\omega}\right)\Big|_{\omega=2\pi k/N}
\qquad
X[k]\,H[k]\;\longleftrightarrow\;x[n]\circledast_N h[n]
```
```math
M=M_1+M_2-1,\quad N=N_1+N_2-1
\qquad
X\!\left(e^{j\omega}\right)=\frac1T\sum_k X_a\!\left(j\Big(\frac{\omega}{T}-k\frac{2\pi}{T}\Big)\right),\quad \omega=\Omega T
```

![Aliasing](python/figures/i02_aliasing.png)

### i03 · Quantizzazione
Quantizzatore uniforme (**+6.02 dB per bit**, misurato da 2 a 8 bit) e quantizzatore ottimo di **Lloyd–Max** per una variabile gaussiana, implementato in modo iterativo.

```math
\Delta=\frac{2A}{2^b},\qquad D=\frac{\Delta^2}{12}=\frac{\sigma_x^2}{2^{2b}}\;\Rightarrow\;\mathrm{SNR}_q\approx 6.02\,b\ \mathrm{dB}
```
Condizioni di ottimo di Lloyd–Max (risolte iterando):
```math
\theta_k=\frac{q_{k-1}+q_k}{2},\qquad q_k=E\{x\mid \theta_k\le x<\theta_{k+1}\}=\frac{\int_{\theta_k}^{\theta_{k+1}}x\,p_x(x)\,dx}{\int_{\theta_k}^{\theta_{k+1}}p_x(x)\,dx}
```

![Quantizzazione](python/figures/i03_quantization.png)

### i04 · Filtraggio spaziale e bordi
Smoothing, sharpening, gradiente di Sobel, Laplaciano; verifica della risposta in frequenza di una maschera 3×3.

```math
y[n_1,n_2]=\sum_{k_1}\sum_{k_2}x[k_1,k_2]\,h[n_1-k_1,n_2-k_2]
\qquad
\nabla^2 f\approx\begin{bmatrix}0&1&0\\1&-4&1\\0&1&0\end{bmatrix}*f
\qquad
x_{sharp}=x-\alpha\,\nabla^2x
```
```math
h=\begin{bmatrix}0&1&0\\1&1&1\\0&1&0\end{bmatrix}\;\Rightarrow\;H\!\left(e^{j\omega_1},e^{j\omega_2}\right)=1+2\cos\omega_1+2\cos\omega_2
```

![Filtri](python/figures/i04_filtering_edges.png)

### i05 · Restauro di immagini
Sfocatura da movimento più rumore; confronto tra filtro inverso, pseudo-inverso, minimi quadrati vincolati e filtro di Wiener.

Modello di osservazione $y=h*x+n$, cioè $Y=HX+N$ in frequenza:
```math
F_{inv}=\frac{1}{H}
\qquad
F_{CLS}=\frac{H^*}{|H|^2+\lambda\,|C|^2}
\qquad
F_{W}=\frac{H^*\,P_{xx}}{|H|^2P_{xx}+P_{nn}}=\frac{H^*}{|H|^2+\mathrm{NSR}}
```
dove $C$ è un passa-alto (Laplaciano) e $\mathrm{NSR}=P_{nn}/P_{xx}$. Dove $|H|\to 0$ il filtro inverso amplifica il rumore, mentre CLS e Wiener restano limitati.
- PSNR: inverso **4.5 dB** (il rumore viene amplificato), Wiener **30.6 dB**.

![Restauro](python/figures/i05_restoration.png)

### i06 · KLT, DCT e codifica JPEG
KLT di un caso 2×2, verifica che la DCT approssima la KLT di un processo di Markov, codificatore JPEG semplificato (DCT 8×8 e quantizzazione con la matrice standard).

```math
\mathbf{K}_x\,\mathbf{b}^{(k)}=\lambda_k\,\mathbf{b}^{(k)},\qquad X_k=\left\langle \mathbf{b}^{(k)},\mathbf{x}\right\rangle,\qquad \varepsilon_M=\sum_{k=M}^{N-1}\lambda_k
```
```math
b^{(k)}_n=\sqrt{\tfrac{2}{N}}\cos\!\left(\frac{\pi k(2n+1)}{2N}\right),\quad b^{(0)}_n=\tfrac{1}{\sqrt N}
\qquad
\tilde X[k_1,k_2]=\mathrm{round}\!\left(\frac{X[k_1,k_2]}{Q[k_1,k_2]}\right)
```
- **13 %** di coefficienti non nulli con **35.6 dB** di PSNR.

![JPEG](python/figures/i06_jpeg.png)

### i07 · Elaborazione di segnali su grafo
Laplaciano di grafo e trasformata di Fourier su grafo (GFT): un segnale liscio concentra il 97 % dell'energia sulle prime tre "frequenze".

```math
\mathbf{L}=\mathbf{D}-\mathbf{A},\qquad \mathbf{L}\,\mathbf{u}_k=\lambda_k\,\mathbf{u}_k,\qquad \mathbf{S}=\mathbf{U}^T\mathbf{s},\quad \mathbf{s}=\mathbf{U}\mathbf{S}
```
```math
\mathbf{s}^T\mathbf{L}\,\mathbf{s}=\sum_{(i,j)\in E}w_{ij}\,(s_i-s_j)^2
```
Gli autovalori $\lambda_k$ fanno la parte delle frequenze: più la variazione del segnale tra nodi collegati è piccola, più l'energia si concentra sui $\lambda_k$ bassi.

![Grafo](python/figures/i07_graph_fourier.png)

---

## Simulink
`simulink/build_models.m` crea i modelli da codice, usando solo la libreria base:

| Modello | Blocchi principali | Cosa mostra |
|---|---|---|
| `sim_aliasing` | Sine Wave → Zero-Order Hold → To Workspace | un tono a 5 kHz campionato a 8 kHz ricompare a 3 kHz |
| `sim_quantizer` | Sine Wave → Quantizer → Sum | rumore di quantizzazione, SNR ≈ 6.02·b + 1.76 dB |
| `sim_fir_lowpass` | Sine Wave ×2 → Sum → Discrete FIR Filter | filtro passa-basso a finestra di Hamming |
| `sim_mti_canceller` | From Workspace → Unit Delay + Sum, Discrete FIR [1 −3 3 −1] | canceller MTI singolo e triplo |

Costruire i modelli da script, invece che a mano, li rende riproducibili e versionabili con git.

## Scelte implementative
- **Immagini di test sintetiche**: sfumatura, bordi netti e una tessitura periodica, così ogni effetto (ringing, aliasing, perdita di dettaglio) è visibile e lo script non dipende da file esterni.
- **Teoria accanto alla misura**: ogni script stampa entrambi i valori; le differenze sono commentate nel codice. Ad esempio, nel SAR la risoluzione a −3 dB vale 0.886·d_a/2.
- **Seed fissi** per la riproducibilità dei risultati Monte Carlo.

## Prossimi passi
- Validare gli script MATLAB e i modelli Simulink su una release recente di MATLAB.
- Aggiungere un rivelatore CFAR (CA-CFAR) alla mappa range-Doppler.
- Estendere i modelli Simulink con DSP System Toolbox (Buffer, FFT, Spectrum Analyzer) e provare una versione in virgola fissa.
