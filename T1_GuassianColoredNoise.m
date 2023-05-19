%
% Comunicaciones Digitales  
% Ruido Gausino de Color
% Autores: Flores Andres, Gomez Franklin, Otavalo David, Zaruma Samantha.
%
%
clc, clear all;
a = 0.9; % parametro para LPF 
L = 50; % Número de muestras utilizadas en el cálculo de covarianza automática
Fs = 1000; % sampling rate
Fc = 10; 
t = 0:1/Fs:2; % vector de tiempo
variance = 1; % Var de white noise

% Generate a dummy signal 
signal = 5*sin(2*pi*Fc*t);

% Gaussian White Noise with media cero and varianza 1
whiteNoise = sqrt(variance) * randn(1, length(signal));

% Calculo auto-correlacion 
[whiteNoiseCorr, lags] = xcorr(whiteNoise, L);

% Calculo auto-covariance 
[whiteNoiseCov, lags] = xcov(whiteNoise, L);

% ruido en el dominio de la frecuencia
NFFT = 2^nextpow2(length(whiteNoise));
whiteNoiseSpectrum = fft(whiteNoise, NFFT) / length(whiteNoise);
f = Fs / 2 * linspace(0, 1, NFFT/2 + 1);

% Ruido Gausino de Color 
x = whiteNoise;

% LPF 1er orden y(n) = a*y(n-1) + (1-a)*x(n)
% Fun. Trans. Fil. Y(Z) = X(Z)*(1-a)/(1-a*Z^-1)
[y, zf] = filter(1 - a, [1, -a], x);
coloredNoise = y;

% auto-correlacion Ruido Gausiano de Color
[coloredNoiseCorr, lags] = xcorr(coloredNoise, L);

% auto-covariance Ruido Gausiano de Color
[coloredNoiseCov, lags] = xcov(coloredNoise, L);

NFFT = 2^nextpow2(length(coloredNoise));
coloredNoiseSpectrum = fft(coloredNoise, NFFT) / length(coloredNoise);
f = Fs / 2 * linspace(0, 1, NFFT/2 + 1);

% Graficas Ruido Blanco
figure(1);
subplot(2, 2, 1);
plot(t, whiteNoise);
title('Ruido Blanco');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 3);
plot(f, 2*abs(whiteNoiseSpectrum(1:NFFT/2 + 1)));
title('PSD Ruido Blanco');
xlabel('Frequency (Hz)');

subplot(2, 2, 4);
stem(f, 2*abs(whiteNoiseSpectrum(1:NFFT/2 + 1)))
title('PSD Ruido Blanco Muestreado');
xlabel('samples');


figure(2);
subplot(2, 2, 1);
plot(lags, whiteNoiseCorr / max(whiteNoiseCorr));
title('Autocorrelacion Ruido Blanco');
xlabel('t');
subplot(2, 2, 2);
stem(lags, whiteNoiseCorr / max(whiteNoiseCorr));
title('Autocorrelacion Ruido Blanco Muestreado');
xlabel('samples');
subplot(2, 2, 3);
plot(lags, whiteNoiseCov / max(whiteNoiseCov));
title('Autocovarianza Ruido Blanco');
xlabel('t');
subplot(2, 2, 4);
stem(lags, whiteNoiseCov / max(whiteNoiseCov));
title('Autocovarianza Ruido Blanco Muestreado');
xlabel('samples');

% Graficas Ruido Guasino de Color
figure(3);
subplot(2, 2, 1);
plot(t, coloredNoise);
title('Ruido Gausiano de color');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2, 2, 3);
plot(f, 2*abs(coloredNoiseSpectrum(1:NFFT/2 + 1)));
title('PSD Ruido Gausiano de Color');
xlabel('Frequency (Hz)');

subplot(2, 2, 4);
stem(f, 2*abs(coloredNoiseSpectrum(1:NFFT/2 + 1)))
title('PSD Ruido Gausiano de Color Muestreado');
xlabel('samples');

figure(4)
subplot(2, 2, 1);
plot(lags, coloredNoiseCorr / max(coloredNoiseCorr));
title(' Autocorrelacion Ruido Gausiano de Color');
xlabel('t');

subplot(2, 2, 2);
stem(lags, coloredNoiseCorr / max(coloredNoiseCorr));
title(' Autocorrelacion Ruido Gausiano de Color Muestreado');
xlabel('samples');

subplot(2, 2, 3);
plot(lags, coloredNoiseCov / max(coloredNoiseCov));
title('AutoCovarianza Ruido Gausiano de Color');
xlabel('t');

subplot(2, 2, 4);
stem(lags, coloredNoiseCov / max(coloredNoiseCov));
title('AutoCovarianza Ruido Gausiano de Color Muestreado');
xlabel('samples');







