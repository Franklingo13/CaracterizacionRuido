clc, clear;
% Parámetros de la simulación
fs = 1000;              % Frecuencia de muestreo (Hz)
dur = 1;                % Duración de la señal (segundos)
fLow = 10;              % Frecuencia de corte baja (Hz)
fHigh = 100;            % Frecuencia de corte alta (Hz)

% Generación de ruido blanco
numSamples = fs * dur;
whiteNoise = randn(1, numSamples);

% Filtro pasabajos
[b, a] = butter(2, fHigh / (fs/2), 'low');
filteredNoise = filter(b, a, whiteNoise);

% Filtro pasaaltos
[b, a] = butter(2, fLow / (fs/2), 'high');
coloredNoise = filter(b, a, filteredNoise);

% Gráfico del ruido Gaussiano de color
figure(1)
time = (0:numSamples-1) / fs;
plot(time, coloredNoise);
xlabel('Tiempo (s)');
ylabel('Amplitud');
title('Ruido Gaussiano de color');

% Reproducir el sonido del ruido Gaussiano de color
%soundsc(coloredNoise, fs);

% Cálculo de la transformada de Fourier del ruido de color
fftColoredNoise = fftshift(fft(coloredNoise));
freq = (-fs/2 : fs/numSamples : fs/2 - fs/numSamples);

% Gráfico en frecuencia del ruido de color
figure(2)
plot(freq, abs(fftColoredNoise));
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro de frecuencia del ruido de color');

% Cálculo de la transformada de Fourier del ruido Gaussiano
fftGaussianNoise = fftshift(fft(whiteNoise));
freq_g = (-fs/2 : fs/numSamples : fs/2 - fs/numSamples);

% Gráfico en frecuencia del ruido de color
figure(3)
plot(freq_g, abs(fftGaussianNoise));
xlabel('Frecuencia (Hz)');
ylabel('Amplitud');
title('Espectro de frecuencia del ruido Gaussiano');
