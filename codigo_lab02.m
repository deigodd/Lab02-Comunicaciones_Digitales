clc; clear; close all;

% Parametros
f0 = 1000;             % Frecuencia fundamental (Hz)
t = linspace(-0.01, 0.01, 1000);  % Tiempo en segundos
f = linspace(-2*f0, 2*f0, 1000);  % Frecuencia en Hz
alphas = [0, 0.25, 0.75, 1];      % Valores de roll-off
colors = ['r', 'g', 'b', 'k'];    % Colores para distinguir

%% --- Graficos de frecuencia e impulso ---
figure;
% Respuesta en frecuencia
subplot(2,1,1);
hold on;
for i = 1:length(alphas)
    alpha = alphas(i);
    f_delta = alpha * f0;
    f1 = f0 - f_delta;
    B = f0 + f_delta;
    He = zeros(size(f));
    for k = 1:length(f)
        fk = abs(f(k));
        if fk < f1
            He(k) = 1;
        elseif fk >= f1 && fk < B
            He(k) = 0.5 *(1 + cos(pi * (fk - f1) / (2 * f_delta)));
        else
            He(k) = 0;
        end
    end
    plot(f, He, 'Color', colors(i), 'LineWidth', 2, 'DisplayName', ['\alpha = ' num2str(alpha)]);
end

% -GRAFICOS DE FRECUENCIA- 
title('Respuesta en Frecuencia para distintos valores de \alpha');
xlabel('Frecuencia (Hz)');
ylabel('Magnitud');
legend('show'); grid on;

% Respuesta al impulso
subplot(2,1,2);
hold on;
impulsos = cell(1,length(alphas));  % Guardamos respuestas para usarlas luego
for i = 1:length(alphas)
    alpha = alphas(i);
    f_delta = alpha * f0;
    he = zeros(size(t));
    for k = 1:length(t)
        if t(k) == 0
            he(k) = 2 * f0;
        elseif abs(4 * f_delta * t(k)) == 1
            he(k) = pi / 2 * sin(2 * pi * f0 * t(k));
        else
            numerator = sin(2 * pi * f0 * t(k));
            denominator = 2 * pi * f0 * t(k);
            cos_part = cos(2 * pi * f_delta * t(k));
            sinc_part = numerator / denominator;
            denominator2 = 1 - (4 * f_delta * t(k))^2;
            he(k) = 2 * f0 * sinc_part * cos_part / denominator2;
        end
    end
    impulsos{i} = he;  
    plot(t, he, 'Color', colors(i), 'LineWidth', 2, 'DisplayName', ['\alpha = ' num2str(alpha)]);
end

% -GRAFICOS DE IMPULTOS- 
title('Respuesta al Impulso para distintos valores de \alpha');
xlabel('Tiempo (s)');
ylabel('Amplitud');
legend('show'); grid on;

%---------------------------------------------------------
% Parametros comunes - Diagrama de ojo
N_bits = 10000;                 % Numero de bits
bit_rate = 1e3;                 % Tasa de bits (bps)
fs = 10*bit_rate;               % Frecuencia de muestreo
alpha_values = [0, 0.25, 0.75, 1];
colors = ['r', 'g', 'b', 'k'];  % Colores para cada ALPHA
span = 10;                      % Duracion del filtro
SNR = 20;                       % Relacion signal-ruido (dB)

% Genera secuencia NRZ-L con ruido
bits = randi([0 1], 1, N_bits);
nrz_signal = 2*bits - 1;        % Mapeo: 0 -> -1, 1 -> 1
upsampled_signal = repelem(nrz_signal, fs/bit_rate);
noisy_signal = awgn(upsampled_signal, SNR, 'measured');
% Itera sobre cada ALPHA y crea una figura para cada uno
% Itera sobre cada ALPHA y crea una figura para cada uno
for idx = 1:length(alpha_values)
    figure;  % Crea una nueva figura para cada ALPHA
    hold on;
    alpha = alpha_values(idx);
    sps = fs/bit_rate;  % Muestras por simbolo
    % Filtro de coseno alzado
    filter_coeff = rcosdesign(alpha, span, sps, 'sqrt');
    % Aplica filtro
    filtered_signal = conv(noisy_signal, filter_coeff, 'same');
    % Extrae segmentos de 2 simbolos para el diagrama de ojo
    segment_length = 2*sps;
    num_segments = floor(length(filtered_signal)/segment_length);

    for k = 1:min(50, num_segments)
        start_idx = (k-1)*segment_length + 1;
        end_idx = k*segment_length;
        segment = filtered_signal(start_idx:end_idx);

        % El tiempo debe cubrir 2 simbolos (cada simbolo dura 1/bit_rate segundos)
        time = linspace(0, 2/bit_rate, segment_length) * 1e3;  % en milisegundos

        % Graficar cada segmento
        plot(time, segment, 'Color', colors(idx),'LineWidth', 0.5);
    end

    % Etiquetas y formato
    title(['Diagrama de Ojo (ALPHA = ', num2str(alpha), ')']);
    xlabel('Tiempo (ms)');  % Etiqueta en milisegundos
    ylabel('Amplitud');
    xlim([0, 2/bit_rate * 1e3]);  % 0 a 2/bit_rate en milisegundos
    grid on;
    hold off;
end
