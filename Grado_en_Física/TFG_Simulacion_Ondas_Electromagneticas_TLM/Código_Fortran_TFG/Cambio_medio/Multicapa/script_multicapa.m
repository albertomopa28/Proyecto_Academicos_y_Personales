% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Cambio_medio\Multicapa'
load salidas_multicapa.txt;
load fteor.txt;
filename = 'datos_multicapa.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Leer el primer número y asignarlo a la variable 'modo'
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo_modo = 'Pulso gaussiano';
elseif modo == 2
    titulo_modo = 'Excitación +1 sólo en línea 6: Modos TM';
else
    titulo_modo = 'Error';
end
fclose(fileID);

iter = salidas_multicapa(:, 1);
t = salidas_multicapa(:, 2);
E1 = salidas_multicapa(:, 3);
H1 = salidas_multicapa(:, 4);
E2 = salidas_multicapa(:, 5);
H2 = salidas_multicapa(:, 6);
E3 = salidas_multicapa(:, 7);
H3 = salidas_multicapa(:, 8);
E4 = salidas_multicapa(:, 9);
H4 = salidas_multicapa(:, 10);

n = length(t);
dt = t(2) - t(1);
c = 2.99792458e8;
dl = 2 * c * dt;

hmin = 10 * dl;
fmaxn = 1.01*2.49999999E+11;
df = 1 / (n * dt);
disp(df);
f = (0:n-1) * df;
fmaxg = fmaxn;
nmf = round(fmaxg / df);

E1f = abs(fft(E1));
E2f = abs(fft(E2));
E3f = abs(fft(E3));
E4f = abs(fft(E4));
H1f = abs(fft(H1));
H2f = abs(fft(H2));
H3f = abs(fft(H3));
H4f = abs(fft(H4));

%Para las frecuencias:
fteor=fteor(:,1);
nt=length(fteor);

figure(1); % Una sola figura
clf;


hold on;
plot(f(2:min(nmf, length(E2f))) * 1e-9, E2f(2:min(nmf, length(E2f))), 'LineWidth', 2);
title('Transformada de Fourier de E2', 'FontSize', 16);
xlabel('Frecuencia (GHz)', 'FontSize', 14);
ylabel('|E(f)|', 'FontSize', 14);
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;
ampteor = ampmin * ones(size(fteor));
ejes(3) = ampmin;
stem(fteor * 1e-9, ampteor, 'r', 'LineWidth', 1.5, 'MarkerSize', 8, 'DisplayName', 'Frecuencias teóricas'); % Cambiado a rojo y ajustado grosor y tamaño
axis(ejes);
legend({'FFT E2', 'Frecuencias teóricas'}, 'FontSize', 12, 'Location', 'best');
set(gca, 'LineWidth', 1.5, 'FontSize', 12); % Cuadro de la gráfica
hold off;

% Encontrar los primeros cinco mínimos de E2f
[~, locs] = findpeaks(-E2f); % Encuentra picos negativos, que corresponden a mínimos
locs = locs(1:5); % Selecciona los primeros cinco mínimos
frecuencias_minimos = f(locs) * 1e-9; % Convierte las frecuencias a GHz

% Guardar los primeros cinco mínimos en un archivo .txt
fileID = fopen('primeros_minimos.txt', 'w');
fprintf(fileID, 'Primeros cinco mínimos de la transformada de Fourier de E2 (en GHz):\n');
fprintf(fileID, '%f\n', frecuencias_minimos);
fclose(fileID);

% Crear figura para visualizar la transformada de Fourier con los mínimos
figure(2);
clf;
plot(f(2:min(nmf, length(E2f))) * 1e-9, E2f(2:min(nmf, length(E2f))), 'LineWidth', 2);
hold on;
plot(f(locs) * 1e-9, E2f(locs), 'ro', 'MarkerSize', 8, 'DisplayName', 'Mínimos encontrados'); % Marcar los mínimos en rojo
title('Transformada de Fourier de E2 con mínimos resaltados', 'FontSize', 16);
xlabel('Frecuencia (GHz)', 'FontSize', 14);
ylabel('|E(f)|', 'FontSize', 14);
legend({'FFT E2', 'Mínimos encontrados'}, 'FontSize', 12, 'Location', 'best');
set(gca, 'LineWidth', 1.5, 'FontSize', 12);
hold off;

% Segundo subplot para la señal H
%subplot(2, 1, 2);
%hold on; % Permitir varias gráficas en la misma figura
%plot(f(5:min(nmf, length(H1f)))*1e-9, H1f(5:min(nmf, length(H1f))));
%plot(f(5:min(nmf, length(H2f)))*1e-9, H2f(5:min(nmf, length(H2f))));
%plot(f(5:min(nmf, length(H3f)))*1e-9, H3f(5:min(nmf, length(H3f))));
%plot(f(5:min(nmf, length(H4f)))*1e-9, H4f(5:min(nmf, length(H4f))));
%hold off;
%title('FFT Señal H');
%xlabel('Frecuencia (GHz)');
%ylabel('mod(H)')
%ejes = axis;
%ampmax = ejes(4);
%ampmin = -0.1 * ampmax;

%ejes(3) = ampmin;
%axis(ejes);

