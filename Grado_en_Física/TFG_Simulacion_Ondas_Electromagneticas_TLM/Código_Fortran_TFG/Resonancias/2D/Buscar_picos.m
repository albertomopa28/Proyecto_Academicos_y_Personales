% Cambiar al directorio donde se encuentra el archivo 'salidas1D.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Resonancias\2D'
load salidas2DTM.txt;

% Cargar los datos
t = salidas2DTM(:,2);
ex = salidas2DTM(:,3);
ey = salidas2DTM(:,4);
ez = salidas2DTM(:,5);

% Calcular el campo eléctrico total E
E = sqrt(ex.^2 + ey.^2 + ez.^2);

% Parámetros de la señal y frecuencias
n = length(t);
dt = t(2) - t(1);
c = 2.99792458e8; % Velocidad de la luz
df = 1 / (n * dt); % Resolución en frecuencia
disp(df);
f = (0:n-1) * df; % Vector de frecuencias

% Limitar frecuencias entre 1 MHz y 50 MHz
f_min = 10e6; % 1 MHz
f_max = 300e6; % 50 MHz
idx_min = find(f >= f_min, 1);
idx_max = find(f <= f_max, 1, 'last');
frecuencias = f(idx_min:idx_max) * 1e-6; % Frecuencias en MHz

% Transformada de Fourier del campo eléctrico
exf = abs(fft(ex));
eyf = abs(fft(ey));
ezf = abs(fft(ez));
ef = sqrt(exf.^2 + eyf.^2 + ezf.^2);
ef = ef(idx_min:idx_max); % Filtrar frecuencias en el rango deseado

% Asegurar que las dimensiones de 'frecuencias' y 'ef' coincidan
if length(frecuencias) ~= length(ef)
    frecuencias = linspace(f_min * 1e-6, f_max * 1e-6, length(ef));
end

% Ajustar la sensibilidad a los picos con 'MinPeakProminence'
% Valor bajo para detectar más picos, ajustar según los datos
[picos, locs] = findpeaks(ef, 'NPeaks', 6, 'MinPeakProminence', 0.01);

% Obtener las frecuencias correspondientes a los picos
frecuencias_picos = frecuencias(locs);

% Guardar las frecuencias de los picos en un archivo .txt
filename = 'frecuencias_picos.txt';
fileID = fopen(filename, 'w');
fprintf(fileID, 'Frecuencias de los 6 primeros picos (MHz):\n');
fprintf(fileID, '%f\n', frecuencias_picos);
fclose(fileID);

% Guardar las frecuencias y amplitudes de los picos en 'picos.txt'
filename_picos = 'picos.txt';
fileID_picos = fopen(filename_picos, 'w');
fprintf(fileID_picos, 'Frecuencia (MHz) - Amplitud:\n');
for i = 1:length(frecuencias_picos)
    fprintf(fileID_picos, '%f MHz - %f\n', frecuencias_picos(i), picos(i));
end
fclose(fileID_picos);

% Graficar el campo eléctrico en función del tiempo
figure(1);
clf;
plot(t, E, 'b-');
title('Campo Eléctrico E en función del tiempo');
xlabel('Tiempo (s)');
ylabel('Amplitud de E');

% Graficar la transformada de Fourier del campo eléctrico y marcar los picos
figure(2);
clf;
plot(frecuencias, ef, 'r-');
hold on;
plot(frecuencias(locs), picos, 'bo', 'MarkerFaceColor', 'b'); % Marcar los picos
title('Transformada de Fourier de E con picos');
xlabel('Frecuencia (MHz)');
ylabel('Amplitud |E(f)|');
legend('FFT de E', 'Picos');
hold off;

