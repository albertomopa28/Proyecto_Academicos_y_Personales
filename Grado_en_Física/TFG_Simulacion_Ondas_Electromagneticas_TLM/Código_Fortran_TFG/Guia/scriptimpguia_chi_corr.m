% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Guia'
load salidaimpguia.txt;
sal = salidaimpguia;

% Constantes
eps0 = 8.854185e-12;
mu0 = 1.256637e-6;
z0 = sqrt(mu0 / eps0);
c = 1 / sqrt(mu0 * eps0);
fc = 256.583740e6; % frecuencia de corte
disp(c);

% Variables a representar
t = sal(:, 2);
E = sal(:, 3);
H = sal(:, 4);

n = length(t);
dt = t(2) - t(1);
nf = 7000;
E = E(1:nf);
H = H(1:nf);
df = 1 / (nf * dt);
f = (0:nf - 1) * df;

% Representamos la alimentación y su transformada de Fourier
load alim.txt;
ei = alim;
eif = abs(fft(ei));
nei = length(eif);
dfi = 1 / (nei * dt);
fi = (0:nei - 1) * dfi;
fmaxi = 700e6;
nmfi = round(fmaxi / dfi);

figure(1);
clf;
plot(f(1:nmfi), eif(1:nmfi));

figure(2);
clf;
plot(ei);

% Graficar los campos E y H 
figure(3);
clf;
plot(E, 'b-'); % Graficar E
hold on;
plot(H * z0, 'r-'); % Graficar H * z0
title('Señal E y H');
xlabel('Iteración');
ylabel('Amplitud E y H');
hold off;

% Calcular la impedancia en función de la frecuencia
Ef = fft(E);
Hf = fft(H);
y = z0 ./ sqrt(1 - (fc ./ f(1:nf)).^2);

% Rango de frecuencias para representar
fmin = 100e6;
fmax = 700e6;
nmfmin = round(fmin / df);
nmfmax = min(round(fmax / df), length(f)); % Para evitar dimensiones fuera del límite
zetan = Ef ./ Hf;

% Verificar y ajustar dimensiones
nmf_range = nmfmin:nmfmax;
zetan_range = zetan(nmf_range);
y_range = y(nmf_range);

% Asegurarse de que las dimensiones de zetan_range y y_range coincidan
min_length = min(length(zetan_range), length(y_range));
zetan_range = zetan_range(1:min_length);
y_range = y_range(1:min_length);

% Figura para el módulo de zetan
figure(4);
clf;
plot(f(nmf_range(1:min_length)) * 1e-6, abs(zetan_range), 'black');
title('Módulo de Z');
xlabel('Frecuencia (MHz)');
ylabel('|Z|');
hold on;
plot(f(nmf_range(1:min_length)) * 1e-6, abs(y_range), 'r');
hold off;

% Figura para la parte real de zetan
figure(5);
plot(f(nmf_range(1:min_length)) * 1e-6, -real(zetan_range), 'black');
title('Parte real de Z');
xlabel('Frecuencia (MHz)');
ylabel('Real(Z)');
hold on;
plot(f(nmf_range(1:min_length)) * 1e-6, real(y_range), 'r');
hold off;

% Figura para la parte imaginaria de zetan
figure(6);
clf;
plot(f(nmf_range(1:min_length)) * 1e-6, imag(zetan_range), 'g');
title('Parte imaginaria de Z');
xlabel('Frecuencia (MHz)');
ylabel('Imag(Z)');
hold on;
plot(f(nmf_range(1:min_length)) * 1e-6, imag(y_range), 'r');
hold off;

% Rango de frecuencias para representar y calcular estadísticas
fmin = 260e6;
fmax = 513.167480e6;
nmfmin = round(fmin / df);
nmfmax = min(round(fmax / df), length(f)); % Para evitar dimensiones fuera del límite
zetan = Ef ./ Hf;

% Verificar y ajustar dimensiones en el nuevo rango
nmf_range = nmfmin:nmfmax;
zetan_range = zetan(nmf_range);
y_range = y(nmf_range);

% Asegurarse de que las dimensiones de zetan_range y y_range coincidan
min_length = min(length(zetan_range), length(y_range));
zetan_range = zetan_range(1:min_length);
y_range = y_range(1:min_length);

% Actualizar N al nuevo rango de frecuencias
N = min_length; % Número de datos en el rango de frecuencias

% Coeficiente de correlación r para el módulo
Z_mean_abs = mean(abs(zetan_range));
y_mean_abs = mean(abs(y_range));
r_abs = sum((abs(zetan_range) - Z_mean_abs) .* (abs(y_range') - y_mean_abs)) / ...
    sqrt(sum((abs(zetan_range) - Z_mean_abs).^2) * sum((abs(y_range') - y_mean_abs).^2));

% Coeficiente de correlación r para la parte real usando -real(zetan)
Z_mean_real = mean(-real(zetan_range));
y_mean_real = mean(real(y_range));
r_real = sum((-real(zetan_range) - Z_mean_real) .* (real(y_range') - y_mean_real)) / ...
    sqrt(sum((-real(zetan_range) - Z_mean_real).^2) * sum((real(y_range') - y_mean_real).^2));


% Coeficiente de correlación r para la parte imaginaria usando Spearman
% Coeficiente de correlación r para la parte real usando imag(zetan)
Z_mean_imag = mean(-imag(zetan_range));
y_mean_imag = mean(imag(y_range));

r_imag = N * sum(imag(zetan_range) .* imag(y_range')) - sum(imag(zetan_range)) * sum(imag(y_range'))/ ...
sqrt((N *sum(imag(zetan_range).^2)-sum(imag(zetan_range))^2)*(N*sum(imag(y_range').^2)-sum(imag(y_range'))^2));

r_imag = sum((-imag(zetan_range) - Z_mean_imag) .* (imag(y_range') - y_mean_imag)) / ...
    sqrt(sum((-imag(zetan_range) - Z_mean_imag).^2) * sum((imag(y_range') - y_mean_imag).^2));



% Mostrar resultados
disp(['Coeficiente de correlación r (módulo): ', num2str(r_abs)]);
disp(['Coeficiente de correlación r (parte real): ', num2str(r_real)]);
disp(['Coeficiente de correlación r (parte imaginaria): ', num2str(r_imag)]);



