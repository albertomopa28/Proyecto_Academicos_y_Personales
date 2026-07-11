% Cambiar al directorio donde se encuentra el archivo 'salidaimpguia.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Guia'
load salidaimpguia.txt;
sal = salidaimpguia;

% Constantes:
eps0 = 8.854185e-12;
mu0 = 1.256637e-6;
z0 = sqrt(mu0 / eps0);
c = 1 ./ sqrt(mu0 * eps0);
fc = 256.583740e6; % frecuencia de corte
disp(c);

% Variables a representar
t = sal(:, 2);
E = sal(:, 3);
H = sal(:, 4);

n = length(t);
dt = t(2) - t(1);
nf = 1000;
E = E(1:nf);
H = H(1:nf);
Erep = sal(:, 3);
Hrep = sal(:, 4);
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

f1 = 256.58; % Frecuencia de corte del modo TE_10
f2 = 513.17; % Frecuencia de corte del modo TE_20

% Gráfica
figure(1);
clf;
plot(f(1:nmfi)* 1e-6, eif(1:nmfi), 'b', 'LineWidth', 1.5); % Señal FFT
hold on;

% Añadir las frecuencias f1 y f2
ymax = max(eif(1:nmfi)); % Obtener el valor máximo para colocar las etiquetas
plot([f1 f1], [0 ymax], 'black--', 'LineWidth', 1); % Línea en f1
plot([f2 f2], [0 ymax], 'black--', 'LineWidth', 1); % Línea en f2

% Etiquetas de las frecuencias
text(f1, ymax * 0.9, 'f_1', 'HorizontalAlignment', 'left', 'Color', 'black', 'FontSize', 25);
text(f2, ymax * 0.9, 'f_2', 'HorizontalAlignment', 'left', 'Color', 'black', 'FontSize', 25);


% Configuración de la gráfica
xlabel('Frecuencia (MHz)', 'FontSize', 20);
ylabel('|E(f)|', 'FontSize', 20);
title('Transformada de Fourier de la alimentación', 'FontSize', 30);
grid on;
hold off;


figure(2);
clf;
plot(ei)




    %Caculamos la impedancia en funcion de la frecuencia y la representamos
    %junto a la impedancia teorica en el rango de frecuencias excitadas:
% Calcular FFT de E y H
Ef = fft(E);
Hf = fft(H);

% Calculamos la función y(x)
y = z0 ./ sqrt(1 - (fc ./ f).^2);

%Rango de frecuencias en el que se representa
fmin=100e6;
fmax = 700e6;
nmfmin=round(fmin / df);
nmfmax = round(fmax / df);
zetan= Ef ./ Hf;


% Primer subplot: Señal E y H
subplot(2, 2, 1);
plot(Erep, 'b-', 'LineWidth', 1.5); % Graficar E
hold on;
plot(Hrep * z0, 'r-', 'LineWidth', 1.5); % Graficar H * z0
title('Señal E y H', 'FontSize', 16);
xlabel('Iteración', 'FontSize', 14);
ylabel('Amplitud (V/m)', 'FontSize', 14);
legend('E', 'Z_0*H', 'FontSize', 12, 'Location', 'best'); % Añadir leyenda
hold off;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Segundo subplot: Módulo de Z
subplot(2, 2, 2);
plot(f(nmfmin:nmfmax) * 1e-6, abs(zetan(nmfmin:nmfmax)), 'blue', 'LineWidth', 1.5);
hold on;
plot(f(nmfmin:nmfmax) * 1e-6, abs(y(nmfmin:nmfmax)), 'r', 'LineWidth', 1.5);
title('Módulo de Z', 'FontSize', 16);
xlabel('Frecuencia (MHz)', 'FontSize', 14);
ylabel('|Z|', 'FontSize', 14);
legend('|Z_{numerico}|', '|Z_{teorico}|', 'FontSize', 12, 'Location', 'best'); % Añadir leyenda
hold off;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Tercer subplot: Parte real de Z
subplot(2, 2, 3);
plot(f(nmfmin:nmfmax) * 1e-6, -real(zetan(nmfmin:nmfmax)), 'blue', 'LineWidth', 1.5);
hold on;
plot(f(nmfmin:nmfmax) * 1e-6, real(y(nmfmin:nmfmax)), 'r', 'LineWidth', 1.5);
title('Parte real de Z', 'FontSize', 16);
xlabel('Frecuencia (MHz)', 'FontSize', 14);
ylabel('Real(Z)', 'FontSize', 14);
legend('Real(Z_{numerico})', 'Real(Z_{teorico})', 'FontSize', 12, 'Location', 'best'); % Añadir leyenda
hold off;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);

% Cuarto subplot: Parte imaginaria de Z
subplot(2, 2, 4);
plot(f(nmfmin:nmfmax) * 1e-6, imag(zetan(nmfmin:nmfmax)), 'blue', 'LineWidth', 1.5);
hold on;
plot(f(nmfmin:nmfmax) * 1e-6, imag(y(nmfmin:nmfmax)), 'r', 'LineWidth', 1.5);
title('Parte imaginaria de Z', 'FontSize', 16);
xlabel('Frecuencia (MHz)', 'FontSize', 14);
ylabel('Imag(Z)', 'FontSize', 14);
legend('Imag(Z_{numerico})', 'Imag(Z_{teorico})', 'FontSize', 12, 'Location', 'best'); % Añadir leyenda
hold off;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);




% Crear una nueva figura para las transformadas de Fourier de E y H
figure(8);
clf;

% Graficar |Ef| y |Hf| en función de la frecuencia
plot(f(1:nf/2) * 1e-6, abs(Ef(1:nf/2)), 'b-', 'LineWidth', 1.5); % Módulo de Ef
hold on;
plot(f(1:nf/2) * 1e-6, abs(z0*Hf(1:nf/2)), 'r-', 'LineWidth', 1.5); % Módulo de Hf
title('Transformadas de Fourier de E y H', 'FontSize', 16);
xlabel('Frecuencia (MHz)', 'FontSize', 14);
ylabel('Amplitud', 'FontSize', 14);
legend('|Ef|', '|Hf|', 'FontSize', 12, 'Location', 'best'); % Añadir leyenda
hold off;
set(gca, 'LineWidth', 1.5, 'FontSize', 12);



