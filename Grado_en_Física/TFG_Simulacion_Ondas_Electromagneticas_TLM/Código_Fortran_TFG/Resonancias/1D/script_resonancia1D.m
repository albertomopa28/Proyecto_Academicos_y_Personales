% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Resonancias\1D'
load salidas1D.txt;
load salidas2puntos.txt;
filename = 'datos1D.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Leer el modo de operación (TE o TM)
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo = 'Excitación +1 sólo en línea 3: Modos TE';
elseif modo == 2
    titulo = 'Excitación +1 sólo en línea 6: Modos TM';
else
    titulo = 'Error';
end
fclose(fileID);

% Definir los campos y el tiempo
it = salidas1D(:,1);
t = salidas1D(:,2);
ex = salidas1D(:,3);
ey = salidas1D(:,4);
ez = salidas1D(:,5);
hx = salidas1D(:,6);
hy = salidas1D(:,7);
hz = salidas1D(:,8);
E = sqrt(ex.^2 + ey.^2 + ez.^2);
H = sqrt(hx.^2 + hy.^2 + hz.^2);

E1 = salidas2puntos(:,3);
H1 = salidas2puntos(:,4);
E2 = salidas2puntos(:,5);
H2 = salidas2puntos(:,6);

% Parámetros de la señal
n = length(t);
dt = t(2) - t(1);
c = 2.99792458e8;
dl = 2 * c * dt;
disp(dl);

hmin = 10 * dl;
fmaxn = c / hmin;
df = 1 / (n * dt);
f = (0:n-1) * df;
fmaxg = 50 * 10^6;
nmf = round(fmaxg / df);

% FFT de los campos
exf = abs(fft(ex));
eyf = abs(fft(ey));
ezf = abs(fft(ez));
ef = sqrt(exf.^2 + eyf.^2 + ezf.^2);

% Cargar las frecuencias teóricas
load freqs1D.txt
fteor = freqs1D(:,4);

% Primera figura con la FFT
figure(1);
clf;
subplot(211);
plot(f(5:nmf) * 1e-6, ef(5:nmf));
xlabel('Frecuencia (MHz)');
ylabel('mod(E)');
title(titulo);
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;
ampteor = ampmin * ones(size(fteor));
ejes(3) = ampmin;
hold on;
stem(fteor * 1e-6, ampteor);
legend('FFT E', 'Frecuencias teóricas'); % Añadir leyenda aquí
axis(ejes);
hold off;


% Segunda figura con E y H en el tiempo
figure(2);
clf;
subplot(2, 1, 1);
plot(t, E, 'b-');
title('Señal E');
xlabel('Tiempo (s)');
ylabel('Amplitud E');

subplot(2, 1, 2);
plot(t, H, 'r-');
title('Señal H');
xlabel('Tiempo (s)');
ylabel('Amplitud H');

e1f = abs(fft(E1));
e2f = abs(fft(E2));

% Tercera figura con subplots para E y su FFT
figure(3);
clf;

% Subplot 1: Campo E en el tiempo
subplot(2, 1, 1);
plot(it, E1, 'b-');
hold on;
plot(it, E2, 'r-');
title('Campo E');
xlabel('Iteración');
ylabel('E (V/m)');
legend('E1', 'E2'); % Leyenda para los campos E1 y E2
hold off;

% Subplot 2: FFT del campo E
subplot(2, 1, 2);
plot(f(5:nmf) * 1e-6, e2f(5:nmf), 'r-');

hold on;
plot(f(5:nmf) * 1e-6, e1f(5:nmf), 'b-');
title('Transformada de Fourier de E');
xlabel('Frecuencia (MHz)');
ylabel('|E(f)|');
legend('FFT E1', 'FFT E2'); % Leyenda para las FFTs de E1 y E2
hold off;

