% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Dispersion'
load salidasdispersion.txt;

filename = 'datosdispersion.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
modo = fscanf(fileID, '%d', 1);
fclose(fileID);

iter = salidasdispersion(:, 1);
t = salidasdispersion(:, 2);
E1 = salidasdispersion(:, 3);
H1 = salidasdispersion(:, 4);
E2 = salidasdispersion(:, 5);
H2 = salidasdispersion(:, 6);
E3 = salidasdispersion(:, 7);
H3 = salidasdispersion(:, 8);
E4 = salidasdispersion(:, 9);
H4 = salidasdispersion(:, 10);

% Seleccionar datos hasta la iteración 400
iter_limit = 600;
idx_limit = find(iter <= iter_limit); % Índice de las iteraciones hasta 400
iter = iter(idx_limit);
E1 = E1(idx_limit);
H1 = H1(idx_limit);
E2 = E2(idx_limit);
H2 = H2(idx_limit);
E3 = E3(idx_limit);
H3 = H3(idx_limit);

n = length(t);
dt = t(2) - t(1);
c = 2.99792458e8;
dl = 2 * c * dt;
disp(dl);

hmin = 10 * dl;
fmaxn = c / hmin;
df = 1 / (n * dt);
f = (0:n-1) * df;
fmaxg = fmaxn;
nmf = round(fmaxg / df);

E1f = abs(fft(E1));
E2f = abs(fft(E2));
E3f = abs(fft(E3));
H1f = abs(fft(H1));
H2f = abs(fft(H2));
H3f = abs(fft(H3));

% Leer una variable numérica desde el teclado
numero = input('Quiere representar el campo electrico (E) o el magnetico (H) (1/2): ');

% Condicional para comprobar el valor introducido
if numero == 1
    X1 = E1;
    X2 = E2;
    X3 = E3;
    titulo = 'Campo E';
elseif numero == 2
    X1 = H1;
    X2 = H2;
    X3 = H3;
    titulo = 'Señal H';
else
    disp('Número no reconocido');
    return; % Salir del script si no se reconoce el número
end

% Graficar los campos E o H en subplots
figure(1);
clf;

plot(iter, X1, 'b-', 'LineWidth', 1.25);
hold on;
plot(iter, X2, 'g-', 'LineWidth', 1.25);
plot(iter, X3, 'r-', 'LineWidth', 1.25);
hold off;

% Títulos y etiquetas con fuente más grande
title(titulo, 'FontSize', 20);
xlabel('Iteración', 'FontSize', 18);
ylabel('E (V/m)', 'FontSize', 18);

% Leyenda con fuente más grande
legend({'E1', 'E2', 'E3'}, 'FontSize', 16, 'Location', 'best');

% Personalización del gráfico
set(gca, 'LineWidth', 1.5, 'FontSize', 16);


