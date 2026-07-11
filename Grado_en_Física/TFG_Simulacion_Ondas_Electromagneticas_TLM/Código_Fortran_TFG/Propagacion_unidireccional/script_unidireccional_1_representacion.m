% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Propagacion_unidireccional'
load salidas_unidireccional.txt;

filename = 'datos_unidireccional.txt';
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

iter = salidas_unidireccional(:, 1);
t = salidas_unidireccional(:, 2);
E1 = salidas_unidireccional(:, 3);
H1 = salidas_unidireccional(:, 4);
E2 = salidas_unidireccional(:, 5);
H2 = salidas_unidireccional(:, 6);

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
H1f = abs(fft(H1));
H2f = abs(fft(H2));


% Leer una variable numérica desde el teclado
numero = input('Quiere representar el campo electrico (E) o el magnetico (H) (1/2): ');

% Condicional para comprobar el valor introducido
if numero == 1
    X1 = E1;
    X2 = E2;
    X1f = E1f;
    X2f = E2f;
    titulo = 'Señal E';
elseif numero == 2
    X1 = H1;
    X2 = H2;
    X1f = H1f;
    X2f = H2f;
    titulo = 'Señal H';
else
    disp('Número no reconocido');
    return; % Salir del script si no se reconoce el número
end

% Graficar los campos E o H en subplots
figure(1);
clf; % Nueva figura

plot(iter, X1, 'b-'); 
hold on; % Permitir que otras curvas se añadan
plot(iter, X2, 'g-');
hold off;
title([titulo, ' - ', titulo_modo]); % Combinar los títulos
xlabel('Iteración'); % Etiqueta del eje x
ylabel('Amplitud'); % Etiqueta del eje y
legend('Punto 1', 'Punto 2'); % Añadir leyenda

figure(2);
clf;
hold on; % Permitir varias gráficas en la misma figura
plot(f(5:min(nmf, length(X1f)))*1e-6, X1f(5:min(nmf, length(X1f))));
plot(f(5:min(nmf, length(X2f)))*1e-6, X2f(5:min(nmf, length(X2f))));
hold off;
xlabel('Frecuencia (MHz)');
ylabel('mod');
title([titulo, ' - Frecuencia']);
legend('Campo 1', 'Campo 2'); % Añadir leyenda

ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ejes(3) = ampmin;
axis(ejes);