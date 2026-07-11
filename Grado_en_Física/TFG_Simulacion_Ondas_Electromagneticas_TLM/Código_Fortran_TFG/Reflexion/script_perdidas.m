% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Perdidas'
load salidasperdidas.txt;

filename = 'datosperdidas.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Leer el primer número y asignarlo a la variable 'modo'
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo_modo = 'Excitación +1 sólo en línea 3: Modos TE';
elseif modo == 2
    titulo_modo = 'Excitación +1 sólo en línea 6: Modos TM';
else
    titulo_modo = 'Error';
end
fclose(fileID);

iter = salidasperdidas(:, 1);
t = salidasperdidas(:, 2);
E1 = salidasperdidas(:, 3);
H1 = salidasperdidas(:, 4);
E2 = salidasperdidas(:, 5);
H2 = salidasperdidas(:, 6);
E3 = salidasperdidas(:, 7);
H3 = salidasperdidas(:, 8);

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
    X1f = E1f;
    X2f = E2f;
    X3f = E3f;
    titulo = 'Señal E';
elseif numero == 2
    X1 = H1;
    X2 = H2;
    X3 = H3;
    X1f = H1f;
    X2f = H2f;
    X3f = H3f;
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
plot(iter, X3, 'r-');
hold off;
title([titulo, ' - ', titulo_modo]); % Combinar los títulos
xlabel('Iteración'); % Etiqueta del eje x
ylabel('Amplitud'); % Etiqueta del eje y
legend('Campo 1', 'Campo 2', 'Campo 3'); % Añadir leyenda

figure(2);
clf;
hold on; % Permitir varias gráficas en la misma figura
plot(f(5:min(nmf, length(X1f)))*1e-6, X1f(5:min(nmf, length(X1f))));
plot(f(5:min(nmf, length(X2f)))*1e-6, X2f(5:min(nmf, length(X2f))));
plot(f(5:min(nmf, length(X3f)))*1e-6, X3f(5:min(nmf, length(X3f))));
hold off;
xlabel('Frecuencia (MHz)');
ylabel('mod');
title([titulo, ' - Frecuencia']);
legend('Campo 1', 'Campo 2', 'Campo 3'); % Añadir leyenda

ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ejes(3) = ampmin;
axis(ejes);
