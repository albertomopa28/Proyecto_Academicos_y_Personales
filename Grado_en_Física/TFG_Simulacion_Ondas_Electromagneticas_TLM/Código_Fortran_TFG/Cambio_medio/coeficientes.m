% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Cambio_medio'
load salidas_cambio_medio.txt;

filename = 'datos_cambio_medio.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo_modo = 'Pulso gaussiano';
elseif modo == 2
    titulo_modo = 'Excitación +1 sólo en línea 6: Modos TM';
else
    titulo_modo = 'Error';
end
fclose(fileID);

iter = salidas_cambio_medio(:, 1);
t = salidas_cambio_medio(:, 2);
E1 = salidas_cambio_medio(:, 3);
H1 = salidas_cambio_medio(:, 4);
E2 = salidas_cambio_medio(:, 5);
H2 = salidas_cambio_medio(:, 6);

% Encontrar todos los picos en E1 y H1, y E2 y H2
[pksE1, locsE1] = findpeaks(E1);
[pksE2, locsE2] = findpeaks(E2);

% Encontrar máximos y mínimos en H1
[pksH1_pos, locsH1_pos] = findpeaks(H1);    % Máximos
[pksH1_neg, locsH1_neg] = findpeaks(-H1);   % Mínimos (como valores positivos)
[pksH2_pos, locsH2_pos] = findpeaks(H2);

% Convertir los mínimos a sus valores negativos originales
pksH1_neg = -pksH1_neg;

% Unir los máximos y mínimos en un solo vector de picos de H1
pksH1 = [pksH1_pos; pksH1_neg];
locsH1 = [locsH1_pos; locsH1_neg];

% Ordenar los picos y sus posiciones para mantener el orden en la señal
[pksE1, idxE1] = sort(pksE1, 'descend');
locsE1 = locsE1(idxE1);
[pksE2, idxE2] = sort(pksE2, 'descend');
locsE2 = locsE2(idxE2);
[pksH1, idxH1] = sort(pksH1, 'descend');
locsH1 = locsH1(idxH1);
[pksH2, idxH2] = sort(pksH2_pos, 'descend');
locsH2 = locsH2_pos(idxH2);

% Obtener el valor máximo absoluto de cada señal
E1_max = max(abs(pksE1));
H1_max = max(abs(pksH1));
E2_max = max(abs(pksE2));
H2_max = max(abs(pksH2));

% Filtrar los picos en función del máximo absoluto de cada señal
pksE1_filt = pksE1(abs(pksE1) > 0.01*E1_max);
locsE1_filt = locsE1(abs(pksE1) > 0.01*E1_max);

pksH1_filt = pksH1(abs(pksH1) > 0.01*H1_max);
locsH1_filt = locsH1(abs(pksH1) > 0.01*H1_max);

pksE2_filt = pksE2(abs(pksE2) > 0.01*E2_max);
locsE2_filt = locsE2(abs(pksE2) > 0.01*E2_max);

pksH2_filt = pksH2(abs(pksH2) > 0.01*H2_max);
locsH2_filt = locsH2(abs(pksH2) > 0.01*H2_max);

% Verificar que se tienen suficientes picos válidos para cada variable
if length(pksE1_filt) >= 2 && length(pksH1_filt) >= 2 && length(pksE2_filt) >= 1 && length(pksH2_filt) >= 1
    % Tomar los primeros tres picos válidos de E1 y H1, y el primero de E2 y H2
    E1_peak = pksE1_filt(1); E2_peak = pksE1_filt(2); E3_peak = pksE2_filt(1);
    H1_peak = pksH1_filt(1); H2_peak = pksH1_filt(2); H3_peak = pksH2_filt(1);

    % Calcular coeficientes de transmisión y reflexión
    Gamma_E = E2_peak / E1_peak;
    T_E = E3_peak / E1_peak;
    Gamma_H = H2_peak / H1_peak;
    T_H = H3_peak / H1_peak;
    Z1=E1_peak/H1_peak;
    Z2=E2_peak/H2_peak;
    Z3=E3_peak/H3_peak;
    % Mostrar los resultados en pantalla
    disp(['Gamma_E = ', num2str(Gamma_E)]);
    disp(['T_E = ', num2str(T_E)]);
    disp(['Gamma_H = ', num2str(Gamma_H)]);
    disp(['T_H = ', num2str(T_H)]);
    disp(['Z1 = ', num2str(Z1)]);
    disp(['Z2 = ', num2str(Z2)]);
    disp(['Z3 = ', num2str(Z3)]);
    % Graficar los picos en la señal original
    figure(1);
    clf;

% Gráfica de Campo E
subplot(2, 1, 1);
plot(iter, E1, 'b-');
hold on;
plot(iter, E2, 'g-');
plot(locsE1_filt(1:2), pksE1_filt(1:2), 'ro', 'MarkerFaceColor', 'r'); % Marcar los dos primeros picos válidos de E1
plot(locsE2_filt(1), pksE2_filt(1), 'mo', 'MarkerFaceColor', 'm'); % Marcar el primer pico válido de E2
title('Campo E');
xlabel('Iteración');
ylabel('E (V/m)');
legend('E1', 'Picos de E1', 'E2', 'Pico de E2');
hold off;

% Gráfica de Campo H
subplot(2, 1, 2);
plot(iter, H1, 'b-');
hold on;
plot(iter, H2, 'g-');
plot(locsH1_filt(1:2), pksH1_filt(1:2), 'ro', 'MarkerFaceColor', 'r'); % Marcar los dos primeros picos válidos de H1
plot(locsH2_filt(1), pksH2_filt(1), 'mo', 'MarkerFaceColor', 'm'); % Marcar el primer pico válido de H2
title('Campo H');
xlabel('Iteración');
ylabel('H (A/m)');
legend('H1', 'Picos de H1', 'H2', 'Pico de H2');
hold off;
else
    disp('No se encontraron suficientes picos válidos en E1, H1, E2 o H2 para el cálculo.');
end


