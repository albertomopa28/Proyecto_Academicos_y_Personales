% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Resonancias\2D'
load salidas2D.txt;
load salidas2DTM.txt;
filename = 'datos2D.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Para saber en que modo estamos leemos el primer número y asignarlo a la variable 'modo'
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo = 'Excitación +1 sólo en línea 3: Modos TE';
elseif modo == 2
    titulo = 'Excitación +1 sólo en línea 6: Modos TM';
else
    titulo = 'Error';
end
fclose(fileID);

%Datos modo TE
tTM=salidas2DTM(:,2);
exTM=salidas2DTM(:,3);
eyTM=salidas2DTM(:,4);
ezTM=salidas2DTM(:,5);
hxTM=salidas2DTM(:,6);
hyTM=salidas2DTM(:,7);
hzTM=salidas2DTM(:,8);
ETM=sqrt(exTM.^2+eyTM.^2+ezTM.^2);
HTM=sqrt(hxTM.^2+hyTM.^2+hzTM.^2);

%Datos modo TM
t=salidas2D(:,2);
ex=salidas2D(:,3);
ey=salidas2D(:,4);
ez=salidas2D(:,5);
hx=salidas2D(:,6);
hy=salidas2D(:,7);
hz=salidas2D(:,8);
E=sqrt(ex.^2+ey.^2+ez.^2);
H=sqrt(hx.^2+hy.^2+hz.^2);

%Variables extraibles de los datos y constantes:
c=2.99792458e8;
n=length(t);
dt=t(2)-t(1);
dl=2*c*dt;
hmin=10*dl;
fmaxn=c/hmin;
df=1./(n*dt);
disp(df);
f=(0:n-1)*df;
fmaxg=fmaxn;
nmf=round(fmaxg/df);

exf=abs(fft(ex));
eyf=abs(fft(ey));
ezf=abs(fft(ez));
hxf=abs(fft(hx));
hyf=abs(fft(hy));
hzf=abs(fft(hz));
ef=sqrt(exf.^2+eyf.^2+ezf.^2);
hf=sqrt(hxf.^2+hyf.^2+hzf.^2);

load freqs2D.txt
fteor=freqs2D(:,4);
nt=length(fteor);

load freqs2DTM.txt
fteorTM=freqs2DTM(:,4);
ntTM=length(fteorTM);

    %Representamos la FFT
figure(1);
clf;
subplot(211);
clear axis;
plot(f(5:nmf)*1e-6,ef(5:nmf));
xlabel('Frecuencia (MHz)');
ylabel('mod(E)')
title(titulo);

ejes=axis;
ampmax=ejes(4);
ampmin=-0.1*ampmax;

ampteor=ampmin*ones(size(fteor));
ejes(3)=ampmin;
hold;
stem(fteor*1e-6,ampteor)
axis(ejes);

subplot(212);

plot(f(5:nmf)*1e-6,hf(5:nmf));
xlabel('Frecuencia (MHz)');
ylabel('mod(H)')
title(titulo);
ejes=axis;
ampmax=ejes(4);
ampmin=-0.1*ampmax;
fmaxima=ejes(2);

ampteor=ampmin*ones(size(fteor));
ejes(3)=ampmin;
hold;
clear axis;
stem(fteor*1e-6,ampteor)
axis(ejes);



    % Graficar los campos E y H en subplots
figure(2);
clf; % Nueva figura

% Primer subplot para la señal E
subplot(2, 1, 1);
plot(t, E, 'b-'); % Graficar E
title('Señal E');
xlabel('Iteración'); % Etiqueta del eje x
ylabel('Amplitud E'); % Etiqueta del eje y

% Segundo subplot para la señal H
subplot(2, 1, 2);
plot(t, H, 'r-'); % Graficar H
title('Señal H');% Etiqueta del eje x
ylabel('Amplitud H')
xlabel('Iteración');  % Etiqueta del eje y




    %FFT modo TE y TM
eTMf=abs(fft(ezTM));
figure(3);
clf;

% Subplot 1: Modos TM
subplot(211);
clear axis;
plot(f(5:nmf)*1e-6, eTMf(5:nmf), 'b-'); % Transformada de Fourier en azul
xlabel('Frecuencia (MHz)');
ylabel('|E(f)|');
title('Transformada de Fourier de E, Modos TM');

ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ampteor = ampmin * ones(size(fteorTM));
ejes(3) = ampmin;
hold on;
stem(fteorTM * 1e-6, ampteor, 'g', 'DisplayName', 'Frecuencias teóricas'); % Frecuencias teóricas en verde
legend('FFT E (Modos TM)', 'Frecuencias teóricas (Modos TM)');
axis(ejes);
hold off;

% Subplot 2: Modos TE
subplot(212);
plot(f(5:nmf) * 1e-6, ef(5:nmf), 'r-'); % Transformada de Fourier en rojo
xlabel('Frecuencia (MHz)');
ylabel('|E(f)|');
title('Transformada de Fourier de E, Modos TE');
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;

ampteor = ampmin * ones(size(fteor));
ejes(3) = ampmin;
hold on;
stem(fteor * 1e-6, ampteor, 'g', 'DisplayName', 'Frecuencias teóricas'); % Frecuencias teóricas en verde
legend('FFT E (Modos TM)', 'Frecuencias teóricas (Modos TM)');
axis(ejes);
hold off;

