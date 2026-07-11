% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Resonancias\3D'
load salidas3D.txt;
filename = 'datos3D.txt';
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

t=salidas3D(:,2);
ex=salidas3D(:,3);
ey=salidas3D(:,4);
ez=salidas3D(:,5);
hx=salidas3D(:,6);
hy=salidas3D(:,7);
hz=salidas3D(:,8);
E=sqrt(ex.^2+ey.^2+ez.^2);
H=sqrt(hx.^2+hy.^2+hz.^2);
n=length(t);
dt=t(2)-t(1);
c=2.99792458e8;
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

load freqs3D.txt
fteor=freqs3D(:,4);
nt=length(fteor);


figure(1);
clf;
plot(f(5:nmf) * 1e-6, ef(5:nmf)); 
xlabel('Frecuencia (MHz)');
ylabel('|E(f)|');
title('Transformada de Fourier de E');
ejes = axis;
ampmax = ejes(4);
ampmin = -0.1 * ampmax;
ampteor = ampmin * ones(size(fteor));
ejes(3) = ampmin;
hold on;
stem(fteor * 1e-6, ampteor, 'g', 'DisplayName', 'Frecuencias teóricas');
axis(ejes);
legend('FFT E', 'Frecuencias teóricas');

% Leyenda para la figura:
% Figura 1: Gráfico de la transformada de Fourier del campo eléctrico E en función de la frecuencia. Se comparan las frecuencias teóricas con las obtenidas por la FFT.


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
title('Señal H');
xlabel('Iteración'); % Etiqueta del eje x
ylabel('Amplitud H'); % Etiqueta del eje y

