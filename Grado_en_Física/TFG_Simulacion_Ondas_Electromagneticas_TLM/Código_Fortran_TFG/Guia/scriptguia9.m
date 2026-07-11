% Cambiar al directorio donde se encuentra el archivo 'salidas.txt'
cd 'C:\Users\usuario\Documents\Cosas_de_Albero\4º_Carrera\TFG\programas\VscodeFortran\Guia'
load salidasguia.txt;
filename = 'datosguia.txt';
fileID = fopen(filename, 'r'); % Abrir el archivo en modo de lectura
% Para saber en que modo estamos leemos el primer n?mero y asignarlo a la variable 'modo'
modo = fscanf(fileID, '%d', 1);
if modo == 1
    titulo = 'Excitacion +1 solo en linea 3: Modos TE';
elseif modo == 2
    titulo = 'Excitacion +1 solo en linea 6: Modos TM';
else
    titulo = 'Excitacion solo en linea 3';
end
fclose(fileID);

t=salidasguia(:,2);
E1=salidasguia(:,3);
E2=salidasguia(:,4);
E3=salidasguia(:,5);
E4=salidasguia(:,6);
H1=salidasguia(:,7);
H2=salidasguia(:,8);
H3=salidasguia(:,9);
H4=salidasguia(:,10);
n=length(t);
dt=t(2)-t(1);
c=2.99792458e8;
dl=2*c*dt;
disp(dl);

hmin=10*dl;
fmaxn=c/hmin;
df=1./(n*dt);
f=(0:n-1)*df;
fmaxg=fmaxn;
nmf=round(fmaxg/df);




% Graficar los campos E y H en subplots
figure(2);
clf; % Nueva figura

% Primer subplot para la se?al E
plot(t, E1, 'b-'); % Graficar 
hold
plot(t, E2, 'r-')
plot(t, E3, 'c-')
plot(t, E4, 'g-')

title('Se?al E');
xlabel('Iteraci?n'); % Etiqueta del eje x
ylabel('Amplitud E'); % Etiqueta del eje y
