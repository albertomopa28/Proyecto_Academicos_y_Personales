% Programa para graficar la función y = c / sqrt(1 - (f/x)^2)

% Definir constantes
eps0=8.854185e-12;
mu0=1.256637e-6;
z0 = sqrt(mu0/eps0);
c=1./sqrt(mu0*eps0);

f = 256.583740e6;  % frecuencia (puedes cambiarlo si es necesario, en Hz)
x = linspace(0.01, 4*f, 1000); % Se evita x = 0 para prevenir división por cero

% Calcular la función y(x)
y = z0 ./ sqrt(1 - (f ./ x).^2);

% Graficar el módulo de y(x)
figure(1);
plot(x, abs(y), 'b', 'LineWidth', 2);
xlabel('x (variable independiente)');
ylabel('|y(x)|');
title('Gráfica del módulo de y = c / sqrt(1 - (f/x)^2)');
grid on;

% Graficar la parte real de y(x)
figure(2);
plot(x, real(y), 'g', 'LineWidth', 2);
xlabel('x (variable independiente)');
ylabel('Real(y(x))');
title('Gráfica de la parte real de y = c / sqrt(1 - (f/x)^2)');
grid on;

% Graficar la parte imaginaria de y(x)
figure(3);
plot(x, imag(y), 'r', 'LineWidth', 2);
xlabel('x (variable independiente)');
ylabel('Imag(y(x))');
title('Gráfica de la parte imaginaria de y = c / sqrt(1 - (f/x)^2)');
grid on;


