#include <stdio.h>
#include <math.h>
#include "gsl_rng.h" // Libreria para generación de números aleatorios
#include <stdlib.h>  // para memoria dinámica
#include <fstream>  //para leer la matriz
#include <string>
#define N 40 // Dimensión de la matriz de los espines
#define mu 1 // número de partones almacenados

gsl_rng *tau;

void inicializarmatrizrandom(int S[N][N]);
void iniciarmatrizdeformada(int S[N][N], int epsilon [mu][N][N],int n);
double calc_DeltaE(int S[N][N], double (*omega)[N][N][N][N], double theta[N][N],int n, int m);
double calcularp(int S[N][N], double T, double DeltaE, int n, int m);
void calcularomega(double (*omega)[N][N][N][N], int epsilon[mu][N][N], double a[mu]);
void calcular_a(int epsilon[mu][N][N], double a[mu]);
void calculartheta(double theta[N][N], double (*omega)[N][N][N][N]);
void calcularsolapamiento(int S[N][N], double a[mu], int epsilon[mu][N][N], double m[mu]);
void leerMatrizDesdeArchivo(int A[mu][N][N], const std::string& nombreArchivo, int a);

int main()
{
    int i, j, k, l, iteraciones, n, m, paso, tiempo;
    double T, random, p, DeltaE;
    int S[N][N];
    double a[mu], theta[N][N], solapamiento[mu];

    double (*omega)[N][N][N][N] = (double(*)[N][N][N][N])malloc(N * N * N * N * sizeof(double));
    int semilla = 1934050348; // semilla con la que vamos a generar los números aleatorios
    FILE *f1, *f2;
    // P configuraciones, la primera coordenada se refiere al patrón y las otras dos a la posición
    int epsilon[mu][N][N];

    extern gsl_rng *tau;               // Puntero al estado del número aleatorio
    tau = gsl_rng_alloc(gsl_rng_taus); // Inicializamos el puntero
    gsl_rng_set(tau, semilla);         // Inicializamos la semilla

    f1 = fopen("hopfield_data.dat", "w");         // Abrimos el archivo en el que vamos a escribir las matrices en el tiempo
    f2 = fopen("hopfield_solapamiento.dat", "w"); // Abrimos el archivo en el que vamos a escribir el solapamiento en el tiempo

    T = 0.03;              // Elegimos la temperatura del sistema entre 0 y 5

    //Primero leemos los distinos patrones y los imprimos en epsilon

    for(i=0; i<mu; i++)
    {
    std::string filename = std::to_string(i) + ".txt";
    leerMatrizDesdeArchivo(epsilon,filename, i);
    }


    // inicializamos la matriz ya sea aleatorio o de patron deformado
    //inicializarmatrizrandom(S);  
    iniciarmatrizdeformada(S,epsilon, 0);  //el tercer numero indica el patron que vamos a deformar

    iteraciones = 100*N*N; // número de saltos que vamos a dar
    paso = N * N;        // cada 'paso' veces se mostrará la nueva matriz


    // Calculamos parámetros constantes como omega, a y theta

    calcular_a(epsilon, a);
    calcularomega(omega, epsilon, a);
    calculartheta(theta, omega);

    for (i = paso; i <= iteraciones; i++)
    {
        // elegimos un punto al azar
        n = gsl_rng_uniform_int(tau, N); // número aleatorio entero [0,N-1]
        m = gsl_rng_uniform_int(tau, N);

        // Calculamos DeltaE para el cálculo de la probabilidad p
        DeltaE = calc_DeltaE(S, omega, theta,n,m);
        // Evaluamos p
        p = calcularp(S, T, DeltaE, n, m);

        random = gsl_rng_uniform(tau); // generamos un número aleatorio entre 0 y 1.

        if (random < p)
        {
            if (S[n][m] == 1)
            {
                S[n][m] = 0;
            }else
            {S[n][m] = 1;}
        }

        if ((i % paso) == 0) // Cada 'paso' veces se muestra el estado de la matriz y el solapamiento
        {
            for (j = 0; j < N; j++)
            {
                for (k = 0; k < (N - 1); k++)
                {
                    fprintf(f1, "%i ,", S[j][k]);
                }
                fprintf(f1, "%i ", S[j][N - 1]); // El último de cada fila va sin coma
                fprintf(f1, "\n");
            }
            fprintf(f1, "\n");

            tiempo = i / paso;

            calcularsolapamiento(S, a, epsilon, solapamiento);
            // Primero imprimimos el tiempo montecarlo
            fprintf(f2, "%i, ",tiempo);
            for (j = 0; j < (mu - 1); j++)
            {
                fprintf(f2, "%lf, ", solapamiento[j]);
            }
            fprintf(f2, "%lf ", solapamiento[j]); // El último de cada fila va sin coma
            fprintf(f2, "\n");
        }
    }
    fclose(f1);
    fclose(f2);
    free(omega); // Liberar la memoria asignada a omega
    return 0;
}

void inicializarmatrizrandom(int S[N][N])
{
    int i, j;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            S[i][j] = gsl_rng_uniform_int(tau, 2);
        }
    }
}

void iniciarmatrizdeformada(int S[N][N], int epsilon [mu][N][N], int n)
{
    int i,j;
    double p;
    for(i=0;i<N;i++)
    {
        for(j=0;j<N;j++)
        {
            p=gsl_rng_uniform(tau);
            if(p<0.2)
            {
                S[i][j]=epsilon[n][i][j];
                if (S[i][j] == 1)
                {
                    S[i][j] = 0;
                }else
                {S[i][j] = 1;}
            }else 
            {
                S[i][j]=epsilon[n][i][j];
            }

        }

    }
}

double calcularp(int S[N][N], double T, double DeltaE, int n, int m)
{
    int i;
    double p;

    // Calculamos el p provisional
    p = exp(-1. * DeltaE / T);

    if (1. < p)
    {
        p = 1.;
    }

    return p;
}

double calc_DeltaE(int S[N][N], double (*omega)[N][N][N][N], double theta[N][N],int n, int m)
{
    // posición de la matriz omega que se está calculando
    size_t i, j, a, b; // Cambio de int a size_t para los índices
    a=n;
    b=m;
    double DeltaE = 0;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            DeltaE = DeltaE-(*omega)[i][j][a][b]*S[i][j];
        }
    }

    //Delta de S(n,m) se puede expresar como 1-2S(n,m) S=0 entonces deltaS 1 y igual para S 1
    DeltaE = (DeltaE+theta[n][m])*(1-2*S[n][m]);




    return DeltaE;
}

void calcularomega(double (*omega)[N][N][N][N], int epsilon[mu][N][N], double a[mu])
{
    int p; // indica el patrón al que nos estamos refiriendo
    // posición de la matriz omega que se está calculando
    size_t i, j, k, l; // Cambio de int a size_t para los índices


    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            for (k = 0; k < N; k++)
            {
                for (l = 0; l < N; l++)
                {

                    if (i == k && j == l)
                    {
                        (*omega)[i][j][k][l] = 0;
                    }
                    else
                    {
                        //iniciamos omega a 0
                        (*omega)[i][j][k][l] = 0;

                        for (p=0; p<mu; p++)
                        {
                            (*omega)[i][j][k][l]=(*omega)[i][j][k][l]+(epsilon[p][i][j]-a[p])*(epsilon[p][k][l]-a[p]);
                        }
                         //Dividimos entre N²
                        (*omega)[i][j][k][l]=(*omega)[i][j][k][l]/(N*N);
                    }
                }
            }
        }
    }

   


}

void calcular_a(int epsilon[mu][N][N], double a[mu])
{
    int i, j, p;

    for (p = 0; p < mu; p++)
    {
        a[p] = 0;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                a[p] = a[p] + epsilon[p][i][j];
            }
        }
        a[p] = a[p] / (N * N);
    }
}

void calculartheta(double theta[N][N], double (*omega)[N][N][N][N])
{
    size_t i, j, k, l;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            theta[i][j] = 0;
            for (k = 0; k < N; k++)
            {
                for (l = 0; l < N; l++)
                {
                    theta[i][j] = theta[i][j] + (*omega)[i][j][k][l];
                }
            }
            theta[i][j] = 0.5 * theta[i][j];
        }
    }
}

void calcularsolapamiento(int S[N][N], double a[mu], int epsilon[mu][N][N], double m[mu])
{
    int i, j, l;
    double aux;
    for (l = 0; l < mu; l++)
    {
        aux = (N * N * a[l])* (1 - a[l]);
        m[l] = 0;
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                m[l] = m[l] + (epsilon[l][i][j] - a[l]) * (S[i][j] - a[l]);
            }
        }
        m[l] = m[l] / aux;
    }
}


void leerMatrizDesdeArchivo(int A[mu][N][N], const std::string& nombreArchivo, int a)
{
    std::ifstream archivo(nombreArchivo);


        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                archivo >> A[a][i][j];
            }
        }

        archivo.close();
}