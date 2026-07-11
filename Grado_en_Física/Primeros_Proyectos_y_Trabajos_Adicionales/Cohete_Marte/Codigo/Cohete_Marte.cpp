
#include <stdio.h>
#include <math.h>

// Ctes con unidades en SI
#define h 0.001                                           // paso en horas
const double resc = 3600;                                // reescalado de una hora
const double D = 2000000;//149.6e9;                                // Reescalado distancia
const double G = (6.67e-11) * resc * resc / (D * D * D); // cte gravitacion univ
#define Mt 5.9736e24                                    // Masa tierra
#define Mm  6.39e23                                      // masa Marte
#define Ms 1.989e30                                     // masa Sol
const double dst = (149.6e9) / D;                        // Dist Sol Tierra
const double dsm = (227.9e9) / D;                        // Dist Sol Marte
const double wt = (1.990987e-7) * resc;                  // frec angular de la tierra respecto a la sol
const double wm = (1.057e-7) * resc;                     // frec angular de Marte respecto a la sol
const double Rt=6.371e6/D ;                                      // Radio terrestre
const double Vesc = 2971.1 * resc / D;                   // velocidad de escape de la Tierra
const double Desfase_m=4.688965-wm*9388.720000;          //wt+theta=phi siendo t el tiempo al que llega el cohete a Marte y phi su coordenada
// Ctes aux para las derivadas
#define Pi 3.141592653589793238462643383279502884197169399
// funciones
void valoresinicales(double *r, double *phi, double *p_r, double *p_phi);
void impulso(double *r, double *phi, double *p_r, double *p_phi);
void orbitamarte(double r, double rm, double phi, double *p_r, double *p_phi);
double derivada_p_phi(double r, double r_t, double r_m, double phi, double t);
double derivada_p_r(double r, double r_t, double r_m, double phi, double p_phi, double t);
void Energia(double *E, double *T, double *V, double r, double r_t, double r_m, double phi, double p_r, double p_phi, double t);

int main()
{

    int i, iteraciones;

    //variables del cohete: posicion y tiempo, r y phi son las coordenadas, p_r y p_ phi los momento y d_ indica derivada
    double r, phi, p_r, p_phi, d_r, d_phi, d_p_r, d_p_phi, t;

    // Distancia a la Tierra y Marte
    double r_t, r_m;

    // Energias
    double E, T, V;

    // definimos variables auxiliares para la posicion que se van a usar para calcular las pendientes en distintos puntos definidos por las
    // variables auxiliares
    double r_aux, phi_aux, p_r_aux, p_phi_aux,ejecutado,ejecutado1;

 // definimos las k que son el paso en cada coordenada, tiene dimension tres pq en cada paso las calculamos en tres puntos distintos
    double kr[4], kphi[4], kp_r[4], kp_phi[4];

    // Archivos en los que vamos a imprimr los resultados
    FILE *f1, *f2, *f3;
    f1 = fopen("Cohete_data.dat", "w");
    f2 = fopen("Energia.dat", "w");
    f3 = fopen("Cohete_respecto_a_Marte.dat", "w");

    // Inicializamos la posicion y los momentos
    valoresinicales(&r, &phi, &p_r, &p_phi);
    // Calculamos tambien los valores iniciales para r_t y r_m
    r_t = sqrt(dst * dst + r * r - 2 * r * dst * cos(phi - wt * t));
    r_m = sqrt(dsm * dsm + r * r - 2 * r * dsm * cos(phi - wm * t));
    // Guardamos los datos iniciales en coordenadas cartesianas
        fprintf(f1, "%lf, %lf \n", r * cos(phi), r * sin(phi));           // cohete
        fprintf(f1, "%lf, %lf \n", dst * cos(wt * t), dst * sin(wt * t)); // Tierra
        fprintf(f1, "%lf, %lf \n", dsm * cos(wm * t+Desfase_m), dsm * sin(wm * t+Desfase_m)); // Marte
        fprintf(f1, "\n");

    // Numero de iteraciones que vamos a realizar
    iteraciones = 20000000;
    i = 0;

    //Para que la funcion impulso se ejecute solo una vez
    ejecutado=0;
    ejecutado1=0;

    // Comenzamos el bucle en el que se calculara iterativamente la posicion respecto al tiempo
    while (i < iteraciones)
    {
        
        // Calculamos el tiempo en cada iteracion
        t = i * h;

        // Calculamos las pendientes para poder calcular las k[0]
        d_r = p_r;
        d_phi = p_phi / (r * r);
        d_p_phi = derivada_p_phi(r, r_t, r_m, phi, t);
        d_p_r = derivada_p_r(r, r_t, r_m, phi, p_phi, t);

        kr[0] = h * d_r;
        kphi[0] = h * d_phi;
        kp_phi[0] = h * d_p_phi;
        kp_r[0] = h * d_p_r;

        //Hacemos lo mismo pero para calcular las k[1] usando las variables auxiliares y desplazandonos k[0]/2 en el espacio y h/2 en t
        r_aux = r + kr[0] / 2;
        r_t = sqrt(dst * dst + r_aux * r_aux - 2 * r_aux * dst * cos(phi - wt * (t + h / 2)));
        r_m = sqrt(dsm * dsm + r * r - 2 * r * dsm * cos(phi - wm * t+Desfase_m));
        phi_aux = phi + kphi[0] / 2;
        p_r_aux = p_r + kp_r[0] / 2;
        p_phi_aux = p_phi + kp_phi[0] / 2;

        d_r = p_r_aux;
        d_phi = p_phi_aux / (r_aux * r_aux);
        d_p_phi = derivada_p_phi(r_aux, r_t, r_m, phi_aux, t + h / 2);
        d_p_r = derivada_p_r(r_aux, r_t, r_m, phi_aux, p_phi_aux, t + h / 2);

        kr[1] = h * d_r;
        kphi[1] = h * d_phi;
        kp_phi[1] = h * d_p_phi;
        kp_r[1] = h * d_p_r;

        // para k[2] nos desplazamos k[1]/2 en el espacio y h/2 en t
        r_aux = r + kr[1] / 2;
        r_t = sqrt(dst * dst + r_aux * r_aux - 2 * r_aux * dst * cos(phi - wt * (t + h / 2)));
        r_m = sqrt(dsm * dsm + r * r - 2 * r * dsm * cos(phi - wm * t+Desfase_m));
        phi_aux = phi + kphi[1] / 2;
        p_r_aux = p_r + kp_r[1] / 2;
        p_phi_aux = p_phi + kp_phi[1] / 2;

        d_r = p_r_aux;
        d_phi = p_phi_aux / (r_aux * r_aux);
        d_p_phi = derivada_p_phi(r_aux, r_t, r_m, phi_aux, t + h / 2);
        d_p_r = derivada_p_r(r_aux, r_t, r_m, phi_aux, p_phi_aux, t + h / 2);

        kr[2] = h * d_r;
        kphi[2] = h * d_phi;
        kp_phi[2] = h * d_p_phi;
        kp_r[2] = h * d_p_r;

        // para k[3] nos desplazamos k[2] en el espacio y h en t
        r_aux = r + kr[2];
        r_t = sqrt(dst * dst + r_aux * r_aux - 2 * r_aux * dst * cos(phi - wt * (t + h / 2)));
        r_m = sqrt(dsm * dsm + r * r - 2 * r * dsm * cos(phi - wm * t+Desfase_m));
        phi_aux = phi + kphi[2];
        p_r_aux = p_r + kp_r[2];
        p_phi_aux = p_phi + kp_phi[2];

        d_r = p_r_aux;
        d_phi = p_phi_aux / (r_aux * r_aux);
        d_p_phi = derivada_p_phi(r_aux, r_t, r_m, phi_aux, t + h);
        d_p_r = derivada_p_r(r_aux, r_t, r_m, phi_aux, p_phi_aux, t + h);

        kr[3] = h * d_r;
        kphi[3] = h * d_phi;
        kp_phi[3] = h * d_p_phi;
        kp_r[3] = h * d_p_r;

        // una vez se han calculado todas las k podemos calcular la posicion para t+h
        r = r + (kr[0] + 2 * kr[1] + 2 * kr[2] + kr[3]) / 6;
        phi = phi + (kphi[0] + 2 * kphi[1] + 2 * kphi[2] + kphi[3]) / 6;
        p_phi = p_phi + (kp_phi[0] + 2 * kp_phi[1] + 2 * kp_phi[2] + kp_phi[3]) / 6;
        p_r = p_r + (kp_r[0] + 2 * kp_r[1] + 2 * kp_r[2] + kp_r[3]) / 6;
        
       if(ejecutado==0)
        {
            if(abs(p_r)<0.07)
            {
                impulso(&r, &phi, &p_r, &p_phi);
                ejecutado=1;
            }
        }



        if (i % 10000 == 0)
        {
            // Guardamos los datos de cada paso pero en coordenas cartesianas
            fprintf(f1, "%lf, %lf \n", r * cos(phi), r * sin(phi));           // cohete
            fprintf(f1, "%lf, %lf \n", dst * cos(wt * t), dst * sin(wt * t)); // Tierra
            fprintf(f1, "%lf, %lf \n", dsm * cos(wm * t+Desfase_m), dsm * sin(wm * t+Desfase_m)); // Marte
            fprintf(f1, "\n");
        }


        // Calculamos las nuevos radios a la tierra y marte
        r_t = sqrt(dst * dst + r * r - 2 * r * dst * cos(phi - wt * t));
        r_m = sqrt(dsm * dsm + r * r - 2 * r * dsm * cos(phi - wm * t+Desfase_m));
      
      if(ejecutado1==0)
        {
        if(t>9388.710000)
        {
             orbitamarte(r,r_m,phi,&p_r,&p_phi);
             ejecutado1=1;
        }
        }



                


        if (i % 10000 == 0)
        {
            if(t>0)
            {
                fprintf(f3, "%lf, %lf \n \n",r*cos(phi)-dsm*cos(wt+Desfase_m),r*sin(phi)-dsm*sin(wt+Desfase_m));
            }
        }

        if (i % 10000 == 0)
        {
            // Calculamos las energias
            Energia(&E, &T, &V, r, r_t, r_m, phi, p_r, p_phi, t);
            fprintf(f2, "%lf   %lf  %lf   %lf\n", t, E, T, V);
        }


        i = i + 1;
    }
}

void valoresinicales(double *r, double *phi, double *p_r, double *p_phi)
{
    double r_orb,v,v_salida,theta;
    //v_salida velocidad con la que  sale el cohete y theta el angulo de salida
    *r = dst + Rt;
    v_salida=sqrt(2*G*Mt/Rt);
    theta=-15*Pi/360;
    v=wt*dst;
    *p_phi=*r*v-v_salida*sin(theta);
    *p_r=v_salida*cos(theta);
}

void impulso(double *r, double *phi, double *p_r, double *p_phi)
{
    double v_hohmann, v,Deltav;
    //p_r=0
    v=*p_phi/(*r);
    v_hohmann=sqrt(2*G*Ms*dsm/(*r*(*r+dsm)));
    Deltav=v_hohmann-v;
    *p_phi=*r*(v+Deltav);
}

void orbitamarte(double r, double rm, double phi, double *p_r, double *p_phi)
{
    double v,Deltav;
    v=*p_phi/r;
    Deltav=wm*dsm+sqrt(G*Mm/rm)-v;
    v=v+Deltav;
    *p_phi=v*r;
}

double derivada_p_phi(double r, double r_t, double r_m, double phi, double t)
{
    // variable a calcular:
    double d_p_phi;
    d_p_phi = -G * ((Mt * r / pow(r_t, 3)) * sin(phi - wt * t) + (Mm * r / pow(r_m, 3)) * sin(phi - wm * t+Desfase_m));
    return d_p_phi;
}

double derivada_p_r(double r, double r_t, double r_m, double phi, double p_phi, double t)
{
    // variable a calcular:
    double d_p_r;

    d_p_r = pow(p_phi, 2) / pow(r, 3) - G * (Ms / pow(r, 2) + (Mt / pow(r_t, 3)) * (r - dst * cos(phi - wt * t)) + (Mm / pow(r_m, 3)) * (r - dsm * cos(phi - wm * t+Desfase_m)));
    return d_p_r;
}


void Energia(double *E, double *T, double *V, double r, double r_t, double r_m, double phi, double p_r, double p_phi, double t)

{
    // Vamos a calcular la energia en las unidades rescaladas, por eso vamos a usar la Gprima rescalada
    // Energia cinetica
    *T = pow(p_r, 2) / 2 + pow(p_phi / (r), 2) / 2;
    // Energia potencial
    *V = -G * Ms / (r)-G * Mt / (r_t)-G * Mm / (r_m);
    *E = *T + *V;
}