#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "fonctions.h"
#include "conio.h"

/* Calcul du delta */
double calcul_delta(double a, double b, double c)
{
    return (b * b) - (4.0 * a * c);
}

/* Calcul une racine réelle */
double racine_reelle(double a, double b, double delta, int numero_racine)
{
    double x;

    if(numero_racine == 1)
    {
        x = (-b - racine_carree(delta)) / (2.0 * a);
    }
    else
    {
        x = (-b + racine_carree(delta)) / (2.0 * a);
    }
    return x;
}

/* Calcul une racine complexe */
double racine_complexe(double a, double b, double delta, char partie)
{
    double z;

    if(partie == 'r')
    {
        z = (-b) / (2.0 * a);
    }
    else
    {
        z = (racine_carree(-delta)) / (2.0 * a);
    }
    return z;
}

/* Résolution d'équation de degré 2 */
void degre2(t_racine *racine, double a, double b, double c, double changementVariable = 0.0, int i = 0)
{
    double x1, x2, delta;

    if(b == 0.0)
    {
        if((-c / a) >= 0.0)
        {
            x1 = racine_carree(-c / a) + changementVariable;
            x2 = -racine_carree(-c / a) + changementVariable;

            racine[i].re = x1;
            racine[i + 1].re = x2;
            racine[i].im = 0.0;
            racine[i + 1].im = 0.0;
            /*
            printf("%c%i = %lf\r\n", variable, numero_racine, x1);
            printf("%c%i = %lf\r\n", variable, numero_racine + 1, x2);
            */
        }
        else
        {
            x1 = -racine_carree(c / a) + changementVariable;
            x2 = racine_carree(c / a) + changementVariable;

            racine[i].re = 0.0;
            racine[i + 1].re = 0.0;
            racine[i].im = x1;
            racine[i + 1].im = x2;
            /*
            printf("%c%i = %lfi\r\n", variable, numero_racine, x1);
            printf("%c%i = %lfi\r\n", variable, numero_racine + 1, x2);
            */
        }
    }
    else if(c == 0.0)
    {
        x1 = degre1(a, b) + changementVariable;
        x2 = 0.0 + changementVariable;

        racine[i].re = x1;
        racine[i + 1].re = x2;
        racine[i].im = 0.0;
        racine[i + 1].im = 0.0;
        /*
        printf("%c%i = %lf\r\n", variable, numero_racine, x1);
        printf("%c%i = %lf\r\n", variable, numero_racine + 1, x2);
        */
    }
    else
    {
        delta = calcul_delta(a, b, c);

        if(delta > 0.0)
        {
            x1 = racine_reelle(a, b, delta, 1) + changementVariable;
            x2 = racine_reelle(a, b, delta, 2) + changementVariable;

            racine[i].re = x1;
            racine[i + 1].re = x2;
            racine[i].im = 0.0;
            racine[i + 1].im = 0.0;
            /*
            printf("%c%i = %lf\r\n", variable, numero_racine, x1);
            printf("%c%i = %lf\r\n", variable, numero_racine + 1, x2);
            */
        }
        else if(delta == 0.0)
        {
            x1 = (-b) / (2.0 * a) + changementVariable;
            x2 = x1;

            racine[i].re = x1;
            racine[i + 1].re = x2;
            racine[i].im = 0.0;
            racine[i + 1].im = 0.0;

            // printf("%c%i = %c%i = %lf\r\n", variable, numero_racine, variable, numero_racine + 1, x1);
        }
        else
        {
            x1 = racine_complexe(a, b, delta, 'r') + changementVariable;
            x2 = racine_complexe(a, b, delta, 'i');

            racine[i].re = x1;
            racine[i + 1].re = x1;

            racine[i].im = -x2;
            racine[i + 1].im = x2;

            /*
            printf("%c%i = ", variable, numero_racine);
            afficheComplexe(x1, -x2);
            printf("%c%i = ", variable, numero_racine + 1);
            afficheComplexe(x1, x2);
            */
        }
    }
    return;
}
