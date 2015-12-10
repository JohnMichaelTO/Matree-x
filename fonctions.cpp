#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "fonctions.h"
#include "conio.h"
#include "puissanceN.h"

/* Affichage d'un complexe */
void afficheComplexe(double re, double im)
{
    if(re != 0.0)
    {
        printf("%g", re);
        if(im > 0.0)
        {
            printf("+%gi", im);
        }
        else if(im < 0.0)
        {
            printf("%gi", im);
        }
    }
    else
    {
        if(im == 0.0)
        {
            printf("0");
        }
        else
        {
            printf("%gi", im);
        }
    }
    return;
}

/* Résolution d'équation de degré 1 */
double degre1(double a, double b)
{
    double x;
    x = -b / a;
    return x;
}

/* Fonction puissance */
double puissance(double x, int n)
{
    int i, imax;
    double sum = 1.0;
    imax = n;
    if(imax < 0)
    {
        imax = -n;
    }

    for(i = 0; i < imax; i++)
    {
        sum = sum * x;
    }

    if(n < 0)
    {
        sum = 1.0 / sum;
    }
    return sum;
}

/* Fonction racine carrée */
double racine_carree(double x)
{
    int i;
    double y = x;

    if(x != 0.0)
    {
        for(i = 0; i < 10; i++)
        {
            x = 0.5 * (x + y / x);
        }
    }
    else
    {
        x = 0.0;
    }

    return x;
}

/* Fonction racine cubique */
/*
double racine_cubique(double x)
{
    int i;
    double y = x;

    if(x != 0.0)
    {
        for(i = 0; i < 10; i++)
        {
            x = (1.0 / 3.0) * (2.0 * x + y / (x * x));
        }
    }
    else
    {
        x = 0.0;
    }

    return x;
}
*/

/* Factoriel */
double factoriel(int n)
{
    int cpt;
    double result = 1.0;

    for(cpt = 1; cpt <= n; cpt++)
    {
        result = result * cpt;
    }

    return result;
}


/* Calcul d'un polynome */
double calculPoly(double polynome[], int degree, double x)
{
    double P = 0.0, S = 0.0;
    int cpt;
    for(cpt = degree; cpt >= 0; cpt--)
    {
        S = polynome[cpt] * puissance(x, cpt);
        P = P + S;
    }
    return P;
}

/* Méthode de Newton pour la recherche d'une racine */
double newton(double polynome[], int degree, double x0)
{
    double polynomeDeriv[4] = {0.0};
    int cpt, p;
    double xi, P = 0.0, P2 = 0.0, Pverif = 0.0, Pverif2 = 0.0;
    xi = x0;
    int n = 0;
    int dicho = 0;

    p = degree - 1;
    for(cpt = p; cpt > 0; cpt--)
    {
        polynomeDeriv[cpt - 1] = polynome[cpt] * cpt;
    }

    do
    {
        n++;
        P = calculPoly(polynome, degree, xi);
        P2 = calculPoly(polynomeDeriv, p, xi);
        if(P2 != 0.0)
        {
            if(P < 0.0 && P2 < 0.0)
            {
                xi = xi + P / P2;
            }
            else
            {
                xi = xi - P / P2;
            }
            Pverif = P * puissance(10.0, 7.0);
            Pverif2 = P * puissance(10.0, 6.0);
            if(Pverif < 0.0)
            {
                Pverif = -Pverif;
            }

            if(Pverif2 < 0.0)
            {
                Pverif2 = -Pverif2;
            }

            if(absolue(xi) > 1000)
            {
                xi = -x0;
            }
        }
        else
        {
            dicho = 1;
        }
    } while((Pverif > 1.0 && Pverif2 >= 0.0) && n < 1000 && P2 != 0.0);

    // Si la méthode de Newton n'a pas réussi à converger vers une racine
    if(dicho == 1)
    {
        xi = dichotomie(polynome, degree);
    }

    return xi;
}

/* Retourner la valeur absolue d'un nombre */
double absolue(double x)
{
    if(x < 0.0)
    {
        x = -x;
    }
    return x;
}

/* Méthode de dichotomie */
double dichotomie(double polynome[], int degree)
{
    double a, b, precision, x = 0.0;

    a = -2;
    b = 2;
    precision = puissance(1.0, -7.0);
    while(absolue(a - b) > precision)
    {
        x = absolue(a + b) / 2.0;
        if(calculPoly(polynome, degree, x) == 0.0)
        {
            return x;
        }
        else if(calculPoly(polynome, degree, a) * calculPoly(polynome, degree, x) < 0.0)
        {
            b = x;
        }
        else
        {
            a = x;
        }
    }
    return x;
}

double cosinus(double x)
{
    double precision = 0.000001;
    int signe = 1;
    int n = -2;
    double cos = 0;
    do
    {
        cos = signe * (puissance(x, 2 + n) / factoriel(2 + n)) + cos;
        signe = -1 * signe;
        n = n + 2;
    } while ((puissance(x, n + 1) / factoriel(n + 1)) >= precision);
    return cos;
}
