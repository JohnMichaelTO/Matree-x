#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"
#include "matrice.h"
#include "fonctions.h"
#include "conio.h"
#define PI 3.14159265358979323846264338327950288419716939937510

/* [DEBUT] RESOLUTION EQUATION DEGRE 3 */
void degre3(t_racine *racine, double a, double b, double c, double d, int i = 0)
{
    // Division par a de tous les coefficients
    b /= a;
    c /= a;
    d /= a;
    a /= a;

    // Définition des variables locales
    double p, q, p2, q2, pose, delta;
    int cpt;
    int degree = 3;

    // Définition des variables servants aux racines u et v
    double u = 0, v = 0, u1 = 0, v1 = 0;
    double u2r = 0, u2i = 0, u3r = 0, u3i = 0, v2r = 0, v2i = 0, v3r = 0, v3i = 0;

    // Définition des variables servants aux variables temporaires (changement de variables)
    double t1 = 0, t2 = 0, t2r = 0, t2i = 0, t3r = 0, t3i = 0, z2r = 0, z3r = 0;

    // Définition des variables stockant le résultat, et le calcul des racines complexes
    double z[5], r, A;

    // Définition de j et j^2 avec les parties réelles et imaginaires
    double jr = -0.5;
    double ji = racine_carree(3.0) / 2.0;
    double jr2 = jr;
    double ji2 = -ji;

    if(b == 0.0)
    {
        pose = 0.0;
        p = c / a;
        q = d / a;
    }
    else
    {
        pose = - b / (3.0 * a);
        p = (-(b * b)/(3.0 * a * a))+(c / a);
        q = (b / (27.0 * a)) * (((2.0 * b * b) / (a * a)) - ((9.0 * c) / a)) + (d / a);
    }

    if(p == 0 && q == 0)
    {
        z[0] = pose;
        racine[i].re = z[0];
        racine[i].im = 0.0;
        i++;
    }
    else if(q == 0)
    {
        z[0] = pose;

        if(-p < 0.0)
        {
            z[1] = -racine_carree(p) + pose;
            z[2] = racine_carree(p) + pose;
        }
        else
        {
            z[1] = -racine_carree(-p) + pose;
            z[2] = racine_carree(-p) + pose;
        }

        for(cpt = i; cpt < degree; cpt++)
        {
            racine[cpt].re = z[cpt];
            racine[cpt].im = 0.0;
        }
    }
    else if(p == 0)
    {
        if(-q < 0.0)
        {
            z[0] = pow(q, (1.0 / 3.0));
            // z[0] = racine_cubique(q);
        }
        else
        {
            z[0] = pow(-q, (1.0 / 3.0));
            // z[0] = racine_cubique(-q);
        }

        z[0] = z[0] + pose;
        racine[i].re = z[i];
        racine[i].im = 0.0;
        i++;
    }
    else
    {
        p2 = -q;
        q2 = -(p * p * p) / 27.0;

        p2 = -p2;

        delta = calcul_delta(1.0, p2, q2);

        if(delta > 0.0)
        {
            u = racine_reelle(1.0, p2, delta, 1);
            v = racine_reelle(1.0, p2, delta, 2);
        }

        if(delta > 0.0)
        {
            if(u < 0.0)
            {
                u1 = -pow(-u, (1.0 / 3.0));
                // u1 = -racine_cubique(-u);
            }
            else
            {
                u1 = pow(u, (1.0 / 3.0));
                // u1 = racine_cubique(u);
            }

            if(v < 0.0)
            {
                v1 = -pow(-v, (1.0 / 3.0));
                // v1 = -racine_cubique(-v);
            }
            else
            {
                v1 = pow(v, (1.0 / 3.0));
                // v1 = racine_cubique(v);
            }

            u2r = u1 * jr;
            u2i = u1 * ji;

            u3r = u1 * jr2;
            u3i = u1 * ji2;

            v2r = v1 * jr2;
            v2i = v1 * ji2;

            v3r = v1 * jr;
            v3i = v1 * ji;

            t1 = u1 + v1;
            t2r = u2r + v2r;
            t2i = u2i + v2i;
            t3r = u3r + v3r;
            t3i = u3i + v3i;
        }

        if(delta > 0.0)
        {

            z[0] = t1 + pose;
            z2r = t2r + pose;
            z3r = t3r + pose;

            // printf("%c1 = %lf\r\n", variable, z[0]);
            racine[i].re = z[i];
            racine[i].im = 0.0;
            i++;

            // printf("%c2 = ", variable);
            racine[i].re = z2r;
            racine[i].im = t2i;
            i++;
            // afficheComplexe(z2r, t2i);

            // printf("%c3 = ", variable);
            racine[i].re = z3r;
            racine[i].im = t3i;
            i++;
            // afficheComplexe(z3r, t3i);
        }
        else if(delta == 0.0)
        {
            t1 = 3.0 * q / p;
            t2 = - (3.0 * q) / (2.0 * p);

            z[0] = t1 + pose;
            z[1] = t2 + pose;
            z[2] = z[1];

            for(cpt = i; cpt < degree; cpt++)
            {
                racine[cpt].re = z[cpt];
                racine[cpt].im = 0.0;
            }

            /*
            if(pose < 0.0)
            {
                printf("%c1 = %lf\r\n", variable, z[0]);
                printf("%c2 = %c3 = %lf", variable, variable, z[1]);
            }
            else if(pose > 0.0)
            {
                printf("%c1 = %lf\r\n", variable, z[0]);
                printf("%c2 = %c3 = %lf", variable, variable, z[1]);
            }
            else
            {
                printf("%c1 = %lf\r\n", variable, z[0]);
                printf("%c2 = %lf", variable, z[1]);
            }
            */
        }
        else
        {
            r = racine_carree(-p / 3.0);
            A = acos((-q / 2.0) * racine_carree(27.0 / (-p * p * p)));

            for (cpt = i; cpt < degree ; cpt++)
            {
                z[cpt] = 2.0 * r * cos((1.0 / 3.0) * A + ((2 * cpt * PI) / 3.0)) + pose;

                if(z[cpt] == -0.0)
                {
                    z[cpt] = 0.0;
                    // printf("%c%i = 0\r\n", variable, cpt);
                }
                else
                {
                    // printf("%c%i = %lf\r\n", variable, cpt, z[cpt]);
                }
                racine[cpt].re = z[cpt];
                racine[cpt].im = 0.0;
            }
        }
    }
}
/* [FIN] RESOLUTION EQUATION DEGRE 3 */
