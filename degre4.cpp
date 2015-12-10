#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "fonctions.h"
#include "conio.h"

/* Résolution d'équation de degré 4 */
void degre4(t_racine *racine, double a, double b, double c, double d, double e, int i = 0)
{
    // Division par a de tous les coefficients
    if(a != 1.0)
    {
        b /= a;
        c /= a;
        d /= a;
        e /= a;
        a /= a;
    }

    double p, q, r, y, pose = 0.0;

    // Calcul des coefficients p, q et r
    if(b == 0.0)
    {
        // Si b = 0, on a p = c/a, q = d/a, r = e/a, or a = 1, donc on a :
        p = c;
        q = d;
        r = e;
    }
    else
    {
        p = (-3.0 * b * b) / (8.0 * a * a) + (c / a);
        // q = pow((b / 2.0), 3.0) / pow(a, 3.0) - ((1.0 / 2.0) * b * c) / (a * a) + (d / a);
        q = puissance((b / 2.0), 3) / puissance(a, 3) - ((1.0 / 2.0) * b * c) / (a * a) + (d / a);
        // r = -3.0 * pow(((b / 4.0) / a), 4.0) + c * (pow((b / 4.0), 2.0) / pow(a, 3.0)) - ((1.0 / 4.0) * b * d) / (a * a) + (e / a);
        r = -3.0 * puissance(((b / 4.0) / a), 4) + c * (puissance((b / 4.0), 2) / puissance(a, 3)) - ((1.0 / 4.0) * b * d) / (a * a) + (e / a);

        // On pose x = z + (-b/4a)
        pose = - b / 4 * a;
    }

    /*************************** Cas où q = 0 : équation bicarrée ***************************/
    /* Pour optimiser le programme, on pourrait faire un cas qui gère une équation bicarrée */
    /*************************** Cas où q = 0 : équation bicarrée ***************************/

    degre3(racine, 8.0, (-4.0 * p), (-8.0 * r), (4.0 * p * r - q * q), i);
    y = racine[0].re;
    /*
    printf("y = %g\n", y);
    system("pause");
    */

    double A, B, C;
    A = 1.0;
    B = racine_carree(2 * y - p);
    C = (q / (2 * racine_carree(2 * y - p)));

    degre2(racine, A, B, y - C, pose, i);
    degre2(racine, A, -B, y + C, pose, i + 2);
}
