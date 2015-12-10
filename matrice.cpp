#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "puissanceN.h"
#include "fonctions.h"
#include "conio.h"
#define EXP 2.71828182845904523536028747135266249775724709369995

/* Saisie d'une matrice */
void saisieMatrice(int tailleMatrice, Matrice *matrice)
{
    int i, X, x, y;
    double saisie;

    printf("\n   Saisir la matrice :\n  ");
    X = wherex();
    x = X;
    y = wherey();
    gotoxy(x, y);
    for (i = 0 ; i < tailleMatrice*tailleMatrice; i++)
    {
        scanf("%lf", &saisie);
        matrice->coefficient[i].re = saisie;
        matrice->coefficient[i].im = 0.0;
        if((i + 1) % tailleMatrice != 0)
        {
            x = x + DEF_X;
            gotoxy(x, y);
        }
        else
        {
            x = X;
            y = y + 2;
            gotoxy(x, y);
        }
    }
}

/* Transposée d'une matrice */
void matTrans(int tailleMatrice, Matrice matrice, Matrice* resultante)
{
    int lig = 0, col = 0;
    t_complexe **coeffTampon, **coeffResultant;

    creatMatrice2D(tailleMatrice, &coeffTampon);
    creatMatrice2D(tailleMatrice, &coeffResultant);

    // Copie 1D->2D
    copie1DtoD2(tailleMatrice, matrice, coeffTampon);

    // Transposition
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            coeffResultant[col][lig].re = coeffTampon[lig][col].re;
            coeffResultant[col][lig].im = coeffTampon[lig][col].im;
        }
    }

    // Copie 2D->1D
    copie2Dto1D(tailleMatrice, coeffResultant, resultante);

    /* Libération de la mémoire */
    freeMatrice2D(tailleMatrice, &coeffTampon);
    freeMatrice2D(tailleMatrice, &coeffResultant);
}

/* Calcul de la trace d'une matrice */
t_complexe matTrace(int tailleMatrice, Matrice matrice)
{
    int i = 0;
    t_complexe trace;
    trace.re = 0.0;
    trace.im = 0.0;

    for (i = 0 ; i < (tailleMatrice * tailleMatrice) ; i += tailleMatrice+1)
    {
        trace.re += matrice.coefficient[i].re;
        trace.im += matrice.coefficient[i].im;
    }
    return trace;
}

/* Inversion de matrice */
int matInverse(int tailleMatrice, Matrice matrice, Matrice* resultante)
{
    // int i, j, k;
    t_complexe det;
    int erreur = 1;
    /*
    t_complexe **matriceA, **matriceI, pivot, coef, temp;
    Matrice matriceIdentite;
    */
    det = matDet(tailleMatrice, matrice);

    if(det.re != 0.0 || det.im != 0.0)
    {
        if(tailleMatrice == 2)
        {
            matrice.coefficient[1].re = -matrice.coefficient[1].re;
            matrice.coefficient[1].im = -matrice.coefficient[1].im;

            matrice.coefficient[2].re = -matrice.coefficient[2].re;
            matrice.coefficient[2].im = -matrice.coefficient[2].im;

            resultante->coefficient[0] = divComplexe(matrice.coefficient[3], det);
            resultante->coefficient[1] = divComplexe(matrice.coefficient[1], det);
            resultante->coefficient[2] = divComplexe(matrice.coefficient[2], det);
            resultante->coefficient[3] = divComplexe(matrice.coefficient[0], det);

            matrice.coefficient[1].re = -matrice.coefficient[1].re;
            matrice.coefficient[1].im = -matrice.coefficient[1].im;

            matrice.coefficient[2].re = -matrice.coefficient[2].re;
            matrice.coefficient[2].im = -matrice.coefficient[2].im;
        }
        else
        {
            matMin(tailleMatrice, matrice, resultante);
            matCom(tailleMatrice, *resultante, resultante);
            matTrans(tailleMatrice, *resultante, resultante);
            det.re = 1.0 / det.re;
            det.im = 1.0 / det.im;
            matMutScalaire(tailleMatrice, *resultante, det, resultante);
        }
        /*else
        {
            creatMatrice2D(tailleMatrice, &matriceA);
            copie1DtoD2(tailleMatrice, matrice, matriceA);

            creatMatrice(tailleMatrice, &matriceIdentite);
            matIdent(tailleMatrice, matriceIdentite, &matriceIdentite);

            creatMatrice2D(tailleMatrice, &matriceI);
            copie1DtoD2(tailleMatrice, matriceIdentite, matriceI);

            for(i = 0; i < tailleMatrice; i++)
            {
            do
            {
                pivot.re = matriceA[i][i].re;
                pivot.im = matriceA[i][i].im;

                if(pivot.re == 0.0 && pivot.im == 0.0)
                {
                    for(j = 0; j < tailleMatrice; j++)
                    {
                        if(i + 1 < tailleMatrice)
                        {
                            swap(&matriceA[i][j], &matriceA[i + 1][j]);
                            swap(&matriceI[i][j], &matriceI[i + 1][j]);
                        }
                    }
                    #if DEBUG
                    printf("Apres echange\n\n");
                    printf("\nMatrice A\n");
                    afficheMatrice2D(tailleMatrice, matriceA);
                    printf("\nMatrice I\n");
                    afficheMatrice2D(tailleMatrice, matriceI);
                    printf("\n");
                    system("pause");
                    #endif
                }
            } while(pivot.re == 0.0 && pivot.im == 0.0);
            #if DEBUG
            printf("---------------\nPivot different de zero trouve : ");
            afficheComplexe(pivot.re, pivot.im);
            printf("\n");
            #endif

            for(j = 0; j < tailleMatrice; j++)
            {
                matriceA[i][j] = divComplexe(matriceA[i][j], pivot);
                matriceI[i][j] = divComplexe(matriceI[i][j], pivot);
            }

            for(k = 0; k < tailleMatrice; k++)
            {
                if(k != i)
                {
                    coef.re = matriceA[k][i].re;
                    coef.im = matriceA[k][i].im;

                    for(j = 0; j < tailleMatrice; j++)
                    {
                        temp = produitComplexe(coef, matriceA[i][j]);
                        matriceA[k][j].re = matriceA[k][j].re - temp.re;
                        matriceA[k][j].im = matriceA[k][j].im - temp.im;

                        temp = produitComplexe(coef, matriceI[i][j]);
                        matriceI[k][j].re = matriceI[k][j].re - temp.re;
                        matriceI[k][j].im = matriceI[k][j].im - temp.im;
                    }
                }
            }

            #if DEBUG
            printf("\nMatrice A\n");
            afficheMatrice2D(tailleMatrice, matriceA);
            printf("\nMatrice I\n");
            afficheMatrice2D(tailleMatrice, matriceI);
            printf("\n");
            system("pause");
            #endif
        }

        copie2Dto1D(tailleMatrice, matriceI, resultante);
        freeMatrice(&matriceIdentite);
        freeMatrice2D(tailleMatrice, &matriceA);
        freeMatrice2D(tailleMatrice, &matriceI);
        */
    }
    else
    {
        erreur = 0; // Matrice non inversible
    }
    return erreur;
}

/* Calcul des coefficients du polynome caractéristique */
void polyCaract(int degree, Matrice matrice, double *Poly)
{
    if(degree == 2)
    {
        Poly[2] = 1.0;
        Poly[1] = -matTrace(degree, matrice).re;
        Poly[0] = matDet(degree, matrice).re;
    }
    else if(degree == 3)
    {
        Poly[3] = 1.0;
        Poly[2] = -matTrace(degree, matrice).re;
        Poly[1] = matrice.coefficient[0].re * matrice.coefficient[4].re + matrice.coefficient[0].re * matrice.coefficient[8].re - matrice.coefficient[1].re * matrice.coefficient[3].re - matrice.coefficient[2].re * matrice.coefficient[6].re + matrice.coefficient[4].re * matrice.coefficient[8].re - matrice.coefficient[5].re * matrice.coefficient[7].re;
        Poly[0] = -matDet(degree, matrice).re;
    }
    else if(degree == 4)
    {
        Poly[4] = 1.0;
        Poly[3] = -matTrace(degree, matrice).re;
        Poly[2] = -matrice.coefficient[1].re * matrice.coefficient[4].re + matrice.coefficient[0].re * matrice.coefficient[5].re - matrice.coefficient[2].re * matrice.coefficient[8].re - matrice.coefficient[6].re * matrice.coefficient[9].re + matrice.coefficient[0].re * matrice.coefficient[10].re + matrice.coefficient[5].re * matrice.coefficient[10].re - matrice.coefficient[3].re * matrice.coefficient[12].re - matrice.coefficient[7].re * matrice.coefficient[13].re - matrice.coefficient[11].re * matrice.coefficient[14].re + matrice.coefficient[0].re * matrice.coefficient[15].re + matrice.coefficient[5].re * matrice.coefficient[15].re + matrice.coefficient[10].re * matrice.coefficient[15].re;
        Poly[1] = matrice.coefficient[2].re * matrice.coefficient[5].re * matrice.coefficient[8].re - matrice.coefficient[1].re * matrice.coefficient[6].re * matrice.coefficient[8].re - matrice.coefficient[2].re * matrice.coefficient[4].re * matrice.coefficient[9].re + matrice.coefficient[0].re * matrice.coefficient[6].re * matrice.coefficient[9].re + matrice.coefficient[1].re * matrice.coefficient[4].re * matrice.coefficient[10].re - matrice.coefficient[0].re * matrice.coefficient[5].re * matrice.coefficient[10].re + matrice.coefficient[3].re * matrice.coefficient[5].re * matrice.coefficient[12].re - matrice.coefficient[1].re * matrice.coefficient[7].re * matrice.coefficient[12].re + matrice.coefficient[3].re * matrice.coefficient[10].re * matrice.coefficient[12].re - matrice.coefficient[2].re * matrice.coefficient[11].re * matrice.coefficient[12].re - matrice.coefficient[3].re * matrice.coefficient[4].re * matrice.coefficient[13].re + matrice.coefficient[0].re * matrice.coefficient[7].re * matrice.coefficient[13].re + matrice.coefficient[7].re * matrice.coefficient[10].re * matrice.coefficient[13].re - matrice.coefficient[6].re * matrice.coefficient[11].re * matrice.coefficient[13].re - matrice.coefficient[3].re * matrice.coefficient[8].re * matrice.coefficient[14].re - matrice.coefficient[7].re * matrice.coefficient[9].re * matrice.coefficient[14].re + matrice.coefficient[0].re * matrice.coefficient[11].re * matrice.coefficient[14].re + matrice.coefficient[5].re * matrice.coefficient[11].re * matrice.coefficient[14].re + matrice.coefficient[1].re * matrice.coefficient[4].re * matrice.coefficient[15].re - matrice.coefficient[0].re * matrice.coefficient[5].re * matrice.coefficient[15].re + matrice.coefficient[2].re * matrice.coefficient[8].re * matrice.coefficient[15].re + matrice.coefficient[6].re * matrice.coefficient[9].re * matrice.coefficient[15].re - matrice.coefficient[0].re * matrice.coefficient[10].re * matrice.coefficient[15].re - matrice.coefficient[5].re * matrice.coefficient[10].re * matrice.coefficient[15].re;
        Poly[0] = matDet(degree, matrice).re;
    }
    else
    {
        // Méthode générale à faire si on a le temps
    }

    affichePoly('x', Poly, degree);
}

/* Affichage du polynome caractéristique */
void affichePoly(char variable, double *poly, int degree)
{
    int n;
    for(n = degree; n > 0; n--)
    {
        if(poly[n] != 0)
        {
            if(poly[n] == 1.0)
            {
                printf("%c", variable);
            }
            else if(poly[n] == -1.0)
            {
                printf("-%c", variable);
            }
            else
            {
                printf("%g%c", poly[n], variable);
            }

            if(n > 1)
            {
                printf("^%i", n);
            }
        }

        if(poly[n - 1] > 0)
        {
            printf("+");
        }
    }

    if(poly[0] == 1.0)
    {
        printf("1");
    }
    else if(poly[0] == -1.0)
    {
        printf("-1");
    }
    else if(poly[0] != 0.0)
    {
        printf("%g", poly[0]);
    }

    return;
}

/* Calcul et Affichage des valeurs propres */
void valeurPropre(int degree, double *poly, t_racine *racines)
{
    int i;
    if(degree == 4)
    {
        degre4(racines, poly[4], poly[3], poly[2], poly[1], poly[0], 0);
    }
    else if(degree == 3)
    {
        degre3(racines, poly[3], poly[2], poly[1], poly[0], 0);
    }
    else if(degree == 2)
    {
        degre2(racines, poly[2], poly[1], poly[0], 0.0, 0);
    }

    for(i = 0; i < degree; i++)
    {
        printf(" x%i = ", i + 1);
        afficheComplexe(racines[i].re, racines[i].im);
        printf("\n");
    }
}

/* Calcul des vecteurs propres */
void vecteurPropre(int tailleMatrice, Matrice matrice, t_complexe inconnue, t_vecteur *vecteur, int nVect = 0)
{
    int i, j, k;
    t_complexe **matriceA, pivot, coef, temp;

    creatMatrice2D(tailleMatrice, &matriceA);
    copie1DtoD2(tailleMatrice, matrice, matriceA);

    for (i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            if(j == i)
            {
                matriceA[i][j].re = inconnue.re - matriceA[i][j].re;
                matriceA[i][j].im = inconnue.im - matriceA[i][j].im;
            }
            else
            {
                matriceA[i][j].re = -matriceA[i][j].re;
                matriceA[i][j].im = -matriceA[i][j].im;
            }
        }
    }

    if(tailleMatrice == 2)
    {
        if(matriceA[0][1].re != 0.0 || matriceA[0][1].im != 0.0)
        {
            matriceA[0][1].re = -matriceA[0][1].re;
            matriceA[0][1].im = -matriceA[0][1].im;
            vecteur[nVect].n[0] = divComplexe(matriceA[0][1], matriceA[0][0]);
        }
        else
        {
            matriceA[1][1].re = -matriceA[1][1].re;
            matriceA[1][1].im = -matriceA[1][1].im;
            vecteur[nVect].n[0] = divComplexe(matriceA[1][1], matriceA[1][0]);
        }
    }
    else
    {
        for(i = 0; i < (tailleMatrice - 1); i++)
        {
            do
            {
                pivot.re = matriceA[i][i].re;
                pivot.im = matriceA[i][i].im;

                if(pivot.re == 0.0 && pivot.im == 0.0)
                {
                    for(j = 0; j < tailleMatrice; j++)
                    {
                        if(i + 1 < tailleMatrice)
                        {
                            swap(&matriceA[i][j], &matriceA[i + 1][j]);
                        }
                    }
                    #if DEBUG
                    printf("Apres echange\n\n");
                    printf("\nMatrice A\n");
                    afficheMatrice2D(tailleMatrice, matriceA);
                    printf("\n");
                    system("pause");
                    #endif
                }
            } while(pivot.re == 0.0 && pivot.im == 0.0);
            #if DEBUG
            printf("---------------\nPivot different de zero trouve : ");
            afficheComplexe(pivot.re, pivot.im);
            printf("\n");
            #endif

            for(j = 0; j < tailleMatrice; j++)
            {
                matriceA[i][j] = divComplexe(matriceA[i][j], pivot);
            }

            for(k = 0; k < tailleMatrice; k++)
            {
                if(k != i)
                {
                    coef.re = matriceA[k][i].re;
                    coef.im = matriceA[k][i].im;

                    for(j = 0; j < tailleMatrice; j++)
                    {
                        temp = produitComplexe(coef, matriceA[i][j]);
                        matriceA[k][j].re = matriceA[k][j].re - temp.re;
                        matriceA[k][j].im = matriceA[k][j].im - temp.im;
                    }
                }
            }
            #if DEBUG
            printf("\nMatrice A\n");
            afficheMatrice2D(tailleMatrice, matriceA);
            printf("\n");
            system("pause");
            #endif
        }

        #if DEBUG
        printf("V = (");
        #endif
        for(i = 0; i < (tailleMatrice - 1); i++)
        {
            vecteur[nVect].n[i].re = -matriceA[i][(tailleMatrice - 1)].re;
            vecteur[nVect].n[i].im = -matriceA[i][(tailleMatrice - 1)].im;

            #if DEBUG
            afficheComplexe(vecteur[nVect].n[i].re, vecteur[nVect].n[i].im);
            printf(", ");
            #endif
        }
        #if DEBUG
        printf("1)\n\n");
        #endif
    }

    vecteur[nVect].n[tailleMatrice - 1].re = 1.0;
    vecteur[nVect].n[tailleMatrice - 1].im = 0.0;

    freeMatrice2D(tailleMatrice, &matriceA);
}

/* Produit de nombres complexes */
t_complexe produitComplexe(t_complexe A, t_complexe B)
{
    t_complexe C;

    C.re = A.re * B.re - A.im * B.im;
    C.im = A.re * B.im + A.im * B.re;

    return C;
}

/* Division de nombres complexes */
t_complexe divComplexe(t_complexe A, t_complexe B)
{
    t_complexe C;
    double temp;
    temp = B.re * B.re + B.im * B.im;

    C.re = (A.re * B.re + A.im * B.im) / temp;
    C.im = (A.im * B.re - A.re * B.im) / temp;
    return C;
}

/* Fonction échange */
void swap(t_complexe *x, t_complexe *y)
{
    t_complexe temp;
    temp = *x;
    *x = *y;
    *y = temp;
}

/* Multiplication de matrice par un scalaire */
void matMutScalaire(int tailleMatrice, Matrice matrice, t_complexe scalaire, Matrice *resultante)
{
    int i;

    // Multiplication
    for(i = 0 ; i < tailleMatrice * tailleMatrice; i++)
    {
        resultante->coefficient[i].re = scalaire.re * matrice.coefficient[i].re;
        resultante->coefficient[i].im = scalaire.im * matrice.coefficient[i].im;
    }
}

/* Addition de matrices */
void matAdd(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice* resultante)
{
    int i = 0;

    for(i = 0 ; i < tailleMatrice*tailleMatrice ; i++)
    {
        resultante->coefficient[i].re = matrice1.coefficient[i].re + matrice2.coefficient[i].re;
        resultante->coefficient[i].im = matrice1.coefficient[i].im + matrice2.coefficient[i].im;
    }
}

/* Soustraction de matrices */
void matDiff(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice* resultante)
{
    int i = 0;

    for(i = 0 ; i < tailleMatrice*tailleMatrice ; i++)
    {
        resultante->coefficient[i].re = matrice1.coefficient[i].re - matrice2.coefficient[i].re;
        resultante->coefficient[i].im = matrice1.coefficient[i].im - matrice2.coefficient[i].im;
    }
}

/* Produit de matrices */
void matMut(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice *resultante)
{
    int lig = 0, col = 0, lig1 = 0, col1 = 0, lig2 = 0, col2 = 0;
    t_complexe **coeffMatrice1, **coeffMatrice2, **coeffResultante, temp;

    creatMatrice2D(tailleMatrice, &coeffMatrice1);
    creatMatrice2D(tailleMatrice, &coeffMatrice2);
    creatMatrice2D(tailleMatrice, &coeffResultante);

    #if DEBUG
    printf("Le produit est :\n\n");
    afficheMatrice(tailleMatrice, matrice1);
    printf("\n\net :\n\n");
    afficheMatrice(tailleMatrice, matrice2);
    printf("\n\n");
    system("pause");
    #endif

    // Copie 1D->2D
    copie1DtoD2(tailleMatrice, matrice1, coeffMatrice1);
    copie1DtoD2(tailleMatrice, matrice2, coeffMatrice2);

    // Multiplication
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            coeffResultante[lig][col].re = 0.0;
            coeffResultante[lig][col].im = 0.0;
            lig1 = lig;
            col1 = 0;
            col2 = col;
            lig2 = 0;
            do
            {
                temp = produitComplexe(coeffMatrice1[lig1][col1], coeffMatrice2[lig2][col2]);

                #if DEBUG
                printf("Produit Complexe :\n%g | %g\n%g | %g\n\n", coeffMatrice1[lig1][col1].re, coeffMatrice1[lig1][col1].im, coeffMatrice2[lig2][col2].re, coeffMatrice2[lig2][col2].im);
                printf("Coef[%d][%d].re = %g\n", lig, col, coeffResultante[lig][col].re);
                printf("Coef[%d][%d].im = %g\n", lig, col, coeffResultante[lig][col].im);

                printf("Temp[%d][%d].re = %g\n", lig, col, temp.re);
                printf("Temp[%d][%d].im = %g\n", lig, col, temp.im);
                #endif

                coeffResultante[lig][col].re += temp.re;
                coeffResultante[lig][col].im += temp.im;

                #if DEBUG
                printf("result[%d][%d].re = %g\n", lig, col, coeffResultante[lig][col].re);
                printf("result[%d][%d].im = %g\n", lig, col, coeffResultante[lig][col].im);
                system("pause");
                #endif

                col1++;
                lig2++;
            } while(col1 < tailleMatrice);
        }
    }

    // Copie 2D->1D
    copie2Dto1D(tailleMatrice, coeffResultante, resultante);

    /* Libération de la mémoire */
    freeMatrice2D(tailleMatrice, &coeffMatrice1);
    freeMatrice2D(tailleMatrice, &coeffMatrice2);
    freeMatrice2D(tailleMatrice, &coeffResultante);
}

/* Matrice à la puissance k */
void matPuissance(int tailleMatrice, Matrice matrice, int n, Matrice *resultante)
{
    int i;
    Matrice matriceA, matriceB, matriceC;

    creatMatrice(tailleMatrice, &matriceA);
    creatMatrice(tailleMatrice, &matriceB);
    creatMatrice(tailleMatrice, &matriceC);

    if(n < 0)
    {
        matInverse(tailleMatrice, matrice, &matriceA);
        n = -n;
    }
    else
    {
        copie(tailleMatrice, matrice, matriceA);
    }
    copie(tailleMatrice, matriceA, matriceB);
    copie(tailleMatrice, matriceB, matriceC);

    if(n != 0)
    {
        for(i = 0; i < (n - 1); i++)
        {
            matMut(tailleMatrice, matriceA, matriceB, &matriceC);
            copie(tailleMatrice, matriceC, matriceB);
        }
        copie(tailleMatrice, matriceC, *resultante);
    }
    else
    {
        matIdent(tailleMatrice, matrice, resultante);
    }

    freeMatrice(&matriceA);
    freeMatrice(&matriceB);
    freeMatrice(&matriceC);
}

/* Exponentielle d'une matrice */
void matExp(int tailleMatrice, Matrice matrice, Matrice *resultante)
{
    int i;

    for(i = 0 ; i < tailleMatrice * tailleMatrice; i++)
    {
        resultante->coefficient[i].re = puissance(EXP, matrice.coefficient[i].re);
        if(matrice.coefficient[i].im != 0.0)
        {
            resultante->coefficient[i].im = puissance(EXP, matrice.coefficient[i].im);
        }
        else
        {
            resultante->coefficient[i].im = 0.0;
        }
    }
}

/* Calcul du mineur d'une matrice */
void matMin(int tailleMatrice, Matrice matrice, Matrice *resultante)
{
    int i, j, k, l, m = 0, n = 0;
    Matrice sousMatrice;
    t_complexe **matrice2D, **sousMatrice2D, **resultante2D;

    if(tailleMatrice == 2)
    {
        resultante->coefficient[0] = matrice.coefficient[3];
        resultante->coefficient[1] = matrice.coefficient[2];
        resultante->coefficient[2] = matrice.coefficient[1];
        resultante->coefficient[3] = matrice.coefficient[0];
    }
    else
    {
        creatMatrice2D(tailleMatrice, &matrice2D);
        copie1DtoD2(tailleMatrice, matrice, matrice2D);

        creatMatrice2D(tailleMatrice, &resultante2D);

        for(i = 0; i < tailleMatrice; i++)
        {
            for(j = 0; j < tailleMatrice; j++)
            {
                m = 0;
                creatMatrice((tailleMatrice - 1), &sousMatrice);

                creatMatrice2D((tailleMatrice - 1), &sousMatrice2D);

                for(k = 0; k < tailleMatrice; k++)
                {
                    n = 0;
                    for(l = 0; l < tailleMatrice; l++)
                    {
                        if(i != k && j != l)
                        {
                            sousMatrice2D[m][n].re = matrice2D[k][l].re;
                            sousMatrice2D[m][n].im = matrice2D[k][l].im;
                            #if DEBUG
                            printf("i = %i\nj = %i\nk = %i\nl = %i\nm = %i\nn = %i\n\n", i, j, k, l, m, n);
                            #endif
                            n++;
                        }
                    }
                    if(i != k)
                    {
                        m++;
                    }
                }
                copie2Dto1D((tailleMatrice - 1), sousMatrice2D, &sousMatrice);

                resultante2D[i][j] = matDet(tailleMatrice - 1, sousMatrice);

                freeMatrice2D((tailleMatrice - 1), &sousMatrice2D);
                freeMatrice(&sousMatrice);
            }
        }
        copie2Dto1D(tailleMatrice, resultante2D, resultante);

        freeMatrice2D(tailleMatrice, &matrice2D);
        freeMatrice2D(tailleMatrice, &resultante2D);
    }
}

/* Pair/Impair dans les coefficients pour générer la comatrice */
void matCom(int tailleMatrice, Matrice matrice, Matrice *resultante)
{
    int i, j;
    t_complexe **M, **R;

    creatMatrice2D(tailleMatrice, &M);
    creatMatrice2D(tailleMatrice, &R);

    copie1DtoD2(tailleMatrice, matrice, M);

    for(i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            R[i][j].re = puissance(-1.0, (i + 1) + (j + 1)) * M[i][j].re;
            R[i][j].im = puissance(-1.0, (i + 1) + (j + 1)) * M[i][j].im;
        }
    }

    copie2Dto1D(tailleMatrice, R, resultante);

    freeMatrice2D(tailleMatrice, &M);
    freeMatrice2D(tailleMatrice, &R);
}

/* Calcul du déterminant d'une matrice */
t_complexe matDet(int tailleMatrice, Matrice matrice)
{
    int i, j, k = 0;
    t_complexe det, coef, temp;
    det.re = 0.0;
    det.im = 0.0;
    coef.re = 0.0;
    coef.im = 0.0;

    Matrice sousMatrice;
    if(tailleMatrice == 2)
    {
        // det.re = (matrice.coefficient[0].re * matrice.coefficient[3].re) - (matrice.coefficient[1].re * matrice.coefficient[2].re); // Old version
        det.re = produitComplexe(matrice.coefficient[0], matrice.coefficient[3]).re - produitComplexe(matrice.coefficient[1], matrice.coefficient[2]).re;
        det.im = produitComplexe(matrice.coefficient[0], matrice.coefficient[3]).im - produitComplexe(matrice.coefficient[1], matrice.coefficient[2]).im;
    }
    else
    {
        for(i = 0; i < tailleMatrice; i++)
        {
            sousMatrice.coefficient = (t_complexe*)calloc((tailleMatrice - 1)*(tailleMatrice - 1), sizeof(t_complexe));
            k = 0;
            for(j = tailleMatrice; j < tailleMatrice * tailleMatrice; j++)
            {
                if(j % tailleMatrice != i)
                {
                    sousMatrice.coefficient[k].re = matrice.coefficient[j].re;
                    sousMatrice.coefficient[k].im = matrice.coefficient[j].im;
                    #if DEBUG
                    printf("\nCopie : Matrice[%d] => SousMatrice[%d] = ", j, k);
                    afficheComplexe(matrice.coefficient[j].re, matrice.coefficient[j].im);
                    #endif
                    k++;
                }
            }

            coef.re = matrice.coefficient[i].re;
            coef.im = matrice.coefficient[i].im;

            if(i % 2 != 0)
            {
                coef.re = -coef.re;
                coef.im = -coef.im;
            }

            #if DEBUG
            printf("\n\n Sous Matrice %d\n\n", i + 1);
            afficheMatrice(tailleMatrice - 1, sousMatrice);
            printf("\n");
            #endif

            temp = produitComplexe(coef, matDet(tailleMatrice - 1, sousMatrice));
            det.re = det.re + temp.re;
            det.im = det.im + temp.im;
            free(sousMatrice.coefficient);
        }
    }
    return det;
}

/* Matrice identité */
void matIdent(int tailleMatrice, Matrice matrice, Matrice *resultante)
{
    int i = 0;

    for (i = 0 ; i < (tailleMatrice * tailleMatrice) ; i += tailleMatrice + 1)
    {
        resultante->coefficient[i].re = 1.0;
        resultante->coefficient[i].im = 0.0;
    }
}

/* Copie d'une matrice 1D vers une matrice 2D */
void copie1DtoD2(int tailleMatrice, Matrice matrice, t_complexe **resultante)
{
    int i = 0, lig = 0, col = 0;
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            resultante[lig][col].re = matrice.coefficient[i].re;
            resultante[lig][col].im = matrice.coefficient[i].im;
            i++;
        }
    }
}

/* Copie d'une matrice 2D vers une matrice 1D */
void copie2Dto1D(int tailleMatrice, t_complexe **matrice, Matrice *resultante)
{
    int i = 0, lig = 0, col = 0;
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            resultante->coefficient[i].re = matrice[lig][col].re;
            resultante->coefficient[i].im = matrice[lig][col].im;
            i++;
        }
    }
}

/* Copie d'une matrice 1D */
void copie(int tailleMatrice, Matrice matrice1, Matrice matrice2)
{
    int i = 0;
    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        matrice2.coefficient[i].re = matrice1.coefficient[i].re;
        matrice2.coefficient[i].im = matrice1.coefficient[i].im;
    }
}

/* Diagonalisation */
void diagonalisation(int tailleMatrice, t_racine *racine, t_vecteur * vecteur, Matrice *resultante)
{
    int i, j, nRacine = 0;
    Matrice D, P, P_1, Temp;
    t_complexe **P_temp;

    creatMatrice2D(tailleMatrice, &P_temp);

    creatMatrice(tailleMatrice, &D);
    creatMatrice(tailleMatrice, &P);
    creatMatrice(tailleMatrice, &P_1);
    creatMatrice(tailleMatrice, &Temp);

    for(i = 0; i < (tailleMatrice * tailleMatrice); i += tailleMatrice + 1)
    {
        D.coefficient[i].re = racine[nRacine].re;
        D.coefficient[i].im = racine[nRacine].im;
        nRacine += 1;
    }

    for(i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            P_temp[j][i] = vecteur[i].n[j];
        }
    }

    copie2Dto1D(tailleMatrice, P_temp, &P);

    printf("Matrice P\n\n");
    afficheMatrice(tailleMatrice, P);
    printf("\n\n");
    verifierWolfram(tailleMatrice, P);

    printf("Matrice D\n\n");
    afficheMatrice(tailleMatrice, D);
    printf("\n\n");
    verifierWolfram(tailleMatrice, D);

    copie(tailleMatrice, P, Temp);
    matInverse(tailleMatrice, Temp, &P_1);

    printf("Matrice P-1\n\n");
    afficheMatrice(tailleMatrice, P_1);
    printf("\n\n");
    verifierWolfram(tailleMatrice, P_1);

    matMut(tailleMatrice, P, D, &Temp);
    printf("Matrice PxD\n\n");
    afficheMatrice(tailleMatrice, Temp);
    printf("\n\n");
    verifierWolfram(tailleMatrice, Temp);

    matMut(tailleMatrice, Temp, P_1, resultante);
    printf("Matrice A = PDP-1\n\n");
    afficheMatrice(tailleMatrice, *resultante);
    printf("\n\n");
    verifierWolfram(tailleMatrice, *resultante);

    freeMatrice(&D);
    freeMatrice(&P);
    freeMatrice(&P_1);
    freeMatrice(&Temp);
    freeMatrice2D(tailleMatrice, &P_temp);
}

/* Création d'une matrice 1D */
void creatMatrice(int tailleMatrice, Matrice *matrice)
{
    matrice->coefficient = (t_complexe*) calloc(tailleMatrice * tailleMatrice, sizeof(t_complexe));
}

/* Libération d'une matrice 1D */
void freeMatrice(Matrice *matrice)
{
    free(matrice->coefficient);
}

/* Création d'une matrice 2D */
void creatMatrice2D(int tailleMatrice, t_complexe ***matrice)
{
    int i;
    (*matrice) = (t_complexe**) calloc(tailleMatrice, sizeof(t_complexe));
    for (i = 0 ; i < tailleMatrice ; i++)
    {
        (*matrice)[i] = (t_complexe*) calloc(tailleMatrice, sizeof(t_complexe));
    }
}

/* Libération d'une matrice 2D */
void freeMatrice2D(int tailleMatrice, t_complexe ***matrice)
{
    int i;
    for (i = 0 ; i < tailleMatrice ; i++)
    {
        free((*matrice)[i]);
    }
    free(*matrice);
}

/* Affichage d'une matrice 1D */
void afficheMatrice(int tailleMatrice, Matrice matrice)
{
    int i, X, x, y;

    printf("   ");
    X = wherex();
    x = X;
    y = wherey();

    for (i = 0 ; i < tailleMatrice * tailleMatrice; i++)
    {
        gotoxy(x, y);
        afficheComplexe(matrice.coefficient[i].re, matrice.coefficient[i].im);
        if((i + 1) % tailleMatrice != 0)
        {
            x = x + DEF_X;
        }
        else
        {
            x = X;
            y = y + 2;
        }
    }
}

/* Affichage d'une matrice 2D */
void afficheMatrice2D(int tailleMatrice, t_complexe **matrice)
{
    int i, j;
    for (i = 0 ; i < tailleMatrice; i++)
    {
        for (j = 0 ; j < tailleMatrice; j++)
        {
            printf("   ");
            afficheComplexe(matrice[i][j].re, matrice[i][j].im);
            printf("   ");
        }
        printf("\n");
    }
}

/* Affichage finale */
void affichageTotal(int tailleMatrice, int pasUneMatrice, int deuxMatrices, int inverse, int mineur, Matrice matrice1, Matrice matrice2, Matrice resultante)
{
    int i, j;
    t_complexe trace, determinant;
    t_vecteur *vecteur;
    double *Poly;
    t_racine *racine;
    Poly = (double*) calloc(tailleMatrice + 1, sizeof(double));
    racine = (t_racine*) calloc(tailleMatrice, sizeof(t_racine));
    vecteur = (t_vecteur*) calloc(tailleMatrice, sizeof(t_vecteur));
    int nbComplexe = 0;

    for(i = 0; i < tailleMatrice; i++)
    {
        vecteur[i].n = (t_complexe*) calloc(tailleMatrice, sizeof(t_complexe));
    }

    for(i = 0; i < tailleMatrice; i++)
    {
        racine[i].im = 0.0;
    }
    system("cls");

    printf("\n   Matrice :\n");
    afficheMatrice(tailleMatrice, matrice1);
    if(deuxMatrices == 1)
    {
        printf("\n\n   Matrice secondaire :\n");
        afficheMatrice(tailleMatrice, matrice2);
    }
    if(mineur == 1)
    {
        matMin(tailleMatrice, matrice1, &resultante);

        printf("\n\n   Mineur :\n");
        afficheMatrice(tailleMatrice, resultante);

        matCom(tailleMatrice, resultante, &resultante);

        printf("\n\n   Comatrice :\n");
        afficheMatrice(tailleMatrice, resultante);
    }
    else if(pasUneMatrice == 0)
    {
        if(inverse == 1)
        {
            printf("\n\n\n\n\n   Matrice resultante :\n");
            afficheMatrice(tailleMatrice, resultante);
        }
        else
        {
            printf("\n\n\n\n\n   La matrice n'est pas inversible.\n");
        }
    }
    else if (pasUneMatrice == 1)
    {
        determinant = matDet(tailleMatrice, matrice1);
        trace = matTrace(tailleMatrice, matrice1);
        printf("\n\n Determinant : ");
        afficheComplexe(determinant.re, determinant.im);
        printf("\n\n Trace : ");
        afficheComplexe(trace.re, trace.im);
        printf("\n\n Polynome caracteristique : ");
        polyCaract(tailleMatrice, matrice1, Poly);
        printf("\n\n Valeurs propres :\n");
        valeurPropre(tailleMatrice, Poly, racine);

        printf("\n\n Vecteurs propres :\n");
        for(i = 0; i < tailleMatrice; i++)
        {
            vecteurPropre(tailleMatrice, matrice1, racine[i], vecteur, i);
        }

        for(i = 0; i < tailleMatrice; i++)
        {
            printf(" V%i = (", i + 1);
            for(j = 0; j < (tailleMatrice - 1 ); j++)
            {
                afficheComplexe(vecteur[i].n[j].re, vecteur[i].n[j].im);
                printf(", ");

                /* Utilisation des complexes */
                if(vecteur[i].n[j].im != 0)
                {
                    nbComplexe = 1;
                }
            }
            afficheComplexe(vecteur[i].n[tailleMatrice - 1].re, vecteur[i].n[tailleMatrice - 1].im);
            printf(")\n");
        }
        printf("\nProchain affichage : Matrices P, D, P-1, PxD\n");
        system("pause");

        if(nbComplexe == 1)
        {
            diagonalisation(tailleMatrice, racine, vecteur, &resultante);
            printf("\n\nOn retrouve la Matrice A :\n");
            afficheMatrice(tailleMatrice, resultante);
        }
        else
        {
            diagonalisationPuissanceN(tailleMatrice, racine, vecteur);
        }
    }
    printf("\n\n\n");

    verifierWolfram(tailleMatrice, matrice1);

    if(deuxMatrices == 1)
    {
        verifierWolfram(tailleMatrice, matrice2);
    }

    printf("Appuyez sur une touche pour revenir au menu");
    getch();
    free(Poly);
    free(racine);
    free(vecteur);
}

/* Vérification facile par un copier coller sur wolfram alpha */
void verifierWolfram(int tailleMatrice, Matrice matrice)
{
    int i;
    printf("Verifier sur wolfram alpha :\n\n");

    printf("\t{{");
    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        afficheComplexe(matrice.coefficient[i].re, matrice.coefficient[i].im);
        if(i < (tailleMatrice * tailleMatrice - 1))
        {
            if((i + 1) % tailleMatrice == 0)
            {
                printf("}, {");
            }
            else
            {
                printf(", ");
            }
        }
    }
    printf("}}\n\n\n");
}
