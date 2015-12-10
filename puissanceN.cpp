#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "puissanceN.h"
#include "fonctions.h"
#include "conio.h"

/* Création d'une matrice de type puissance N */
void creatMatriceN(int tailleMatrice, t_matriceN *matrice, int nb = 0)
{
    int i;

    matrice->mat = (t_puissanceN*) calloc((tailleMatrice * tailleMatrice), sizeof(t_puissanceN));

    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        matrice->mat[i].Case = (t_case*) calloc(5, sizeof(t_case));
    }
}

/* Libération de matrice de type puissance N */
void freeMatriceN(int tailleMatrice, t_matriceN *matrice)
{
    int i;
    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        free(matrice->mat[i].Case);
    }
    free(matrice->mat);
}

/* Création d'une matrice 2D de type puissance N */
void creatMatrice2DN(int tailleMatrice, t_puissanceN ***matrice, int nb = 0)
{
    int i, j;
    (*matrice) = (t_puissanceN**) calloc(tailleMatrice, sizeof(t_puissanceN));
    for (i = 0 ; i < tailleMatrice; i++)
    {
        (*matrice)[i] = (t_puissanceN*) calloc(tailleMatrice, sizeof(t_puissanceN));
    }

    for(i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            (*matrice)[i][j].Case = (t_case*) calloc(5, sizeof(t_case));
        }
    }
}

/* Libération de matrice 2D de type puissance N */
void freeMatrice2DN(int tailleMatrice, t_puissanceN ***matrice)
{
    int i, j;

    for(i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            free((*matrice)[i][j].Case);
        }
    }

    for (i = 0 ; i < tailleMatrice; i++)
    {
        free((*matrice)[i]);
    }
    free(*matrice);
}

/* Ajouter un coefficient + un coefficient à la puissance n */
void AddPuissanceN(int tailleMatrice, int nb, t_matriceN *matrice, int i, double coef, double n)
{
    matrice->mat[i].Case[nb].coef = coef;
    matrice->mat[i].Case[nb].n = n;
    matrice->mat[i].nb += 1;
}

/* Diagonalisation et matrice à la puissance N */
void diagonalisationPuissanceN(int tailleMatrice, t_racine *racine, t_vecteur * vecteur)
{
    int i, j, nRacine = 0;
    Matrice Pmat, Pmat_1;
    t_complexe **P_temp;
    t_matriceN P, P_1, D, resultante, Temp;
    t_puissanceN **coeffMatrice1, **coeffMatrice2, **coeffResultante, Temp2;

    // Création : Matrice 2D à la puissance N
    creatMatrice2DN(tailleMatrice, &coeffMatrice1);
    creatMatrice2DN(tailleMatrice, &coeffMatrice2);
    creatMatrice2DN(tailleMatrice, &coeffResultante);

    Temp2.Case = (t_case*) calloc(5, sizeof(t_case));

    // Création : Matrice à la puissance N
    creatMatriceN(tailleMatrice, &P, 0);
    creatMatriceN(tailleMatrice, &P_1, 0);
    creatMatriceN(tailleMatrice, &D, 0);
    creatMatriceN(tailleMatrice, &resultante, 0);
    creatMatriceN(tailleMatrice, &Temp, 0);

    // Création : Matrice 2D et 1D
    creatMatrice2D(tailleMatrice, &P_temp);
    creatMatrice(tailleMatrice, &Pmat);
    creatMatrice(tailleMatrice, &Pmat_1);

    // Place les valeurs propres dans la matrice D
    for(i = 0; i < (tailleMatrice * tailleMatrice); i += tailleMatrice + 1)
    {
        AddPuissanceN(tailleMatrice, 0, &D, i, 1.0, racine[nRacine].re);
        nRacine += 1;
    }

    // Place les vecteurs propres dans la matrice P
    for(i = 0; i < tailleMatrice; i++)
    {
        for(j = 0; j < tailleMatrice; j++)
        {
            P_temp[j][i] = vecteur[i].n[j];
        }
    }

    copie2Dto1D(tailleMatrice, P_temp, &Pmat);

    printf("\n\nMatrice P\n\n");
    afficheMatrice(tailleMatrice, Pmat);
    printf("\n\n");

    printf("Matrice D\n\n");
    afficheMatricePuissanceN(tailleMatrice, D);
    printf("\n\n");

    // Calcul de la matrice P-1
    matInverse(tailleMatrice, Pmat, &Pmat_1);

    printf("Matrice P-1\n\n");
    afficheMatrice(tailleMatrice, Pmat_1);
    printf("\n\n");

    copieToN(tailleMatrice, Pmat, &P);
    copieToN(tailleMatrice, Pmat_1, &P_1);

    // Calcul de la matrice PxD
    matMutPuissanceN(tailleMatrice, coeffMatrice1, coeffMatrice2, coeffResultante, Temp2, P, D, &Temp);

    printf("Matrice PxD\n\n");
    afficheMatricePuissanceN(tailleMatrice, Temp);
    printf("\n\n");

    // Calcul de la matrice A^n = P D^n P-1
    matMutPuissanceN(tailleMatrice, coeffMatrice1, coeffMatrice2, coeffResultante, Temp2, Temp, P_1, &resultante);

    printf("Prochain affichage : Matrice A^n. Appuyez sur une touche pour continuer..");
    getch();

    // Redimentionnement pour l'affichage de la matrice puissance N
    //system("mode con cols=100 lines=25");

    printf("Matrice A^n\n\n");
    afficheMatricePuissanceN(tailleMatrice, resultante);
    printf("\n\n");

    /* Libération de la mémoire */
    free(Temp2.Case);

    freeMatrice2D(tailleMatrice, &P_temp);

    freeMatrice(&Pmat);
    freeMatrice(&Pmat_1);

    freeMatrice2DN(tailleMatrice, &coeffMatrice1);
    freeMatrice2DN(tailleMatrice, &coeffMatrice2);
    freeMatrice2DN(tailleMatrice, &coeffResultante);

    freeMatriceN(tailleMatrice, &P_1);
    freeMatriceN(tailleMatrice, &resultante);
    freeMatriceN(tailleMatrice, &P);
    freeMatriceN(tailleMatrice, &D);
    freeMatriceN(tailleMatrice, &Temp);
}

/* Copie une matrice 1D dans une matrice de type puissance N */
void copieToN(int tailleMatrice, Matrice Pmat, t_matriceN *P)
{
    int i;

    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        AddPuissanceN(tailleMatrice, 0, P, i, Pmat.coefficient[i].re, 1.0);
    }
}

/* Copie d'une matrice 1D de type puissance N à une matrice 2D de type puissance N */
void copie1DtoD2N(int tailleMatrice, t_matriceN matrice, t_puissanceN **resultante)
{
    int i = 0, j, lig = 0, col = 0;
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            for(j = 0; j < 5; j++)
            {
                resultante[lig][col].Case[j].coef = matrice.mat[i].Case[j].coef;
                resultante[lig][col].Case[j].n = matrice.mat[i].Case[j].n;
            }
            resultante[lig][col].nb = matrice.mat[i].nb;
            i++;
        }
    }
}

/* Copie d'une matrice 2D de type puissance N à une matrice 1D de type puissance N */
void copie2Dto1DN(int tailleMatrice, t_puissanceN **matrice, t_matriceN *resultante)
{
    int i = 0, j, lig = 0, col = 0;
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            for(j = 0; j < 5; j++)
            {
                resultante->mat[i].Case[j].coef = matrice[lig][col].Case[j].coef;
                resultante->mat[i].Case[j].n = matrice[lig][col].Case[j].n;
            }
            resultante->mat[i].nb = matrice[lig][col].nb;
            i++;
        }
    }
}

/* Copie d'une matrice de type puissance N vers une autre */
void copieN(int tailleMatrice, t_matriceN matrice1, t_matriceN *matrice2)
{
    int i = 0, j;
    for(i = 0; i < (tailleMatrice * tailleMatrice); i++)
    {
        for(j = 0; j < 5; j++)
        {
            matrice2->mat[i].Case[j].coef = matrice1.mat[i].Case[j].coef;
            matrice2->mat[i].Case[j].n = matrice1.mat[i].Case[j].n;
        }
        matrice2->mat[i].nb = matrice1.mat[i].nb;
    }
}

/* Addition d'une "case" puissance N */
void AddCase(t_puissanceN A, t_puissanceN B, t_puissanceN *C)
{
    int i, nb = 0;

    for(i = 0; i < A.nb; i++)
    {
        C->Case[nb] = A.Case[i];
        nb += 1;
    }

    for(i = 0; i < B.nb; i++)
    {
        C->Case[nb] = B.Case[i];
        nb += 1;
    }
    C->nb = nb;
}

/* Produit de deux cases puissance N */
t_case ProduitDeuxCase(t_case A, t_case B)
{
    t_case C;

    C.coef = A.coef * B.coef;
    C.n = A.n * B.n;
    return C;
}

/* Produit de coefficients d'une matrice puissance N */
void ProduitCase(t_puissanceN A, t_puissanceN B, t_puissanceN *C)
{
    int i, j, nbF = 0;

    for(i = 0; i < A.nb; i++)
    {
        for(j = 0; j < B.nb; j++)
        {
            #if DEBUG
            printf("A[%i].coef = %g\n", i, A.Case[i].coef);
            printf("A[%i].n = %g\n\n", i, A.Case[i].n);
            printf("B[%i].coef = %g\n", j, B.Case[i].coef);
            printf("B[%i].n = %g\n\n", j, B.Case[i].n);
            #endif
            C->Case[nbF] = ProduitDeuxCase(A.Case[i], B.Case[j]);
            #if DEBUG
            printf("C[%i].coef = %g\n", nbF, C->Case[nbF].coef);
            printf("C[%i].n = %g\n\n", nbF, C->Case[nbF].n);
            system("pause");
            #endif
            nbF += 1;
        }
    }
    C->nb = nbF;
}

/* Produit de matrice à la puissance N */
void matMutPuissanceN(int tailleMatrice, t_puissanceN **coeffMatrice1, t_puissanceN **coeffMatrice2, t_puissanceN **coeffResultante, t_puissanceN temp, t_matriceN matrice1, t_matriceN matrice2, t_matriceN *resultante)
{
    int lig = 0, col = 0, lig1 = 0, col1 = 0, lig2 = 0, col2 = 0;

    #if DEBUG
    printf("Le produit est :\n\n");
    afficheMatricePuissanceN(tailleMatrice, matrice1);
    printf("\n\net :\n\n");
    afficheMatricePuissanceN(tailleMatrice, matrice2);
    printf("\n\n");
    system("pause");
    #endif

    copie1DtoD2N(tailleMatrice, matrice1, coeffMatrice1);
    copie1DtoD2N(tailleMatrice, matrice2, coeffMatrice2);

    // Multiplication
    for(lig = 0 ; lig < tailleMatrice ; lig++)
    {
        for(col = 0 ; col < tailleMatrice ; col++)
        {
            coeffResultante[lig][col].Case[0].coef = 1.0;
            coeffResultante[lig][col].Case[0].n = 1.0;
            coeffResultante[lig][col].nb = 0;

            lig1 = lig;
            col1 = 0;
            col2 = col;
            lig2 = 0;
            do
            {
                ProduitCase(coeffMatrice1[lig1][col1], coeffMatrice2[lig2][col2], &temp);

                #if DEBUG
                printf("Coef[%d][%d].re = %g\n", lig, col, coeffResultante[lig][col].re);
                printf("Coef[%d][%d].im = %g\n", lig, col, coeffResultante[lig][col].im);

                printf("Temp[%d][%d].re = %g\n", lig, col, temp.re);
                printf("Temp[%d][%d].im = %g\n", lig, col, temp.im);
                #endif

                AddCase(coeffResultante[lig][col], temp, &coeffResultante[lig][col]);

                #if DEBUG/*
                printf("result[%d][%d].re = %g\n", lig, col, coeffResultante[lig][col].re);
                printf("result[%d][%d].im = %g\n", lig, col, coeffResultante[lig][col].im);
                system("pause");*/
                #endif

                col1++;
                lig2++;
            } while(col1 < tailleMatrice);
        }
    }

    copie2Dto1DN(tailleMatrice, coeffResultante, resultante);
}

/* Affichage d'une matrice de type puissance N */
void afficheMatricePuissanceN(int tailleMatrice, t_matriceN matrice)
{
    int i, j, X, x, y;

    printf("   ");
    X = wherex();
    x = X;
    y = wherey();

    for (i = 0 ; i < tailleMatrice * tailleMatrice; i++)
    {
        gotoxy(x, y);

        affichePuissanceN(matrice.mat[i].Case[0].coef, matrice.mat[i].Case[0].n, 1);
        for(j = 1; j < matrice.mat[i].nb; j++)
        {
            affichePuissanceN(matrice.mat[i].Case[j].coef, matrice.mat[i].Case[j].n, 0);
        }

        if((i + 1) % tailleMatrice != 0)
        {
            x = x + 45;
        }
        else
        {
            x = X;
            y = y + 2;
        }
    }
}

/* Affichage optimisée d'un coefficient d'une matrice puissance N */
void affichePuissanceN(double coef, double n, int debut)
{
    if((coef == 0.0 || n == 0.0) && debut == 1)
    {
        printf("0");
    }
    else
    {
        if(debut == 0 && coef > 0.0)
        {
            printf("+");
        }

        if(coef == 1.0 && n != 1.0)
        {
            if(n > 0.0)
            {
                printf("%g^n", n);
            }
            else
            {
                printf("(%g)^n", n);
            }
        }
        else
        {
            printf("%g", coef);

            if(n != 1.0)
            {
                if(n > 0.0)
                {
                    printf("x%g^n", n);
                }
                else if(n < 0.0)
                {
                    printf("x(%g)^n", n);
                }
            }
        }
    }
    return;
}
