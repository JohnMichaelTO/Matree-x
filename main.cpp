#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "main.h"
#include "matrice.h"
#include "puissanceN.h"
#include "fonctions.h"
#include "conio.h"

#define SYMB 42 // Symbole utilisé pour les listes
#define NOMBRE_OPTIONS 3
#define NOMBRE_OPERATIONS 10

using namespace std;

int main()
{
    int selection = 0, run = 1;

    do
    {
        system("cls");
        printf("\n   %c Matree'x %c\n\n\n %c 0 - Quitter\n\n %c 1 - Taille 2\n %c 2 - Taille 3\n %c 3 - Taille 4\n\n\n", SYMB, SYMB, SYMB, SYMB, SYMB, SYMB);
        do
        {
            printf("\n  Selection : ");
            scanf("%d",&selection);
            if (selection < 0 || selection > NOMBRE_OPTIONS) // Erreur dans la saisie.
            {
                printf(" /!\\ Erreur : Ceci n'est pas un item selectionnable. /!\\ \n");
            }
            else if (selection != 0) // Sous-programme valide sélectionné : lancement.
            {
                printf(" Execution de la fonction desiree.\n\n");
            }
            else // Arrêt
            {
                run = 0;
            }
        } while (selection < 0 || selection > NOMBRE_OPTIONS);

        system("cls");
        switch (selection) // Vérification et lancement du sous-programme approprié.
        {
            case 1:
                menuMatrice(2);
                break;
            case 2:
                menuMatrice(3);
                break;
            case 3:
                menuMatrice(4);
                break;
        }
    } while(run != 0);
    return 0;
}

void menuMatrice (int tailleMatrice)
{
    int selection = 0, run = 0, pasUneMatrice = 0, deuxMatrices = 0, k, inverse = 1, mineur = 0;
    t_complexe scalaire;
    scalaire.im = 0.0;
    Matrice matrice1, matrice2, resultante;

    do
    {
        system("cls");
        printf("\n   %c Matrice de taille %d %c\n\n\n %c 0 - Retour\n\n\n   Operations sur une matrice\n\n %c 1 - Informations Diverses\n %c 2 - Transposee\n %c 3 - Inverse\n %c 4 - Mineur et Comatrice\n %c 5 - Multiplication par un scalaire\n %c 6 - Puissance k\n %c 7 - Exponentielle de matrice\n\n\n   Operations sur deux matrices\n\n %c 8 - Additionner\n %c 9 - Soustraire\n %c 10 - Multiplier\n\n\n", SYMB, tailleMatrice, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB, SYMB);
        do
        {
            printf("\n  Selection : ");
            scanf("%d", &selection);
            if (selection < 0 || selection > NOMBRE_OPERATIONS) // Erreur dans la saisie.
            {
                printf(" /!\\ Erreur : Ceci n'est pas un item selectionnable. /!\\ \n");
            }
            else if (selection != 0) // Sous-programme valide sélectionné : lancement.
            {
                printf(" Execution de la fonction desiree.\n\n");
            }
            else // Arrêt
            {
                run = 0;
            }
        } while (selection < 0 || selection > NOMBRE_OPERATIONS);
        system("cls");

        if(selection > 0)
        {
            // Création des matrices
            creatMatrice(tailleMatrice, &matrice1);
            creatMatrice(tailleMatrice, &matrice2);
            creatMatrice(tailleMatrice, &resultante);

            saisieMatrice(tailleMatrice, &matrice1); // Saisie de la première matrice

            if(selection > 7) // Si besoin, saisie de la deuxième matrice
            {
                system("cls");
                saisieMatrice(tailleMatrice, &matrice2);
                deuxMatrices = 1;
            }
            system("cls");

            switch (selection) // Vérification et lancement du sous-programme approprié.
            {
                case 1:
                    pasUneMatrice = 1;
                    break;
                case 2:
                    matTrans(tailleMatrice, matrice1, &resultante);
                    break;
                case 3:
                    inverse = matInverse(tailleMatrice, matrice1, &resultante);
                    break;
                case 4:
                    mineur = 1;
                    break;
                case 5:
                    // Saisie d'un scalaire réel
                    printf("Saisir le scalaire voulu : ");
                    scanf("%lf", &scalaire.re);
                    matMutScalaire(tailleMatrice, matrice1, scalaire, &resultante);
                    break;
                case 6:
                    // Saisie d'une puissance k
                    printf("\n Entrez une puissance k != 0 : ");
                    scanf("%d", &k);
                    matPuissance(tailleMatrice, matrice1, k, &resultante);
                    break;
                case 7:
                    matExp(tailleMatrice, matrice1, &resultante);
                    break;
                case 8:
                    matAdd(tailleMatrice, matrice1, matrice2, &resultante);
                    break;
                case 9:
                    matDiff(tailleMatrice, matrice1, matrice2, &resultante);
                    break;
                case 10:
                    matMut(tailleMatrice, matrice1, matrice2, &resultante);
                    break;
            }
            // Affichage
            affichageTotal(tailleMatrice, pasUneMatrice, deuxMatrices, inverse, mineur, matrice1, matrice2, resultante);
            // Libération mémoire
            freeMatrice(&matrice1);
            freeMatrice(&matrice2);
            freeMatrice(&resultante);
        }
    } while(run != 0);
}
