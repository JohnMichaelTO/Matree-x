#define DEF_X 18

/* Structure d'une racine */
typedef struct s_racine
{
    double re;
    double im;
} t_racine;

/* Structure d'un complexe */
typedef t_racine t_complexe;

/* Structure d'une matrice 1D */
typedef struct s_matrice
{
    t_complexe *coefficient;
} Matrice;

/* Structure d'un vecteur propre */
typedef struct s_vecteur
{
    t_complexe *n;
} t_vecteur;

/* Saisie d'une matrice */
void saisieMatrice(int tailleMatrice, Matrice *matrice);

/* Transposée d'une matrice */
void matTrans(int tailleMatrice, Matrice matrice, Matrice* resultante);

/* Calcul de la trace d'une matrice */
t_complexe matTrace(int tailleMatrice, Matrice matrice);

/* Inversion de matrice */
int matInverse(int tailleMatrice, Matrice matrice, Matrice* resultante);

/* Calcul des coefficients du polynome caractéristique */
void polyCaract(int degree, Matrice matrice, double *Poly);

/* Affichage du polynome caractéristique */
void affichePoly(char variable, double *poly, int degree);

/* Calcul et Affichage des valeurs propres */
void valeurPropre(int degree, double *poly, t_racine *racines);

/* Calcul des vecteurs propres */
void vecteurPropre(int tailleMatrice, Matrice matrice, t_complexe inconnue, t_vecteur *vecteur, int nVect);

/* Produit de nombres complexes */
t_complexe produitComplexe(t_complexe A, t_complexe B);

/* Division de nombres complexes */
t_complexe divComplexe(t_complexe A, t_complexe B);

/* Fonction échange */
void swap(t_complexe *x, t_complexe *y);

/* Multiplication de matrice par un scalaire */
void matMutScalaire(int tailleMatrice, Matrice matrice, t_complexe scalaire, Matrice *resultante);

/* Addition de matrices */
void matAdd(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice* resultante);

/* Soustraction de matrices */
void matDiff(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice* resultante);

/* Produit de matrices */
void matMut(int tailleMatrice, Matrice matrice1, Matrice matrice2, Matrice *resultante);

/* Matrice à la puissance k */
void matPuissance(int tailleMatrice, Matrice matrice, int n, Matrice *resultante);

/* Exponentielle d'une matrice */
void matExp(int tailleMatrice, Matrice matrice, Matrice *resultante);

/* Calcul du mineur d'une matrice */
void matMin(int tailleMatrice, Matrice matrice, Matrice *resultante);

/* Pair/Impair dans les coefficients pour générer la comatrice */
void matCom(int tailleMatrice, Matrice matrice, Matrice *resultante);

/* Calcul du déterminant d'une matrice */
t_complexe matDet(int tailleMatrice, Matrice matrice);

/* Matrice identité */
void matIdent(int tailleMatrice, Matrice matrice, Matrice *resultante);

/* Copie d'une matrice 1D vers une matrice 2D */
void copie1DtoD2(int tailleMatrice, Matrice matrice, t_complexe **resultante);

/* Copie d'une matrice 2D vers une matrice 1D */
void copie2Dto1D(int tailleMatrice, t_complexe **matrice, Matrice *resultante);

/* Copie d'une matrice 1D */
void copie(int tailleMatrice, Matrice matrice1, Matrice matrice2);

/* Diagonalisation */
void diagonalisation(int tailleMatrice, t_racine *racine, t_vecteur * vecteur, Matrice *resultante);

/* Création d'une matrice 1D */
void creatMatrice(int tailleMatrice, Matrice *matrice);

/* Libération d'une matrice 1D */
void freeMatrice(Matrice *matrice);

/* Création d'une matrice 2D */
void creatMatrice2D(int tailleMatrice, t_complexe ***matrice);

/* Libération d'une matrice 2D */
void freeMatrice2D(int tailleMatrice, t_complexe ***matrice);

/* Affichage d'une matrice 1D */
void afficheMatrice(int tailleMatrice, Matrice matrice);

/* Affichage d'une matrice 2D */
void afficheMatrice2D(int tailleMatrice, t_complexe **matrice);

/* Affichage finale */
void affichageTotal(int tailleMatrice, int pasUneMatrice, int deuxMatrices, int inverse, int mineur, Matrice matrice1, Matrice matrice2, Matrice resultante);

/* Vérification facile par un copier coller sur wolfram alpha */
void verifierWolfram(int tailleMatrice, Matrice matrice);
