/* Contenu dans un coefficient d'une matrice puissance N */
typedef struct s_case
{
    double coef;
    double n;
} t_case;

/* Coefficient d'une matrice puissance N */
typedef struct s_puissanceN
{
    t_case *Case;
    int nb;
} t_puissanceN;

/* Matrice puissance N */
typedef struct s_matriceN
{
    t_puissanceN *mat;
} t_matriceN;

/* Création d'une matrice de type puissance N */
void creatMatriceN(int tailleMatrice, t_matriceN *matrice, int nb);

/* Libération de matrice de type puissance N */
void freeMatriceN(int tailleMatrice, t_matriceN *matrice);

/* Création d'une matrice 2D de type puissance N */
void creatMatrice2DN(int tailleMatrice, t_puissanceN ***matrice, int nb);

/* Libération de matrice 2D de type puissance N */
void freeMatrice2DN(int tailleMatrice, t_puissanceN ***matrice);

/* Ajouter un coefficient + un coefficient à la puissance n */
void AddPuissanceN(int tailleMatrice, int nb, t_matriceN *matrice, int i, double coef, double n);

/* Diagonalisation et matrice à la puissance N */
void diagonalisationPuissanceN(int tailleMatrice, t_racine *racine, t_vecteur * vecteur);

/* Copie une matrice 1D dans une matrice de type puissance N */
void copieToN(int tailleMatrice, Matrice Pmat, t_matriceN *P);

/* Copie d'une matrice 1D de type puissance N à une matrice 2D de type puissance N */
void copie1DtoD2N(int tailleMatrice, t_matriceN matrice, t_puissanceN **resultante);

/* Copie d'une matrice 2D de type puissance N à une matrice 1D de type puissance N */
void copie2Dto1DN(int tailleMatrice, t_puissanceN **matrice, t_matriceN *resultante);

/* Copie d'une matrice de type puissance N vers une autre */
void copieN(int tailleMatrice, t_matriceN matrice1, t_matriceN *matrice2);

/* Addition d'une "case" puissance N */
void AddCase(t_puissanceN A, t_puissanceN B, t_puissanceN *C);

/* Produit de deux cases puissance N */
t_case ProduitDeuxCase(t_case A, t_case B);

/* Produit de coefficients d'une matrice puissance N */
void ProduitCase(t_puissanceN A, t_puissanceN B, t_puissanceN *C);

/* Produit de matrice à la puissance N */
void matMutPuissanceN(int tailleMatrice, t_puissanceN **coeffMatrice1, t_puissanceN **coeffMatrice2, t_puissanceN **coeffResultante, t_puissanceN temp, t_matriceN matrice1, t_matriceN matrice2, t_matriceN *resultante);

/* Affichage d'une matrice de type puissance N */
void afficheMatricePuissanceN(int tailleMatrice, t_matriceN matrice);

/* Affichage optimisée d'un coefficient d'une matrice puissance N */
void affichePuissanceN(double coef, double n, int debut);
