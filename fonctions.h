/* Affichage d'un complexe */
void afficheComplexe(double re, double im);

/* Résolution d'équation de degré 1 */
double degre1(double a, double b);

/* Fonction puissance */
double puissance(double x, int n);

/* Fonction racine carrée */
double racine_carree(double x);

/* Fonction racine cubique */
/*
double racine_cubique(double x);
*/

/* Factoriel */
double factoriel(int n);

/* Calcul d'un polynome */
double calculPoly(double polynome[], int degree, double x);

/* Méthode de Newton pour la recherche d'une racine */
double newton(double polynome[], int degree, double x0);

/* Retourner la valeur absolue d'un nombre */
double absolue(double x);

/* Méthode de dichotomie */
double dichotomie(double polynome[], int degree);

double cosinus(double x);

/* Calcul du delta */
double calcul_delta(double a, double b, double c);

/* Calcul une racine réelle */
double racine_reelle(double a, double b, double delta, int numero_racine);

/* Calcul une racine complexe */
double racine_complexe(double a, double b, double delta, char partie);

/* Résolution d'équation de degré 2 */
void degre2(t_racine *racine, double a, double b, double c, double changementVariable, int i);

/* Résolution d'équation de degré 3 */
void degre3(t_racine *racine, double a, double b, double c, double d, int i);

/* Résolution d'équation de degré 4 */
void degre4(t_racine *racine, double a, double b, double c, double d, double e, int i);
