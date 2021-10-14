#include <stdio.h>
#include <assert.h>

/** \brief find_roots() is needed to find roots from giving parameters
 *
 * \param 'a' 'b' 'c' param. form equation 'ax^2+bx+c=0'
 * \param x1 and x2 solutions of equation if they exist
 * \return returns number of roots
 *
 */
int find_roots (double a, double b, double c, double *x1, double *x2);

int inputs (double *a, double *b, double *c);
double discriminant(double a, double b, double c);
int mod_more_0(double a);
int solve_linear(double k, double b, double* x1);
int output_roots (int num_of_roots, double x1, double x2);
int test_project ();
int check_input (double* a);

/**< const for comparing with 0  */
const double epsilon = 1e-6;/// константа для сранения вводимых чисел с нулем

/**< enum is needed to denote the numbers as number of roots */
enum roots
{
    ZERO_ROOTS = 0,
    ONE_ROOT   = 1,
    TWO_ROOTS  = 2,
    INF_ROOTS  = 3
};

int main()
{
    test_project();

    double a = 0, b = 0, c = 0;
    inputs (&a, &b, &c);

    double x1 = 0, x2 = 0;

    int num_of_roots = find_roots (a, b, c, &x1, &x2);
    output_roots(num_of_roots, x1, x2);

    return 0;
}

/** \brief find_roots() is needed to find roots from giving parameters
 *
 * \param 'a' 'b' 'c' param. form equation 'ax^2+bx+c=0'
 * \param x1 and x2 solutions of equation if they exist
 * \return number of roots
 *
 */
int find_roots (double a, double b, double c, double* x1, double* x2)
{
    assert (x1 != x2);
    if (!mod_more_0(a))
    {
        return solve_linear (b, c, &x1);
    }

    double Discr = discriminant(a, b, c);

    if (!mod_more_0(Discr))
    {
        *x1 = -b / (2*a);
        return ONE_ROOT;
    }
    else if (Discr < 0)
    {
        return ZERO_ROOTS;
    }
    else if (Discr > 0)
    {
        if (!mod_more_0(c))
        {
            *x2 = 0;
            return solve_linear(a, b, &x1) + 1;
        }
        else
        {
            double Sq_Di = sqrt(Discr);
            *x1 = ( -b + Sq_Di ) / (2*a);
            *x2 = ( -b - Sq_Di ) / (2*a);
            return TWO_ROOTS;
        }
    }
}

/** \brief inputs() writing words for undestanding function of programm
 *
 * \param a b c param. of equation to count from console
 *
 */
int inputs (double* a, double* b, double* c)
{
    printf("Enter coefficient a b c from your equation ax^2+bx+c=0\n"
           "a:");
    check_input (a);
    printf("b:");
    check_input (b);
    printf("c:");
    check_input (c);
    return 0;
}

/** \brief checking param. counting from console not to write a letter or symbols instead of number
 *
 * \param a param. is coefficient from equation
 * \return 0 anyway
 *
 */
int check_input (double* a)
{
   while (scanf("%lg", a) == 0)
   {
       printf("pls enter right number\n");
       while (getchar() != '\n') {;}
       printf(":");
   }

   while (getchar() != '\n') {;}

   return 0;
}

/** \brief counting discriminant while solving equation
 *
 * \param a b c from equation
 * \return meaning of discriminant
 *
 */
double discriminant (double a, double b, double c)
{
    return b*b - 4*a*c;
}

/** \brief solving linear equation kx+b=0
 *
 * \param k b from new equation
 * \param x1 is solution of this equation
 * \return number of roots in this equation
 *
 */
int solve_linear (double k, double b, double* x1)
{
    if (!mod_more_0(k))
    {
        if (!mod_more_0(b))
            return INF_ROOTS;
        else
            return ZERO_ROOTS;
    }
    else
    {
        *x1 = -b / k;
        return ONE_ROOT;
    }
}

/** \brief comparing with 0
 *
 * \param a is param which we are comparing
 * \return 0 if 'a' equal to 0, else 1
 *
 */
int mod_more_0 (double a)
{
    if (fabs(a) <= epsilon)
        return 0;
    else
        return 1;
}

/** \brief printing solutions of equation
 *
 * \param num_of_roots
 * \param x1, x2 are wrote on console
 *
 */
int output_roots (int num_of_roots, double x1, double x2)
{
    switch (num_of_roots)
    {
    case ZERO_ROOTS:
        printf("No solution\n");
        break;

    case ONE_ROOT:
        printf("x = %lg\n", x1);
        break;

    case TWO_ROOTS:
        printf("x1 = %lg\n"
               "x2 = %lg\n", x1, x2);
        break;

    case INF_ROOTS:
        printf("many solution\n");
        break;

    default:
        break;
    }
    return 0;
}

/** \brief testing our project to find out errors
 * \return num of errors
 */
int test_project ()
{
    double x1 = 0, x2 = 0;
    int ERRORS = 0;
    double solut[] = {INF_ROOTS, ONE_ROOT, ZERO_ROOTS, TWO_ROOTS, TWO_ROOTS};
    struct koef{
    double a; double b; double c;
    };
    struct koef* k = (struct koef*) calloc(5, sizeof(*k));
    k[0].a = 0, k[0].b = 0, k[0].c = 0;
    k[1].a = 1, k[1].b = 2, k[1].c = 1;
    k[2].a = 0, k[2].b = 0, k[2].c = 1;
    k[3].a = 1, k[3].b = 5, k[3].c = 4;
    k[4].a = 1, k[4].b = 3, k[4].c = 2;
    for (int i = 0; i <= 4; i++)
    {
        if (find_roots(k[i].a, k[i].b, k[i].c,  &x1, &x2) != solut[i])
        {
            ERRORS++;
        }
    }
    free(k);
    printf("%d ERRORS\n", ERRORS);
    return 0;
}
