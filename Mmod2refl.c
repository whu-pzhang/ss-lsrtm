/* Convert acoustic impedance to normal-indicent reflectivity.
 */

/**
 *      author: pzhang
 *      date: 2017-09-28
 */

#include <rsf.h>

int main(int argc, char *argv[])
{
    sf_init(argc, argv);

    sf_file Fvel = sf_input("in");
    /*< input velocity file >*/
    sf_file Fden = sf_input("den");
    /*< input density file >*/
    sf_file Fref = sf_output("out");
    /*< output reflectivity file >*/

    int n1, n2;
    if (!sf_histint(Fvel, "n1", &n1)) sf_error("No n1= in input");
    /*< fast dissmension of input >*/
    n2 = sf_leftsize(Fvel, 1);

    float *vtr = malloc(n1 * sizeof(float));
    float *sig = malloc(n1 * sizeof(float));

    int n11, n22;
    if (!sf_histint(Fden, "n1", &n11)) sf_error("No n1= in density");
    n22 = sf_leftsize(Fden, 1);

    if (n11 != n1 || n22 != n2) sf_error("Dissmension not match!");

    float *dtr = malloc(n11 * sizeof(float));

    for (int i2 = 0; i2 < n2; i2++) {
        sf_floatread(vtr, n1, Fvel);  // read in velocity
        sf_floatread(dtr, n11, Fden); // read in density
        float imp1 = dtr[0] * vtr[0];
        for (int i1 = 0; i1 < n1 - 1; i1++) {
            float imp2 = dtr[i1 + 1] * vtr[i1 + 1];
            sig[i1] = (imp2 - imp1) / (imp2 + imp1 + FLT_EPSILON);
            imp1 = imp2;
        }

        sig[n1 - 1] = 0.;

        sf_floatwrite(sig, n1, Fref);
    }

    free(vtr);
    free(dtr);
    free(sig);

    return 0;
} /* main */
