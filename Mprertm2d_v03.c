/* pre-stack depth RTM migration and its adjoint
*
*   with source illumination
*/
#include <rsf.h>

#include "prertm2d_v03.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

int main(int argc, char *argv[])
{
    bool adj, verb;
    //bool illum;
    int nx, nz, nt, ns, nr, n1, n2, nss;
    float dx, dz, dt, ox, oz;

    // Initialize Madagascar
    sf_init(argc, argv);

// initialize OpenMP support
#ifdef _OPENMP
    omp_init();
    double wall_clock_time_s = omp_get_wtime();
#else
    double wall_clock_time_s = (double)clock() / CLOCKS_PER_SEC;
#endif

    if (!sf_getbool("verb", &verb)) verb = false;
    if (!sf_getbool("adj", &adj)) adj = true;
    //if (!sf_getbool("illum", &illum)) illum = true;

    // setup I/O files
    sf_file Finp = sf_input("in");
    sf_file Fout = sf_output("out");

    sf_file Fvel = sf_input("vel");  // migration velocity
    sf_file Fwlt = sf_input("wlt");  // wavelet
    sf_file Fsou = sf_input("sou");  // sources coordinates
    sf_file Frec = sf_input("rec");  // receiver coordinates

    // velocity dimensions
    if (!sf_histint(Fvel, "n1", &nz)) sf_error("No n1= in velocity");
    if (!sf_histint(Fvel, "n2", &nx)) sf_error("No n2= in velocity");
    if (!sf_histfloat(Fvel, "d1", &dz)) sf_error("No d1= in velocity");
    if (!sf_histfloat(Fvel, "d2", &dx)) sf_error("No d2= in velocity");
    sf_axis az = sf_iaxa(Fvel, 1);
    sf_axis ax = sf_iaxa(Fvel, 2);
    ox = sf_o(ax);
    oz = sf_o(az);

    // source and reciever coordinates
    if (!sf_histint(Fsou, "n2", &ns)) sf_error("No n2= in source");
    /* ns: total shot numbers() */
    sf_axis as = sf_iaxa(Fsou, 2);
    if (!sf_getint("nss", &nss)) nss = ns;
    /* nss: supergather number in simulatneous-source or
           shot number(ns) in conventional acquisition */
    if (!sf_histint(Frec, "n2", &nr)) sf_error("No n2= in reciever");
    sf_axis ar = sf_iaxa(Frec, 2);

    if (adj) {  // migration
        if (!sf_histint(Finp, "n1", &nt)) sf_error("No n1= in data");
        if (!sf_histfloat(Finp, "d1", &dt)) sf_error("No d1= in data");
        if (!sf_histint(Finp, "n2", &nr)) sf_error("No n2= in data");

        // setup output file headers
        sf_oaxa(Fout, az, 1);
        sf_oaxa(Fout, ax, 2);
        sf_putint(Fout, "n3", 1);
    } else {  // modeling
        if (!sf_histint(Finp, "n1", &n1) || n1 != nz)
            sf_error("Need n1= %d in mod", nz);
        if (!sf_histint(Finp, "n2", &n2) || n2 != nx)
            sf_error("Need n2= %d in mod", nx);
        if (!sf_histint(Fwlt, "n2", &nt)) sf_error("No n2= in wavelet");
        if (!sf_histfloat(Fwlt, "d2", &dt)) sf_error("No d2= in wavelet");

        sf_putint(Fout, "n1", nt);
        sf_putfloat(Fout, "d1", dt);
        sf_putfloat(Fout, "o1", 0.0);
        sf_putstring(Fout, "label1", "Time");
        sf_putstring(Fout, "unit1", "s");

        sf_oaxa(Fout, ar, 2);
        sf_putint(Fout, "n3", nss);
        sf_putfloat(Fout, "d3", sf_d(as));
        sf_putfloat(Fout, "o3", sf_o(as));
    }

    int nb;
    if (!sf_getint("nb", &nb)) nb = 30;

    float **vv = sf_floatalloc2(nz, nx);
    float *wlt = sf_floatalloc(nt);
    pt2d *src2d = pt2dalloc1(ns);
    pt2d *rec2d = pt2dalloc1(nr);

    float *dat = sf_floatalloc(nt * nr * nss);
    float *mod = sf_floatalloc(nz * nx);
    if (adj)
        sf_floatread(dat, nt * nr * nss, Finp);
    else
        sf_floatread(mod, nz * nx, Finp);

    // read velocity
    sf_floatread(vv[0], nz * nx, Fvel);
    // read wavelet
    sf_floatread(wlt, nt, Fwlt);
    // 读入炮点位置
    pt2dread1(Fsou, src2d, ns, 2);  // 只读取坐标信息（x,z），不读值value
    // 读入接收点位置
    pt2dread1(Frec, rec2d, nr, 2);

    prertm2d_init(verb, nx, nz, nb, nt, ns, nss, nr, dx, dz, dt, ox, oz, vv,
                  wlt, src2d, rec2d);

    prertm2d_loop(adj, false, nz * nx, nt * nr * nss, mod, dat);

    prertm2d_close();

    if (adj)  // migration
        // output image
        sf_floatwrite(mod, nz * nx, Fout);
    else  // modeling
        // output data
        sf_floatwrite(dat, nt * nr * nss, Fout);

    // clean up
    free(*vv);
    free(vv);
    free(wlt);
    free(src2d);
    free(rec2d);
    free(dat);
    free(mod);
    sf_close();

#ifdef _OPENMP
    double wall_clock_time_e = omp_get_wtime();
#else
    double wall_clock_time_e = (double)clock() / CLOCKS_PER_SEC;
#endif
    if (verb) {
        fprintf(stderr, "\nElapsed time: %lf s\n",
                wall_clock_time_e - wall_clock_time_s);
    }

    return EXIT_SUCCESS;
}
