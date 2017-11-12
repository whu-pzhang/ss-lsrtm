/* pre-stack depth RTM migration and its adjoint(single-shot)
 *
 *   with source illumination
 */
#include <rsf.h>
#include "prertm2d_v04.h"
#include "dbg.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

int main(int argc, char *argv[])
{
    bool adj, verb, illum;
    // bool illum;
    int nx, nz, nt, nr, n1, n2;
    float dx, dz, dt, dr, or, ox, oz, ot;

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
    if (!sf_getbool("illum", &illum)) illum = false;

    // setup I/O files
    sf_file Finp = sf_input("in");
    sf_file Fout = sf_output("out");

    sf_file Fvel = sf_input("vel"); // migration velocity
    sf_file Fwlt = sf_input("wlt"); // wavelet
    sf_file Fsou = sf_input("sou"); // sources coordinates
    sf_file Frec = sf_input("rec"); // receiver coordinates

    // velocity dimensions
    if (!sf_histint(Fvel, "n1", &nz)) sf_error("No n1= in velocity");
    if (!sf_histint(Fvel, "n2", &nx)) sf_error("No n2= in velocity");
    if (!sf_histfloat(Fvel, "d1", &dz)) sf_error("No d1= in velocity");
    if (!sf_histfloat(Fvel, "d2", &dx)) sf_error("No d2= in velocity");
    if (!sf_histfloat(Fvel, "o1", &oz)) sf_error("No o1= in velocity");
    if (!sf_histfloat(Fvel, "o2", &ox)) sf_error("No o2= in velocity");

    sf_axis az = sf_iaxa(Fvel, 1);
    sf_axis ax = sf_iaxa(Fvel, 2);

    if (!sf_histint(Frec, "n2", &nr)) sf_error("No n2= in reciever");
    if (!sf_histfloat(Frec, "d2", &dr)) sf_error("No d2= in reciever");
    if (!sf_histfloat(Frec, "o2", & or)) sf_error("No o2= in reciever");

    sf_axis ar = sf_iaxa(Frec, 2);

    if (adj) { // migration
        if (!sf_histint(Finp, "n1", &nr)) sf_error("No n1= in data");
        if (!sf_histfloat(Finp, "d1", &dr)) sf_error("No d1= in data");
        if (!sf_histint(Finp, "n2", &nt)) sf_error("No n2= in data");
        if (!sf_histfloat(Finp, "d2", &dt)) sf_error("No d1= in data");
        if (!sf_histfloat(Finp, "o2", &ot)) sf_error("No o1= in data");

        // setup output file headers
        sf_oaxa(Fout, az, 1);
        sf_oaxa(Fout, ax, 2);
    } else { // modeling
        if (!sf_histint(Finp, "n1", &n1) || n1 != nz)
            sf_error("Need n1= %d in mod", nz);
        if (!sf_histint(Finp, "n2", &n2) || n2 != nx)
            sf_error("Need n2= %d in mod", nx);
        if (!sf_histint(Fwlt, "n2", &nt)) sf_error("No n2= in wavelet");
        if (!sf_histfloat(Fwlt, "d2", &dt)) sf_error("No d2= in wavelet");
        if (!sf_histfloat(Fwlt, "o2", &ot)) sf_error("No o2= in wavelet");

        sf_putint(Fout, "n2", nt);
        sf_putfloat(Fout, "d2", dt);
        sf_putfloat(Fout, "o2", 0.0);
        sf_putstring(Fout, "label2", "Time");
        sf_putstring(Fout, "unit2", "s");

        sf_oaxa(Fout, ar, 1);
    }

    int nb;
    if (!sf_getint("nb", &nb)) nb = 30;

    float **vv = sf_floatalloc2(nz, nx);
    float *wlt = sf_floatalloc(nt);
    pt2d *src2d = pt2dalloc1(1);
    pt2d *rec2d = pt2dalloc1(nr);

    float *dat = sf_floatalloc(nt * nr);
    float *mod = sf_floatalloc(nz * nx);
    if (adj)
        sf_floatread(dat, nt * nr, Finp);
    else
        sf_floatread(mod, nz * nx, Finp);

    // read velocity
    sf_floatread(vv[0], nz * nx, Fvel);
    // read wavelet
    sf_floatread(wlt, nt, Fwlt);
    // 读入炮点位置
    pt2dread1(Fsou, src2d, 1, 2); // 只读取坐标信息（x,z），不读值value
    // 读入接收点位置
    pt2dread1(Frec, rec2d, nr, 2);

    if (verb) {
        log_info("nx=%7d\tdx=%7.3f\tox=%7.3f", nx, dx, ox);
        log_info("nz=%7d\tdz=%7.3f\toz=%7.3f", nz, dz, oz);
        log_info("nt=%7d\tdt=%7.3f\tot=%7.3f", nt, dt, ot);
        log_info("nr=%7d\tdr=%7.3f\tor=%7.3f", nr, dr, or);
        log_info("source position: x=%7.3f\tz=%7.3f", src2d->x, src2d->z);
    }

    prertm2d_init(verb, illum, nx, nz, nb, nt, nr, dx, dz, dt, ox, oz, vv, wlt,
                  src2d, rec2d);
    //
    prertm2d_loop(adj, false, nz * nx, nt * nr, mod, dat);
    //
    prertm2d_close();

    if (adj) { // migration
        // output image
        sf_floatwrite(mod, nz * nx, Fout);
    } else { // modeling
        // output data
        sf_floatwrite(dat, nt * nr, Fout);
    }

    // clean up
    free(*vv);
    free(vv);
    free(wlt);
    free(src2d);
    free(rec2d);
    free(dat);
    free(mod);

#ifdef _OPENMP
    double wall_clock_time_e = omp_get_wtime();
#else
    double wall_clock_time_e = (double)clock() / CLOCKS_PER_SEC;
#endif
    // if (verb) {
    fprintf(stderr, "\nElapsed time: %lf s\n",
            wall_clock_time_e - wall_clock_time_s);
    // }

    return EXIT_SUCCESS;
}
