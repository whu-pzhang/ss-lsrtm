/* least-square pre-stack depth RTM migration with structure-enhancing constrained
 */
#include <rsf.h>
#include <rsfpwd.h>

#include "lsrtmsr.h"

#ifdef _OPENMP
#include <omp.h>
#else
#include <time.h>
#endif

const float TOLERANCE = 1.e-12;

int main(int argc, char *argv[])
{
    bool adj, verb;
    int nx, nz, nt, ns, nr, nss, niter;
    float dx, dz, dt, ox, oz;
    float eps;
    int rect1, rect2, order, radius;

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
    if (!sf_getint("niter", &niter)) niter = 5;
    if (!sf_getfloat("eps", &eps)) eps = 0.01;
    if (!sf_getint("order", &order)) order=1;
    if (!sf_getint("rect1", &rect1)) rect1=5;
    if (!sf_getint("rect2", &rect2)) rect2=5;
    if (!sf_getint("radius", &radius)) radius=10;

    // setup I/O files
    sf_file Finp = sf_input("in");
    sf_file Fout = sf_output("out");

    sf_file Ferr = sf_output("error");
    //sf_file Fres = sf_output("res"); // residual r = Lm - d

    sf_file Fvel = sf_input("vel");  // migration velocity
    sf_file Fwlt = sf_input("wlt");  // wavelet
    sf_file Fsou = sf_input("sou");  // sources coordinates
    sf_file Frec = sf_input("rec");  // receiver coordinates

    sf_file Fdip = sf_input("dip");  // structure dip file

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

    if (!sf_histint(Finp, "n1", &nt)) sf_error("No n1= in data");
    if (!sf_histfloat(Finp, "d1", &dt)) sf_error("No d1= in data");
    if (!sf_histint(Finp, "n2", &nr)) sf_error("No n2= in data");

    // setup output file headers
    sf_oaxa(Fout, az, 1);
    sf_oaxa(Fout, ax, 2);
    sf_putint(Fout, "n3", 1);
    // error
    sf_putint(Ferr, "n1", niter);
    sf_putfloat(Ferr, "d1", 1);
    sf_putfloat(Ferr, "o1", 1);
    sf_putint(Ferr, "n2", 1);
    sf_putint(Ferr, "n3", 1);

    int nb;
    if (!sf_getint("nb", &nb)) nb = 30;

    float **vv = sf_floatalloc2(nz, nx);
    float *wlt = sf_floatalloc(nt);
    pt2d *src2d = pt2dalloc1(ns);
    pt2d *rec2d = pt2dalloc1(nr);

    int nm = nz * nx;
    int nd = nt * nr * nss;
    float *dat = sf_floatalloc(nd);
    float *mod = sf_floatalloc(nm);
    float **slope = sf_floatalloc2(nz, nx);

    // float *pp = sf_floatalloc(nm);

    float *err = sf_floatalloc(niter);

    // 读入数据
    sf_floatread(dat, nd, Finp);
    memset(mod, 0.0, nm * sizeof(float));

    // read velocity
    sf_floatread(vv[0], nz * nx, Fvel);
    // read wavelet
    sf_floatread(wlt, nt, Fwlt);
    // 读入炮点位置
    pt2dread1(Fsou, src2d, ns, 2);  // 只读取坐标信息（x,z），不读值value
    // 读入接收点位置
    pt2dread1(Frec, rec2d, nr, 2);


    // read dip
    sf_floatread(slope[0], nm, Fdip);

    prertm2d_init(verb, nx, nz, nb, nt, ns, nss, nr, dx, dz, dt, ox, oz, vv,
                  wlt, src2d, rec2d);

    // prertm2d_loop(adj, false, nz * nx, nt * nr * nss, mod, dat);
    //  least-squares migration
    // sf_solver(prertm2d_loop, sf_cgstep, nz * nx, nt * nr * nss, mod, dat,
    // niter,
    //           "verb", verb, "end");
    // 初始化precondition算子
    // weight_init(wgt);   // Fm=d, m=Sp;
    // sf_triangle2_init(3,20,nz,nx,2);
    // pwdsl_init(nz, nx, order, rect1, rect2, eps);
    // pwdsl_set(wgt);
    pwsmooth_init(radius, nz, nx, order, eps);
    pwsmooth_set(slope);


    sf_solver_prec(prertm2d_loop, sf_cgstep, pwsmooth_lop,
                   nm, nm, nd, mod, dat,
                   niter, eps, "verb", verb, "err", err, "end");
    sf_cgstep_close();

    // conjugate gradient with shaping filter(regularization) to the input data
    // sf_conjgrad_init(nm, nm, nd, nd, eps, 1.e-6, true, false);
    // sf_conjgrad(NULL, prertm2d_loop, pwsmooth_lop, pp, mod, dat, niter);
    // sf_conjgrad_close();


    prertm2d_close();
    pwsmooth_close();

    // output image
    sf_floatwrite(mod, nz * nx, Fout);
    // output residual
    sf_floatwrite(err, niter, Ferr);

    // clean up
    free(*vv);
    free(vv);
    free(wlt);
    free(src2d);
    free(rec2d);
    free(dat);
    free(mod);
    free(*slope);
    free(slope);
    // free(pp);
    free(err);
    free(ar); free(as);
    free(ax); free(az);

#ifdef _OPENMP
    double wall_clock_time_e = omp_get_wtime();
#else
    double wall_clock_time_e = (double)clock() / CLOCKS_PER_SEC;
#endif
    if (verb) {
        fprintf(stderr, "\nElapsed time: %lf s\n",
                wall_clock_time_e - wall_clock_time_s);
    }

    exit(EXIT_SUCCESS);
}
