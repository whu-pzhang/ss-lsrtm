#include <rsf.h>
/*^*/
#include <assert.h>
// #include "dbg.h"

#include "prertm2d_v04.h"
// #include "laplac2.h"

#ifndef EPS
#define EPS FLT_EPSILON
#endif

static float c0, c11, c12, c21, c22;
static bool verb, illum;
static int nx, nz, nb, nxpad, nzpad, nt;
static int nr;
static float dx, dz, ox, oz;
static float *bndr;
static float **padvv, **sp0, **sp1, **gp0, **gp1, *rwbndr, *wlt;
static pt2d *src2d, *rec2d;

static void step_forward(bool adj, float **u0, float **u1);

void prertm2d_init(bool verb_, bool illum_, int nx_, int nz_, int nb_, int nt_,
                   int nr_, float dx_, float dz_, float dt_, float ox_,
                   float oz_, float **vv, float *wlt_, pt2d *src2d_,
                   pt2d *rec2d_)
/*< init >*/
{
    verb = verb_;
    illum = illum_;
    nx = nx_;
    nz = nz_;
    nb = nb_;
    nt = nt_;
    nr = nr_;
    dx = dx_;
    dz = dz_;
    ox = ox_;
    oz = oz_;

    nxpad = nx + 2 * nb;
    nzpad = nz + 2 * nb;

    // set laplacian coefficients
    float idz = 1.0 / (dz * dz);
    float idx = 1.0 / (dx * dx);

    c11 = 4.0 * idz / 3.0;
    c12 = -idz / 12.0;
    c21 = 4.0 * idx / 3.0;
    c22 = -idx / 12.0;
    c0 = -2.0 * (c11 + c12 + c21 + c22);

    wlt = wlt_;
    src2d = src2d_;
    rec2d = rec2d_;

    bndr = (float *)malloc(nb * sizeof(float));
    padvv = sf_floatalloc2(nzpad, nxpad);
    sp0 = sf_floatalloc2(nzpad, nxpad);
    sp1 = sf_floatalloc2(nzpad, nxpad);
    gp0 = sf_floatalloc2(nzpad, nxpad);
    gp1 = sf_floatalloc2(nzpad, nxpad);
    rwbndr = (float *)malloc((nt * 4 * (nx + nz)) * sizeof(float));

    /* initialized sponge ABC coefficients */
    for (int ib = 0; ib < nb; ib++) {
        float t = 0.015f * (nb - 1 - ib);
        bndr[ib] = expf(-t * t);
    }

    // precompute vp^2 * dt*2
    for (int ix = 0; ix < nx * nz; ++ix)
        *(vv[0] + ix) *= *(vv[0] + ix) * dt_ * dt_; // vel=vel^2*dt^2
    expand2d(padvv, vv);

    // log_info("nx=%7d\tdx=%7.3f\tox=%7.3f", nx, dx_, ox);
    //
    // log_info("nt=%7d", nt);
    // log_info("nr=%7d", nr);
    // log_info("source position: x=%7.3f\tz=%7.3f", src2d->x, src2d->z);
}

void prertm2d_close()
/*< free allocated memory >*/
{
    free(*padvv);
    free(padvv);
    free(*sp0);
    free(sp0);
    free(*sp1);
    free(sp1);
    free(*gp0);
    free(gp0);
    free(*gp1);
    free(gp1);
    free(rwbndr);
    free(bndr);
    // free(abc->bzl);
    // free(abc->bzh);
    // free(abc->bxl);
    // free(abc->bxh);
    // free(abc);
    // free(spo->w);
    // free(spo);
}

void prertm2d_loop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< linear operator >*/
{
    int it, ix, iz;
    float **ptr = NULL;
    sf_adjnull(adj, add, nm, nd, mod, dat);

    float **Iss = sf_floatalloc2(nzpad, nxpad);
    float **Isg = sf_floatalloc2(nzpad, nxpad);

    // source loop

    // initialize is-th source wavefield
    memset(sp0[0], 0, nzpad * nxpad * sizeof(float));
    memset(sp1[0], 0, nzpad * nxpad * sizeof(float));
    memset(gp0[0], 0, nzpad * nxpad * sizeof(float));
    memset(gp1[0], 0, nzpad * nxpad * sizeof(float));

    if (adj) { /* migration m = Lt d */
        for (it = 0; it < nt; it++) {
            if (verb) sf_warning("Migration: %d;", it + 1);
            // sp0: U@t+1 and U@t-1
            // sp1: U@t
            add_source(adj, src2d, sp1, wlt[it]);
            step_forward(adj, sp0, sp1);
            // abcone2d_apply(sp0, sp1, 2, abc);
            apply_sponge2d(sp0);
            apply_sponge2d(sp1);
            ptr = sp0;
            sp0 = sp1;
            sp1 = ptr;
            boundary_rw(false, sp0, &rwbndr[it * 4 * (nx + nz)]);
        }

        for (it = nt - 1; it >= 0; --it) {
            //  reverse time order Img[] += Ps[] * Pg[];
            if (verb) sf_warning("Migration: %d;", it + 1);
            // reconstruct forward wavefield
            boundary_rw(true, sp0, &rwbndr[it * 4 * (nx + nz)]);
            ptr = sp0;
            sp0 = sp1;
            sp1 = ptr;
            step_forward(adj, sp0, sp1);
            add_source(false, src2d, sp1, wlt[it]);
            // backward wavefield
            for (int ir = 0; ir < nr; ++ir) {
                int rx = (int)((rec2d[ir].x - ox) / dx) + nb;
                int rz = (int)((rec2d[ir].z - oz) / dz) + nb;
                gp1[rx][rz] += dat[it * nr + ir];
            }
            //
            step_forward(adj, gp0, gp1);
            // abcone2d_apply(gp0, gp1, 2, abc);
            apply_sponge2d(gp0);
            apply_sponge2d(gp1);
            ptr = gp0;
            gp0 = gp1;
            gp1 = ptr;

            // cross-correlation image condition with source illumination
            // normalization
            for (ix = 0; ix < nx; ++ix) {
                for (iz = 0; iz < nz; ++iz) {
                    Iss[ix][iz] +=
                        sp1[ix + nb][iz + nb] * sp1[ix + nb][iz + nb];
                    Isg[ix][iz] +=
                        gp1[ix + nb][iz + nb] * sp1[ix + nb][iz + nb];
                }
            }
        } // end of it
        for (ix = 0; ix < nx; ++ix) {
            for (iz = 0; iz < nz; ++iz) {
                // apply IC
                if (illum) { // with source illumination
                    mod[iz + ix * nz] = Isg[ix][iz] / (Iss[ix][iz] + EPS);
                } else { // without source illumination
                    mod[iz + ix * nz] = Isg[ix][iz];
                }
            }
        }

        // laplacian filter
        // laplac2_init(nz, nx);
        // laplac2_lop(false, false, nm, nm, pp, mod);

        // sf_warning(".");

    } else { // Born modeling/demigration: m = L d
        for (it = 0; it < nt; ++it) {
            // forward time order: Pg[] += Ps[] * Img[];
            if (verb) sf_warning("Modeling: %d;", it + 1);

            add_source(true, src2d, sp1, wlt[it]);
            step_forward(true, sp0, sp1);
            apply_sponge2d(sp0);
            apply_sponge2d(sp1);
            ptr = sp0;
            sp0 = sp1;
            sp1 = ptr;
            for (ix = 0; ix < nx; ++ix)
                for (iz = 0; iz < nz; ++iz)
                    gp1[ix + nb][iz + nb] +=
                        mod[ix * nz + iz] * sp0[ix + nb][iz + nb];

            step_forward(adj, gp0, gp1);
            apply_sponge2d(gp0);
            apply_sponge2d(gp1);
            ptr = gp0;
            gp0 = gp1;
            gp1 = ptr;

            for (int ir = 0; ir < nr; ++ir) {
                int rx = (int)((rec2d[ir].x - ox) / dx) + nb;
                int rz = (int)((rec2d[ir].z - oz) / dz) + nb;
                dat[it * nr + ir] += gp1[rx][rz];
            }
        } // end of it
        // sf_warning(".");
    } // end of if
    free(*Iss);
    free(Iss);
    free(*Isg);
    free(Isg);
}

void expand2d(float **padvv, float **vv)
/*< expand domain of 'a' to 'b': sources(a)-->destination(b) >*/
{
    int ix, iz;

#ifdef _OPENMP
#pragma omp parallel for private(ix, iz) shared(nx, nz, padvv, vv, nb)
#endif
    for (ix = 0; ix < nx; ++ix) {
        for (iz = 0; iz < nz; ++iz) {
            padvv[ix + nb][iz + nb] = vv[ix][iz]; // centered area
        }
    }

    // for (ix = nb; ix < nx + nb; ++ix) {
    //     for (iz = 0; iz < nb; ++iz) {
    //         padvv[ix][iz] = padvv[ix][nb];                    // top
    //         padvv[ix][iz + nz + nb] = padvv[ix][nz + nb - 1]; // bottom
    //     }
    // }
    for (ix = 0; ix < nxpad; ix++) {
        for (iz = 0; iz < nb; iz++) {
            padvv[ix][iz] = padvv[ix][nb];
            padvv[ix][nzpad - iz - 1] = padvv[ix][nzpad - nb - 1];
        }
    }

    for (ix = 0; ix < nb; ++ix) {
        for (iz = 0; iz < nzpad; ++iz) {
            padvv[ix][iz] = padvv[nb][iz];                         // left
            padvv[nxpad - ix - 1][iz] = padvv[nxpad - nb - 1][iz]; // right
        }
    }
}

void add_source(bool add, pt2d *src2d, float **u, float wlt)
/*< add/substract seismic sources in grid >*/
{
    int sx = (int)((src2d->x - ox) / dx) + nb;
    int sz = (int)((src2d->z - oz) / dz) + nb;
    if (add) { // add sources
        u[sx][sz] -= wlt;
    } else { // substract sources
        u[sx][sz] += wlt;
    }
}

static void step_forward(bool adj, float **u0, float **u1)
{
    int ix, iz;

    if (adj) {
#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix, iz)                         \
    shared(nxpad, nzpad, u0, u1, padvv, c0, c11, c12, c21, c22)
#endif
        for (ix = 2; ix < nxpad - 2; ix++)
            for (iz = 2; iz < nzpad - 2; iz++) {
                u0[ix][iz] =
                    2. * u1[ix][iz] - u0[ix][iz] +
                    padvv[ix][iz] * (c0 * u1[ix][iz] +
                                     c11 * (u1[ix][iz - 1] + u1[ix][iz + 1]) +
                                     c12 * (u1[ix][iz - 2] + u1[ix][iz + 2]) +
                                     c21 * (u1[ix - 1][iz] + u1[ix + 1][iz]) +
                                     c22 * (u1[ix - 2][iz] + u1[ix + 2][iz]));
            }
    } else {
#ifdef _OPENMP
#pragma omp parallel for default(none) private(ix, iz)                         \
    shared(nxpad, nzpad, u0, u1, padvv, c0, c11, c12, c21, c22)
#endif
        for (ix = 2; ix < nxpad - 2; ix++)
            for (iz = 2; iz < nzpad - 2; iz++) {
                u0[ix][iz] = 2. * u1[ix][iz] - u0[ix][iz] +
                             c0 * padvv[ix][iz] * u1[ix][iz] +
                             c11 * (padvv[ix][iz - 1] * u1[ix][iz - 1] +
                                    padvv[ix][iz + 1] * u1[ix][iz + 1]) +
                             c12 * (padvv[ix][iz - 2] * u1[ix][iz - 2] +
                                    padvv[ix][iz + 2] * u1[ix][iz + 2]) +
                             c21 * (padvv[ix - 1][iz] * u1[ix - 1][iz] +
                                    padvv[ix + 1][iz] * u1[ix + 1][iz]) +
                             c22 * (padvv[ix - 2][iz] * u1[ix - 2][iz] +
                                    padvv[ix + 2][iz] * u1[ix + 2][iz]);
            }
    }
}

void apply_sponge2d(float **u)
/*< apply sponge (Gaussian taper) absorbing boundary condition
L=Gaussian taper ABC; L=L*, L is self-adjoint operator >*/
{
    int ib, ibx, ibz, ix, iz;
    float w;

    for (ib = 0; ib < nb; ++ib) {
        w = bndr[ib];
        ibz = nzpad - ib - 1;
        for (ix = 0; ix < nxpad; ix++) {
            u[ix][ib] *= w;  /*    top sponge */
            u[ix][ibz] *= w; /* bottom sponge */
        }

        ibx = nxpad - ib - 1;
        for (iz = 0; iz < nzpad; iz++) {
            u[ib][iz] *= w;  /*   left sponge */
            u[ibx][iz] *= w; /*  right sponge */
        }
    }
}

void boundary_rw(bool read, float **p, float *boundary)
/*< read/write using effective boundary saving strategy;
    if read=true, read the boundary out; else save/write the boundary >*/
{
    int ix, iz;

    if (read) { // read the boundary out
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz) shared(nx, nz, nb, p, boundary)
#endif
        for (ix = 0; ix < nx; ix++)
            for (iz = 0; iz < 2; iz++) {
                p[ix + nb][iz - 2 + nb] = boundary[iz + 4 * ix];
                p[ix + nb][iz + nz + nb] = boundary[iz + 2 + 4 * ix];
            }
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz) shared(nx, nz, nb, p, boundary)
#endif
        for (iz = 0; iz < nz; iz++)
            for (ix = 0; ix < 2; ix++) {
                p[ix - 2 + nb][iz + nb] = boundary[4 * nx + iz + nz * ix];
                p[ix + nx + nb][iz + nb] =
                    boundary[4 * nx + iz + nz * (ix + 2)];
            }
    } else { // save/write boundary
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz) shared(nx, nz, nb, p, boundary)
#endif
        for (ix = 0; ix < nx; ix++)
            for (iz = 0; iz < 2; iz++) {
                boundary[iz + 4 * ix] = p[ix + nb][iz - 2 + nb];
                boundary[iz + 2 + 4 * ix] = p[ix + nb][iz + nz + nb];
            }
#ifdef _OPENMP
#pragma omp parallel for private(ix, iz) shared(nx, nz, nb, p, boundary)
#endif
        for (iz = 0; iz < nz; iz++)
            for (ix = 0; ix < 2; ix++) {
                boundary[4 * nx + iz + nz * ix] = p[ix - 2 + nb][iz + nb];
                boundary[4 * nx + iz + nz * (ix + 2)] =
                    p[ix + nx + nb][iz + nb];
            }
    }
}
