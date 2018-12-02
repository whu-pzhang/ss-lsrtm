#include <rsf.h>

#include "prertm2d_v03.h"

#ifndef EPS
#define EPS FLT_EPSILON
#endif

typedef struct {
    bool free;
    float *bzl;
    float *bzh;
    float *bxl;
    float *bxh;
} abcone2d;


typedef struct {
    float *w;
} sponge;


static float c0, c11, c12, c21, c22;
static bool verb;
static int nx, nz, nb, nxpad, nzpad, nt;
static int ns, nss, sns, nr;
static float dx, dz, ox, oz;
static float **padvv, **su0, **su1, **gu0, **gu1, *rwbndr, *wlt;
static pt2d *src2d, *rec2d;
static abcone2d *abc;
static sponge *spo;

static void expand2d(float **padvv, float **vv);
static abcone2d *abcone2d_make(float dx, float dz, float dt, float **padvv);
static sponge *sponge_make(int nb);
static void step_forward(bool adj, float **u0, float **u1);
static void inject_source(bool adj, pt2d *src2d, float **u, int sns, float wlt);
static void abcone2d_apply(float **uo, float **um, int nop, abcone2d *abc);
static void sponge2d_apply(float **u0, float **u1, sponge *spo);
static void boundary_rw(bool read, float **p, float *boundary);

void prertm2d_init(bool verb_, int nx_, int nz_, int nb_, int nt_, int ns_,
                   int nss_, int nr_, float dx_, float dz_, float dt_,
                   float ox_, float oz_, float **vv, float *wlt_, pt2d *src2d_,
                   pt2d *rec2d_)
/*< init >*/
{
    verb = verb_;
    nx = nx_;
    nz = nz_;
    nb = nb_;
    nt = nt_;
    ns = ns_;
    nss = nss_;  // number of ss
    sns = ns / nss;
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

    padvv = sf_floatalloc2(nzpad, nxpad);
    su0 = sf_floatalloc2(nzpad, nxpad);
    su1 = sf_floatalloc2(nzpad, nxpad);
    gu0 = sf_floatalloc2(nzpad, nxpad);
    gu1 = sf_floatalloc2(nzpad, nxpad);
    rwbndr = sf_floatalloc(nt * 4 * (nx + nz));

    // precompute vp^2 * dt*2
    for (int ix = 0; ix < nx * nz; ++ix)
        *(vv[0] + ix) *= *(vv[0] + ix) * dt_ * dt_;  // vel=vel^2*dt^2
    expand2d(padvv, vv);

    // absorbing boundary
    abc = abcone2d_make(dx, dz, dt_, padvv);
    spo = sponge_make(nb);
}

void prertm2d_close()
/*< free allocated memory >*/
{
    free(*padvv);
    free(padvv);
    free(*su0);
    free(su0);
    free(*su1);
    free(su1);
    free(*gu0);
    free(gu0);
    free(*gu1);
    free(gu1);
    free(rwbndr);
    free(abc->bzl);
    free(abc->bzh);
    free(abc->bxl);
    free(abc->bxh);
    free(abc);
    free(spo->w);
    free(spo);
}

void prertm2d_loop(bool adj, bool add, int nm, int nd, float *mod, float *dat)
/*< linear operator >*/
{
    int it, ix, iz;
    float **ptr = NULL;
    sf_adjnull(adj, add, nm, nd, mod, dat);

    float **Iss = sf_floatalloc2(nzpad, nxpad);
    float **Isg = sf_floatalloc2(nzpad, nxpad);

    for (int iss = 0; iss < nss; ++iss) {  // loop over shot(supergather)
        if (verb) sf_warning("%d", iss + 1);
        // initialize is-th source wavefield
        memset(su0[0], 0, nzpad * nxpad * sizeof(float));
        memset(su1[0], 0, nzpad * nxpad * sizeof(float));
        memset(gu0[0], 0, nzpad * nxpad * sizeof(float));
        memset(gu1[0], 0, nzpad * nxpad * sizeof(float));

        if (adj) {  // migration
            for (it = 0; it < nt; ++it) {
                if (verb) sf_warning("Migration: %d;", it + 1);
                // su0: U@t+1 and U@t-1
                // su1: U@t
                inject_source(adj, &src2d[iss * sns], su1, sns, wlt[it]);
                step_forward(adj, su0, su1);

                abcone2d_apply(su0, su1, 2, abc);
                sponge2d_apply(su0, su1, spo);
                ptr = su0;
                su0 = su1;
                su1 = ptr;
                boundary_rw(false, su0, &rwbndr[it * 4 * (nx + nz)]);
            }

            for (it = nt - 1; it >= 0; --it) {
                if (verb) sf_warning("Migration: %d;", it + 1);
                // reconstruct forward wavefield
                boundary_rw(true, su0, &rwbndr[it * 4 * (nx + nz)]);
                ptr = su0;
                su0 = su1;
                su1 = ptr;
                step_forward(adj, su0, su1);
                inject_source(false, &src2d[iss * sns], su1, sns, wlt[it]);
                // backward wavefield
                // inject data
                for (int ir = 0; ir < nr; ++ir) {
                    int rx = (int)((rec2d[ir].x - ox) / dx) + nb;
                    int rz = (int)((rec2d[ir].z - oz) / dz) + nb;
                    gu1[rx][rz] += dat[iss * nr * nt + ir * nt + it];
                }
                //
                step_forward(adj, gu0, gu1);

                abcone2d_apply(gu0, gu1, 2, abc);
                sponge2d_apply(gu0, gu1, spo);

                ptr = gu0;
                gu0 = gu1;
                gu1 = ptr;

                // cross-correlation image condition
                // with source illumination normalization
                for (ix = 0; ix < nx; ++ix) {
                    for (iz = 0; iz < nz; ++iz) {
                        Iss[ix][iz] +=
                            su1[ix + nb][iz + nb] * su1[ix + nb][iz + nb];
                        Isg[ix][iz] +=
                            gu1[ix + nb][iz + nb] * su1[ix + nb][iz + nb];
                    }
                }
            } // end of it
            for (ix = 0; ix < nx; ++ix) {
                for (iz = 0; iz < nz; ++iz) {
                    // apply IC
                    mod[iz + ix * nz] =
                        Isg[ix][iz] / (Iss[ix][iz] + EPS);
                }
            }
            sf_warning(".");
        } else {  // modeling
            for (it = 0; it < nt; ++it) {
                if (verb) sf_warning("Modeling: %d;", it + 1);
                // source wavefield
                inject_source(true, &src2d[iss * sns], su1, sns, wlt[it]);
                step_forward(true, su0, su1);

                abcone2d_apply(su0, su1, 2, abc);
                sponge2d_apply(su0, su1, spo);
                ptr = su0;
                su0 = su1;
                su1 = ptr;

                // 计算二次震源,即 reflector
                for (ix = 0; ix < nx; ++ix)
                    for (iz = 0; iz < nz; ++iz)
                        gu1[ix + nb][iz + nb] +=
                            mod[ix * nz + iz] * su1[ix + nb][iz + nb];

                // reflector progragation
                step_forward(adj, gu0, gu1);

                // absorbing boundary
                abcone2d_apply(gu0, gu1, 2, abc);
                sponge2d_apply(gu0, gu1, spo);

                ptr = gu0;
                gu0 = gu1;
                gu1 = ptr;

                // record data
                for (int ir = 0; ir < nr; ++ir) {
                    int rx = (int)((rec2d[ir].x - ox) / dx) + nb;
                    int rz = (int)((rec2d[ir].z - oz) / dz) + nb;
                    dat[iss * nr * nt + ir * nt + it] += gu1[rx][rz];
                }
            }  // end of it
            sf_warning(".");
        }  // end of if
    }
    free(*Iss); free(Iss);
    free(*Isg); free(Isg);
}

static void expand2d(float **padvv, float **vv)
{
    int ix, iz;

#ifdef _OPENMP
    #pragma omp parallel for default(none) private(ix, iz) shared(nx, nz, padvv,   \
    vv, nb)
#endif
    for (ix = 0; ix < nx; ++ix) {
    for (iz = 0; iz < nz; ++iz) {
            padvv[ix + nb][iz + nb] = vv[ix][iz];  //
        }
    }

    for (ix = nb; ix < nx + nb; ++ix) {
    for (iz = 0; iz < nb; ++iz) {
            padvv[ix][iz] = padvv[ix][nb];  // top
            padvv[ix][iz + nz + nb] = padvv[ix][nz + nb - 1];  // bottom
        }
    }

    for (ix = 0; ix < nb; ++ix) {
    for (iz = 0; iz < nzpad; ++iz) {
            padvv[ix][iz] = padvv[nb][iz];  // left
            padvv[ix + nx + nb][iz] = padvv[nx + nb - 1][iz];  // right
        }
    }
}

static abcone2d *abcone2d_make(float dx, float dz, float dt, float **padvv)
{
    abcone2d *abc;
    int iz, ix;
    float d;

    abc = (abcone2d *)sf_alloc(1, sizeof(abcone2d));

    abc->bzl = sf_floatalloc(nxpad);
    abc->bzh = sf_floatalloc(nxpad);
    abc->bxl = sf_floatalloc(nzpad);
    abc->bxh = sf_floatalloc(nzpad);

    for (ix = 0; ix < nxpad; ++ix) {
        d = padvv[ix][nb] * dt / dz;
        abc->bzl[ix] = (1 - d) / (1 + d);  // top
        d = padvv[ix][nzpad - nb - 1] * dt / dz;
        abc->bzh[ix] = (1 - d) / (1 + d);  // bottom
    }
    for (iz = 0; iz < nzpad; ++iz) {
        d = padvv[nb][iz] * dt / dx;
        abc->bxl[iz] = (1 - d) / (1 + d);  // left
        d = padvv[nxpad - nb - 1][iz] * dt / dx;
        abc->bxh[iz] = (1 - d) / (1 + d);  // right
    }

    return abc;
}

static sponge *sponge_make(int nb)
{
    sponge *spo;
    int ib;
    float sb, fb;

    spo = (sponge *)sf_alloc(1, sizeof(sponge));
    spo->w = sf_floatalloc(nb);
    sb = 4.0 * nb;
    for (ib = 0; ib < nb; ++ib) {
        fb = ib / (sqrt(2.0) * sb);
        spo->w[ib] = exp(-fb * fb);
    }
    return spo;
}

static void inject_source(bool adj, pt2d *src2d, float **u, int sns, float wlt)
{
    for (int is = 0; is < sns; ++is) {
        int sx = (int)((src2d[is].x - ox) / dx) + nb;
        int sz = (int)((src2d[is].z - oz) / dz) + nb;
        if (adj) {
            u[sx][sz] -= wlt;
        } else {
            u[sx][sz] += wlt;
        }
    }
}

static void step_forward(bool adj, float **u0, float **u1)
{
    int ix, iz;

    if (adj) {
#ifdef _OPENMP
        #pragma omp parallel for default(none) private(ix, iz) shared(                 \
        nxpad, nzpad, u0, u1, padvv, c0, c11, c12, c21, c22)
#endif
        for (ix = 2; ix < nxpad - 2; ix++)
    for (iz = 2; iz < nzpad - 2; iz++) {
        u0[ix][iz] = 2. * u1[ix][iz] - u0[ix][iz] +
                         padvv[ix][iz] * (c0 * u1[ix][iz] +
                                          c11 * (u1[ix][iz - 1] + u1[ix][iz + 1]) +
                                          c12 * (u1[ix][iz - 2] + u1[ix][iz + 2]) +
                                          c21 * (u1[ix - 1][iz] + u1[ix + 1][iz]) +
                                          c22 * (u1[ix - 2][iz] + u1[ix + 2][iz]));
        }
    } else {
#ifdef _OPENMP
        #pragma omp parallel for default(none) private(ix, iz) shared(                 \
        nxpad, nzpad, u0, u1, padvv, c0, c11, c12, c21, c22)
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

static void abcone2d_apply(float **uo, float **um, int nop, abcone2d *abc)
{
    int iz, ix, iop;

    for (ix = 0; ix < nxpad; ++ix) {
        for (iop = 0; iop < nop; ++iop) {

            /* top BC */
            if (!abc->free) { /* not free surface, apply ABC */
                iz = nop - iop;
                uo[ix][iz] = um[ix][iz + 1] +
                             (um[ix][iz] - uo[ix][iz + 1]) * abc->bzl[ix];
            }

            /* bottom BC */
            iz = nzpad - nop + iop - 1;
            uo[ix][iz] =
                um[ix][iz - 1] + (um[ix][iz] - uo[ix][iz - 1]) * abc->bzh[ix];
        }
    }

    for (iz = 0; iz < nzpad; iz++) {
        for (iop = 0; iop < nop; iop++) {

            /* left BC */
            ix = nop - iop;
            uo[ix][iz] =
                um[ix + 1][iz] + (um[ix][iz] - uo[ix + 1][iz]) * abc->bxl[iz];

            /* right BC */
            ix = nxpad - nop + iop - 1;
            uo[ix][iz] =
                um[ix - 1][iz] + (um[ix][iz] - uo[ix - 1][iz]) * abc->bxh[iz];
        }
    }
}

static void sponge2d_apply(float **u0, float **u1, sponge *spo)
{
    float w;

    for (int ib = 0; ib < nb; ++ib) {
        w = spo->w[nb - ib - 1];

        int ibz = nzpad - ib - 1;
        for (int ix = 0; ix < nxpad; ix++) {
            u0[ix][ib] *= w; /*    top sponge */
            u1[ix][ib] *= w;
            u0[ix][ibz] *= w; /* bottom sponge */
            u1[ix][ibz] *= w;
        }
    }

    for (int ib = 0; ib < nb; ++ib) {
        w = spo->w[nb - ib - 1];

        int ibx = nxpad - ib - 1;
        for (int iz = 0; iz < nzpad; iz++) {
            u0[ib][iz] *= w; /*   left sponge */
            u1[ib][iz] *= w;
            u0[ibx][iz] *= w; /*  right sponge */
            u1[ibx][iz] *= w;
        }
    }
}

static void boundary_rw(bool read, float **p, float *boundary)
{
    int ix, iz;

    if (read) {

        for (ix = 0; ix < nx; ix++)
            for (iz = 0; iz < 2; iz++) {
                p[ix + nb][iz - 2 + nb] = boundary[iz + 4 * ix];
                p[ix + nb][iz + nz + nb] = boundary[iz + 2 + 4 * ix];
            }

        for (iz = 0; iz < nz; iz++)
            for (ix = 0; ix < 2; ix++) {
                p[ix - 2 + nb][iz + nb] = boundary[4 * nx + iz + nz * ix];
                p[ix + nx + nb][iz + nb] =
                    boundary[4 * nx + iz + nz * (ix + 2)];
            }
    } else {

        for (ix = 0; ix < nx; ix++)
            for (iz = 0; iz < 2; iz++) {
                boundary[iz + 4 * ix] = p[ix + nb][iz - 2 + nb];
                boundary[iz + 2 + 4 * ix] = p[ix + nb][iz + nz + nb];
            }

        for (iz = 0; iz < nz; iz++)
            for (ix = 0; ix < 2; ix++) {
                boundary[4 * nx + iz + nz * ix] = p[ix - 2 + nb][iz + nb];
                boundary[4 * nx + iz + nz * (ix + 2)] =
                    p[ix + nx + nb][iz + nb];
            }
    }
}
