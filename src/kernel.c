/* Copyright (C) 2013-2015 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <math.h>

// naive implementation of the kernel function
////////////////////////////////////////////////////////////////////////////////

SEXP exponential_kernel(SEXP x, SEXP y, SEXP l, SEXP var)
{
        R_len_t i, j, k;
        R_len_t nx;
        R_len_t mx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rans, r;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if (length(dim) != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];
        mx = INTEGER(dim)[1];

        dim = getAttrib(y, R_DimSymbol);
        if (length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        r = 0.0;
                        for (k = 0; k < mx; k++) {
                                r += (rx[i + nx*k] - ry[j + ny*k])*(rx[i + nx*k] - ry[j + ny*k]);
                        }
                        rans[i + nx*j] =
                                (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*r);
                }
        }
        UNPROTECT(1);

        return(ans);
}

SEXP exponential_kernel_gradient(SEXP x, SEXP y, SEXP l, SEXP var, SEXP li)
{
        R_len_t i, j, k;
        R_len_t nx;
        R_len_t mx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        int    *ri   = INTEGER(li);
        double *rans, r;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if (length(dim) != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];
        mx = INTEGER(dim)[1];

        dim = getAttrib(y, R_DimSymbol);
        if (length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        if (*ri == 1) {
                for(i = 0; i < nx; i++) {
                        for(j = 0; j < ny; j++) {
                                r = 0.0;
                                for (k = 0; k < mx; k++) {
                                        r += (rx[i + nx*k] - ry[j + ny*k])*(rx[i + nx*k] - ry[j + ny*k]);
                                }
                                rans[i + nx*j] =
                                        (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*r);
                                rans[i + nx*j] =
                                        rans[i + nx*j]*r/((*rl)*(*rl)*(*rl));
                        }
                }
        }
        else if (*ri == 2) {
                for(i = 0; i < nx; i++) {
                        for(j = 0; j < ny; j++) {
                                r = 0.0;
                                for (k = 0; k < mx; k++) {
                                        r += (rx[i + nx*k] - ry[j + ny*k])*(rx[i + nx*k] - ry[j + ny*k]);
                                }
                                rans[i + nx*j] =
                                        2*sqrt(*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*r);
                        }
                }
        }
        else {
                error("invalid argument i");
        }
        UNPROTECT(1);

        return(ans);
}

////////////////////////////////////////////////////////////////////////////////

SEXP ornstein_uhlenbeck_kernel(SEXP x, SEXP y, SEXP l)
{
        R_len_t i, j, k;
        R_len_t nx;
        R_len_t mx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rans, r;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if (length(dim) != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];
        mx = INTEGER(dim)[1];

        dim = getAttrib(y, R_DimSymbol);
        if (length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        r = 0.0;
                        for (k = 0; k < mx; k++) {
                                r += (rx[i + nx*k] - ry[j + ny*k])*(rx[i + nx*k] - ry[j + ny*k]);
                        }
                        rans[i + nx*j] = exp(-sqrt(r)/(*rl));
                }
        }
        UNPROTECT(1);

        return(ans);
}

// Matern kernel
////////////////////////////////////////////////////////////////////////////////

SEXP matern_kernel(SEXP x, SEXP y, SEXP l, SEXP nu)
{
        R_len_t i, j, k;
        R_len_t nx;
        R_len_t mx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rnu  = REAL(nu);
        double *rans, r, t, c;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if (length(dim) != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];
        mx = INTEGER(dim)[1];

        dim = getAttrib(y, R_DimSymbol);
        if (length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(nu) != 1) {
                error("nu is not a scalar");
        }

        /* required constant */
        c = pow(2.0, 1.0-(*rnu))/gammafn(*rnu);
        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        r = 0.0;
                        for (k = 0; k < mx; k++) {
                                r += sqrt((rx[i + nx*k] - ry[j + ny*k])*(rx[i + nx*k] - ry[j + ny*k]));
                        }
                        if (sqrt(r) < 1.0e-8) {
                                rans[i + nx*j] = 1.0;
                        }
                        else {
                                t = sqrt(2.0*(*rnu))/(*rl)*r;
                                rans[i + nx*j] = c*pow(t, (*rnu))*bessel_k(t, (*rnu), 1);
                        }
                }
        }
        UNPROTECT(1);

        return(ans);
}
