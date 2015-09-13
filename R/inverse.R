# Copyright (C) 2015 Philipp Benner
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

cholesky.update <- function(L11, L11inv, A12, A22) {
    # derived from the blocked cholesky algorithm
    if (!is.matrix(A12)) {
        A12 <- as.matrix(A12)
    }
    if (!is.matrix(A22)) {
        A22 <- as.matrix(A22)
    }
    if (is.null(L11) && is.null(L11inv) && is.null(A12)) {
        return (chol(A22))
    }
    L12 <- t(L11inv)%*%A12
    L22 <- chol(A22 - t(L12)%*%L12)
    M   <- matrix(0.0, nrow(L11)+nrow(A22),ncol(L11)+ncol(A22))
    M[1:nrow(L11),           1:ncol(L11)          ] <- L11
    M[1:nrow(L11),           (ncol(L11)+1):ncol(M)] <- L12
    M[(nrow(L11)+1):nrow(M), (ncol(L11)+1):ncol(M)] <- L22
    M
}

block.inversion <- function(A, B, C, D, Ainv, Dinv) {
    if (is.null(A) && is.null(B) && is.null(C)) {
        return (Dinv)
    }
    if (!is.matrix(A)) {
        A <- as.matrix(A)
    }
    if (!is.matrix(D)) {
        D <- as.matrix(D)
    }
    if (!is.matrix(B)) {
        B <- matrix(B, nrow(A), ncol(D))
    }
    if (!is.matrix(C)) {
        C <- matrix(C, nrow(D), ncol(A))
    }
    Sa <- D - C%*%Ainv%*%B
    Sd <- A - B%*%Dinv%*%C
    Sainv <- solve(Sa)
    Sdinv <- solve(Sd)
    M <- matrix(0.0, nrow(A)+nrow(C), ncol(A)+ncol(B))
    M[1:nrow(A),           1:ncol(A)          ] <- Sdinv
    M[(nrow(A)+1):nrow(M), (ncol(A)+1):ncol(M)] <- Sainv
    M[1:nrow(A),           (ncol(A)+1):ncol(M)] <- -Ainv%*%B%*%Sainv
    M[(nrow(A)+1):nrow(M), 1:ncol(A)          ] <- -Dinv%*%C%*%Sdinv
    M
}

cholesky.inverse.update <- function(L11, L11inv, A12, A22) {
    if (is.null(L11) && is.null(L11inv) && is.null(A12)) {
        K    <- chol(A22)
        Kinv <- solve(K)
    }
    else {
        K    <- cholesky.update(L11, L11inv, A12, A22)
        Kinv <- block.inversion(L11,
                                K[1:nrow(L11),(ncol(L11)+1):ncol(K)],
                                K[(nrow(L11)+1):nrow(K),1:ncol(L11)],
                                K[(nrow(L11)+1):nrow(K),(ncol(L11)+1):ncol(K)],
                                L11inv,
                                solve(K[(nrow(L11)+1):nrow(K),(ncol(L11)+1):ncol(K)]))
    }
    list(K=K, Kinv=Kinv)
}

# ------------------------------------------------------------------------------
if (FALSE) {
    M <- matrix(c(18, 22,   54,  42,
                  22, 70,   86,  62,
                  54, 86,  174, 134,
                  42, 62,  134, 106), 4,4)
    A <- M[1:2,1:2]

    L    <- chol(A)
    Linv <- solve(L)

    # chol(M) == cholesky.update(L, Linv, M[1:2,3:4], M[3:4,3:4])
    K    <- chol(M)
    Kinv <- solve(K)

    # compute Kinv with the block inversion formula
    block.inversion(K[1:2,1:2], K[1:2,3:4], K[3:4,1:2], K[3:4,3:4],
                    solve(K[1:2,1:2]), solve(K[3:4,3:4]))

    # combining both functions to compute K and Kinv directly
    # from L and Linv given the missing blocks of M...
    cholesky.inverse.update(L, Linv, M[1:2,3:4], M[3:4,3:4])
}
