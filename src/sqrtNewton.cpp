#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


// Compute the square root of a matrix using Newton's method
// does arma has a faster implementation? (yes)
arma::mat sqrtmNewt_dense(arma::mat C, arma::mat sqrt0, double errTol) {
    arma::mat X = sqrt0;
    arma::mat X_new;
//    arma::mat XC;
    double err = std::numeric_limits<double>::infinity();
    while (err > errTol) {
        X_new = 0.5*(X + solve(X, C, solve_opts::fast));
//        X_new = 0.5*(X + XC);
        err   = norm(X_new - X, "fro");
        X     = X_new;
    }
    return X_new;
}


arma::mat sqrtmNewt_sp(arma::mat C, arma::sp_mat sqrt0, double errTol) {
    arma::sp_mat X = sqrt0;
    arma::mat X_new(X);
//    arma::mat XC;
    double err = std::numeric_limits<double>::infinity();
    while (err > errTol) {
        X_new =  0.5*(X_new + spsolve(X, C, "lapack"));
//        X_new = + XC);
        err   = norm(X_new - X, "fro");
        X     = sp_mat(X_new);
    }
    return X_new;
}


arma::mat sqrtmNewt(arma::mat C, SEXP sqrt0, double errTol = 1e-3) {
    if (Rf_isS4(sqrt0)) {
        if (Rf_inherits(sqrt0, "sparseMatrix")) {
            //TODO: make this faster!
            return sqrtmNewt_sp(C, as<arma::sp_mat>(sqrt0), errTol);
        } ;
        stop("matrix has unknown class") ;
    } else {
        return sqrtmNewt_dense(C, as<arma::mat>(sqrt0), errTol);
    }
}



arma::mat solveCpp(arma::mat C, arma::mat X) {
    arma::mat out;
    out = solve(C, X);
    return out;
}
