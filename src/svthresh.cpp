#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;


//soft thresholding operator of symmetric matrix
arma::sp_mat SOFTTHRESH(arma::mat Sig, float lambda, bool shrinkDiag = true) {
    arma::mat A = symmatu(Sig);
    int n = A.n_cols;
    arma::sp_mat S(sign(A));
    arma::mat    M(abs(A)-lambda);
    M.elem(find(M < 0)).zeros();
    if (!shrinkDiag)
      M.diag() = A.diag();
    return M % S;
}


List softSVT(arma::mat M, float tau=0, int k=0) {
    arma::mat U, V;
    arma::vec d;

    svd(U, d, V, M);
    if (k != 0)
      tau = d(k); //*1.01;
    arma::vec tmpd = d - tau;
    tmpd.elem(find(tmpd < 0)).zeros();
    M = U*diagmat(tmpd)*V.t();
    return List::create(Named("M")   = M,
                        Named("tau") = tau,
                        Named("d")   = tmpd);
}
