#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
 
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::mat sqrtmNewt2(arma::mat& C, arma::mat& sqrt0, const double& errTol) {
    arma::mat X(sqrt0);
    arma::mat X_new;
//    arma::mat XC;
    double err = std::numeric_limits<double>::infinity();
    while (err > errTol) {
        X_new = 0.5*(X + solve(X, C, solve_opts::fast));
        err   = norm(X_new - X, "fro");
        X     = X_new;
    }
    return X_new;
}

// [[Rcpp::export]]
arma::sp_mat SOFTTHRESH2(arma::mat& Sig, const float& lambda, const bool& shrinkDiag = true) {
    arma::mat A = symmatu(Sig);
    int n = A.n_cols;
    arma::sp_mat S(sign(A));
    arma::mat    M(abs(A)-lambda);
    M.elem(find(M < 0)).zeros();
    if (!shrinkDiag)
      M.diag() = A.diag();
    return M % S;
}

// [[Rcpp::export]]
void svdPowSym(arma::mat& U, arma::vec& d, arma::mat& V, arma::mat& A, int k, int q) {
    int l = k;
    int m = A.n_rows;
    int n = A.n_cols;
    arma::mat P = randu<mat>(m,k-1);
    arma::mat Q, R;
    qr(Q, R, (P.t()*A).t());
    for (int j = 0; j<q; j++) {
        qr(P, R, A*Q);
        qr(Q, R, (P.t()*A).t());
    }
    svd_econ(U, d, V, A*Q);
//    d = d.subvec(0,k-1);
//    U = U.cols(0,k-1);
}

// [[Rcpp::export]]
arma::mat softSVT3(arma::mat& M, int k, double beta=0) {
    arma::mat U, V;
    arma::vec d;
    if (beta != 0)
      k = M.n_cols;
    svdPowSym(U, d, V, M, k, k+1);
    if (beta == 0)
      beta = d(k);
    arma::vec tmpd = d - beta;
    tmpd.elem(find(tmpd < 0)).zeros();
    return U*diagmat(tmpd)*V;
}


// [[Rcpp::export]]
arma::mat softSVT2(arma::mat& M, int k, double beta=0) {
    arma::mat U, V;
    arma::vec d;
    svd_econ(U, d, V, M);
    if (beta == 0)
      beta = d(k);
    arma::vec tmpd = d - beta;
    tmpd.elem(find(tmpd < 0)).zeros();
    return U*diagmat(tmpd)*V.t();
}


// [[Rcpp::export]]
List ADMM(const arma::mat& SigmaO, const double& lambda, arma::mat& I, arma::mat& Lambda,
          arma::mat& Y, double beta=0, int r=0, double mu=0, const double& eta=4/5, 
          const double& muf=1e-6, const int& maxiter=500, const double& newtol=1e-5,
          const double& tol=1e-5, const double& over_relax_par=1.6, bool shrinkDiag=true) {

    int n = SigmaO.n_rows;
    if (mu == 0) mu = n;
    arma::vec r_norm(maxiter), eps_pri(maxiter);
    arma::mat Lambda1 = Lambda.cols(0  , n-1  );
    arma::mat Lambda2 = Lambda.cols(n  , n*2-1);
    arma::mat Lambda3 = Lambda.cols(n*2, n*3-1);
    arma::mat R = Y.cols(0,n-1), L = Y.cols(n*2,n*3-1);
    arma::sp_mat S(Y.cols(n,n*2-1));
    arma::mat RY(R), LY(L), SY(S);
    arma::mat RO(size(R)), SO(size(S)), LO(size(L));
    arma::mat RA(size(R)), SA(size(S)), LA(size(L)), TA(size(L));
    arma::mat X(size(Y));
//    arma::mat Y_old(Y);
    arma::mat K, K2, KI;
    int iter=0;
    for (iter; iter < maxiter; iter++) {
      // update X = (R,S,L)
//      if (iter > 0) {
        RA = RY + mu*Lambda1;
        SA = SY + mu*Lambda2;
        LA = LY + mu*Lambda3;
//      }
      // update R
      K  = mu*SigmaO - RA;
      K2 = K*K + 4*mu*I;
//      KI = sqrtmNewt2(K2, I, newtol);
      arma::sqrtmat_sympd(KI, K2);
      R  = (KI-K)/2;

      // update S
      S = SOFTTHRESH2(SA, lambda*mu, shrinkDiag);
      if (iter == 0) {
        L = LA;
      } else {
        L = softSVT2(LA, r, mu*beta);
      }

      X.cols(0  , n-1  ) = R;
      X.cols(n  , n*2-1) = S;
      X.cols(n*2, n*3-1) = L;

      RO = over_relax_par*R + (1-over_relax_par)*RY;
      SO = over_relax_par*S + (1-over_relax_par)*SY;
      LO = over_relax_par*L + (1-over_relax_par)*LY;

      RA = RO - mu*Lambda1;
      SA = SO - mu*Lambda2;
      LA = LO - mu*Lambda3;
      TA = (RA-SA+LA)/3;
      RY = RA - TA;
      SY = SA + TA;
      LY = LA - TA;

     // update Lambda
      Lambda1 = Lambda1 - (RO-RY)/mu;
      Lambda2 = Lambda2 - (SO-SY)/mu;
      Lambda3 = Lambda3 - (LO-LY)/mu;

      Y.cols(0  ,n-1  ) = RY;
      Y.cols(n  ,n*2-1) = SY;
      Y.cols(n*2,n*3-1) = LY;

      r_norm(iter)  = norm(X - Y, "fro");
      eps_pri(iter) = sqrt(3*n*n)*tol + tol*std::max(norm(X,"fro"), norm(Y,"fro"));

      if (r_norm(iter) < eps_pri(iter)) {
         break;
      }
      mu = std::max(mu*eta, muf);

    }

    Lambda.cols(0  ,n-1  ) = Lambda1;
    Lambda.cols(n  ,n*2-1) = Lambda2;
    Lambda.cols(n*2,n*3-1) = Lambda3;

    return List::create(_["R"] = R,  _["S"] = S, _["L"] = L,
                        _["Y"] = Y,
                        _["Lambda"] = Lambda,
//                        _["RA"] = RA, _["SA"] = SA, _["LA"] = LA,
                        _["iter"] = iter+1, _["beta"] = beta, _["mu"]   = mu,
                        _["history"] = List::create(
                          _["r_norm"]  =  r_norm,
                          _["eps_pri"] = eps_pri)
                        );
}

