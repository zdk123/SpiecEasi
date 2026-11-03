#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// LAPACK symbols declared directly (with hidden Fortran string lengths) to avoid
// clashing with Armadillo's own LAPACK prototypes.
extern "C" void dsyevr_(const char* jobz, const char* range, const char* uplo,
                        const int* n, double* a, const int* lda,
                        const double* vl, const double* vu,
                        const int* il, const int* iu, const double* abstol,
                        int* m, double* w, double* z, const int* ldz, int* isuppz,
                        double* work, const int* lwork, int* iwork, const int* liwork,
                        int* info, size_t, size_t, size_t);
extern "C" double dlamch_(const char* cmach, size_t);


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



arma::sp_mat SOFTTHRESH2(arma::mat& Sig, const double& lambda, const bool& shrinkDiag = true) {
    arma::mat A = symmatu(Sig);
//    int n = A.n_cols;
    arma::sp_mat S(sign(A));
    arma::mat    M(abs(A)-lambda);
    M.elem(find(M < 0)).zeros();
    if (!shrinkDiag)
      M.diag() = A.diag();
    return M % S;
}

// Dense form of SOFTTHRESH2, for use inside the ADMM loop (no sparse round trip).
arma::mat SOFTTHRESHD(const arma::mat& Sig, const double& lambda, const bool& shrinkDiag = true) {
    arma::mat A = symmatu(Sig);
    arma::mat M(abs(A)-lambda);
    M.elem(find(M < 0)).zeros();
    if (!shrinkDiag)
      M.diag() = A.diag();
    return M % sign(A);
}


void svdPowSym(arma::mat& U, arma::vec& d, arma::mat& V, arma::mat& A, int k, int q) {
//    int l = k;
    int m = A.n_rows;
//    int n = A.n_cols;
    arma::mat P = randu<mat>(m,k-1);
    arma::mat Q, R;
    qr(Q, R, (P.t()*A).t());
    for (int j = 0; j<q; j++) {
        qr(P, R, A*Q);
        qr(Q, R, (P.t()*A).t());
    }
    svd_econ(U, d, V, A*Q);
}

List svdPowSym2(arma::mat& M, int k) {
    arma::mat U, V;
    arma::vec d;
    svdPowSym(U, d, V, M, k, k+1);
    return List::create(_["d"] = d, _["U"] = U, _["V"] = V) ;
}


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

arma::mat eigPSDThresh2(arma::mat& M, int k, double beta=0) {
    arma::mat A = symmatu(M);
    A = 0.5 * (A + A.t());

    arma::vec d;
    arma::mat Q;
    eig_sym(d, Q, A);

    if (beta == 0) {
      if (k <= 0) {
        beta = std::numeric_limits<double>::infinity();
      } else if (k < static_cast<int>(d.n_elem)) {
        beta = std::max(d(d.n_elem - k - 1), 0.0);
      } else {
        beta = 0.0;
      }
    }

    arma::vec tmpd = d - beta;
    tmpd.elem(find(tmpd < 0)).zeros();
    arma::mat L = Q * diagmat(tmpd) * Q.t();
    return 0.5 * (L + L.t());
}

// Same result as eigPSDThresh2, but asks LAPACK (dsyevr) only for the eigenpairs
// that survive the threshold: the top k+1 when beta == 0, else those above beta.
arma::mat eigPSDThreshPartial(arma::mat& M, int k, double beta=0) {
    int n = M.n_rows;
    if (beta == 0 && k <= 0)
      return arma::zeros<arma::mat>(n, n);
    if (beta == 0 && k >= n)
      return eigPSDThresh2(M, k, beta);

    arma::mat A = symmatu(M);
    A = 0.5 * (A + A.t());

    const char* range = (beta == 0) ? "I" : "V";
    int il = n - k, iu = n, ncol = (beta == 0) ? k + 1 : n;
    double vl = beta, vu = std::numeric_limits<double>::max();
    double abstol = 2 * dlamch_("S", 1);
    int m = 0, info = 0, lwork = -1, liwork = -1, iwq = 0;
    double wq = 0;
    arma::vec w(n);
    arma::mat Z(n, ncol);
    std::vector<int> isuppz(2 * ncol);

    dsyevr_("V", range, "U", &n, A.memptr(), &n, &vl, &vu, &il, &iu, &abstol, &m,
            w.memptr(), Z.memptr(), &n, isuppz.data(), &wq, &lwork, &iwq, &liwork, &info, 1, 1, 1);
    lwork = static_cast<int>(wq); liwork = iwq;
    std::vector<double> work(lwork);
    std::vector<int> iwork(liwork);
    dsyevr_("V", range, "U", &n, A.memptr(), &n, &vl, &vu, &il, &iu, &abstol, &m,
            w.memptr(), Z.memptr(), &n, isuppz.data(), work.data(), &lwork, iwork.data(), &liwork, &info, 1, 1, 1);
    if (info != 0)
      return eigPSDThresh2(M, k, beta);
    if (m == 0)
      return arma::zeros<arma::mat>(n, n);

    arma::vec d = w.head(m);
    arma::mat Q = Z.cols(0, m - 1);
    if (beta == 0)
      beta = std::max(d(0), 0.0);  // d(0) is the (k+1)-th largest eigenvalue
    arma::vec tmpd = d - beta;
    tmpd.elem(find(tmpd < 0)).zeros();
    arma::mat L = Q * diagmat(tmpd) * Q.t();
    return 0.5 * (L + L.t());
}

// R-update prox, (sqrt(K^2 + 4 mu I) - K)/2, as originally computed.
arma::mat Rprox2(const arma::mat& K, double mu, const arma::mat& I) {
    arma::mat KI, K2 = K*K + 4*mu*I;
    arma::sqrtmat_sympd(KI, K2);
    return (KI-K)/2;
}

// Same prox from one eigendecomposition of K (K and K^2 + 4 mu I share
// eigenvectors). Each eigenvalue uses the cancellation-free branch:
// 2 mu / (sqrt(d^2 + 4 mu) + d) for d >= 0, (sqrt(d^2 + 4 mu) - d)/2 for d < 0.
arma::mat RproxEig(const arma::mat& K, double mu) {
    arma::mat A = 0.5 * (K + K.t());
    arma::vec d;
    arma::mat Q;
    eig_sym(d, Q, A);
    arma::vec s = sqrt(d % d + 4*mu);
    arma::vec t(d.n_elem);
    for (arma::uword i = 0; i < d.n_elem; i++)
      t(i) = (d(i) >= 0) ? 2*mu / (s(i) + d(i)) : (s(i) - d(i)) / 2;
    arma::mat R = Q * diagmat(t) * Q.t();
    return 0.5 * (R + R.t());
}

// [[Rcpp::export(.Rupdate_old)]]
arma::mat RupdateOldR(arma::mat K, double mu) {
    return Rprox2(K, mu, arma::eye(K.n_rows, K.n_rows));
}

// [[Rcpp::export(.Rupdate_new)]]
arma::mat RupdateNewR(arma::mat K, double mu) {
    return RproxEig(K, mu);
}

// [[Rcpp::export(.eigPSDThresh2)]]
arma::mat eigPSDThresh2R(arma::mat M, int k, double beta=0) {
    return eigPSDThresh2(M, k, beta);
}

// [[Rcpp::export(.eigPSDThreshPartial)]]
arma::mat eigPSDThreshPartialR(arma::mat M, int k, double beta=0) {
    return eigPSDThreshPartial(M, k, beta);
}


List SVD2(arma::mat& M) {
    arma::mat U, V;
    arma::vec d;
    svd_econ(U, d, V, M);
    return List::create(_["d"] = d, _["U"] = U, _["V"] = V) ;
}

// Rcpp::List SVD3(arma::mat& A, int k=5) {
//
//   Rcpp::Environment base("package:irlba");
//   Rcpp::Function f = base["irlba"];
//   return f(A, k);
// }

// // [[Rcpp::export]]
// DL_FUNC SVD4(SEXP xp, arma::mat& A, int k=5) {
//   DL_FUNC IRLB = reinterpret_cast<DL_FUNC>( R_ExternalPtrAddr(xp) ) ;
//
// //   int maxit=1000;
// //   double tol=1e-5;
// //   double eps2=1e-12;
// //   int SP;
// //   int RESTART = 0;
// //   // SEXP RV = Rcpp::Nullable<>();
// //   // SEXP RW = Rcpp::Nullable<>();
// //   // SEXP RS = Rcpp::Nullable<>();
// //   // SEXP SCALE = Rcpp::Nullable<>();
// //   // SEXP SHIFT = Rcpp::Nullable<>();
// //   // SEXP CENTER = Rcpp::Nullable<>();
// //   Rcpp::Nullable<arma::mat> RV = R_NilValue;
// //   Rcpp::Nullable<arma::vec> RW = R_NilValue;
// //   Rcpp::Nullable<arma::vec> RS = R_NilValue;
// //   Rcpp::Nullable<float> SCALE = R_NilValue;
// //   Rcpp::Nullable<float> SHIFT = R_NilValue;
// //   Rcpp::Nullable<float> CENTER = R_NilValue;
// //
// //   double svtol = tol;
// // //  RESTART <- 0L
// // //RV <- RW <- RS <- NULL
// //
// //   int n = A.n_rows;
// //   int work = k + 7;
// //
// //   arma::vec v = randn(n);
//
//   return IRLB;
//   //A, k, v, work, maxit, tol, eps2, SP, RESTART, RV, RW, RS, SCALE, SHIFT, CENTER, svtol);
//
// }


// // [[Rcpp::export]]
// arma::mat softSVT4(arma::mat& M, int k, double beta=0) {
//   List sout = SVD3(M, k);
//   arma::mat U = as<arma::mat>(sout["u"]);
//   arma::mat V = as<arma::mat>(sout["v"]);
//   arma::vec d = as<arma::vec>(sout["d"]);
//   if (beta == 0)
//     beta = d(k-1);
//   arma::vec tmpd = d - beta;
//   return U*diagmat(tmpd)*V.t();
// }

// ---- Damped Anderson acceleration of the ADMM fixed-point map -------------------
// Clean-room implementation after Henderson & Varadhan (2019), "Damped Anderson
// acceleration with restarts and monotonicity control for accelerating EM and EM-like
// algorithms", JCGS 28:834-846, written from the paper's description.
//
// The ADMM updates keep LY = SY - RY and Lambda1 = -Lambda2 = Lambda3, and all blocks
// are symmetric, so the state is the upper triangles of (RY, SY, Lambda1).
// Each step solves the damped least squares problem
//   gamma = argmin ||f - dF gamma||^2 + lam ||gamma||^2,
// with lam set so ||gamma(lam)||^2 = delta ||gamma(0)||^2, delta = 1/(1 + alpha^(kappa - s)):
// heavy damping at first, relaxing as accepted steps s accumulate. The history is a
// sliding window of m steps. A step whose fixed-point residual grows is rejected: the
// plain ADMM step is taken instead, the history is cleared and s steps back.

static arma::vec packUpper(const arma::mat& A, const arma::mat& B, const arma::mat& C) {
    arma::uvec idx = arma::trimatu_ind(arma::size(A));
    return arma::join_cols(A.elem(idx), arma::join_cols(B.elem(idx), C.elem(idx)));
}

static void unpackUpper(const arma::vec& x, arma::mat& A, arma::mat& B, arma::mat& C) {
    arma::uvec idx = arma::trimatu_ind(arma::size(A));
    arma::uword q = idx.n_elem;
    A.zeros(); B.zeros(); C.zeros();
    A.elem(idx) = x.subvec(0, q-1);
    B.elem(idx) = x.subvec(q, 2*q-1);
    C.elem(idx) = x.subvec(2*q, 3*q-1);
    A = symmatu(A); B = symmatu(B); C = symmatu(C);
}

// Tikhonov parameter giving ||gamma(lam)||^2 = delta ||gamma(0)||^2, by bisection on log(lam).
static double dampingLambda(const arma::vec& d, const arma::vec& uf, double delta) {
    double g0 = arma::accu(arma::square(uf / d));
    if (delta >= 1 || g0 == 0) return 0;
    double target = delta * g0;
    auto h = [&](double loglam) {
      double lam = std::exp(loglam);
      return arma::accu(arma::square(d % uf / (arma::square(d) + lam))) - target;
    };
    double scale = std::log(d.max() * d.max());
    double lo = scale - 40, hi = scale + 40;
    if (h(lo) < 0) return 0;
    for (int i = 0; i < 100 && hi - lo > 1e-10; i++) {
      double mid = 0.5 * (lo + hi);
      (h(mid) > 0 ? lo : hi) = mid;
    }
    return std::exp(0.5 * (lo + hi));
}

// Runs the accelerated iteration; step(it, thresholdL) performs one ADMM iteration on
// the state and returns true once ADMM's own stopping rule is met. Returns the index of
// the last step, so the caller reports iter+1 ADMM steps (rejected steps included).
template <typename Step>
int andersonDriver(Step& step, arma::mat& RY, arma::mat& SY, arma::mat& LY,
                   arma::mat& Lambda1, arma::mat& Lambda2, arma::mat& Lambda3,
                   int m, int maxiter, bool threshold_first, double& s) {
    const double kappa = 25, alpha = 1.2;
    int it = 0;
    auto set_state = [&](const arma::vec& x) {
      unpackUpper(x, RY, SY, Lambda1);
      LY = SY - RY; Lambda2 = -Lambda1; Lambda3 = Lambda1;
    };
    // G(x): one ADMM step from x; false once maxiter steps are used or ADMM has converged
    bool done = false;
    auto G = [&](const arma::vec& x) {
      set_state(x);
      done = step(it, it > 0 || threshold_first) || it + 1 >= maxiter;
      it++;
      return packUpper(RY, SY, Lambda1);
    };

    arma::vec x = packUpper(RY, SY, Lambda1);
    arma::vec gx = G(x);
    if (done) return it - 1;
    arma::vec f = gx - x;
    double rAcc = arma::norm(f);
    arma::mat dX, dF;
    // s: accepted-step counter driving the damping schedule; starts from the caller's value

    while (true) {
      arma::vec xn = gx;
      if (dF.n_cols > 0) {
        arma::mat U, V;
        arma::vec d;
        arma::svd_econ(U, d, V, dF);
        arma::uvec keep = arma::find(d > 1e-14 * d.max());
        U = U.cols(keep); V = V.cols(keep); d = d.elem(keep);
        arma::vec uf = U.t() * f;
        double lam = dampingLambda(d, uf, 1 / (1 + std::pow(alpha, kappa - s)));
        arma::vec gamma = V * (d % uf / (arma::square(d) + lam));
        xn = gx - (dX + dF) * gamma;
      }
      arma::vec gxn = G(xn);
      if (done) return it - 1;
      arma::vec fn = gxn - xn;
      double rn = arma::norm(fn);

      if (dF.n_cols > 0 && rn > rAcc) {
        // reject: plain step from x, clear the history, relax the damping less
        xn = gx;
        gxn = G(xn);
        if (done) return it - 1;
        fn = gxn - xn;
        rn = arma::norm(fn);
        dX.reset(); dF.reset();
        s = std::max(s - 1, 0.0);
      } else {
        dX = arma::join_rows(dX, xn - x);
        dF = arma::join_rows(dF, fn - f);
        if (static_cast<int>(dX.n_cols) > m) { dX.shed_col(0); dF.shed_col(0); }
        s += 1;
      }
      x = xn; gx = gxn; f = fn; rAcc = rn;
    }
}

//' @noRd
// [[Rcpp::export]]
List ADMM(const arma::mat& SigmaO, const double& lambda, arma::mat& I,
          arma::mat& Lambda, arma::mat& Y, double beta=0, int r=0, double mu=0,
          const double& eta=4/5, const double& muf=1e-6, const int&
          maxiter=500, const double& newtol=1e-5, const double& tol=1e-5,
          const double& over_relax_par=1.6, bool shrinkDiag=true,
          bool partial_eig=false, bool fast_R=false, bool balance_mu=false,
          bool threshold_first=false, int anderson=0,
          double anderson_s0=0) {

    int n = SigmaO.n_rows;
    if (mu == 0) mu = n;
    if (anderson > 0 && (balance_mu || eta != 1))
      Rcpp::stop("anderson acceleration needs a fixed mu (eta = 1, balance.mu = FALSE)");
    // balance_mu: stop on primal AND dual residuals (Boyd et al. 3.3) and adapt mu
    // by residual balancing (3.4.1), frozen after balance_iters so the changes stay summable.
    const int balance_iters = 50;
    const double balance_ratio = 10, balance_factor = 2;
    arma::vec r_norm(maxiter), eps_pri(maxiter), s_norm(maxiter), eps_dual(maxiter);
    arma::mat Lambda1 = Lambda.cols(0  , n-1  );
    arma::mat Lambda2 = Lambda.cols(n  , n*2-1);
    arma::mat Lambda3 = Lambda.cols(n*2, n*3-1);
    arma::mat R = Y.cols(0,n-1), L = Y.cols(n*2,n*3-1);
    arma::mat S = Y.cols(n,n*2-1);
    arma::mat RY(R), LY(L), SY(S);
    arma::mat RO(size(R)), SO(size(S)), LO(size(L));
    arma::mat RA(size(R)), SA(size(S)), LA(size(L)), TA(size(L));
    arma::mat RYp, SYp, LYp;
    // squared Frobenius norm summed over the three blocks
    auto sq = [](const arma::mat& M) { double v = norm(M, "fro"); return v*v; };
//    arma::mat Y_old(Y);
    arma::mat K, K2, KI;

    // One ADMM iteration on the state (RY, SY, LY, Lambda1..3); records the
    // primal residual in slot `it` of the history.
    auto step = [&](int it, bool thresholdL) {
      // update X = (R,S,L)
        RA = RY + mu*Lambda1;
        SA = SY + mu*Lambda2;
        LA = LY + mu*Lambda3;
      // update R
      K  = mu*SigmaO - RA;
      if (fast_R) {
        R = RproxEig(K, mu);
      } else {
        K2 = K*K + 4*mu*I;
//        KI = sqrtmNewt2(K2, I, newtol);
        arma::sqrtmat_sympd(KI, K2);
        R  = (KI-K)/2;
      }

      // update S
      S = SOFTTHRESHD(SA, lambda*mu, shrinkDiag);
      if (!thresholdL) {
        L = LA;
      } else {
        L = partial_eig ? eigPSDThreshPartial(LA, r, mu*beta)
                        : eigPSDThresh2(LA, r, mu*beta);
//        L = softSVT4(LA, r+1, mu*beta);
      }

      RO = over_relax_par*R + (1-over_relax_par)*RY;
      SO = over_relax_par*S + (1-over_relax_par)*SY;
      LO = over_relax_par*L + (1-over_relax_par)*LY;

      if (balance_mu) {
        RYp = RY; SYp = SY; LYp = LY;
      }
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

      r_norm(it)  = std::sqrt(sq(R - RY) + sq(S - SY) + sq(L - LY));
      eps_pri(it) = sqrt(3*n*n)*tol + tol*std::sqrt(std::max(sq(R) + sq(S) + sq(L),
                                                             sq(RY) + sq(SY) + sq(LY)));
      return r_norm(it) < eps_pri(it);
    };

    int iter;
    double anderson_s = anderson_s0;
    if (anderson > 0) {
      iter = andersonDriver(step, RY, SY, LY, Lambda1, Lambda2, Lambda3, anderson,
                            maxiter, threshold_first, anderson_s);
    } else
    for (iter = 0; iter < maxiter; iter++) {
      bool converged = step(iter, iter > 0 || threshold_first);

      if (!balance_mu) {
        if (converged) {
           break;
        }
        mu = std::max(mu*eta, muf);
        continue;
      }

      s_norm(iter)   = std::sqrt(sq(RY - RYp) + sq(SY - SYp) + sq(LY - LYp)) / mu;
      eps_dual(iter) = sqrt(3*n*n)*tol + tol*std::sqrt(sq(Lambda1) + sq(Lambda2) + sq(Lambda3));
      if (r_norm(iter) < eps_pri(iter) && s_norm(iter) < eps_dual(iter)) {
         break;
      }
      // mu is 1/rho: a large primal residual calls for a larger rho, i.e. a smaller mu
      if (iter < balance_iters) {
        if (r_norm(iter) > balance_ratio * s_norm(iter))
          mu = std::max(mu / balance_factor, muf);
        else if (s_norm(iter) > balance_ratio * r_norm(iter))
          mu = mu * balance_factor;
      }

    }

    Lambda.cols(0  , n-1  ) = Lambda1;
    Lambda.cols(n  , n*2-1) = Lambda2;
    Lambda.cols(n*2, n*3-1) = Lambda3;
    Y.cols(0  , n-1  ) = RY;
    Y.cols(n  , n*2-1) = SY;
    Y.cols(n*2, n*3-1) = LY;

    return List::create(_["R"] = R,  _["S"] = arma::sp_mat(S), _["L"] = L,
                        _["Y"] = Y,
                        _["Lambda"] = Lambda,
//                        _["RA"] = RA, _["SA"] = SA, _["LA"] = LA,
                        _["iter"] = iter+1, _["beta"] = beta, _["mu"]   = mu,
                        _["anderson_s"] = anderson_s,
                        _["history"] = List::create(
                          _["r_norm"]  =  r_norm,
                          _["eps_pri"] = eps_pri,
                          _["s_norm"]  = s_norm,
                          _["eps_dual"] = eps_dual)
                        );
}
