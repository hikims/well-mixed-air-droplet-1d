#pragma once

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "params1.h"

inline double clamp_value(double x, double lo, double hi)
{
    return x < lo ? lo : (x > hi ? hi : x);
}

constexpr double PI_CONST = 3.141592653589793238462643383279502884;

#ifdef USE_REAL_TBB
#include <oneapi/tbb.h>
namespace ptbb = oneapi::tbb;
#else
#include <cstddef>
namespace ptbb {
template<class T>
class blocked_range {
public:
    blocked_range(T b, T e) : b_(b), e_(e) {}
    T begin() const { return b_; }
    T end() const { return e_; }
private:
    T b_, e_;
};

template<class R, class B>
void parallel_for(const R& r, const B& b) { b(r); }

template<class R, class V, class F, class J>
V parallel_reduce(const R& r, V init, const F& f, const J&) { return f(r, init); }

class global_control {
public:
    enum parameter { max_allowed_parallelism };
    global_control(parameter, std::size_t) {}
};
}
#endif

#define UBEGIN (0)
#define UEND (UBEGIN + Nx)

#define SBEGIN UEND
#define SEND (SBEGIN + Nx)

#define CBEGIN SEND
#define CEND (CBEGIN + Nx)

#define WBEGIN CEND
#define WEND (WBEGIN + 1)
#define IWBAR (WBEGIN)

#define XCBEGIN WEND
#define XCEND (XCBEGIN + N_BS)

#define NV XCEND

#define IU(ix) (UBEGIN + (ix))
#define IS(ix) (SBEGIN + (ix))
#define IC(ix) (CBEGIN + (ix))
#define IXC(i) (XCBEGIN + (i))

#define U(ix) Y[IU(ix)]
#define S(ix) Y[IS(ix)]
#define C(ix) Y[IC(ix)]
#define WBAR Y[IWBAR]
#define XC(i) Y[IXC(i)]

inline double wrap_periodic(double x)
{
    x = std::fmod(x, Lx);
    if (x < 0.0) x += Lx;
    return x;
}

inline double pdiff(double p1, double p2)
{
    double d = p1 - p2;

    while (d >= 0.5 * Lx) d -= Lx;
    while (d < -0.5 * Lx) d += Lx;

    return d;
}

inline double periodic_absdist(double p1, double p2)
{
    return std::abs(pdiff(p1, p2));
}

template<class T>
T Ga(T u, T s)
{
    using std::pow;

    return 0.5 *
           (
               (pow(aa, 2) / (pow(aa, 2) + pow(u, 2)))
               +
               (bb / (bb + s))
           )
           *
           (G0 - G1)
           +
           G1;
}

template<class T>
T f(T x)
{
    return alpha * x / (beta + x);
}

template<class T>
T f1(T x)
{
    return 1.0 - (x / C0);
}

template<class T>
T f2(T x)
{
    return 1.0 - (x / W0);
}

template<class T>
T f3(T x)
{
    return 1.0 - (x / U0);
}

template<class T, class T2>
auto tanh_switch(T x, T2 at)
{
    using std::tanh;
    return 0.5 * (tanh((x - at) / taneps) + 1.0);
}

template<class T, class T2, class T3>
auto sigmoid(T x, T2 center, T3 radius)
{
    return tanh_switch(x, center - radius)
           *
           tanh_switch(-x, -(center + radius));
}

inline double interp_periodic_field(
    const double* Y,
    int field_begin,
    double x
){
    x = wrap_periodic(x);

    double q = x * inv_dx;
    int i0f = static_cast<int>(std::floor(q));
    double ratio = q - static_cast<double>(i0f);

    int i0 = ((i0f % Nx) + Nx) % Nx;
    int i1 = (i0 + 1) % Nx;

    return Y[field_begin + i0] * (1.0 - ratio)
           +
           Y[field_begin + i1] * ratio;
}

inline double covered_filter(int ix, const double* Y)
{
    double x = dx * ix;
    double filter = 0.0;

    for (int i = 0; i < N_BS; ++i) {
        double center = XC(i);

        while (center - x >= 0.5 * Lx) center -= Lx;
        while (center - x < -0.5 * Lx) center += Lx;

        filter = std::max(filter, sigmoid(x, center, rad));
    }

    return clamp_value(filter, 0.0, 1.0);
}

inline double G_uncovered(int ix, const double* Y)
{
    double filter = covered_filter(ix, Y);
    return U(ix) * (1.0 - filter);
}

inline double F_supply(int ix, const double* Y)
{
    double x = dx * ix;
    double ret = 0.0;

    for (int i = 0; i < N_BS; ++i) {
        double center = XC(i);

        while (center - x >= 0.5 * Lx) center -= Lx;
        while (center - x < -0.5 * Lx) center += Lx;

        ret += K0 * C1 * f3(U(ix)) * sigmoid(x, center, rad);
    }

    return ret;
}

inline double droplet_velocity_one(const double* Y, int i)
{
    double ub = interp_periodic_field(Y, UBEGIN, XC(i) - rad);
    double sb = interp_periodic_field(Y, SBEGIN, XC(i) - rad);

    double ue = interp_periodic_field(Y, UBEGIN, XC(i) + rad);
    double se = interp_periodic_field(Y, SBEGIN, XC(i) + rad);

    double fst = fir_co * (Ga(ue, se) - Ga(ub, sb));

    double frep = 0.0;

    if (N_BS > 1) {
        double r1 = periodic_absdist(XC((i - 1 + N_BS) % N_BS), XC(i));
        double e1 = (r1 <= 2.0 * rad)
                    ? (Xi_0 * (2.0 * rad - r1)) / rad
                    : 0.0;

        double r2 = periodic_absdist(XC(i), XC((i + 1) % N_BS));
        double e2 = (r2 <= 2.0 * rad)
                    ? (Xi_0 * (2.0 * rad - r2)) / rad
                    : 0.0;

        frep = sec_co * (-e1 + e2);
    }

    return fst + frep;
}

class TridiagonalLU {
public:
    TridiagonalLU() = default;

    TridiagonalLU(
        int n,
        double lower,
        const std::vector<double>& diag,
        double upper
    ){
        factorize(n, lower, diag, upper);
    }

    void factorize(
        int n,
        double lower,
        const std::vector<double>& diag,
        double upper
    ){
        n_ = n;
        lower_ = lower;
        upper_ = upper;

        cprime_.assign(n_, 0.0);
        invden_.assign(n_, 0.0);

        double den = diag[0];

        if (std::abs(den) < 1e-300) {
            throw std::runtime_error("zero pivot");
        }

        invden_[0] = 1.0 / den;
        cprime_[0] = upper_ * invden_[0];

        for (int i = 1; i < n_; ++i) {
            den = diag[i] - lower_ * cprime_[i - 1];

            if (std::abs(den) < 1e-300) {
                throw std::runtime_error("zero pivot");
            }

            invden_[i] = 1.0 / den;
            cprime_[i] = (i < n_ - 1) ? upper_ * invden_[i] : 0.0;
        }
    }

    void solve(
        const std::vector<double>& rhs,
        std::vector<double>& x
    ) const {
        x.assign(n_, 0.0);

        std::vector<double> dp(n_, 0.0);

        dp[0] = rhs[0] * invden_[0];

        for (int i = 1; i < n_; ++i) {
            dp[i] = (rhs[i] - lower_ * dp[i - 1]) * invden_[i];
        }

        x[n_ - 1] = dp[n_ - 1];

        for (int i = n_ - 2; i >= 0; --i) {
            x[i] = dp[i] - cprime_[i] * x[i + 1];
        }
    }

private:
    int n_ = 0;
    double lower_ = 0.0;
    double upper_ = 0.0;

    std::vector<double> cprime_;
    std::vector<double> invden_;
};

class PeriodicDiffusionLU {
public:
    PeriodicDiffusionLU() = default;

    PeriodicDiffusionLU(
        int n,
        double diffusion,
        double dt,
        double h
    ){
        setup(n, diffusion, dt, h);
    }

    void setup(
        int n,
        double diffusion,
        double dt,
        double h
    ){
        n_ = n;
        lambda_ = dt * diffusion / (h * h);

        if (lambda_ == 0.0) {
            identity_ = true;
            return;
        }

        identity_ = false;

        double diag = 1.0 + 2.0 * lambda_;
        double off = -lambda_;

        alpha_ = off;
        beta_ = off;
        gamma_ = -diag;

        std::vector<double> bdiag(n_, diag);

        bdiag[0] = diag - gamma_;
        bdiag[n_ - 1] = diag - alpha_ * beta_ / gamma_;

        trilu_.factorize(n_, off, bdiag, off);

        std::vector<double> uvec(n_, 0.0);
        uvec[0] = gamma_;
        uvec[n_ - 1] = alpha_;

        trilu_.solve(uvec, z_);
    }

    void solve(
        const std::vector<double>& rhs,
        std::vector<double>& x
    ) const {
        if (identity_) {
            x = rhs;
            return;
        }

        trilu_.solve(rhs, x);

        double num = x[0] + beta_ * x[n_ - 1] / gamma_;
        double den = 1.0 + z_[0] + beta_ * z_[n_ - 1] / gamma_;
        double fact = num / den;

        for (int i = 0; i < n_; ++i) {
            x[i] -= fact * z_[i];
        }
    }

private:
    int n_ = 0;
    double lambda_ = 0.0;
    double alpha_ = 0.0;
    double beta_ = 0.0;
    double gamma_ = 0.0;

    bool identity_ = false;

    TridiagonalLU trilu_;
    std::vector<double> z_;
};

struct PhiPsi {
    double phi = 0.0;
    double psi = 0.0;
};

struct SolverPack {
    PeriodicDiffusionLU lu_u;
    PeriodicDiffusionLU lu_s;
    PeriodicDiffusionLU lu_c;

    std::vector<double> rhs_u;
    std::vector<double> rhs_s;
    std::vector<double> rhs_c;

    std::vector<double> sol_u;
    std::vector<double> sol_s;
    std::vector<double> sol_c;

    explicit SolverPack(double dt)
        :
        lu_u(Nx, du, dt, dx),
        lu_s(Nx, ds, dt, dx),
        lu_c(Nx, dc, dt, dx),
        rhs_u(Nx, 0.0),
        rhs_s(Nx, 0.0),
        rhs_c(Nx, 0.0),
        sol_u(Nx, 0.0),
        sol_s(Nx, 0.0),
        sol_c(Nx, 0.0)
    {}
};

inline void initialize_state(
    std::vector<double>& y,
    double s0,
    unsigned seed,
    bool centered
){
    y.assign(NV, 0.0);

    std::mt19937 mt(seed);
    std::uniform_real_distribution<double> ur(0.0, 1.0);

    double* Y = y.data();

    for (int ix = 0; ix < Nx; ++ix) {
        U(ix) = U0 * 0.001 * ur(mt);
        S(ix) = s0;
        C(ix) = 0.0;
    }

    WBAR = 0.0;

    for (int i = 0; i < N_BS; ++i) {
        XC(i) = (centered && N_BS == 1)
                ? 0.5 * Lx
                : wrap_periodic(
                    Lx * ((i / static_cast<double>(N_BS)) + 0.5)
                    +
                    1.0 * (ur(mt) - 0.5)
                  );
    }
}

inline void load_lastdata_or_initialize(
    std::vector<double>& y,
    const std::string& lastdata,
    double s0,
    unsigned seed,
    bool centered
){
    if (!lastdata.empty()
        &&
        lastdata != "none"
        &&
        lastdata != "NULL"
        &&
        lastdata != "-"
    ){
        y.assign(NV, 0.0);

        std::ifstream ifs(lastdata, std::ios::binary);

        if (!ifs) {
            throw std::runtime_error("cannot open lastdata: " + lastdata);
        }

        ifs.read(reinterpret_cast<char*>(y.data()), sizeof(double) * y.size());

        if (!ifs) {
            throw std::runtime_error("failed read lastdata: " + lastdata);
        }
    } else {
        initialize_state(y, s0, seed, centered);
    }
}

inline void save_lastdata(
    const std::string& file,
    const std::vector<double>& y
){
    std::ofstream ofs(file, std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(y.data()), sizeof(double) * y.size());
}

inline PhiPsi compute_phi_psi(const double* Y)
{
    PhiPsi ret = ptbb::parallel_reduce(
        ptbb::blocked_range<int>(0, Nx),
        PhiPsi{},
        [&](const ptbb::blocked_range<int>& r, PhiPsi local) {
            for (int ix = r.begin(); ix != r.end(); ++ix) {
                local.phi += G_uncovered(ix, Y);
                local.psi += f3(U(ix));
            }
            return local;
        },
        [](PhiPsi a, PhiPsi b) {
            return PhiPsi{a.phi + b.phi, a.psi + b.psi};
        }
    );

    ret.phi /= static_cast<double>(Nx);
    ret.psi /= static_cast<double>(Nx);

    return ret;
}

inline double update_wbar_exact(
    double W,
    const PhiPsi& pp,
    double kw,
    double dt
){
    double A = kw + Ku * pp.phi / W0 + Ka * pp.psi;
    double B = Ku * pp.phi;

    if (std::abs(A) < 1e-14) {
        W += dt * B;
    } else {
        double Wss = B / A;
        W = Wss + (W - Wss) * std::exp(-A * dt);
    }

    if (!std::isfinite(W)) {
        W = 0.0;
    }

    return clamp_value(W, 0.0, 1.5 * W0);
}

inline void wrap_droplets(std::vector<double>& y)
{
    double* Y = y.data();

    for (int i = 0; i < N_BS; ++i) {
        XC(i) = wrap_periodic(XC(i));
    }
}

inline void imex_lu_step(
    std::vector<double>& y,
    SolverPack& solvers,
    double dt,
    double s0,
    double kw,
    bool move_droplets
){
    double* Y = y.data();

    ptbb::parallel_for(
        ptbb::blocked_range<int>(0, Nx),
        [&](const ptbb::blocked_range<int>& r) {
            for (int ix = r.begin(); ix != r.end(); ++ix) {
                double uu = std::max(0.0, U(ix));
                double ss = std::max(0.0, S(ix));
                double cc = clamp_value(C(ix), 0.0, 0.999 * C0);
                double ww = clamp_value(WBAR, 0.0, 1.5 * W0);

                double reac = f1(cc) * f(ss) * uu;

                double Ru =
                    -Kus * reac
                    -
                    Ku * f2(ww) * G_uncovered(ix, Y)
                    +
                    Ka * f3(uu) * ww
                    +
                    F_supply(ix, Y);

                double Rs =
                    -Kus * reac
                    +
                    D1COEF * (s0 - ss);

                double Rc =
                    Kus * reac
                    -
                    D2COEF * cc;

                solvers.rhs_u[ix] = uu + dt * Ru;
                solvers.rhs_s[ix] = ss + dt * Rs;
                solvers.rhs_c[ix] = cc + dt * Rc;
            }
        }
    );

    solvers.lu_u.solve(solvers.rhs_u, solvers.sol_u);
    solvers.lu_s.solve(solvers.rhs_s, solvers.sol_s);
    solvers.lu_c.solve(solvers.rhs_c, solvers.sol_c);

    ptbb::parallel_for(
        ptbb::blocked_range<int>(0, Nx),
        [&](const ptbb::blocked_range<int>& r) {
            for (int ix = r.begin(); ix != r.end(); ++ix) {
                U(ix) = std::max(0.0, solvers.sol_u[ix]);
                S(ix) = std::max(0.0, solvers.sol_s[ix]);
                C(ix) = clamp_value(solvers.sol_c[ix], 0.0, 0.999 * C0);
            }
        }
    );

    PhiPsi pp = compute_phi_psi(Y);
    WBAR = update_wbar_exact(WBAR, pp, kw, dt);

    if (move_droplets) {
        std::vector<double> vel(N_BS);

        for (int i = 0; i < N_BS; ++i) {
            vel[i] = droplet_velocity_one(Y, i);
        }

        for (int i = 0; i < N_BS; ++i) {
            XC(i) += dt * vel[i];
        }

        wrap_droplets(y);
    }
}

inline void ensure_dir(const std::string& path)
{
    std::string cmd = std::string("mkdir -p ") + path;
    std::system(cmd.c_str());
}

inline void save_snapshot(
    const std::string& path,
    int it,
    const std::vector<double>& y,
    double s0
){
    ensure_dir(path);

    const double* Y = y.data();

    std::string up = path + "/u" + std::to_string(it);
    std::string sp = path + "/s" + std::to_string(it);
    std::string cp = path + "/c" + std::to_string(it);
    std::string wp = path + "/w" + std::to_string(it);
    std::string xp = path + "/x" + std::to_string(it);

    /*
        Save numbered snapshot files.
    */
    {
        std::ofstream ou(up);
        std::ofstream os(sp);
        std::ofstream oc(cp);
        std::ofstream ow(wp);
        std::ofstream ox(xp);

        ou << std::setprecision(15);
        os << std::setprecision(15);
        oc << std::setprecision(15);
        ow << std::setprecision(15);
        ox << std::setprecision(15);

        for (int ix = 0; ix < Nx; ++ix) {
            const double x = ix * dx;

            ou << x << " " << U(ix) / U0 << "\n";
            os << x << " " << S(ix) / s0 << "\n";
            oc << x << " " << C(ix) / C0 << "\n";
            ow << x << " " << WBAR / W0 << "\n";
        }

        for (int i = 0; i < N_BS; ++i) {
            ox << XC(i) << " " << 0.0 << "\n";
        }
    }

    /*
        Save latest snapshot files directly.
        These files are used by plotting scripts:
            u, s, c, w, x
    */
    {
        std::ofstream ou(path + "/u");
        std::ofstream os(path + "/s");
        std::ofstream oc(path + "/c");
        std::ofstream ow(path + "/w");
        std::ofstream ox(path + "/x");

        ou << std::setprecision(15);
        os << std::setprecision(15);
        oc << std::setprecision(15);
        ow << std::setprecision(15);
        ox << std::setprecision(15);

        for (int ix = 0; ix < Nx; ++ix) {
            const double x = ix * dx;

            ou << x << " " << U(ix) / U0 << "\n";
            os << x << " " << S(ix) / s0 << "\n";
            oc << x << " " << C(ix) / C0 << "\n";
            ow << x << " " << WBAR / W0 << "\n";
        }

        for (int i = 0; i < N_BS; ++i) {
            ox << XC(i) << " " << 0.0 << "\n";
        }
    }
}

inline void save_profile_csv(
    const std::string& file,
    const std::vector<double>& y,
    double
){
    const double* Y = y.data();

    std::ofstream ofs(file);
    ofs << std::setprecision(15);

    ofs << "x,u,s,c,Wbar,xc\n";

    for (int ix = 0; ix < Nx; ++ix) {
        ofs
            << ix * dx << ","
            << U(ix) << ","
            << S(ix) << ","
            << C(ix) << ","
            << WBAR << ","
            << XC(0) << "\n";
    }
}

inline void open_xc_files(
    std::vector<std::ofstream>& ofs,
    std::vector<double>& prev,
    const std::vector<double>& y,
    const std::string& path,
    const std::string& suffix
){
    ensure_dir(path);

    const double* Y = y.data();

    ofs.resize(N_BS);
    prev.resize(N_BS);

    for (int i = 0; i < N_BS; ++i) {
        ofs[i].open(path + "/xc" + std::to_string(i) + suffix);
        ofs[i] << std::setprecision(15);
        prev[i] = XC(i);
    }
}

inline void append_xc(
    std::vector<std::ofstream>& ofs,
    std::vector<double>& prev,
    const std::vector<double>& y,
    double t
){
    const double* Y = y.data();

    for (int i = 0; i < N_BS; ++i) {
        if (std::abs(prev[i] - XC(i)) > Lx * 0.5) {
            ofs[i] << "\n";
        }

        ofs[i] << t << " " << XC(i) << "\n";
        prev[i] = XC(i);
    }
}

inline void run_cover_segment(
    const std::string& lastdata,
    const std::string& output_lastdata,
    const std::string& path,
    const std::string& suffix,
    double s0,
    double kw,
    double t_end,
    double dt,
    int save_every,
    unsigned seed
){
    ensure_dir(path);

    std::vector<double> y;
    load_lastdata_or_initialize(y, lastdata, s0, seed, false);

    SolverPack solvers(dt);

    std::vector<std::ofstream> ofsv;
    std::vector<double> prev;

    open_xc_files(ofsv, prev, y, path, suffix);

    std::ofstream time_ofs(path + "/time");

    int nsteps = static_cast<int>(std::ceil(t_end / dt));
    int it = 0;

    for (int step = 0; step <= nsteps; ++step) {
        double t = step * dt;

        if (step % save_every == 0) {
            save_snapshot(path, it, y, s0);
            append_xc(ofsv, prev, y, t);
            time_ofs << std::setprecision(15) << t << "\n";

            std::cerr
                << "t:" << t
                << " dt:" << dt
                << " Wbar:" << y[IWBAR]
                << std::endl;

            ++it;
        }

        if (step == nsteps) {
            break;
        }

        imex_lu_step(y, solvers, dt, s0, kw, true);
    }

    save_lastdata(output_lastdata, y);
}

inline void run_amplitude(
    double s0,
    double kw,
    double t_end,
    double dt,
    int save_every,
    unsigned seed,
    const std::string& path
){
    ensure_dir(path);

    std::vector<double> y;
    initialize_state(y, s0, seed, false);

    SolverPack solvers(dt);

    std::ofstream amp(path + "/amplitude", std::ios::app);
    std::ofstream xco(path + "/xc" + std::to_string(s0));

    amp << std::setprecision(15);
    xco << std::setprecision(15);

    int nsteps = static_cast<int>(std::ceil(t_end / dt));
    int nstart = static_cast<int>(0.5 * nsteps);

    double amin = DBL_MAX;
    double amax = -DBL_MAX;

    double prev = y[IXC(0)];

    for (int step = 0; step <= nsteps; ++step) {
        double t = step * dt;

        const double* Y = y.data();

        if (step >= nstart) {
            amin = std::min(amin, XC(0));
            amax = std::max(amax, XC(0));
        }

        if (step % save_every == 0) {
            if (std::abs(prev - XC(0)) > Lx * 0.5) {
                xco << "\n";
            }

            xco << t << " " << XC(0) << "\n";
            prev = XC(0);
        }

        if (step == nsteps) {
            break;
        }

        imex_lu_step(y, solvers, dt, s0, kw, true);
    }

    amp << s0 << " " << (amax - amin) << "\n";
}

inline void run_rest_perturb(
    double s0,
    double kw,
    double t_relax,
    double t_free,
    double dt,
    int save_every,
    double eps_shift,
    unsigned seed,
    const std::string& path
){
    ensure_dir(path);

    std::vector<double> y;
    initialize_state(y, s0, seed, true);

    SolverPack solvers(dt);

    int nrelax = static_cast<int>(std::ceil(t_relax / dt));

    for (int step = 0; step < nrelax; ++step) {
        imex_lu_step(y, solvers, dt, s0, kw, false);
    }

    save_lastdata(path + "/lastdata_rest", y);
    save_snapshot(path, 0, y, s0);
    save_profile_csv(path + "/rest_state_profile.csv", y, s0);

    double* Y = y.data();

    double x_rest = XC(0);
    XC(0) = wrap_periodic(XC(0) + eps_shift);

    std::ofstream traj(path + "/rest_perturb_trajectory.csv");
    traj << std::setprecision(15);
    traj << "t,delta_x,theta,Wbar,velocity,xc\n";

    int nfree = static_cast<int>(std::ceil(t_free / dt));
    int it = 1;

    for (int step = 0; step <= nfree; ++step) {
        double t = step * dt;

        Y = y.data();

        double delta = pdiff(XC(0), x_rest);
        double theta = 2.0 * PI_CONST * XC(0) / Lx;
        double vel = droplet_velocity_one(Y, 0);

        if (step % save_every == 0) {
            traj
                << t << ","
                << delta << ","
                << theta << ","
                << WBAR << ","
                << vel << ","
                << XC(0) << "\n";

            save_snapshot(path, it, y, s0);
            ++it;
        }

        if (step == nfree) {
            break;
        }

        imex_lu_step(y, solvers, dt, s0, kw, true);
    }

    save_lastdata(path + "/lastdata_free", y);
    save_profile_csv(path + "/final_state_profile.csv", y, s0);
}