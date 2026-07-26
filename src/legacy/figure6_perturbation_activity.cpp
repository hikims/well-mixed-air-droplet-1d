#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Fig6Case {
    double s0 = 0.0;
    std::string folder;
};

inline std::string s0_folder_name(double s0)
{
    std::ostringstream ss;
    ss << "s0_" << std::fixed << std::setprecision(3) << s0;
    return ss.str();
}

inline void initialize_clean_centered_rest_state(
    std::vector<double>& y,
    double s0
){
    y.assign(NV, 0.0);

    double* Y = y.data();

    for (int ix = 0; ix < Nx; ++ix) {
        U(ix) = 0.0;
        S(ix) = s0;
        C(ix) = 0.0;
    }

    WBAR = 0.0;

    for (int i = 0; i < N_BS; ++i) {
        if (N_BS == 1) {
            XC(i) = 0.5 * Lx;
        } else {
            XC(i) = wrap_periodic(Lx * ((i / static_cast<double>(N_BS)) + 0.5));
        }
    }
}

inline void run_one_perturbation_case(
    const Fig6Case& cas,
    const std::string& out_root,
    double kw,
    double t_relax,
    double t_free,
    double dt,
    int save_every,
    double eps_shift
){
    const std::string out_dir = out_root + "/" + cas.folder;
    ensure_dir(out_dir);

    std::vector<double> y;
    initialize_clean_centered_rest_state(y, cas.s0);

    SolverPack solvers(dt);

    /*
        Step 1. Fixed-droplet relaxation to symmetric rest state.
    */
    const int nrelax = static_cast<int>(std::ceil(t_relax / dt));

    for (int step = 0; step < nrelax; ++step) {
        imex_lu_step(y, solvers, dt, cas.s0, kw, false);

        if (step % (save_every * 500) == 0) {
            std::cerr
                << "relax s0=" << cas.s0
                << " t=" << step * dt
                << " Wbar=" << y[IWBAR]
                << std::endl;
        }
    }

    save_snapshot(out_dir + "/rest_profile", 0, y, cas.s0);
    save_lastdata(out_dir + "/lastdata_rest", y);

    /*
        Step 2. Position perturbation.
    */
    double* Y = y.data();
    const double x_rest = XC(0);
    XC(0) = wrap_periodic(XC(0) + eps_shift);

    /*
        Step 3. Free evolution and record delta x.
    */
    std::ofstream traj(out_dir + "/perturbation.csv");
    traj << std::setprecision(15);
    traj
        << "t,"
        << "delta_x,"
        << "abs_delta_x,"
        << "theta,"
        << "Wbar,"
        << "velocity,"
        << "xc"
        << "\n";

    const int nfree = static_cast<int>(std::ceil(t_free / dt));

    for (int step = 0; step <= nfree; ++step) {
        const double t = step * dt;

        Y = y.data();

        if (step % save_every == 0) {
            const double delta = pdiff(XC(0), x_rest);
            const double theta = 2.0 * PI_CONST * XC(0) / Lx;
            const double vel = droplet_velocity_one(Y, 0);

            traj
                << t << ","
                << delta << ","
                << std::abs(delta) << ","
                << theta << ","
                << WBAR << ","
                << vel << ","
                << XC(0)
                << "\n";
        }

        if (step == nfree) {
            break;
        }

        imex_lu_step(y, solvers, dt, cas.s0, kw, true);
    }

    save_snapshot(out_dir + "/final_profile", 0, y, cas.s0);
    save_lastdata(out_dir + "/lastdata_free", y);

    std::cerr
        << "done perturbation s0=" << cas.s0
        << " folder=" << out_dir
        << std::endl;
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/figure6_lu out_dir smin smax ds0 t_relax t_free dt save_every threads eps

        Example:
            ./bin/figure6_lu figure6_data 0.085 0.120 0.001 80 80 0.0005 20 8 0.001
    */

    std::string out_root = "figure6_data";

    double smin = 0.085;
    double smax = 0.120;
    double ds0 = 0.001;

    double t_relax = 80.0;
    double t_free = 80.0;

    double dt = default_dt;
    int save_every = 20;
    int threads = 8;

    double eps_shift = 0.001;

    if (argc >= 2) out_root = argv[1];
    if (argc >= 3) smin = std::atof(argv[2]);
    if (argc >= 4) smax = std::atof(argv[3]);
    if (argc >= 5) ds0 = std::atof(argv[4]);
    if (argc >= 6) t_relax = std::atof(argv[5]);
    if (argc >= 7) t_free = std::atof(argv[6]);
    if (argc >= 8) dt = std::atof(argv[7]);
    if (argc >= 9) save_every = std::atoi(argv[8]);
    if (argc >= 10) threads = std::atoi(argv[9]);
    if (argc >= 11) eps_shift = std::atof(argv[10]);

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    ensure_dir(out_root);

    std::vector<Fig6Case> cases;

    for (double s = smin; s <= smax + 1.0e-12; s += ds0) {
        Fig6Case cas;
        cas.s0 = s;
        cas.folder = s0_folder_name(s);
        cases.push_back(cas);
    }

    std::ofstream summary(out_root + "/perturbation_cases.csv");
    summary << std::setprecision(15);
    summary << "s0,folder,eps_shift,t_relax,t_free,dt,save_every\n";

    for (const auto& cas : cases) {
        summary
            << cas.s0 << ","
            << cas.folder << ","
            << eps_shift << ","
            << t_relax << ","
            << t_free << ","
            << dt << ","
            << save_every << "\n";
    }

    std::cerr << "Figure 6 perturbation sweep\n";
    std::cerr << "out_root=" << out_root
              << " smin=" << smin
              << " smax=" << smax
              << " ds0=" << ds0
              << " t_relax=" << t_relax
              << " t_free=" << t_free
              << " dt=" << dt
              << " save_every=" << save_every
              << " eps=" << eps_shift
              << " cases=" << cases.size()
              << std::endl;

    /*
        Parallel over s0 cases. Each case still uses TBB inside grid loops.
    */
    ptbb::parallel_for(
        ptbb::blocked_range<int>(0, static_cast<int>(cases.size())),
        [&](const ptbb::blocked_range<int>& range) {
            for (int k = range.begin(); k != range.end(); ++k) {
                run_one_perturbation_case(
                    cases[k],
                    out_root,
                    Kw_open,
                    t_relax,
                    t_free,
                    dt,
                    save_every,
                    eps_shift
                );
            }
        }
    );

    std::cerr << "Done. Data saved to " << out_root << std::endl;

    return 0;
}
