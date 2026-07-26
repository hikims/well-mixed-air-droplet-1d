#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct PerturbationCase {
    double s0 = 0.0;
    std::string folder;
};

inline std::string s0_folder_name(double s0)
{
    std::ostringstream ss;
    // Six decimals are required because this validation uses ds0 = 0.000125.
    ss << "s0_" << std::fixed << std::setprecision(6) << s0;
    return ss.str();
}

inline double grid_aligned_single_droplet_center()
{
    // For Nx=401 and dx=Lx/Nx, Lx/2 lies halfway between two grid points.
    // Place the rest-state center exactly at grid index Nx/2 instead.
    return wrap_periodic(static_cast<double>(Nx / 2) * dx);
}

inline void initialize_clean_grid_aligned_rest_state(
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
            XC(i) = grid_aligned_single_droplet_center();
        } else {
            // This validation is designed for N_BS=1. Retain a deterministic
            // periodic placement for completeness in multi-droplet builds.
            XC(i) = wrap_periodic(
                Lx * ((i / static_cast<double>(N_BS)) + 0.5)
            );
        }
    }
}

inline void run_one_perturbation_case(
    const PerturbationCase& cas,
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
    initialize_clean_grid_aligned_rest_state(y, cas.s0);
    SolverPack solvers(dt);

    // Step 1: fixed-droplet relaxation to the symmetric stationary profile.
    const int nrelax = static_cast<int>(std::ceil(t_relax / dt));
    for (int step = 0; step < nrelax; ++step) {
        imex_lu_step(y, solvers, dt, cas.s0, kw, false);

        if (step % (save_every * 2000) == 0) {
            std::cerr
                << "relax s0=" << cas.s0
                << " t=" << step * dt
                << " Wbar=" << y[IWBAR]
                << std::endl;
        }
    }

    save_snapshot(out_dir + "/rest_profile", 0, y, cas.s0);
    save_lastdata(out_dir + "/lastdata_rest", y);

    // Step 2: perturb only the droplet position, leaving the fields unchanged.
    double* Y = y.data();
    const double x_rest = XC(0);
    XC(0) = wrap_periodic(x_rest + eps_shift);

    // Step 3: release the droplet and record the free response.
    std::ofstream traj(out_dir + "/perturbation.csv");
    if (!traj) {
        throw std::runtime_error("failed to open perturbation.csv");
    }
    traj << std::setprecision(15);
    traj
        << "t,"
        << "delta_x,"
        << "abs_delta_x,"
        << "theta,"
        << "Wbar,"
        << "velocity,"
        << "xc,"
        << "x_rest"
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
                << XC(0) << ","
                << x_rest
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
        << " eps=" << eps_shift
        << " folder=" << out_dir
        << std::endl;
}

int main(int argc, char** argv)
{
    // Usage:
    //   build/figure6_grid_aligned OUT_DIR SMIN SMAX DS0
    //       T_RELAX T_FREE DT SAVE_EVERY THREADS EPS_SHIFT
    if (argc != 11) {
        std::cerr
            << "Usage:\n"
            << "  " << argv[0]
            << " OUT_DIR SMIN SMAX DS0 T_RELAX T_FREE DT"
            << " SAVE_EVERY THREADS EPS_SHIFT\n";
        return 2;
    }

    const std::string out_root = argv[1];
    const double smin = std::atof(argv[2]);
    const double smax = std::atof(argv[3]);
    const double ds0 = std::atof(argv[4]);
    const double t_relax = std::atof(argv[5]);
    const double t_free = std::atof(argv[6]);
    const double dt = std::atof(argv[7]);
    const int save_every = std::atoi(argv[8]);
    const int threads = std::atoi(argv[9]);
    const double eps_shift = std::atof(argv[10]);

    if (!(smax >= smin) || !(ds0 > 0.0) || !(dt > 0.0) ||
        !(t_relax >= 0.0) || !(t_free > 0.0) || save_every <= 0 ||
        threads <= 0) {
        std::cerr << "Invalid numerical arguments.\n";
        return 2;
    }

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    ensure_dir(out_root);

    std::vector<PerturbationCase> cases;
    for (double s = smin; s <= smax + 1.0e-12; s += ds0) {
        cases.push_back(PerturbationCase{s, s0_folder_name(s)});
    }

    std::ofstream summary(out_root + "/perturbation_cases.csv");
    if (!summary) {
        std::cerr << "Failed to open perturbation_cases.csv\n";
        return 1;
    }
    summary << std::setprecision(15);
    summary
        << "s0,folder,eps_shift,t_relax,t_free,dt,save_every,"
        << "center_mode,x_rest,Nx,dx\n";

    const double x_rest = grid_aligned_single_droplet_center();
    for (const auto& cas : cases) {
        summary
            << cas.s0 << ","
            << cas.folder << ","
            << eps_shift << ","
            << t_relax << ","
            << t_free << ","
            << dt << ","
            << save_every << ","
            << "grid_aligned,"
            << x_rest << ","
            << Nx << ","
            << dx << "\n";
    }

    std::cerr
        << "Grid-aligned perturbation validation\n"
        << "out_root=" << out_root
        << " smin=" << smin
        << " smax=" << smax
        << " ds0=" << ds0
        << " t_relax=" << t_relax
        << " t_free=" << t_free
        << " dt=" << dt
        << " save_every=" << save_every
        << " eps=" << eps_shift
        << " x_rest=" << x_rest
        << " Nx=" << Nx
        << " dx=" << dx
        << " cases=" << cases.size()
        << " threads=" << threads
        << std::endl;

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
