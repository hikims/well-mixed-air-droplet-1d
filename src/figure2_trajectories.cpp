#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Fig2Case {
    double s0 = 0.0;
    std::string label;
    std::string folder;
};

inline std::string s0_folder_name(double s0)
{
    std::ostringstream ss;
    ss << "s0_" << std::fixed << std::setprecision(3) << s0;
    return ss.str();
}

inline void run_representative_case(
    const Fig2Case& cas,
    const std::string& out_root,
    double kw,
    double t_end,
    double dt,
    int save_every,
    unsigned seed
){
    const std::string out_dir = out_root + "/" + cas.folder;
    ensure_dir(out_dir);

    std::vector<double> y;
    initialize_state(y, cas.s0, seed, false);

    SolverPack solvers(dt);

    const int nsteps = static_cast<int>(std::ceil(t_end / dt));

    double x_unwrapped = y[IXC(0)];

    std::ofstream traj(out_dir + "/trajectory.csv");
    traj << std::setprecision(15);
    traj << "t,x_wrapped,x_unwrapped,theta_wrapped,theta_unwrapped,velocity,Wbar\n";

    std::ofstream info(out_dir + "/info.txt");
    info << "s0 " << std::setprecision(15) << cas.s0 << "\n";
    info << "label " << cas.label << "\n";
    info << "t_end " << t_end << "\n";
    info << "dt " << dt << "\n";
    info << "save_every " << save_every << "\n";

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;
        double* Y = y.data();

        if (step % save_every == 0) {
            const double xw = XC(0);
            const double theta_wrapped = 2.0 * PI_CONST * xw / Lx;
            const double theta_unwrapped = 2.0 * PI_CONST * x_unwrapped / Lx;
            const double vel = droplet_velocity_one(Y, 0);

            traj
                << t << ","
                << xw << ","
                << x_unwrapped << ","
                << theta_wrapped << ","
                << theta_unwrapped << ","
                << vel << ","
                << WBAR << "\n";
        }

        if (step == nsteps) {
            break;
        }

        const double old_wrapped = XC(0);

        imex_lu_step(y, solvers, dt, cas.s0, kw, true);

        const double new_wrapped = y[IXC(0)];
        const double dx_move = pdiff(new_wrapped, old_wrapped);

        x_unwrapped += dx_move;

        if (step % (save_every * 100) == 0) {
            std::cerr
                << "case " << cas.label
                << " s0=" << cas.s0
                << " t=" << t
                << " x=" << new_wrapped
                << " Wbar=" << y[IWBAR]
                << std::endl;
        }
    }

    /*
        Save final profile in the same format as opendata/u, s, c, w.
        This is useful if you later want representative profile figures.
    */
    save_snapshot(out_dir, 0, y, cas.s0);
    save_lastdata(out_dir + "/lastdata", y);
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/figure2_lu out_dir t_end dt save_every threads [s0_1 s0_2 s0_3 ...]

        Default:
            s0 = 0.08, 0.125, 0.18
    */

    std::string out_root = "figure2_data";
    double t_end = 120.0;
    double dt = default_dt;
    int save_every = 20;
    int threads = 8;

    if (argc >= 2) out_root = argv[1];
    if (argc >= 3) t_end = std::atof(argv[2]);
    if (argc >= 4) dt = std::atof(argv[3]);
    if (argc >= 5) save_every = std::atoi(argv[4]);
    if (argc >= 6) threads = std::atoi(argv[5]);

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    std::vector<Fig2Case> cases;

    if (argc > 6) {
        for (int i = 6; i < argc; ++i) {
            Fig2Case cas;
            cas.s0 = std::atof(argv[i]);
            cas.label = std::string("case ") + std::to_string(i - 5);
            cas.folder = s0_folder_name(cas.s0);
            cases.push_back(cas);
        }
    } else {
        cases.push_back(Fig2Case{0.080, "stationary",  s0_folder_name(0.080)});
        cases.push_back(Fig2Case{0.125, "oscillatory", s0_folder_name(0.125)});
        cases.push_back(Fig2Case{0.180, "rotating",    s0_folder_name(0.180)});
    }

    ensure_dir(out_root);

    std::ofstream summary(out_root + "/representative_cases.csv");
    summary << "s0,label,folder\n";

    for (const auto& cas : cases) {
        summary
            << std::setprecision(15)
            << cas.s0 << ","
            << cas.label << ","
            << cas.folder << "\n";
    }

    std::cerr << "Figure 2 representative simulations\n";
    std::cerr << "out_root=" << out_root
              << " t_end=" << t_end
              << " dt=" << dt
              << " save_every=" << save_every
              << " cases=" << cases.size()
              << std::endl;

    /*
        Representative cases are only three by default.
        Sequential execution is clearer and avoids heavy simultaneous output.
        The internal grid loops still use TBB.
    */
    for (std::size_t k = 0; k < cases.size(); ++k) {
        const unsigned seed = 777u + static_cast<unsigned>(k);
        run_representative_case(
            cases[k],
            out_root,
            Kw_open,
            t_end,
            dt,
            save_every,
            seed
        );
    }

    std::cerr << "Done. Data saved to " << out_root << std::endl;

    return 0;
}
