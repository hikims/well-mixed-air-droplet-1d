#include "model_lu_tbb.hpp"

#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct PhaseResult {
    double s0 = 0.0;
    double amplitude_theta = 0.0;
    double period = 0.0;
    double net_theta = 0.0;
    double mean_abs_velocity = 0.0;
    std::string state = "unknown";
};

inline PhaseResult run_one_phase_case(
    double s0,
    double kw,
    double t_end,
    double dt,
    unsigned seed
){
    std::vector<double> y;
    initialize_state(y, s0, seed, false);

    SolverPack solvers(dt);

    const int nsteps = static_cast<int>(std::ceil(t_end / dt));
    const int nstart = static_cast<int>(0.5 * nsteps);

    double prev_wrapped = 0.0;
    double x_unwrap = 0.0;
    double x_start = 0.0;
    double x_end = 0.0;

    double xmin = DBL_MAX;
    double xmax = -DBL_MAX;

    double sum_abs_v = 0.0;
    int count_v = 0;

    double prev_v = 0.0;
    bool have_prev_v = false;
    std::vector<double> sign_change_times;

    bool measuring = false;

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;
        double* Y = y.data();

        if (step == nstart) {
            prev_wrapped = XC(0);
            x_unwrap = XC(0);
            x_start = x_unwrap;
            x_end = x_unwrap;

            xmin = x_unwrap;
            xmax = x_unwrap;

            sum_abs_v = 0.0;
            count_v = 0;

            prev_v = 0.0;
            have_prev_v = false;
            sign_change_times.clear();

            measuring = true;
        }

        if (measuring && step > nstart) {
            /*
                Convert periodic wrapped position into an unwrapped trajectory.
                pdiff(XC(0), prev_wrapped) gives the shortest periodic displacement.
            */
            const double dx_move = pdiff(XC(0), prev_wrapped);
            x_unwrap += dx_move;
            prev_wrapped = XC(0);

            xmin = std::min(xmin, x_unwrap);
            xmax = std::max(xmax, x_unwrap);
            x_end = x_unwrap;

            const double v = droplet_velocity_one(Y, 0);
            sum_abs_v += std::abs(v);
            ++count_v;

            if (have_prev_v) {
                if ((prev_v < 0.0 && v >= 0.0) || (prev_v > 0.0 && v <= 0.0)) {
                    sign_change_times.push_back(t);
                }
            }

            prev_v = v;
            have_prev_v = true;
        }

        if (step == nsteps) {
            break;
        }

        imex_lu_step(y, solvers, dt, s0, kw, true);
    }

    PhaseResult r;
    r.s0 = s0;

    const double amp_x = xmax - xmin;
    const double net_x = std::abs(x_end - x_start);

    r.amplitude_theta = 2.0 * PI_CONST * amp_x / Lx;
    r.net_theta = 2.0 * PI_CONST * net_x / Lx;

    r.mean_abs_velocity =
        count_v > 0 ? sum_abs_v / static_cast<double>(count_v) : 0.0;

    /*
        Period estimate:
        For oscillatory motion, velocity sign changes occur approximately
        twice per period.
    */
    if (sign_change_times.size() >= 4) {
        double sum_half_period = 0.0;
        int cnt = 0;

        for (std::size_t k = 1; k < sign_change_times.size(); ++k) {
            sum_half_period += sign_change_times[k] - sign_change_times[k - 1];
            ++cnt;
        }

        if (cnt > 0) {
            r.period = 2.0 * sum_half_period / static_cast<double>(cnt);
        }
    }

    /*
        Classification thresholds.
        Adjust these if the boundary looks too noisy.
    */
    const double stationary_threshold = 0.05;          // in theta
    const double rotating_net_threshold = 1.6 * PI_CONST;

    if (r.net_theta > rotating_net_threshold || net_x > 0.8 * Lx) {
        r.state = "rotating";
        r.amplitude_theta = 2.0 * PI_CONST;
    } else if (r.amplitude_theta < stationary_threshold) {
        r.state = "stationary";
        r.amplitude_theta = 0.0;
        r.period = 0.0;
    } else {
        r.state = "oscillatory";
    }

    return r;
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/phase_sweep_lu out_csv smin smax ds0 t_end dt threads

        Example:
            ./bin/phase_sweep_lu results/phase_diagram_long.csv 0.02 0.40 0.005 120 0.0005 8
    */

    std::string out_csv = "results/phase_diagram.csv";

    double smin = 0.02;
    double smax = 0.40;
    double ds0 = 0.005;

    double t_end = 120.0;
    double dt = default_dt;

    int threads = 8;

    if (argc >= 2) out_csv = argv[1];
    if (argc >= 3) smin = std::atof(argv[2]);
    if (argc >= 4) smax = std::atof(argv[3]);
    if (argc >= 5) ds0 = std::atof(argv[4]);
    if (argc >= 6) t_end = std::atof(argv[5]);
    if (argc >= 7) dt = std::atof(argv[6]);
    if (argc >= 8) threads = std::atoi(argv[7]);

    ensure_dir("results");

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    std::vector<double> svals;
    for (double s = smin; s <= smax + 1.0e-12; s += ds0) {
        svals.push_back(s);
    }

    std::vector<PhaseResult> results(svals.size());

    std::cerr << "Phase sweep started\n";
    std::cerr << "out_csv=" << out_csv
              << " smin=" << smin
              << " smax=" << smax
              << " ds=" << ds0
              << " t_end=" << t_end
              << " dt=" << dt
              << " cases=" << svals.size()
              << " threads=" << threads
              << std::endl;

    ptbb::parallel_for(
        ptbb::blocked_range<int>(0, static_cast<int>(svals.size())),
        [&](const ptbb::blocked_range<int>& range) {
            for (int k = range.begin(); k != range.end(); ++k) {
                const double s0 = svals[k];
                const unsigned seed = 1234u + static_cast<unsigned>(k);

                results[k] = run_one_phase_case(
                    s0,
                    Kw_open,
                    t_end,
                    dt,
                    seed
                );

                std::cerr
                    << "done s0=" << s0
                    << " amp=" << results[k].amplitude_theta
                    << " net=" << results[k].net_theta
                    << " period=" << results[k].period
                    << " state=" << results[k].state
                    << std::endl;
            }
        }
    );

    std::ofstream ofs(out_csv);
    ofs << std::setprecision(15);

    ofs << "s0,amplitude_theta,period,net_theta,mean_abs_velocity,state\n";

    for (const auto& r : results) {
        ofs
            << r.s0 << ","
            << r.amplitude_theta << ","
            << r.period << ","
            << r.net_theta << ","
            << r.mean_abs_velocity << ","
            << r.state << "\n";
    }

    std::cerr << "saved: " << out_csv << std::endl;

    return 0;
}
