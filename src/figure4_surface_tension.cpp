#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct Fig4Case {
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

inline double gamma_u_part(double u)
{
    return 0.5 * ((aa * aa) / (aa * aa + u * u)) * (G0 - G1);
}

inline double gamma_s_part(double s)
{
    return 0.5 * (bb / (bb + s)) * (G0 - G1);
}

inline double gamma_total_from_parts(double u, double s)
{
    return gamma_u_part(u) + gamma_s_part(s) + G1;
}

inline void write_surface_tension_state(
    std::ofstream& ofs,
    const std::vector<double>& y,
    double t,
    double x_unwrapped
){
    const double* Y = y.data();

    const double xw = XC(0);
    const double theta_wrapped = 2.0 * PI_CONST * xw / Lx;
    const double theta_unwrapped = 2.0 * PI_CONST * x_unwrapped / Lx;

    const double x_left = XC(0) - rad;
    const double x_right = XC(0) + rad;

    const double u_left = interp_periodic_field(Y, UBEGIN, x_left);
    const double s_left = interp_periodic_field(Y, SBEGIN, x_left);
    const double u_right = interp_periodic_field(Y, UBEGIN, x_right);
    const double s_right = interp_periodic_field(Y, SBEGIN, x_right);

    const double g_left = gamma_u_part(u_left);
    const double g_right = gamma_u_part(u_right);

    const double h_left = gamma_s_part(s_left);
    const double h_right = gamma_s_part(s_right);

    const double gamma_left = gamma_total_from_parts(u_left, s_left);
    const double gamma_right = gamma_total_from_parts(u_right, s_right);

    const double delta_g = g_right - g_left;
    const double delta_h = h_right - h_left;
    const double delta_gamma = gamma_right - gamma_left;

    const double velocity = droplet_velocity_one(Y, 0);

    ofs
        << t << ","
        << xw << ","
        << x_unwrapped << ","
        << theta_wrapped << ","
        << theta_unwrapped << ","
        << velocity << ","
        << WBAR << ","
        << u_left << ","
        << u_right << ","
        << s_left << ","
        << s_right << ","
        << g_left << ","
        << g_right << ","
        << h_left << ","
        << h_right << ","
        << gamma_left << ","
        << gamma_right << ","
        << delta_g << ","
        << delta_h << ","
        << delta_gamma << ","
        << delta_g / G0 << ","
        << delta_h / G0 << ","
        << delta_gamma / G0
        << "\n";
}

inline void run_surface_tension_case(
    const Fig4Case& cas,
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

    std::ofstream ofs(out_dir + "/surface_tension.csv");
    ofs << std::setprecision(15);
    ofs
        << "t,"
        << "x_wrapped,"
        << "x_unwrapped,"
        << "theta_wrapped,"
        << "theta_unwrapped,"
        << "velocity,"
        << "Wbar,"
        << "u_left,"
        << "u_right,"
        << "s_left,"
        << "s_right,"
        << "g_left,"
        << "g_right,"
        << "h_left,"
        << "h_right,"
        << "gamma_left,"
        << "gamma_right,"
        << "delta_g,"
        << "delta_h,"
        << "delta_gamma,"
        << "delta_g_norm,"
        << "delta_h_norm,"
        << "delta_gamma_norm"
        << "\n";

    std::ofstream info(out_dir + "/info.txt");
    info << "s0 " << std::setprecision(15) << cas.s0 << "\n";
    info << "label " << cas.label << "\n";
    info << "t_end " << t_end << "\n";
    info << "dt " << dt << "\n";
    info << "save_every " << save_every << "\n";

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;

        if (step % save_every == 0) {
            write_surface_tension_state(ofs, y, t, x_unwrapped);
        }

        if (step == nsteps) {
            break;
        }

        const double old_wrapped = y[IXC(0)];

        imex_lu_step(y, solvers, dt, cas.s0, kw, true);

        const double new_wrapped = y[IXC(0)];
        const double dx_move = pdiff(new_wrapped, old_wrapped);
        x_unwrapped += dx_move;

        if (step % (save_every * 200) == 0) {
            std::cerr
                << "case " << cas.label
                << " s0=" << cas.s0
                << " t=" << t
                << " x=" << new_wrapped
                << " Wbar=" << y[IWBAR]
                << std::endl;
        }
    }

    save_snapshot(out_dir, 0, y, cas.s0);
    save_lastdata(out_dir + "/lastdata", y);
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/figure4_lu out_dir t_end dt save_every threads [s0_1 s0_2 s0_3 ...]

        Default:
            s0 = 0.08, 0.125, 0.18
    */

    std::string out_root = "figure4_data";
    double t_end = 200.0;
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

    std::vector<Fig4Case> cases;

    if (argc > 6) {
        for (int i = 6; i < argc; ++i) {
            Fig4Case cas;
            cas.s0 = std::atof(argv[i]);
            cas.label = std::string("case ") + std::to_string(i - 5);
            cas.folder = s0_folder_name(cas.s0);
            cases.push_back(cas);
        }
    } else {
        cases.push_back(Fig4Case{0.080, "stationary",  s0_folder_name(0.080)});
        cases.push_back(Fig4Case{0.125, "oscillatory", s0_folder_name(0.125)});
        cases.push_back(Fig4Case{0.180, "rotating",    s0_folder_name(0.180)});
    }

    ensure_dir(out_root);

    std::ofstream summary(out_root + "/surface_tension_cases.csv");
    summary << "s0,label,folder\n";

    for (const auto& cas : cases) {
        summary
            << std::setprecision(15)
            << cas.s0 << ","
            << cas.label << ","
            << cas.folder << "\n";
    }

    std::cerr << "Figure 4 surface-tension simulations\n";
    std::cerr << "out_root=" << out_root
              << " t_end=" << t_end
              << " dt=" << dt
              << " save_every=" << save_every
              << " cases=" << cases.size()
              << std::endl;

    for (std::size_t k = 0; k < cases.size(); ++k) {
        const unsigned seed = 987u + static_cast<unsigned>(k);

        run_surface_tension_case(
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
