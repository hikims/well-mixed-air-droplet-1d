#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

struct SegmentInfo {
    std::string name;
    double kw;
    double duration;
};

inline void write_cover_state(
    std::ofstream& ofs,
    const std::vector<double>& y,
    double t_global,
    double t_local,
    const std::string& segment_name,
    double x_unwrapped
){
    const double* Y = y.data();

    const double xw = XC(0);
    const double theta_wrapped = 2.0 * PI_CONST * xw / Lx;
    const double theta_unwrapped = 2.0 * PI_CONST * x_unwrapped / Lx;
    const double velocity = droplet_velocity_one(Y, 0);

    ofs
        << t_global << ","
        << t_local << ","
        << segment_name << ","
        << xw << ","
        << x_unwrapped << ","
        << theta_wrapped << ","
        << theta_unwrapped << ","
        << velocity << ","
        << WBAR << "\n";
}

inline void run_cover_segment_for_figure5(
    std::vector<double>& y,
    SolverPack& solvers,
    std::ofstream& traj,
    const SegmentInfo& seg,
    const std::string& out_root,
    double s0,
    double dt,
    int save_every,
    double& t_global,
    double& x_unwrapped
){
    const int nsteps = static_cast<int>(std::ceil(seg.duration / dt));

    std::cerr
        << "segment=" << seg.name
        << " duration=" << seg.duration
        << " kw=" << seg.kw
        << " t_global_start=" << t_global
        << std::endl;

    for (int step = 0; step <= nsteps; ++step) {
        const double t_local = step * dt;

        if (step % save_every == 0) {
            write_cover_state(
                traj,
                y,
                t_global + t_local,
                t_local,
                seg.name,
                x_unwrapped
            );
        }

        if (step == nsteps) {
            break;
        }

        const double old_wrapped = y[IXC(0)];

        imex_lu_step(y, solvers, dt, s0, seg.kw, true);

        const double new_wrapped = y[IXC(0)];
        const double dx_move = pdiff(new_wrapped, old_wrapped);
        x_unwrapped += dx_move;

        if (step % (save_every * 200) == 0) {
            std::cerr
                << "  t_local=" << t_local
                << " t_global=" << t_global + t_local
                << " x=" << new_wrapped
                << " Wbar=" << y[IWBAR]
                << std::endl;
        }
    }

    /*
        Save the profile at the end of this segment.
    */
    const std::string profile_dir = out_root + "/profile_" + seg.name + "_end";
    save_snapshot(profile_dir, 0, y, s0);

    /*
        Save binary lastdata for possible continuation.
    */
    save_lastdata(out_root + "/lastdata_" + seg.name, y);

    t_global += seg.duration;
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/figure5_cover_lu out_dir s0 open_time close_time reopen_time dt save_every threads

        Example:
            ./bin/figure5_cover_lu figure5_data 0.125 40 40 40 0.0005 20 8
    */

    std::string out_root = "figure5_data";

    double s0 = 0.125;
    double t_open = 40.0;
    double t_close = 40.0;
    double t_reopen = 40.0;

    double dt = default_dt;
    int save_every = 20;
    int threads = 8;

    if (argc >= 2) out_root = argv[1];
    if (argc >= 3) s0 = std::atof(argv[2]);
    if (argc >= 4) t_open = std::atof(argv[3]);
    if (argc >= 5) t_close = std::atof(argv[4]);
    if (argc >= 6) t_reopen = std::atof(argv[5]);
    if (argc >= 7) dt = std::atof(argv[6]);
    if (argc >= 8) save_every = std::atoi(argv[7]);
    if (argc >= 9) threads = std::atoi(argv[8]);

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    ensure_dir(out_root);

    std::vector<double> y;
    initialize_state(y, s0, 555u, false);

    SolverPack solvers(dt);

    std::ofstream traj(out_root + "/cover_protocol.csv");
    traj << std::setprecision(15);
    traj
        << "t_global,"
        << "t_local,"
        << "segment,"
        << "x_wrapped,"
        << "x_unwrapped,"
        << "theta_wrapped,"
        << "theta_unwrapped,"
        << "velocity,"
        << "Wbar"
        << "\n";

    std::ofstream info(out_root + "/cover_protocol_info.txt");
    info << std::setprecision(15);
    info << "s0 " << s0 << "\n";
    info << "open_time " << t_open << "\n";
    info << "close_time " << t_close << "\n";
    info << "reopen_time " << t_reopen << "\n";
    info << "dt " << dt << "\n";
    info << "save_every " << save_every << "\n";
    info << "Kw_open " << Kw_open << "\n";
    info << "Kw_close " << Kw_close << "\n";

    std::vector<SegmentInfo> segments;
    segments.push_back(SegmentInfo{"open", Kw_open, t_open});
    segments.push_back(SegmentInfo{"closed", Kw_close, t_close});
    segments.push_back(SegmentInfo{"reopened", Kw_open, t_reopen});

    double t_global = 0.0;
    double x_unwrapped = y[IXC(0)];

    for (const auto& seg : segments) {
        run_cover_segment_for_figure5(
            y,
            solvers,
            traj,
            seg,
            out_root,
            s0,
            dt,
            save_every,
            t_global,
            x_unwrapped
        );
    }

    save_lastdata(out_root + "/lastdata_final", y);

    std::cerr << "Figure 5 cover protocol data saved to " << out_root << std::endl;

    return 0;
}
