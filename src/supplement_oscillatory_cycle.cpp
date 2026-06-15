#include "model_lu_tbb.hpp"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

inline std::string snap_folder_name(int snap_id)
{
    std::ostringstream ss;
    ss << "snap_" << std::setw(5) << std::setfill('0') << snap_id;
    return ss.str();
}

inline void write_snapshot_metadata(
    std::ofstream& ofs,
    int snap_id,
    const std::string& folder,
    const std::vector<double>& y,
    double t,
    double x_unwrapped
){
    const double* Y = y.data();

    const double xw = XC(0);
    const double theta_wrapped = 2.0 * PI_CONST * xw / Lx;
    const double theta_unwrapped = 2.0 * PI_CONST * x_unwrapped / Lx;
    const double velocity = droplet_velocity_one(Y, 0);

    ofs
        << snap_id << ","
        << t << ","
        << folder << ","
        << xw << ","
        << x_unwrapped << ","
        << theta_wrapped << ","
        << theta_unwrapped << ","
        << velocity << ","
        << WBAR
        << "\n";
}

inline void write_trajectory_state(
    std::ofstream& ofs,
    const std::vector<double>& y,
    double t,
    double x_unwrapped
){
    const double* Y = y.data();

    const double xw = XC(0);
    const double theta_wrapped = 2.0 * PI_CONST * xw / Lx;
    const double theta_unwrapped = 2.0 * PI_CONST * x_unwrapped / Lx;
    const double velocity = droplet_velocity_one(Y, 0);

    ofs
        << t << ","
        << xw << ","
        << x_unwrapped << ","
        << theta_wrapped << ","
        << theta_unwrapped << ","
        << velocity << ","
        << WBAR
        << "\n";
}

int main(int argc, char** argv)
{
    /*
        Usage:
            ./bin/figure7_lu out_dir s0 t_end final_window dt save_every threads

        Example:
            ./bin/figure7_lu figure7_data 0.125 200 2.0 0.0005 20 8
    */

    std::string out_root = "figure7_data";

    double s0 = 0.125;
    double t_end = 200.0;
    double final_window = 2.0;

    double dt = default_dt;
    int save_every = 20;
    int threads = 8;

    if (argc >= 2) out_root = argv[1];
    if (argc >= 3) s0 = std::atof(argv[2]);
    if (argc >= 4) t_end = std::atof(argv[3]);
    if (argc >= 5) final_window = std::atof(argv[4]);
    if (argc >= 6) dt = std::atof(argv[5]);
    if (argc >= 7) save_every = std::atoi(argv[6]);
    if (argc >= 8) threads = std::atoi(argv[7]);

#ifdef USE_REAL_TBB
    ptbb::global_control control(
        ptbb::global_control::max_allowed_parallelism,
        static_cast<std::size_t>(threads)
    );
#endif

    ensure_dir(out_root);
    ensure_dir(out_root + "/snapshots");

    std::vector<double> y;
    initialize_state(y, s0, 2468u, false);

    SolverPack solvers(dt);

    const int nsteps = static_cast<int>(std::ceil(t_end / dt));
    const double t_snap_start = std::max(0.0, t_end - final_window);

    double x_unwrapped = y[IXC(0)];

    std::ofstream traj(out_root + "/trajectory.csv");
    traj << std::setprecision(15);
    traj
        << "t,"
        << "x_wrapped,"
        << "x_unwrapped,"
        << "theta_wrapped,"
        << "theta_unwrapped,"
        << "velocity,"
        << "Wbar"
        << "\n";

    std::ofstream meta(out_root + "/snapshot_index.csv");
    meta << std::setprecision(15);
    meta
        << "snap_id,"
        << "t,"
        << "folder,"
        << "x_wrapped,"
        << "x_unwrapped,"
        << "theta_wrapped,"
        << "theta_unwrapped,"
        << "velocity,"
        << "Wbar"
        << "\n";

    std::ofstream info(out_root + "/info.txt");
    info << std::setprecision(15);
    info << "s0 " << s0 << "\n";
    info << "t_end " << t_end << "\n";
    info << "final_window " << final_window << "\n";
    info << "dt " << dt << "\n";
    info << "save_every " << save_every << "\n";
    info << "Kw_open " << Kw_open << "\n";

    int snap_id = 0;

    std::cerr
        << "Figure 7 oscillatory-cycle simulation\n"
        << "out_root=" << out_root
        << " s0=" << s0
        << " t_end=" << t_end
        << " final_window=" << final_window
        << " dt=" << dt
        << " save_every=" << save_every
        << std::endl;

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;

        if (step % save_every == 0) {
            write_trajectory_state(traj, y, t, x_unwrapped);

            if (t >= t_snap_start) {
                const std::string folder = snap_folder_name(snap_id);
                const std::string snap_dir = out_root + "/snapshots/" + folder;

                save_snapshot(snap_dir, 0, y, s0);
                write_snapshot_metadata(meta, snap_id, folder, y, t, x_unwrapped);

                ++snap_id;
            }
        }

        if (step == nsteps) {
            break;
        }

        const double old_wrapped = y[IXC(0)];

        imex_lu_step(y, solvers, dt, s0, Kw_open, true);

        const double new_wrapped = y[IXC(0)];
        const double dx_move = pdiff(new_wrapped, old_wrapped);

        x_unwrapped += dx_move;

        if (step % (save_every * 500) == 0) {
            std::cerr
                << "t=" << t
                << " x=" << new_wrapped
                << " Wbar=" << y[IWBAR]
                << " snaps=" << snap_id
                << std::endl;
        }
    }

    save_lastdata(out_root + "/lastdata_final", y);

    std::cerr
        << "Done. Saved " << snap_id
        << " snapshots to " << out_root << "/snapshots"
        << std::endl;

    return 0;
}
