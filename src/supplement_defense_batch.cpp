#include "model_lu_tbb.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cerrno>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct OnlineLateStats {
    double y_min = std::numeric_limits<double>::infinity();
    double y_max = -std::numeric_limits<double>::infinity();
    double v2_sum = 0.0;
    double abs_v_sum = 0.0;
    double first_unwrapped = 0.0;
    double last_unwrapped = 0.0;
    bool has_first = false;
    long long n = 0;

    void add(double y, double v, double unwrapped_theta) {
        y_min = std::min(y_min, y);
        y_max = std::max(y_max, y);
        v2_sum += v * v;
        abs_v_sum += std::abs(v);
        if (!has_first) {
            first_unwrapped = unwrapped_theta;
            has_first = true;
        }
        last_unwrapped = unwrapped_theta;
        ++n;
    }

    double amplitude() const {
        if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
        return y_max - y_min;
    }

    double vrms() const {
        if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
        return std::sqrt(v2_sum / static_cast<double>(n));
    }

    double mean_abs_v() const {
        if (n <= 0) return std::numeric_limits<double>::quiet_NaN();
        return abs_v_sum / static_cast<double>(n);
    }

    double wraps() const {
        if (!has_first) return std::numeric_limits<double>::quiet_NaN();
        return std::abs(last_unwrapped - first_unwrapped) / (2.0 * M_PI);
    }
};

inline double wrap_x(double x) {
    double y = std::fmod(x, Lx);
    if (y < 0.0) y += Lx;
    return y;
}

inline double theta_wrapped_from_x(double x) {
    double th = 2.0 * M_PI * wrap_x(x) / Lx;
    th = std::fmod(th, 2.0 * M_PI);
    if (th < 0.0) th += 2.0 * M_PI;
    return th;
}

inline double theta_unwrapped_from_x(double x_unwrapped) {
    return 2.0 * M_PI * x_unwrapped / Lx;
}

inline std::string classify_state(double amp_theta, double wraps_late, double vrms_late) {
    // Conservative thresholds for nondimensional simulations.
    // The paper only needs qualitative regime classification.
    if (std::isfinite(wraps_late) && wraps_late > 0.75) return "rotating";
    if (std::isfinite(amp_theta) && std::isfinite(vrms_late) && amp_theta < 1.0e-3 && vrms_late < 1.0e-4) return "stationary";
    return "oscillatory";
}

inline bool path_exists(const std::string& path) {
    struct stat st;
    return ::stat(path.c_str(), &st) == 0;
}

inline bool file_is_empty_or_missing(const std::string& file) {
    struct stat st;
    if (::stat(file.c_str(), &st) != 0) return true;
    return st.st_size == 0;
}

inline bool is_directory(const std::string& path) {
    struct stat st;
    if (::stat(path.c_str(), &st) != 0) return false;
    return S_ISDIR(st.st_mode);
}

inline void create_directories_compat(const std::string& dir) {
    if (dir.empty()) return;
    if (is_directory(dir)) return;

    std::size_t start = 0;

    // Preserve an absolute-path leading slash.
    if (!dir.empty() && dir[0] == '/') {
        start = 1;
    }

    while (start <= dir.size()) {
        std::size_t pos = dir.find('/', start);
        std::string part;

        if (pos == std::string::npos) {
            part = dir;
            start = dir.size() + 1;
        } else {
            part = dir.substr(0, pos);
            start = pos + 1;
        }

        if (part.empty() || part == "/") continue;

        if (!is_directory(part)) {
            if (::mkdir(part.c_str(), 0755) != 0 && errno != EEXIST) {
                std::cerr << "mkdir failed for " << part << ": " << std::strerror(errno) << std::endl;
                std::exit(1);
            }
        }
    }
}

inline std::string parent_path_compat(const std::string& file) {
    std::size_t pos = file.find_last_of('/');
    if (pos == std::string::npos) return std::string();
    if (pos == 0) return std::string("/");
    return file.substr(0, pos);
}

inline void ensure_parent(const std::string& file) {
    const std::string parent = parent_path_compat(file);
    if (!parent.empty()) create_directories_compat(parent);
}

inline void append_header_if_new(const std::string& file, const std::string& header) {
    if (file_is_empty_or_missing(file)) {
        ensure_parent(file);
        std::ofstream ofs(file, std::ios::app);
        ofs << header << "\n";
    }
}

struct OpenMetrics {
    std::string suite;
    double s0 = 0.0;
    double dt = 0.0;
    double t_end = 0.0;
    double late_window = 0.0;
    double kw = 0.0;
    unsigned seed = 0;
    double amp_theta = 0.0;
    double vrms_late = 0.0;
    double mean_abs_v = 0.0;
    double wraps_late = 0.0;
    double theta_drift_late = 0.0;
    double wbar_end = 0.0;
    std::string state;
};

OpenMetrics run_open_case(const std::string& suite, double s0, double dt, double t_end, double late_window, double kw, unsigned seed) {
    std::vector<double> y;
    initialize_state(y, s0, seed, false);
    SolverPack solvers(dt);

    const int nsteps = static_cast<int>(std::ceil(t_end / dt));
    double x_unwrapped = y[IXC(0)];
    OnlineLateStats stats;

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;
        const double* Y = y.data();
        const double xw = y[IXC(0)];
        const double theta_w = theta_wrapped_from_x(xw);
        const double theta_u = theta_unwrapped_from_x(x_unwrapped);
        const double v = droplet_velocity_one(Y, 0);

        if (t >= t_end - late_window) {
            stats.add(theta_w, v, theta_u);
        }

        if (step == nsteps) break;

        const double old_x = y[IXC(0)];
        imex_lu_step(y, solvers, dt, s0, kw, true);
        const double new_x = y[IXC(0)];
        x_unwrapped += pdiff(new_x, old_x);
    }

    OpenMetrics m;
    m.suite = suite;
    m.s0 = s0;
    m.dt = dt;
    m.t_end = t_end;
    m.late_window = late_window;
    m.kw = kw;
    m.seed = seed;
    m.amp_theta = stats.amplitude();
    m.vrms_late = stats.vrms();
    m.mean_abs_v = stats.mean_abs_v();
    m.wraps_late = stats.wraps();
    m.theta_drift_late = stats.last_unwrapped - stats.first_unwrapped;
    m.wbar_end = y[IWBAR];
    m.state = classify_state(m.amp_theta, m.wraps_late, m.vrms_late);

    return m;
}

void write_open_metrics(const std::string& out_dir, const OpenMetrics& m) {
    const std::string file = out_dir + "/open_case_metrics.csv";
    append_header_if_new(file, "suite,Nx,dx,dt,s0,t_end,late_window,kw,seed,amp_theta,vrms_late,mean_abs_v,wraps_late,theta_drift_late,wbar_end,state");
    std::ofstream ofs(file, std::ios::app);
    ofs << std::setprecision(15)
        << m.suite << ","
        << Nx << ","
        << dx << ","
        << m.dt << ","
        << m.s0 << ","
        << m.t_end << ","
        << m.late_window << ","
        << m.kw << ","
        << m.seed << ","
        << m.amp_theta << ","
        << m.vrms_late << ","
        << m.mean_abs_v << ","
        << m.wraps_late << ","
        << m.theta_drift_late << ","
        << m.wbar_end << ","
        << m.state << "\n";
}

struct PerturbMetrics {
    std::string suite;
    double s0 = 0.0;
    double eps = 0.0;
    double dt = 0.0;
    double t_relax = 0.0;
    double t_free = 0.0;
    double late_window = 0.0;
    double kw = 0.0;
    unsigned seed = 0;
    double late_delta_amp = 0.0;
    double late_vrms = 0.0;
    double max_abs_delta = 0.0;
    double max_abs_v = 0.0;
    double x_rest = 0.0;
    double x_final = 0.0;
};

PerturbMetrics run_perturb_case(const std::string& suite, double s0, double eps, double dt, double t_relax, double t_free, double late_window, double kw, unsigned seed) {
    std::vector<double> y;
    initialize_state(y, s0, seed, false);
    SolverPack solvers(dt);

    const int nrelax = static_cast<int>(std::ceil(t_relax / dt));
    for (int step = 0; step < nrelax; ++step) {
        imex_lu_step(y, solvers, dt, s0, kw, false); // fixed droplet, relax chemistry
    }

    const double x_rest = y[IXC(0)];
    y[IXC(0)] = wrap_x(x_rest + eps);

    OnlineLateStats stats;
    double max_abs_delta = 0.0;
    double max_abs_v = 0.0;

    const int nfree = static_cast<int>(std::ceil(t_free / dt));
    for (int step = 0; step <= nfree; ++step) {
        const double t = step * dt;
        const double* Y = y.data();
        const double delta = pdiff(y[IXC(0)], x_rest);
        const double v = droplet_velocity_one(Y, 0);
        const double dummy_theta = theta_unwrapped_from_x(delta);

        max_abs_delta = std::max(max_abs_delta, std::abs(delta));
        max_abs_v = std::max(max_abs_v, std::abs(v));

        if (t >= t_free - late_window) {
            stats.add(delta, v, dummy_theta);
        }

        if (step == nfree) break;
        imex_lu_step(y, solvers, dt, s0, kw, true);
    }

    PerturbMetrics m;
    m.suite = suite;
    m.s0 = s0;
    m.eps = eps;
    m.dt = dt;
    m.t_relax = t_relax;
    m.t_free = t_free;
    m.late_window = late_window;
    m.kw = kw;
    m.seed = seed;
    m.late_delta_amp = stats.amplitude();
    m.late_vrms = stats.vrms();
    m.max_abs_delta = max_abs_delta;
    m.max_abs_v = max_abs_v;
    m.x_rest = x_rest;
    m.x_final = y[IXC(0)];
    return m;
}

void write_perturb_metrics(const std::string& out_dir, const PerturbMetrics& m) {
    const std::string file = out_dir + "/perturbation_metrics.csv";
    append_header_if_new(file, "suite,Nx,dx,dt,s0,eps,t_relax,t_free,late_window,kw,seed,late_delta_amp,late_vrms,max_abs_delta,max_abs_v,x_rest,x_final");
    std::ofstream ofs(file, std::ios::app);
    ofs << std::setprecision(15)
        << m.suite << ","
        << Nx << ","
        << dx << ","
        << m.dt << ","
        << m.s0 << ","
        << m.eps << ","
        << m.t_relax << ","
        << m.t_free << ","
        << m.late_window << ","
        << m.kw << ","
        << m.seed << ","
        << m.late_delta_amp << ","
        << m.late_vrms << ","
        << m.max_abs_delta << ","
        << m.max_abs_v << ","
        << m.x_rest << ","
        << m.x_final << "\n";
}

struct SegmentStats {
    OnlineLateStats activity;
    double wbar_end = 0.0;
    double x_end = 0.0;
};

SegmentStats run_cover_segment(std::vector<double>& y, SolverPack& solvers, double s0, double dt, double t_segment, double late_window, double kw) {
    SegmentStats seg;
    const int nsteps = static_cast<int>(std::ceil(t_segment / dt));
    double x_unwrapped = y[IXC(0)];

    for (int step = 0; step <= nsteps; ++step) {
        const double t = step * dt;
        const double* Y = y.data();
        const double theta = theta_wrapped_from_x(y[IXC(0)]);
        const double theta_u = theta_unwrapped_from_x(x_unwrapped);
        const double v = droplet_velocity_one(Y, 0);

        if (t >= t_segment - late_window) {
            seg.activity.add(theta, v, theta_u);
        }

        if (step == nsteps) break;
        const double old_x = y[IXC(0)];
        imex_lu_step(y, solvers, dt, s0, kw, true);
        const double new_x = y[IXC(0)];
        x_unwrapped += pdiff(new_x, old_x);
    }

    seg.wbar_end = y[IWBAR];
    seg.x_end = y[IXC(0)];
    return seg;
}

void run_cover_case(const std::string& out_dir, const std::string& suite, double s0, double dt, double t_open, double t_closed, double t_reopen, double kw_open, double kw_close, unsigned seed) {
    std::vector<double> y;
    initialize_state(y, s0, seed, false);
    SolverPack solvers(dt);

    const double late_open = std::min(5.0, 0.5 * t_open);
    const double late_closed = std::min(5.0, 0.5 * t_closed);
    const double late_reopen = std::min(5.0, 0.5 * t_reopen);

    SegmentStats a = run_cover_segment(y, solvers, s0, dt, t_open, late_open, kw_open);
    SegmentStats b = run_cover_segment(y, solvers, s0, dt, t_closed, late_closed, kw_close);
    SegmentStats c = run_cover_segment(y, solvers, s0, dt, t_reopen, late_reopen, kw_open);

    const std::string file = out_dir + "/cover_sensitivity.csv";
    append_header_if_new(file, "suite,Nx,dx,dt,s0,kw_open,kw_close,t_open,t_closed,t_reopen,seed,open_amp,open_vrms,open_wbar_end,closed_amp,closed_vrms,closed_wbar_end,reopen_amp,reopen_vrms,reopen_wbar_end");
    std::ofstream ofs(file, std::ios::app);
    ofs << std::setprecision(15)
        << suite << ","
        << Nx << ","
        << dx << ","
        << dt << ","
        << s0 << ","
        << kw_open << ","
        << kw_close << ","
        << t_open << ","
        << t_closed << ","
        << t_reopen << ","
        << seed << ","
        << a.activity.amplitude() << ","
        << a.activity.vrms() << ","
        << a.wbar_end << ","
        << b.activity.amplitude() << ","
        << b.activity.vrms() << ","
        << b.wbar_end << ","
        << c.activity.amplitude() << ","
        << c.activity.vrms() << ","
        << c.wbar_end << "\n";
}

std::vector<double> parse_doubles(int start, int argc, char** argv) {
    std::vector<double> out;
    for (int i = start; i < argc; ++i) out.push_back(std::atof(argv[i]));
    return out;
}

void usage() {
    std::cerr
        << "Usage:\n"
        << "  defense_batch_lu open OUT_DIR SUITE dt t_end late_window kw s0...\n"
        << "  defense_batch_lu perturb OUT_DIR SUITE dt t_relax t_free late_window kw eps s0...\n"
        << "  defense_batch_lu cover OUT_DIR SUITE dt t_open t_closed t_reopen kw_open kw_close s0\n";
}

int main(int argc, char** argv) {
    if (argc < 2) {
        usage();
        return 1;
    }

#ifdef USE_REAL_TBB
    // Respect the usual TBB_NUM_THREADS environment variable if set.
    const char* env_threads = std::getenv("TBB_NUM_THREADS");
    if (env_threads) {
        int nthreads = std::max(1, std::atoi(env_threads));
        static ptbb::global_control control(ptbb::global_control::max_allowed_parallelism, static_cast<std::size_t>(nthreads));
    }
#endif

    const std::string mode = argv[1];

    if (mode == "open") {
        if (argc < 9) {
            usage();
            return 1;
        }
        const std::string out_dir = argv[2];
        const std::string suite = argv[3];
        const double dt = std::atof(argv[4]);
        const double t_end = std::atof(argv[5]);
        const double late_window = std::atof(argv[6]);
        const double kw = std::atof(argv[7]);
        create_directories_compat(out_dir);
        auto s0s = parse_doubles(8, argc, argv);
        unsigned seed = 1001u;
        for (double s0 : s0s) {
            std::cerr << "[open] suite=" << suite << " Nx=" << Nx << " dt=" << dt << " s0=" << s0 << std::endl;
            OpenMetrics m = run_open_case(suite, s0, dt, t_end, late_window, kw, seed++);
            write_open_metrics(out_dir, m);
        }
        return 0;
    }

    if (mode == "perturb") {
        if (argc < 11) {
            usage();
            return 1;
        }
        const std::string out_dir = argv[2];
        const std::string suite = argv[3];
        const double dt = std::atof(argv[4]);
        const double t_relax = std::atof(argv[5]);
        const double t_free = std::atof(argv[6]);
        const double late_window = std::atof(argv[7]);
        const double kw = std::atof(argv[8]);
        const double eps = std::atof(argv[9]);
        create_directories_compat(out_dir);
        auto s0s = parse_doubles(10, argc, argv);
        unsigned seed = 2001u;
        for (double s0 : s0s) {
            std::cerr << "[perturb] suite=" << suite << " Nx=" << Nx << " dt=" << dt << " eps=" << eps << " s0=" << s0 << std::endl;
            PerturbMetrics m = run_perturb_case(suite, s0, eps, dt, t_relax, t_free, late_window, kw, seed++);
            write_perturb_metrics(out_dir, m);
        }
        return 0;
    }

    if (mode == "cover") {
        if (argc != 12) {
            usage();
            return 1;
        }
        const std::string out_dir = argv[2];
        const std::string suite = argv[3];
        const double dt = std::atof(argv[4]);
        const double t_open = std::atof(argv[5]);
        const double t_closed = std::atof(argv[6]);
        const double t_reopen = std::atof(argv[7]);
        const double kw_open = std::atof(argv[8]);
        const double kw_close = std::atof(argv[9]);
        const double s0 = std::atof(argv[10]);
        const unsigned seed = static_cast<unsigned>(std::atoi(argv[11]));
        create_directories_compat(out_dir);
        std::cerr << "[cover] suite=" << suite << " Nx=" << Nx << " dt=" << dt << " kw_open=" << kw_open << " s0=" << s0 << std::endl;
        run_cover_case(out_dir, suite, s0, dt, t_open, t_closed, t_reopen, kw_open, kw_close, seed);
        return 0;
    }

    usage();
    return 1;
}
