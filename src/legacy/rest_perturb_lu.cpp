#include <cstdlib>
#include <string>
#include "model_lu_tbb.hpp"
int main(int argc,char** argv){
    double s0=(argc>=2)?std::atof(argv[1]):S0; std::string path=(argc>=3)?argv[2]:"restdata"; double dt=(argc>=4)?std::atof(argv[3]):default_dt; double trelax=(argc>=5)?std::atof(argv[4]):rest_relax_time; double tfree=(argc>=6)?std::atof(argv[5]):rest_free_time; double eps=(argc>=7)?std::atof(argv[6]):eps_shift_default; int save_every=(argc>=8)?std::atoi(argv[7]):default_save_every; int threads=(argc>=9)?std::atoi(argv[8]):0;
    if(threads>0){ ptbb::global_control c(ptbb::global_control::max_allowed_parallelism,threads); run_rest_perturb(s0,Kw_open,trelax,tfree,dt,save_every,eps,21u,path); }
    else run_rest_perturb(s0,Kw_open,trelax,tfree,dt,save_every,eps,21u,path);
}
