#include <cstdlib>
#include <string>
#include "model_lu_tbb.hpp"
int main(int argc,char** argv){
    double s0=(argc>=2)?std::atof(argv[1]):S0; std::string path=(argc>=3)?argv[2]:"data"; double dt=(argc>=4)?std::atof(argv[3]):default_dt; double tend=(argc>=5)?std::atof(argv[4]):open_time; int save_every=(argc>=6)?std::atoi(argv[5]):default_save_every; int threads=(argc>=7)?std::atoi(argv[6]):0;
    if(threads>0){ ptbb::global_control c(ptbb::global_control::max_allowed_parallelism,threads); run_amplitude(s0,Kw_open,tend,dt,save_every,11u,path); }
    else run_amplitude(s0,Kw_open,tend,dt,save_every,11u,path);
}
