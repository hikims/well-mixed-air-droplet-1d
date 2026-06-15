#include <cstdlib>
#include <string>
#include "model_lu_tbb.hpp"
int main(int argc,char** argv){
    std::string last=(argc>=2)?argv[1]:""; double s0=(argc>=3)?std::atof(argv[2]):S0; std::string path=(argc>=4)?argv[3]:"opendata"; double dt=(argc>=5)?std::atof(argv[4]):default_dt; double tend=(argc>=6)?std::atof(argv[5]):open_time; int save_every=(argc>=7)?std::atoi(argv[6]):default_save_every; int threads=(argc>=8)?std::atoi(argv[7]):0;
    if(threads>0){ ptbb::global_control c(ptbb::global_control::max_allowed_parallelism,threads); run_cover_segment(last,"lastdata_open",path,"_open",s0,Kw_open,tend,dt,save_every,1u); }
    else run_cover_segment(last,"lastdata_open",path,"_open",s0,Kw_open,tend,dt,save_every,1u);
}
