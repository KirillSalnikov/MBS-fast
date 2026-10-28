#define main publication_original_main
#include "fullauto_native_publication_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    Require(argc==3||argc==4,"OUTPUT WORKER [COUNT]");const int count=argc==4?std::stoi(argv[3]):64;const std::string root=argv[1];
    for(bool grouped:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);
        SetupPublication(e,root+(grouped?"/grouped":"/interleaved"),argv[2],false);
        e.balancedCorrection=false;e.depthSharing=false;e.workQueuePolicy="on";
        e.base.insert(e.base.end(),{"--trace-cache-mib","256"});
        e.theta={170,171,172,173,176,177,179,180};
        const double start=Seconds();
        if(grouped)e.Prefetch({{count,0,8,0,e.theta}},e.seeds.data(),2);
        std::ostringstream proof;proof<<std::setprecision(17);
        for(size_t row=0;row<e.theta.size();row+=2){
            const auto jobs=e.DirectStage(count,e.seeds.data(),8,8,{e.theta[row],e.theta[row+1]},true);
            for(const auto& job:jobs)for(const auto& m:job.matrix)for(double v:m)proof<<v<<' ';
        }
        Require(e.nativeCommitted==32,"packet coalescing changed the native work budget");
        Atomic(e.root+"/proof.txt",proof.str());
        Atomic(e.root+"/timing.json","{\"seconds\":"+Number(Seconds()-start)+",\"worker_seconds\":"+Number(e.nativeWorkerSeconds)+",\"jobs\":"+Number(e.nativeCommitted)+"}");
    }
    Require(Read(root+"/interleaved/proof.txt")==Read(root+"/grouped/proof.txt"),"packet ordering changed a Mueller value");
    std::cout<<"PASS identical all16, sampling and native job count\n";return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
