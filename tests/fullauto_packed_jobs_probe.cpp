#define main unused_lazy_restore_main
#include "fullauto_lazy_restore_probe.cpp"
#undef main
int main(int argc,char** argv){try{
    if(argc==4&&std::string(argv[3])=="timing"){
        Samples reference;double raw=0,packed=0;
        for(int repeat=0;repeat<2;++repeat)for(bool packing:{false,true}){
            ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,std::string(argv[1])+"/trial"+Number(repeat)+"_packed"+Number(packing),argv[2],false);
            e.theta={165,170};e.phi=4;e.coherent=false;e.base.push_back("--incoherent");e.ConfigureIndex(true);e.packJobs=packing;e.workerPolicy="persistent";e.workQueuePolicy="off";
            Samples values;const double begin=Seconds();
            for(int chunk=0;chunk<8;++chunk){const auto jobs=e.Stage(512,e.seeds.data(),8,6,e.theta,true,false,chunk*512);const auto data=e.Data(jobs);values.insert(values.end(),data.begin(),data.end());}
            const double elapsed=Seconds()-begin;if(reference.empty())reference=values;else Require(values==reference,"packing performance fixture changed all16");
            (packing?packed:raw)+=elapsed;std::cout<<"TIMING packed="<<packing<<" repeat="<<repeat<<" native_jobs="<<e.cache.size()<<" seconds="<<elapsed<<" inode_release="<<e.packedInodesReleased<<'\n';
        }
        std::cout<<"TIMING packed_over_raw="<<packed/raw<<" raw_seconds="<<raw<<" packed_seconds="<<packed<<" all16_exact=true\n";return 0;
    }
    Require(argc==3,"OUTPUT BINARY");
    for(bool controls:{false,true})for(bool coherent:{false,true}){
        ArgPP args;RunConfig config;Engine e(args,config);SetupMemory(e,std::string(argv[1])+"/controls"+Number(controls)+"_coherent"+Number(coherent),argv[2],controls);
        e.coherent=coherent;if(!coherent)e.base.push_back("--incoherent");e.ConfigureIndex(true);e.workerPolicy="persistent";e.workQueuePolicy="off";e.initial=128;
        const auto stage=e.SharedStage(128,e.seeds.data(),8,6,4,e.theta);const auto high=e.Data(stage.first),low=e.Data(stage.second);
        std::vector<IntegralScalar> original;for(const auto& job:stage.first){IntegralReader reader;original.push_back(reader.Resolve(job.folder,6));}
        size_t released=0;for(const auto& item:e.cache)released+=PackedJobStore::Pack(item.second.folder);
        Require(released>e.cache.size()*3,"packing did not release job inodes");
        DIR* directory=opendir((e.root+"/jobs").c_str());size_t remaining=0;while(auto* entry=readdir(directory))if(entry->d_name[0]!='.')++remaining;closedir(directory);Require(remaining==0,"completed physical jobs retained");
        // Simulate a crash after durable pack commit, before the discovery
        // journal append. Full records must be found without retracing.
        Atomic(e.root+"/completed_jobs.journal","");
        setenv("MBS_FULLAUTO_PACK_JOBS","1",1);
        Engine resumed(args,config);SetupMemory(resumed,e.root,argv[2],controls);resumed.coherent=coherent;if(!coherent)resumed.base.push_back("--incoherent");resumed.ConfigureIndex(false);unsetenv("MBS_FULLAUTO_PACK_JOBS");
        resumed.initial=128;resumed.payloads.SetLimit(1);resumed.forbidNativeWork=true;resumed.WarmCache();Require(resumed.cache.size()==e.cache.size(),"packed recovery lost jobs");
        const auto recovered=resumed.SharedStage(128,resumed.seeds.data(),8,6,4,resumed.theta);
        Require(resumed.Data(recovered.first)==high&&resumed.Data(recovered.second)==low,"packed resume changed any matrix element");
        for(size_t seed=0;seed<8;++seed){IntegralReader reader;const auto value=reader.Resolve(recovered.first[seed].folder,6);Require(value.area==original[seed].area&&value.ext==original[seed].ext,"packed integral changed");}
        // A new derived angular view is automatically packed too. It may read
        // native and derived sources out of several existing archive shards.
        const std::vector<double> subset={0,170,180};const auto view=resumed.SharedStage(128,resumed.seeds.data(),8,6,4,subset);
        const auto values=resumed.Data(view.first);for(size_t seed=0;seed<8;++seed)for(size_t row=0;row<3;++row)for(size_t c=0;c<16;++c)
            Require(std::fabs(values[seed][row][c]-high[seed][row*2][c])<=1e-12*std::fabs(high[seed][row*2][0]),"new packed view changed matrix");
        Require(resumed.packedInodesReleased>0,"new derived jobs not packed");
        std::cout<<"PASS packed real PO: coherent="<<coherent<<" controls="<<controls<<" cached_jobs="<<e.cache.size()<<" released_inodes="<<released<<" all16/paired-depth/integrals exact; empty-journal resume with no native work; new packed views\n";
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
