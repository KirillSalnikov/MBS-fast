// Verify that a recovery-only change of planned counts preserves the optical
// estimator definition, supports, saved matrices, phi and execution settings.
#include "../src/FullAuto.cpp"
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==4,"expected old checkpoint root, new checkpoint root, theta file");
    std::ifstream grid(argv[3]);std::vector<double> theta;double value;
    while(grid>>value)theta.push_back(value);Require(theta.size()>1,"invalid theta grid");
    ArgPP args;RunConfig config;Engine a(args,config),b(args,config);
    AdaptiveState x,y;
    for(auto pair:{std::make_pair(&a,&x),std::make_pair(&b,&y)}){
        auto& e=*pair.first;auto& s=*pair.second;e.root=pair.first==&a?argv[1]:argv[2];
        auto text=Read(e.root+"/adaptive_checkpoint.tsv");std::istringstream meta(text.substr(text.find('\n')+1));meta>>e.signature;
        e.theta=theta;e.minPhi=1;e.maxPhi=4096;e.maxCount=1<<29;
        s.counts=Budgets(3,std::vector<int>(theta.size(),256));s.previous.resize(theta.size());
        s.valid.assign(theta.size(),false);s.coarse.assign(theta.size(),false);
        Require(RestoreAdaptive(e,s,false),"checkpoint failed native parser");
    }
    Require(a.signature==b.signature&&x.round==y.round&&x.streak==y.streak&&x.auditSerial==y.auditSerial
        &&x.knots==y.knots&&x.previous==y.previous&&x.valid==y.valid&&x.coarse==y.coarse
        &&a.thetaPhi==b.thetaPhi&&a.phi==b.phi&&a.controls==b.controls&&a.weights==b.weights
        &&a.polarWeights==b.polarWeights&&a.weightTheta==b.weightTheta&&a.mirror==b.mirror
        &&a.kernel.block==b.kernel.block&&a.kernel.warp==b.kernel.warp&&a.kernel.pipeline==b.kernel.pipeline
        &&a.depthSharing==b.depthSharing&&a.depthVerified==b.depthVerified&&a.local==b.local,
        "recovery changed fields other than pending counts");
    Require(x.counts[1]==y.counts[1]&&x.counts[2]==y.counts[2],"correction counts changed");
    size_t changed=0;for(size_t id=0;id<theta.size();++id){Require(y.counts[0][id]<=x.counts[0][id],"recovery increased base work");changed+=x.counts[0][id]!=y.counts[0][id];}
    Require(changed>0,"no deferred base work tested");
    std::cout<<"native checkpoint parser: "<<changed<<" pending base targets deferred; all other state exact\n";return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
