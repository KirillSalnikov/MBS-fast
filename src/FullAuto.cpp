#include "FullAuto.h"
#include "ArgPP.h"
#include "CliOptions.h"
#include "RunConfig.h"
#include "cuda/GpuSupport.h"
#include "RuntimeInfo.h"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <cerrno>
#include <csignal>
#include <sys/file.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>

namespace FullAuto {
namespace {
const double t95 = 2.3646242510102993;
const unsigned production[8] = {11,23,37,53,71,89,107,131};
const unsigned training[8] = {1009,1013,1019,1021,1031,1033,1039,1049};
const unsigned validation[8] = {2003,2011,2017,2027,2029,2039,2053,2063};
const unsigned mirrorValidation[8] = {4003,4007,4013,4019,4021,4027,4049,4051};
const unsigned angularValidation[8] = {6007,6011,6029,6037,6043,6047,6053,6067};
volatile sig_atomic_t activeWorkers[8] = {};
void StopWorkers(int signal) {
    for(int i=0;i<8;++i)if(activeWorkers[i]>0)kill(activeWorkers[i],SIGTERM);
    _exit(128+signal);
}
struct WorkerCleanup {
    ~WorkerCleanup() {
        for(int i=0;i<8;++i)if(activeWorkers[i]>0)kill(activeWorkers[i],SIGTERM);
        for(int i=0;i<8;++i)if(activeWorkers[i]>0){int status;while(waitpid(activeWorkers[i],&status,0)<0&&errno==EINTR){}activeWorkers[i]=0;}
    }
};
void Require(bool ok,const std::string& message) { if(!ok) throw std::runtime_error("fullauto: "+message); }
double Seconds() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
int Power2(double n) {
    Require(std::isfinite(n) && n>=0,"nonfinite sample allocation");
    int v=1;while(v<n && v<(1<<30))v*=2;return v;
}
std::string Number(double x) { std::ostringstream s;s<<std::setprecision(17)<<x;return s.str(); }
std::string Hash(const std::string& s) {
    uint64_t h=14695981039346656037ULL;for(unsigned char c:s){h^=c;h*=1099511628211ULL;}
    std::ostringstream o;o<<std::hex<<h;return o.str();
}
std::string Read(const std::string& p) {
    std::ifstream f(p.c_str(),std::ios::binary);Require(bool(f),"cannot read "+p);
    return std::string(std::istreambuf_iterator<char>(f),std::istreambuf_iterator<char>());
}
bool Exists(const std::string& p) { struct stat s;return stat(p.c_str(),&s)==0; }
void Directory(const std::string& p) {
    if(p.empty() || p=="/" || Exists(p))return;
    size_t slash=p.find_last_of('/');if(slash!=std::string::npos)Directory(p.substr(0,slash));
    Require(mkdir(p.c_str(),0755)==0 || errno==EEXIST,"cannot create "+p);
}
void Atomic(const std::string& p,const std::string& text) {
    const std::string tmp=p+".tmp";std::ofstream f(tmp.c_str());Require(bool(f),"cannot write "+tmp);
    f<<text;f.close();Require(bool(f),"write failed "+tmp);Require(rename(tmp.c_str(),p.c_str())==0,"rename failed "+p);
}
std::string Quoted(const std::string& x) {
    std::ostringstream o;o<<'"';for(unsigned char c:x){if(c=='"'||c=='\\')o<<'\\'<<c;
        else if(c=='\n')o<<"\\n";else if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else o<<c;}o<<'"';return o.str();
}
struct Lock {
    int fd;explicit Lock(const std::string& p):fd(open(p.c_str(),O_CREAT|O_RDWR|O_CLOEXEC,0644)) {
        if(fd<0)throw std::runtime_error("fullauto: cannot open controller lock");
        if(flock(fd,LOCK_EX|LOCK_NB)!=0){close(fd);throw std::runtime_error("fullauto: another controller owns this output");}
    }
    ~Lock(){close(fd);}
};
Matrix Mean(const Samples& s) {
    Require(s.size()==8 && !s[0].empty(),"need eight nonempty independent scrambles");
    Matrix m(s[0].size());for(const auto& a:s){Require(a.size()==m.size(),"angular row mismatch");
        for(size_t i=0;i<m.size();++i)for(int j=0;j<16;++j){Require(std::isfinite(a[i][j]),"nonfinite matrix");m[i][j]+=a[i][j]/8.;}}
    return m;
}
Matrix Variance(const Samples& s) {
    Matrix m=Mean(s),v(m.size());for(const auto& a:s)for(size_t i=0;i<m.size();++i)for(int j=0;j<16;++j)v[i][j]+=std::pow(a[i][j]-m[i][j],2)/7.;return v;
}
Samples Subtract(Samples a,const Samples& b) {
    Require(a.size()==b.size(),"seed mismatch");for(size_t s=0;s<a.size();++s){Require(a[s].size()==b[s].size(),"row mismatch");
        for(size_t i=0;i<a[s].size();++i)for(int j=0;j<16;++j)a[s][i][j]-=b[s][i][j];}return a;
}
Samples Combine(const std::vector<Samples>& levels) {
    Require(!levels.empty(),"no multilevel samples");Samples a=levels[0];
    for(size_t l=1;l<levels.size();++l)for(size_t s=0;s<a.size();++s)for(size_t i=0;i<a[s].size();++i)
        for(int j=0;j<16;++j)a[s][i][j]+=levels[l][s][i][j];
    return a;
}
std::vector<double> ReadGrid(const std::string& p) {
    std::istringstream f(Read(p));std::vector<double> x;std::string line;
    while(std::getline(f,line)){size_t c=line.find('#');if(c!=std::string::npos)line.erase(c);
        std::replace(line.begin(),line.end(),',',' ');std::istringstream row(line);double v;
        while(row>>v){Require(std::isfinite(v)&&v>=0&&v<=180,"theta must be finite in [0,180]");x.push_back(v);}
        Require(row.eof(),"nonnumeric theta row");}
    Require(x.size()>=2,"theta interpolation needs at least two distinct angles");
    for(size_t i=1;i<x.size();++i)Require(x[i]>x[i-1],"theta must be strictly increasing");
    return x;
}
void WriteGrid(const std::string& p,const std::vector<double>& x) {
    std::ostringstream f;for(double v:x)f<<Number(v)<<'\n';Atomic(p,f.str());
}
void WriteMatrix(const std::string& p,const std::vector<double>& x,const Matrix& y) {
    std::ostringstream f;f<<"ScAngle 2pi*dcos M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44\n"<<std::setprecision(17);
    for(size_t i=0;i<x.size();++i){double lower=i?(.5*(x[i-1]+x[i])):x[i];
        double upper=i+1<x.size()?(.5*(x[i+1]+x[i])):x[i];
        f<<x[i]<<' '<<2*3.14159265358979323846*(std::cos(lower*3.14159265358979323846/180)-std::cos(upper*3.14159265358979323846/180));
        for(double v:y[i])f<<' '<<v;
        f<<'\n';}Atomic(p,f.str());
}
}

// Not-a-knot cubic spline. Its coefficients are linear in the observations,
// so the same interpolation can propagate each independent seed estimate.
Matrix CubicInterpolate(const std::vector<double>& x,const Matrix& y,const std::vector<double>& query) {
    const size_t n=x.size();Require(n>=2&&y.size()==n,"invalid spline dimensions");
    std::vector<double> h(n-1);for(size_t i=0;i+1<n;++i){h[i]=x[i+1]-x[i];Require(h[i]>0,"nonincreasing spline coordinates");}
    Matrix second(n);
    for(int j=0;j<16;++j){
        if(n==3){double v=2*((y[2][j]-y[1][j])/h[1]-(y[1][j]-y[0][j])/h[0])/(h[0]+h[1]);for(auto& r:second)r[j]=v;}
        if(n>=4){size_t k=n-2;std::vector<double> low(k),diag(k),up(k),rhs(k);
            for(size_t i=1;i+1<n;++i){size_t r=i-1;low[r]=h[i-1];diag[r]=2*(h[i-1]+h[i]);up[r]=h[i];
                rhs[r]=6*((y[i+1][j]-y[i][j])/h[i]-(y[i][j]-y[i-1][j])/h[i-1]);}
            diag[0]+=h[0]*(h[0]+h[1])/h[1];up[0]-=h[0]*h[0]/h[1];low[0]=0;
            diag[k-1]+=h[n-2]*(h[n-3]+h[n-2])/h[n-3];low[k-1]-=h[n-2]*h[n-2]/h[n-3];up[k-1]=0;
            for(size_t i=1;i<k;++i){double f=low[i]/diag[i-1];diag[i]-=f*up[i-1];rhs[i]-=f*rhs[i-1];}
            for(size_t ii=k;ii-->0;){second[ii+1][j]=(rhs[ii]-(ii+1<k?up[ii]*second[ii+2][j]:0))/diag[ii];}
            second[0][j]=((h[0]+h[1])*second[1][j]-h[0]*second[2][j])/h[1];
            second[n-1][j]=((h[n-3]+h[n-2])*second[n-2][j]-h[n-2]*second[n-3][j])/h[n-3];}
    }
    Matrix result(query.size());for(size_t r=0;r<query.size();++r){double q=query[r];Require(q>=x.front()-1e-12&&q<=x.back()+1e-12,"spline extrapolation forbidden");
        size_t i=std::upper_bound(x.begin(),x.end(),q)-x.begin();i=i?i-1:0;i=std::min(i,n-2);
        double a=(x[i+1]-q)/h[i],b=(q-x[i])/h[i];
        for(int j=0;j<16;++j)result[r][j]=a*y[i][j]+b*y[i+1][j]+((a*a*a-a)*second[i][j]+(b*b*b-b)*second[i+1][j])*h[i]*h[i]/6.;}
    return result;
}
Statistics Estimate(const Samples& s) {
    Statistics r;r.mean=Mean(s);Matrix v=Variance(s);r.interval.resize(v.size());
    for(size_t i=0;i<v.size();++i){Require(r.mean[i][0]>0,"mean M11 must be positive");
        for(int j=0;j<16;++j){r.interval[i][j]=t95*std::sqrt(v[i][j]/8)/r.mean[i][0];r.maxInterval=std::max(r.maxInterval,r.interval[i][j]);}}
    return r;
}
double Difference(const Matrix& a,const Matrix& b,const Matrix& reference) {
    Require(a.size()==b.size()&&a.size()==reference.size(),"difference dimensions differ");double e=0;
    for(size_t i=0;i<a.size();++i){Require(reference[i][0]>0,"invalid difference scale");for(int j=0;j<16;++j)e=std::max(e,std::fabs(a[i][j]-b[i][j])/reference[i][0]);}return e;
}
std::vector<int> Allocate(const std::vector<Samples>& levels,const std::vector<double>& costs,
        const std::vector<int>& counts,const Matrix& reference,double eps) {
    Require(levels.size()==costs.size()&&levels.size()==counts.size()&&eps>0,"invalid allocation dimensions");
    std::vector<Matrix> v;for(size_t l=0;l<levels.size();++l){Require(costs[l]>0&&std::isfinite(costs[l]),"invalid native cost");v.push_back(Variance(levels[l]));}
    std::vector<double> maximum(levels.size(),1.);
    for(size_t i=0;i<reference.size();++i){Require(reference[i][0]>0,"invalid allocation M11");
        for(int j=0;j<16;++j){double root=0;for(size_t l=0;l<v.size();++l)root+=std::sqrt(v[l][i][j]*counts[l]*costs[l]);
            double target=8*std::pow(eps*reference[i][0]/t95,2);
            for(size_t l=0;l<v.size();++l)maximum[l]=std::max(maximum[l],std::sqrt(v[l][i][j]*counts[l]/costs[l])*root/target);}}
    std::vector<int> result;for(double n:maximum)result.push_back(Power2(n));return result;
}

namespace {
struct Kernel {bool pipeline=false;int block=64;std::string warp="auto";};
struct Job {
    std::string folder,data,summary,identity,weightHash;
    Matrix matrix;std::vector<double> angles,inputAngles;std::vector<std::array<double,5>> controls;
    double seconds=0,controlRefinement=0,controlSetup=0;pid_t pid=-1;
    int count=0,depth=0,phi=0,block=0;unsigned seed=0;
    bool pipeline=false,mirrored=false,paired=false,assembled=false,hasControls=false;std::string warp;
};
void ValidateStoredJob(const std::string& folder,std::set<std::string>& visited){
    if(!visited.insert(folder).second)return;
    std::istringstream metadata(Read(folder+"/complete.txt"));std::string id,dataHash,summaryHash,sourceHash;double seconds=0;
    Require(bool(metadata>>id>>seconds>>dataHash>>summaryHash)&&seconds>0,"invalid cached metadata");metadata>>sourceHash;
    Require(Hash(Read(folder+"/result/result.dat"))==dataHash&&(summaryHash=="none"||Hash(Read(folder+"/result/result_analytic_facets.tsv"))==summaryHash),"cached job changed; do not reuse corrupted files");
    if(Exists(folder+"/angular_sources.tsv")){
        const std::string sources=Read(folder+"/angular_sources.tsv");Require(Hash(sources)==sourceHash,"angular cache provenance changed");
        std::istringstream rows(sources);std::string row;while(std::getline(rows,row)){
            size_t tab=row.find('\t');Require(tab!=std::string::npos,"invalid angular cache source");
            const std::string source=row.substr(0,tab);Require(Hash(Read(source+"/complete.txt"))==row.substr(tab+1),"angular cache source metadata changed");
            ValidateStoredJob(source,visited);
        }
    }
}
struct Engine {
    const ArgPP& args;const RunConfig& config;std::string exe,root,signature,binaryHash,weightFile;
    std::vector<std::string> base,gpus;std::vector<double> theta,weightTheta;
    std::vector<std::array<double,2>> weights;std::vector<unsigned> seeds;
    std::map<std::string,Job> cache;Kernel kernel;int phi=4,target=12,pilot=32768,kernelCount=8192;
    int initial=131072,minimum=8192,maxCount=67108864,roundLimit=32,startTheta=33;
    int maxTheta=2147483647,stablePasses=2,minPhi=1,maxPhi=2147483647;
    int workerThreads=4;
    double eps=.01,interpolation=.0025,statistical=.0075;bool controls=true;
    bool mirror=false;std::string mirrorPolicy="auto";double mirrorError=0;
    std::set<std::string> verifiedMirrorGrids;
    bool angularCache=true,sharedMeans=false;size_t assembledJobs=0,reusedRows=0,newRows=0;
    std::string controlPolicy="auto";int thetaAuditCount=0;
    Engine(const ArgPP& a,const RunConfig& c):args(a),config(c){}
    int Int(const std::string& key,int fallback){return args.IsCatched(key)?args.GetIntValue(key):fallback;}
    void Init(int argc,const char* argv[]) {
        Require(config.method==RunMethod::PhysicalOptics&&!config.useFft,"requires coherent direct PO");
        for(const char* key:{"incoh","analytic_mean_reference","haar_mirror_audit","multigrid","multikeq","multikeq_list","gpu_trace","orientfile_dedup","analytic_facet_samples","analytic_control_weights","analytic_backscatter","analytic_azimuth_gaussian","owen_avg","checkpoint"})
            Require(!args.IsCatched(key),std::string("unsupported modifier --")+key+"; run one coherent size without external controls");
        const RuntimeEnvironmentSnapshot environment=QueryRuntimeEnvironmentSnapshot();
        Require(environment.experimental.empty() || args.IsCatched("allow_experimental_environment"),"unset experimental MBS_* variables or acknowledge --allow-experimental-environment");
        for(const char* name:{"OMPI_COMM_WORLD_SIZE","PMI_SIZE","MV2_COMM_WORLD_SIZE"}){const char* v=getenv(name);Require(!v||std::atoi(v)<=1,"use one MPI rank; controller launches its own GPU workers");}
        char path[4096];ssize_t n=readlink("/proc/self/exe",path,sizeof(path)-1);Require(n>0,"cannot locate executable");path[n]=0;exe=path;binaryHash=Hash(Read(exe));
        root=args.IsCatched("o")?args.GetStringValue("o"):"fullauto_result";
        Require(root.find('%')==std::string::npos,"output placeholders are not supported; use an explicit directory");Directory(root);
        char* full=realpath(root.c_str(),nullptr);Require(full!=nullptr,"cannot resolve output directory");root=full;free(full);
        eps=args.GetDoubleValue("fullauto");interpolation=eps*.25;statistical=eps-interpolation;
        angularCache=!args.IsCatched("fullauto_angular_cache")||args.GetStringValue("fullauto_angular_cache")!="off";
        mirrorPolicy=args.IsCatched("mirror_gamma")?"on":(args.IsCatched("fullauto_mirror")?args.GetStringValue("fullauto_mirror"):"auto");
        pilot=Int("fullauto_pilot",32768);kernelCount=Int("fullauto_kernel",8192);initial=Int("fullauto_initial",131072);
        thetaAuditCount=Int("fullauto_theta_audit",(!args.IsCatched("grid")&&!args.IsCatched("tgrid"))?std::min(pilot,8192):0);
        Require(thetaAuditCount>=0&&thetaAuditCount<=2147483646,"invalid theta audit count");
        minimum=Int("fullauto_correction",8192);maxCount=Int("maxorient",67108864);roundLimit=Int("fullauto_rounds",32);startTheta=Int("fullauto_theta_start",33);
        for(int v:{pilot,kernelCount,initial,minimum,maxCount})Require(v>0&&v<2147483647,"counts must fit positive signed32bit");
        Require(initial<=maxCount&&minimum<=maxCount&&roundLimit>0&&startTheta>=2,"invalid count/round/theta limits");
        if(config.adaptive.loadedFromFile){
            pilot=Int("fullauto_pilot",std::max(config.adaptive.minPilotOrientations,std::min(32768,config.adaptive.maxPilotOrientations)));
            initial=Int("fullauto_initial",config.adaptive.minOrientations);
            if(!args.IsCatched("maxorient")&&config.adaptive.maxOrientations>0)maxCount=config.adaptive.maxOrientations;
            roundLimit=Int("fullauto_rounds",config.adaptive.maxJointSweeps);
            startTheta=Int("fullauto_theta_start",config.adaptive.minThetaPoints);maxTheta=config.adaptive.maxThetaPoints;
            minPhi=config.adaptive.minPhiPoints;maxPhi=config.adaptive.maxPhiPoints;
        }
        if(args.IsCatched("max_theta_points"))maxTheta=config.adaptive.maxThetaPoints;
        if(args.IsCatched("max_phi_points"))maxPhi=config.adaptive.maxPhiPoints;
        stablePasses=std::max(2,config.adaptive.stablePasses);
        initial=Power2(initial);minimum=Power2(minimum);Require(initial<=maxCount&&minimum<=maxCount,"rounded counts exceed max-orientations");
        for(int v:{pilot,kernelCount,initial,minimum,maxCount})Require(v>0&&v<2147483647,"invalid configured counts");
        Require(roundLimit>0&&startTheta>=2&&maxTheta>=2&&minPhi>0&&maxPhi>=minPhi,"invalid configured limits");
        target=args.IsCatched("n")?config.maxReflections:12;Require(target>0,"target depth must be positive");
        controlPolicy=args.IsCatched("fullauto_controls")?args.GetStringValue("fullauto_controls"):(args.IsCatched("analytic_facet_average")?"on":"auto");
        controls=controlPolicy!="off"&&config.refractiveReal>=1.;
        Require(controlPolicy!="on"||controls,"analytic facet controls require refractive REAL>=1");
        seeds.assign(production,production+8);if(args.IsCatched("owen_seeds")){Require(args.GetArgNumber("owen_seeds")==8,"need exactly eight production seeds");
            for(int i=0;i<8;++i){int s=args.GetIntValue("owen_seeds",i);Require(s>=0,"negative seed");seeds[i]=s;}}
        std::set<unsigned> unique(seeds.begin(),seeds.end());Require(unique.size()==8,"production seeds must be distinct");
        for(unsigned s:seeds)for(int i=0;i<8;++i)Require(s!=training[i]&&s!=validation[i]&&s!=4001,"production/calibration seed overlap");
        for(unsigned s:seeds){Require(s!=3001&&s!=5003,"production/geometry/mean-provider seed overlap");for(unsigned v:mirrorValidation)Require(s!=v,"production/mirror verification seed overlap");for(unsigned v:angularValidation)Require(s!=v,"production/theta pilot seed overlap");}
        if(args.IsCatched("tgrid"))theta=ReadGrid(args.GetStringValue("tgrid"));
        else {double lo=0,hi=180;int intervals=720;
            if(args.IsCatched("fullauto_theta_range")){lo=args.GetDoubleValue("fullauto_theta_range",0);hi=args.GetDoubleValue("fullauto_theta_range",1);}
            if(!args.IsCatched("grid")){
                // A default 0.25-degree reference grid would miss narrow
                // forward lobes of large particles. Bound phase variation by
                // a geometric diameter before the interpolation search.
                double extent=0,volume=0;
                if(args.IsCatched("p")){
                    Require(args.GetIntValue("p")==1,"supply --scattering-grid or --theta-grid-file for this built-in geometry");
                    double height=args.GetDoubleValue("p",1),diameter=args.GetDoubleValue("p",2);
                    extent=std::hypot(height,diameter);volume=3*std::sqrt(3.)*diameter*diameter*height/8.;
                }else{
                    std::istringstream input(Read(args.GetStringValue("pf")));std::string line;int header=0;
                    std::array<double,3> loBox={{1e300,1e300,1e300}},hiBox={{-1e300,-1e300,-1e300}};
                    std::vector<std::array<double,3>> facet;size_t vertices=0;
                    auto finish=[&](){for(size_t i=1;i+1<facet.size();++i){const auto& a=facet[0];const auto& b=facet[i];const auto& c=facet[i+1];
                        volume+=(a[0]*(b[1]*c[2]-b[2]*c[1])+a[1]*(b[2]*c[0]-b[0]*c[2])+a[2]*(b[0]*c[1]-b[1]*c[0]))/6.;}facet.clear();};
                    while(std::getline(input,line)){size_t comment=line.find('#');if(comment!=std::string::npos)line.erase(comment);
                        if(line.find_first_not_of(" \t\r")==std::string::npos){if(header>=3)finish();continue;}
                        if(header<3){++header;continue;}std::istringstream row(line);std::array<double,3> p;std::string extra;
                        Require(bool(row>>p[0]>>p[1]>>p[2])&&!(row>>extra),"unrecognized particle coordinates; supply an explicit theta grid");
                        for(int j=0;j<3;++j){Require(std::isfinite(p[j]),"nonfinite particle coordinate");loBox[j]=std::min(loBox[j],p[j]);hiBox[j]=std::max(hiBox[j],p[j]);}facet.push_back(p);++vertices;
                    }finish();Require(vertices>=4,"cannot bound particle coordinates");
                    extent=std::hypot(std::hypot(hiBox[0]-loBox[0],hiBox[1]-loBox[1]),hiBox[2]-loBox[2]);volume=std::fabs(volume);
                    if(args.IsCatched("rs"))extent=std::sqrt(3.)*args.GetDoubleValue("rs");
                }
                if(args.IsCatched("k_eq")){Require(volume>0,"cannot determine k-eq scaling; supply an explicit theta grid");
                    double radius=config.wavelengthUm*args.GetDoubleValue("k_eq")/(2*3.14159265358979323846);
                    extent*=radius/std::cbrt(3*volume/(4*3.14159265358979323846));}
                Require(extent>0&&std::isfinite(extent),"invalid particle diameter bound");
                double required=std::ceil((hi-lo)*3.14159265358979323846/180*16*extent/config.wavelengthUm);
                Require(required<=1000000,"default reference grid exceeds one million points; supply a narrower explicit theta grid");
                intervals=std::max(32,int(required));
            }
            if(args.IsCatched("grid")){unsigned count=args.GetArgNumber("grid");lo=count==3?180-args.GetDoubleValue("grid"):args.GetDoubleValue("grid");hi=count==3?180:args.GetDoubleValue("grid",1);intervals=args.GetIntValue("grid",count-1);}
            Require(hi>lo&&intervals>0,"fullauto needs a nonzero theta range");for(int i=0;i<=intervals;++i)theta.push_back(lo+(hi-lo)*i/intervals);}
        // The generated reference grid already bounds the angular search.
        // Prepare its analytic means once, including large grids: per-row
        // storage is linear and the radial kernel depends on the maximum
        // angle, not on the row count. Every refined subset then uses the
        // same means, without repeated cold preparations or proxy jobs.
        sharedMeans=controls;
        if(config.useGpu){int count=VisibleGpuDeviceCount();Require(count>0,"no visible CUDA GPU");
            std::vector<std::string> visible;const char* env=getenv("CUDA_VISIBLE_DEVICES");if(env){std::istringstream f(env);std::string v;while(std::getline(f,v,','))visible.push_back(v);}
            std::string list=args.IsCatched("fullauto_gpus")?args.GetStringValue("fullauto_gpus"):"0";std::istringstream f(list);std::string v;std::set<int> used;
            while(std::getline(f,v,',')){size_t pos=0;int id=std::stoi(v,&pos);Require(pos==v.size()&&id>=0&&id<count&&used.insert(id).second,"invalid/duplicate visible GPU index");gpus.push_back(visible.empty()?v:visible.at(id));}
        Require(!gpus.empty(),"empty GPU list");}else{Require(!args.IsCatched("fullauto_gpus"),"GPU list requires CUDA backend");gpus.push_back("");}
        workerThreads=config.threads>0?config.threads:std::max(1,std::min(16,QueryRuntimeResourceSnapshot().physicalCores/int(gpus.size())));
        std::set<std::string> drop={"fullauto","autofull","o","n","grid","tgrid","nphi","haar_alpha","mirror_gamma","analytic_facet_average","analytic_shadow_control","analytic_mean_cache","analytic_control_weights","close","profile_phases","orientation_pipeline","maxorient","max_theta_points","max_phi_points","adaptive_config","owen_seeds","convergence_passes","threads"};
        std::ostringstream fingerprint;fingerprint<<binaryHash<<'\n';for(int i=1;i<argc;++i)fingerprint<<argv[i]<<'\n';
        // Canonical reconstruction avoids shell quoting and option aliases.
        std::set<std::string> seen;for(const auto& spec:GetCliOptionSpecs())if(args.IsCatched(spec.key)&&seen.insert(spec.key).second){
            if(spec.key=="pf"||spec.key=="tgrid"||spec.key=="adaptive_config")fingerprint<<Hash(Read(args.GetStringValue(spec.key)))<<'\n';
            if(drop.count(spec.key)||spec.key.compare(0,9,"fullauto_")==0)continue;
            base.push_back("--"+spec.canonical);for(unsigned i=0;i<args.GetArgNumber(spec.key);++i)base.push_back(args.GetStringValue(spec.key,i));}
        // Result-affecting environment forms part of resume/cache identity.
        fingerprint<<FormatRuntimeEnvironmentReport(environment);
        char host[256]={};Require(gethostname(host,sizeof(host)-1)==0,"cannot determine host identity");
        fingerprint<<"host="<<host<<"\nthreads="<<workerThreads<<'\n';
        const char* visible=getenv("CUDA_VISIBLE_DEVICES");fingerprint<<"visible="<<(visible?visible:"")<<'\n';
        for(const auto& gpu:gpus)fingerprint<<"gpu="<<gpu<<'\n';
        signature=Hash(fingerprint.str());
    }
    void State(const std::string& status,int round,const std::vector<int>& counts,size_t nodes,double ci,double ie,int streak) {
        std::ostringstream s;s<<std::setprecision(17)<<"{\n\"status\":"<<Quoted(status)<<",\"configuration_fnv1a64\":"<<Quoted(signature)
            <<",\"controller_pid\":"<<getpid()<<",\"round\":"<<round<<",\"counts_per_seed\":[";
        for(size_t i=0;i<counts.size();++i)s<<(i?",":"")<<counts[i];
        s<<"],\"seeds\":[";
        for(size_t i=0;i<seeds.size();++i)s<<(i?",":"")<<seeds[i];
        s<<"],\"phi\":"<<phi<<",\"target_depth\":"<<target<<",\"requested_theta\":"<<theta.size()<<",\"evaluated_theta\":"<<nodes
            <<",\"statistical_budget\":"<<statistical<<",\"interpolation_budget\":"<<interpolation<<",\"max_pointwise95_scaled_M11\":"<<ci
            <<",\"max_guard_residual95_scaled_M11\":"<<ie<<",\"pass_streak\":"<<streak<<",\"pipeline\":"<<(kernel.pipeline?"true":"false")
            <<",\"block\":"<<kernel.block<<",\"warp\":"<<Quoted(kernel.warp)
            <<",\"mirror_gamma\":"<<(mirror?"true":"false")<<",\"mirror_policy\":"<<Quoted(mirrorPolicy)
            <<",\"analytic_controls\":"<<(controls?"true":"false")<<",\"analytic_control_policy\":"<<Quoted(controlPolicy)<<",\"shared_analytic_reference\":"<<(sharedMeans?"true":"false")<<",\"angular_cache\":"<<(angularCache?"true":"false")<<",\"assembled_jobs_this_invocation\":"<<assembledJobs<<",\"reused_angular_rows_this_invocation\":"<<reusedRows<<",\"computed_new_angular_rows_this_invocation\":"<<newRows
            <<",\"mirror_budget\":"<<(mirror?eps/256:0)<<",\"max_mirror_residual95_scaled_M11\":"<<mirrorError
            <<",\"scope\":\"eight-seed pointwise intervals; interpolation guards on requested grid; mirror paired pilot when enabled; fixed target physics; no global simultaneous confidence claim\"}\n";
        Atomic(root+"/fullauto_status.json",s.str());std::cout<<"fullauto "<<status<<" round="<<round<<" theta="<<nodes<<" ci="<<ci<<" interpolation="<<ie<<" streak="<<streak<<std::endl;
    }
    Job Load(Job j) {
        std::istringstream f(Read(j.data));std::string line;std::getline(f,line);while(std::getline(f,line)){if(line.empty())continue;std::istringstream r(line);double thetaValue,area;Require(bool(r>>thetaValue>>area)&&std::isfinite(thetaValue)&&std::isfinite(area),"invalid worker angle");j.angles.push_back(thetaValue);Mueller m;
            for(double& v:m)Require(bool(r>>v)&&std::isfinite(v),"invalid worker matrix "+j.data);
            std::string extra;Require(!(r>>extra),"extra worker columns");j.matrix.push_back(m);}
        if(controls){std::istringstream s(Read(j.summary));std::getline(s,line);while(std::getline(s,line)){if(line.empty())continue;std::istringstream r(line);double angle;std::array<double,5> a;
            Require(bool(r>>angle),"invalid controls");for(double& v:a)Require(bool(r>>v)&&std::isfinite(v),"invalid control summary");j.controls.push_back(a);
            double residual,hybrid,refinement,setup,wr,ws;
            Require(bool(r>>residual>>hybrid>>refinement>>setup>>wr>>ws)&&std::isfinite(refinement)&&std::isfinite(setup),"invalid control diagnostic fields");
            j.controlRefinement=refinement;j.controlSetup=setup;
        }Require(j.controls.size()==j.matrix.size(),"incomplete controls");}
        const std::string log=Read(j.folder+"/stdout.log");for(const char* label:{"Hard tree-limit hits:","Orientations still incomplete:"}){
            size_t pos=log.find(label);Require(pos!=std::string::npos,"missing tree diagnostics");std::istringstream r(log.substr(pos+std::string(label).size()));int n=-1;r>>n;Require(n==0,"worker tree limit or incomplete orientation");}
        return j;
    }
    std::string WriteWeights(const std::vector<double>& x) {
        if(!controls||weights.empty())return "";
        std::ostringstream f;f<<"theta_deg reflection_weight shadow_weight\n"<<std::setprecision(17);
        for(double q:x){size_t k=std::upper_bound(weightTheta.begin(),weightTheta.end(),q)-weightTheta.begin();k=k?k-1:0;k=std::min(k,weightTheta.size()-2);
            double t=(q-weightTheta[k])/(weightTheta[k+1]-weightTheta[k]);f<<q;for(int j=0;j<2;++j)f<<' '<<((1-t)*weights[k][j]+t*weights[k+1][j]);f<<'\n';}
        std::string path=root+"/weights_"+Hash(f.str())+".tsv";Atomic(path,f.str());return path;
    }
    void SetSampling(Job& j,int n,unsigned seed,int depth,const std::vector<double>& x,bool paired){
        j.count=n;j.seed=seed;j.depth=depth;j.phi=phi;j.block=kernel.block;j.pipeline=kernel.pipeline;j.warp=kernel.warp;
        j.mirrored=mirror&&!paired;j.paired=paired;j.inputAngles=x;j.hasControls=controls;
    }
    bool SameSampling(const Job& j,int n,unsigned seed,int depth)const{
        return j.count==n&&j.seed==seed&&j.depth==depth&&j.phi==phi&&j.block==kernel.block
            &&j.pipeline==kernel.pipeline&&j.warp==kernel.warp&&j.mirrored==mirror&&!j.paired&&j.hasControls==controls;
    }
    std::vector<Job> DirectStage(int n,unsigned const* seedSet,size_t ns,int depth,const std::vector<double>& x,bool weighted=true,bool pairedMirror=false) {
        Require(Hash(Read(exe))==binaryHash,"binary changed during run");std::ostringstream angles;for(double q:x)angles<<Number(q)<<'\n';std::string gridHash=Hash(angles.str());
        std::string grid=root+"/theta_"+gridHash+".csv";Atomic(grid,angles.str());std::string wf=weighted?WriteWeights(x):"";
        std::vector<Job> result(ns);std::vector<std::vector<std::string>> commands(ns);std::vector<size_t> todo;
        for(size_t s=0;s<ns;++s){std::vector<std::string> command=base;
            for(const std::string& v:std::vector<std::string>{"--max-reflections",std::to_string(depth),"--sobol-seed",std::to_string(n),std::to_string(seedSet[s]),"--haar-alpha","--theta-grid-file",grid,"--phi-points",std::to_string(phi),"--threads",std::to_string(workerThreads),"--close"})command.push_back(v);
            if(pairedMirror)command.push_back("--haar-mirror-audit");else if(mirror)command.push_back("--mirror-gamma");
            if(controls){command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(sharedMeans?root+"/means_shared.cache":root+"/mean_"+gridHash+".cache");if(sharedMeans){command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");}}
            if(kernel.pipeline)command.push_back("--orientation-pipeline");
            if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}
            std::ostringstream identity;identity<<signature<<'\n'<<kernel.block<<'\n'<<kernel.warp<<'\n';for(const auto& v:command)identity<<v<<'\n';if(!wf.empty())identity<<Read(wf);
            Job j;SetSampling(j,n,seedSet[s],depth,x,pairedMirror);j.weightHash=Hash(wf.empty()?"":Read(wf));j.identity=Hash(identity.str());j.folder=root+"/jobs/j"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";Directory(j.folder);
            auto correctAngles=[&](const Job& job){Require(job.angles.size()==x.size(),"wrong theta count");for(size_t r=0;r<x.size();++r)Require(std::fabs(job.angles[r]-x[r])<1e-10,"worker theta grid differs");};
            if(cache.count(j.identity)){result[s]=cache.at(j.identity);correctAngles(result[s]);continue;}
            const std::string meta=j.folder+"/complete.txt";
            if(Exists(meta)){std::istringstream m(Read(meta));std::string id,dataHash,summaryHash;m>>id>>j.seconds>>dataHash>>summaryHash;
                Require(id==j.identity&&j.seconds>0&&std::isfinite(j.seconds)&&Hash(Read(j.data))==dataHash&&(!controls||Hash(Read(j.summary))==summaryHash),"cached job changed; do not reuse corrupted files");
                j=Load(j);correctAngles(j);cache[j.identity]=j;result[s]=j;continue;}
            if(Exists(j.folder+"/worker.pid")){int pid=std::atoi(Read(j.folder+"/worker.pid").c_str());Require(pid<=0||kill(pid,0)!=0,"an incomplete cached worker is still running; avoid duplicate launch");}
            if(Exists(j.folder+"/result")){std::string saved=j.folder+"/result_incomplete_"+std::to_string(getpid())+"_"+std::to_string(s);
                Require(rename((j.folder+"/result").c_str(),saved.c_str())==0,"cannot preserve incomplete worker output");}
            command.push_back("--output");command.push_back(j.folder+"/result");commands[s]=command;result[s]=j;todo.push_back(s);
        }
        WorkerCleanup cleanup;
        const size_t slots=std::min(gpus.size(),size_t(8));
        std::vector<int> slotSeed(slots,-1);std::vector<double> begun(ns);
        size_t next=0,completed=0;
        auto launch=[&](size_t slot,size_t s){Job& j=result[s];begun[s]=Seconds();std::string text=exe+"\n";for(const auto& v:commands[s])text+=v+"\n";Atomic(j.folder+"/command.args",text);
                std::vector<char*> argv;argv.push_back(const_cast<char*>(exe.c_str()));for(auto& v:commands[s])argv.push_back(const_cast<char*>(v.c_str()));argv.push_back(nullptr);
                j.pid=fork();Require(j.pid>=0,"cannot fork worker");
                if(j.pid==0){signal(SIGTERM,SIG_DFL);signal(SIGINT,SIG_DFL);int log=open((j.folder+"/stdout.log").c_str(),O_WRONLY|O_CREAT|O_TRUNC,0644);if(log<0)_exit(126);dup2(log,1);dup2(log,2);close(log);
                    if(config.useGpu){setenv("CUDA_VISIBLE_DEVICES",gpus[slot].c_str(),1);unsetenv("MBS_GPU_DEVICE");}
                    setenv("MBS_GPU_BLOCK",std::to_string(kernel.block).c_str(),1);
                    if(kernel.warp=="auto")unsetenv("MBS_GPU_WARP_BEAMS");else setenv("MBS_GPU_WARP_BEAMS",kernel.warp=="warp"?"1":"0",1);
                    execv(exe.c_str(),argv.data());_exit(127);}
                activeWorkers[slot]=j.pid;slotSeed[slot]=int(s);Atomic(j.folder+"/worker.pid",std::to_string(j.pid)+"\n");};
        for(size_t slot=0;slot<slots&&next<todo.size();++slot)launch(slot,todo[next++]);
        while(completed<todo.size()){
                int status=0;pid_t waited;do{waited=waitpid(-1,&status,0);}while(waited<0&&errno==EINTR);
                Require(waited>0,"cannot reap worker");
                size_t slot=0;while(slot<slots&&activeWorkers[slot]!=waited)++slot;Require(slot<slots,"unknown worker PID");
                size_t s=size_t(slotSeed[slot]);Job& j=result[s];j.pid=-1;activeWorkers[slot]=0;slotSeed[slot]=-1;j.seconds=Seconds()-begun[s];std::remove((j.folder+"/worker.pid").c_str());
                Require(WIFEXITED(status)&&WEXITSTATUS(status)==0,"worker failed; see "+j.folder+"/stdout.log");
                j=Load(j);Require(j.matrix.size()==x.size(),"worker returned wrong theta count");for(size_t r=0;r<x.size();++r)Require(std::fabs(j.angles[r]-x[r])<1e-10,"worker returned wrong theta coordinates");
                Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+Hash(Read(j.data))+" "+(controls?Hash(Read(j.summary)):"none")+"\n");cache[j.identity]=j;++completed;
                if(next<todo.size())launch(slot,todo[next++]);
        }
        return result;
    }
    double ForecastFullGridCost(int n,unsigned seed,int depth,size_t points,double setup,double conservative)const{
        double sx=0,sy=0,sxx=0,sxy=0;size_t count=0;
        for(const auto& entry:cache){const Job& j=entry.second;
            if(j.assembled||!SameSampling(j,n,seed,depth))continue;
            const double x=j.inputAngles.size(),y=std::max(1e-9,j.seconds-j.controlSetup);
            sx+=x;sy+=y;sxx+=x*x;sxy+=x*y;++count;
        }
        const double denominator=count*sxx-sx*sx;
        if(count<2||denominator<=0)return conservative;
        const double slope=std::max(0.,(count*sxy-sx*sy)/denominator);
        const double intercept=std::max(0.,(sy-slope*sx)/count);
        return std::max(1e-9,std::min(conservative,setup+intercept+slope*points));
    }
    std::vector<Job> Stage(int n,unsigned const* seedSet,size_t ns,int depth,const std::vector<double>& x,bool weighted=true,bool pairedMirror=false){
        if(!angularCache||pairedMirror)return DirectStage(n,seedSet,ns,depth,x,weighted,pairedMirror);
        Require(Hash(Read(exe))==binaryHash,"binary changed during run");
        const std::string wf=weighted?WriteWeights(x):"",weightHash=Hash(wf.empty()?"":Read(wf));
        typedef std::pair<const Job*,size_t> SourceRow;std::vector<std::map<double,SourceRow>> rows(ns);
        std::set<double> requested(x.begin(),x.end());bool everySeedHasRows=true,exact=true;
        for(size_t s=0;s<ns;++s){bool exactSeed=false;
            for(const auto& cached:cache){const Job& j=cached.second;if(!SameSampling(j,n,seedSet[s],depth))continue;
                if(!j.assembled&&j.inputAngles==x&&j.weightHash==weightHash)exactSeed=true;
                for(size_t r=0;r<j.inputAngles.size();++r)if(requested.count(j.inputAngles[r]))rows[s].insert({j.inputAngles[r],{&j,r}});
            }
            exact=exact&&exactSeed;everySeedHasRows=everySeedHasRows&&!rows[s].empty();
        }
        // Keep all eight GPU workers available when a complete new sample
        // budget is needed. Only the theta-refinement path assembles rows.
        if(!everySeedHasRows||exact)return DirectStage(n,seedSet,ns,depth,x,weighted,pairedMirror);
        std::ostringstream coordinates;for(double q:x)coordinates<<Number(q)<<'\n';const std::string gridHash=Hash(coordinates.str());
        const std::string grid=root+"/theta_"+gridHash+".csv";Atomic(grid,coordinates.str());
        std::vector<Job> result(ns);bool allCached=true;
        for(size_t s=0;s<ns;++s){Job j;SetSampling(j,n,seedSet[s],depth,x,false);j.weightHash=weightHash;j.assembled=true;
            j.identity=Hash(signature+"\nangular assembly\n"+Number(n)+" "+Number(seedSet[s])+" "+Number(depth)+" "+Number(phi)+" "+Number(kernel.block)+" "+Number(kernel.pipeline)+" "+kernel.warp+" "+Number(mirror)+" "+Number(controls)+" "+Number(sharedMeans)+"\n"+gridHash+"\n"+weightHash);
            j.folder=root+"/jobs/a"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";Directory(j.folder);
            if(cache.count(j.identity)){result[s]=cache.at(j.identity);continue;}
            if(Exists(j.folder+"/complete.txt")){
                std::set<std::string> visited;ValidateStoredJob(j.folder,visited);
                std::istringstream meta(Read(j.folder+"/complete.txt"));std::string id;meta>>id>>j.seconds;Require(id==j.identity,"angular cache identity changed");
                j=Load(j);Require(j.matrix.size()==x.size(),"wrong assembled theta count");cache[j.identity]=j;result[s]=j;
            }else{allCached=false;result[s]=j;}
        }
        if(allCached)return result;
        std::vector<double> missing;
        for(double q:x){bool found=true;for(const auto& r:rows)found=found&&r.count(q);if(!found)missing.push_back(q);}
        // Native theta files need two endpoints. A single missing point uses
        // one existing row as an anchor; its cached estimate remains valid.
        if(missing.size()==1){for(double q:x)if(q!=missing[0]){missing.push_back(q);break;}std::sort(missing.begin(),missing.end());}
        std::vector<Job> fresh;if(!missing.empty())fresh=DirectStage(n,seedSet,ns,depth,missing,weighted);
        Job means;if(controls&&!sharedMeans){unsigned seed=5003;means=DirectStage(1,&seed,1,1,x,false)[0];}
        std::vector<std::array<double,2>> coefficients(x.size(),{{1.,1.}});
        if(!wf.empty()){std::istringstream file(Read(wf));std::string header;std::getline(file,header);
            for(size_t r=0;r<x.size();++r){double angle;Require(bool(file>>angle>>coefficients[r][0]>>coefficients[r][1]),"invalid assembled control weights");}}
        for(size_t s=0;s<ns;++s){Job& j=result[s];if(cache.count(j.identity))continue;
            std::set<std::string> sources;size_t reused=0;
            j.angles=x;j.matrix.resize(x.size());if(controls)j.controls.resize(x.size());
            for(size_t r=0;r<x.size();++r){const Job* source;size_t index;auto existing=rows[s].find(x[r]);
                if(existing!=rows[s].end()){source=existing->second.first;index=existing->second.second;++reused;}
                else{source=&fresh[s];index=std::lower_bound(missing.begin(),missing.end(),x[r])-missing.begin();}
                sources.insert(source->folder);j.matrix[r]=source->matrix.at(index);
                if(controls){j.controls[r]=source->controls.at(index);
                    if(!sharedMeans){j.controls[r][1]=means.controls[r][1];j.controls[r][2]=means.controls[r][2];}
                    else{j.controlRefinement=std::max(j.controlRefinement,source->controlRefinement);j.controlSetup=source->controlSetup;}
                    const auto& c=j.controls[r];const auto& b=coefficients[r];j.matrix[r][0]=c[0]+b[0]*(c[1]-c[3])+b[1]*(c[2]-c[4]);}
            }
            if(controls&&!sharedMeans){sources.insert(means.folder);j.controlRefinement=means.controlRefinement;j.controlSetup=means.controlSetup;}
            std::set<std::string> visited;std::ostringstream dependencies;
            for(const std::string& source:sources){ValidateStoredJob(source,visited);dependencies<<source<<'\t'<<Hash(Read(source+"/complete.txt"))<<'\n';
                for(const auto& c:cache)if(c.second.folder==source){j.seconds+=c.second.seconds;break;}}
            Require(j.seconds>0,"missing angular source cost");const double representedSeconds=j.seconds;
            j.seconds=ForecastFullGridCost(n,seedSet[s],depth,x.size(),j.controlSetup,j.seconds);Directory(j.folder+"/result");WriteMatrix(j.data,x,j.matrix);
            if(controls){std::ostringstream table;table<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n"<<std::setprecision(17);
                for(size_t r=0;r<x.size();++r){const auto& c=j.controls[r];const auto& b=coefficients[r];table<<x[r];for(double v:c)table<<'\t'<<v;
                    table<<'\t'<<c[0]-b[0]*c[3]-b[1]*c[4]<<'\t'<<j.matrix[r][0]<<'\t'<<j.controlRefinement<<'\t'<<j.controlSetup<<'\t'<<b[0]<<'\t'<<b[1]<<'\n';}Atomic(j.summary,table.str());}
            std::vector<std::string> command=base;
            for(const std::string& v:std::vector<std::string>{"--max-reflections",std::to_string(depth),"--sobol-seed",std::to_string(n),std::to_string(seedSet[s]),"--haar-alpha","--theta-grid-file",grid,"--phi-points",std::to_string(phi),"--threads",std::to_string(workerThreads),"--close"})command.push_back(v);
            if(mirror)command.push_back("--mirror-gamma");
            if(kernel.pipeline)command.push_back("--orientation-pipeline");
            if(controls){command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(sharedMeans?root+"/means_shared.cache":root+"/mean_"+gridHash+".cache");if(sharedMeans){command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");}}
            if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}command.push_back("--output");command.push_back(j.folder+"/result");
            std::string recipe=exe+"\n";for(const auto& v:command)recipe+=v+"\n";Atomic(j.folder+"/command.args",recipe);
            Atomic(j.folder+"/stdout.log","Controller angular assembly; command.args is a reproducible full-grid recipe, not a launched worker.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
            Atomic(j.folder+"/angular_sources.tsv",dependencies.str());
            std::ostringstream provenance;provenance<<"{\"kind\":\"angular_assembly\",\"requested_rows\":"<<x.size()<<",\"reused_rows\":"<<reused<<",\"new_worker_rows\":"<<missing.size()<<",\"cost_scope\":\"full-grid forecast from matching native jobs; not assembly wall time\",\"represented_native_seconds\":"<<Number(representedSeconds)<<",\"forecast_full_grid_seconds\":"<<Number(j.seconds)<<",\"M11_rebased_to_authoritative_full_grid_means\":"<<(controls?"true":"false")<<"}\n";
            Atomic(j.folder+"/angular_reuse.json",provenance.str());
            Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+Hash(Read(j.data))+" "+(controls?Hash(Read(j.summary)):"none")+" "+Hash(dependencies.str())+"\n");
            cache[j.identity]=j;++assembledJobs;reusedRows+=reused;newRows+=missing.size();
            std::cout<<"fullauto angular cache: seed="<<seedSet[s]<<" N="<<n<<" depth="<<depth<<" reused="<<reused<<" new="<<missing.size()<<std::endl;
        }
        return result;
    }
    bool DisableUnsupportedControls(const std::exception& error){
        if(!controls||controlPolicy!="auto")return false;
        const std::string message=error.what();size_t start=message.find("see ");
        const std::string log=start==std::string::npos?message:Read(message.substr(start+4));
        if(log.find("analytic facet kernel exceeds radial-table budget")==std::string::npos
           &&log.find("analytic facet mean did not converge")==std::string::npos
           &&log.find("analytic facet mean is nonfinite or negative")==std::string::npos
           &&log.find("mean M11 must be positive")==std::string::npos)return false;
        controls=false;sharedMeans=false;weights.clear();SetMirror(false);verifiedMirrorGrids.clear();
        Atomic(root+"/analytic_control_fallback.json","{\"enabled\":false,\"reason\":\"analytic mean budget/convergence limit; full sampled PO retained\"}\n");
        std::cout<<"fullauto: analytic mean unavailable; retaining complete sampled PO without facet controls"<<std::endl;return true;
    }
    Samples Data(const std::vector<Job>& jobs){Samples s;for(const auto& j:jobs)s.push_back(j.matrix);return s;}
    double Cost(const std::vector<Job>& j,int n){double s=0;for(const auto& v:j)s+=v.seconds;return s/j.size()/n;}
    void SetMirror(bool value){mirror=value;mirrorError=0;statistical=eps-interpolation-(mirror?eps/256:0);}
    void PrepareMirror(const std::vector<double>& x){
        if(mirrorPolicy=="off")return;
        SetMirror(true);unsigned seed=3001;
        try{Stage(16,&seed,1,std::max(1,target-4),x,false);}
        catch(const std::exception& error){
            const std::string message=error.what();size_t start=message.find("see ");
            bool geometryFailure=start!=std::string::npos&&Read(message.substr(start+4)).find("geometry has no verified body xz reflection")!=std::string::npos;
            if(mirrorPolicy=="on"||!geometryFailure)throw;
            SetMirror(false);Atomic(root+"/mirror_verification.json","{\"enabled\":false,\"geometry_verified\":false,\"reason\":\"no body xz reflection; full gamma selected\"}\n");
            std::cout<<"fullauto mirror: geometry rejected, using full gamma"<<std::endl;
        }
    }
    bool VerifyMirror(const std::vector<double>& x){
        if(!mirror)return true;
        std::ostringstream coordinates;for(double q:x)coordinates<<Number(q)<<'\n';const std::string key=Hash(coordinates.str());
        if(verifiedMirrorGrids.count(key))return true;
        const int n=std::min(pilot,256);double worst=0,mirrorSeconds=0,pairedSeconds=0;std::ostringstream records;bool first=true;
        std::vector<int> depths;for(int d:{std::max(1,target-4),std::max(1,target-2),target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
        auto finalPairs=Stage(n,mirrorValidation,8,target,x,true,true);Matrix reference=Mean(Data(finalPairs));
        bool positive=true;for(const auto& row:reference)positive=positive&&row[0]>0;
        if(!positive){Atomic(root+"/mirror_verification.json","{\"enabled\":false,\"geometry_verified\":true,\"paired_verified\":false,\"reason\":\"nonpositive corrected mirror pilot M11; full gamma selected\"}\n");
            if(mirrorPolicy=="on")throw std::runtime_error("fullauto: mirror pilot M11 is not positive; use --fullauto-mirror off");
            return false;}
        for(int depth:depths){auto reduced=Stage(n,mirrorValidation,8,depth,x);auto paired=depth==target?finalPairs:Stage(n,mirrorValidation,8,depth,x,true,true);
            Samples delta=Subtract(Data(reduced),Data(paired));Matrix mean=Mean(delta),variance=Variance(delta);double maximum=0;
            for(size_t r=0;r<x.size();++r)for(int j=0;j<16;++j)maximum=std::max(maximum,(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/reference[r][0]);
            worst=std::max(worst,maximum);for(const auto& j:reduced)mirrorSeconds+=j.seconds;for(const auto& j:paired)pairedSeconds+=j.seconds;
            records<<(first?"":",")<<"{\"depth\":"<<depth<<",\"max_residual_plus95_scaled_M11\":"<<Number(maximum)<<"}";first=false;
        }
        const double bound=(2*depths.size()-1)*worst;mirrorError=std::max(mirrorError,bound);
        bool pass=bound<=eps/256;std::ostringstream report;
        report<<"{\"enabled\":"<<(pass?"true":"false")<<",\"geometry_verified\":true,\"paired_verified\":"<<(pass?"true":"false")
              <<",\"base_count_per_seed\":"<<n<<",\"phi\":"<<phi<<",\"theta_count\":"<<x.size()
              <<",\"seeds\":[4003,4007,4013,4019,4021,4027,4049,4051],\"levels\":["<<records.str()<<"]"
              <<",\"mlmc_residual_bound_plus95_scaled_M11\":"<<Number(bound)<<",\"budget\":"<<Number(eps/256)
              <<",\"reduced_native_seconds\":"<<Number(mirrorSeconds)<<",\"paired_native_seconds\":"<<Number(pairedSeconds)<<"}\n";
        Atomic(root+"/mirror_verification_"+key+".json",report.str());Atomic(root+"/mirror_verification.json",report.str());
        std::cout<<"fullauto mirror: bound="<<bound<<" budget="<<eps/256<<" verified="<<pass<<std::endl;
        if(pass)verifiedMirrorGrids.insert(key);
        else if(mirrorPolicy=="on")throw std::runtime_error("fullauto: paired Haar mirror verification failed; use --fullauto-mirror off");
        return pass;
    }
    std::vector<size_t> AuditThetaSupport(std::vector<size_t> knots){
        if(!thetaAuditCount)return knots;
        std::cout<<"fullauto: independent dense theta pilot N="<<thetaAuditCount<<" rows="<<theta.size()<<std::endl;
        auto jobs=DirectStage(thetaAuditCount,angularValidation,8,target,theta);Samples data=Data(jobs);Matrix reference=Mean(data);
        for(const auto& r:reference)Require(r[0]>0,"mean M11 must be positive in theta pilot");
        double maximum=0;int refinements=0;
        for(;;){Samples support(8),predicted(8);std::vector<double> kx;
            for(size_t k:knots){kx.push_back(theta[k]);for(int seed=0;seed<8;++seed)support[seed].push_back(data[seed][k]);}
            for(int seed=0;seed<8;++seed)predicted[seed]=CubicInterpolate(kx,support[seed],theta);
            Samples delta=Subtract(predicted,data);Matrix mean=Mean(delta),variance=Variance(delta);std::vector<double> errors(theta.size());maximum=0;
            for(size_t r=0;r<theta.size();++r)for(int j=0;j<16;++j){double v=(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/reference[r][0];errors[r]=std::max(errors[r],v);maximum=std::max(maximum,v);}
            if(maximum<=interpolation*.5||knots.size()==theta.size())break;
            std::set<size_t> next(knots.begin(),knots.end());for(size_t i=1;i<knots.size();++i){size_t worst=knots[i-1];for(size_t r=knots[i-1]+1;r<knots[i];++r)if(errors[r]>errors[worst])worst=r;
                if(errors[worst]>interpolation*.5)next.insert(worst);}
            Require(next.size()>knots.size(),"theta audit cannot refine support");knots.assign(next.begin(),next.end());++refinements;
        }
        std::ostringstream report;report<<"{\"count_per_seed\":"<<thetaAuditCount<<",\"reference_rows\":"<<theta.size()<<",\"support_rows\":"<<knots.size()<<",\"refinements\":"<<refinements<<",\"max_paired_residual95_scaled_M11\":"<<Number(maximum)<<",\"pilot_budget\":"<<Number(interpolation*.5)<<",\"independent_production_seeds\":true}\n";
        Atomic(root+"/dense_theta_pilot.json",report.str());return knots;
    }
    void Calibrate(const std::vector<double>& x) {
        weightTheta=x;weights.clear();phi=4;Matrix ref;double fastest=std::numeric_limits<double>::max();Kernel winner;
        if(config.useGpu){std::vector<Kernel> candidates(5);candidates[1].pipeline=candidates[2].pipeline=candidates[3].pipeline=candidates[4].pipeline=true;candidates[2].block=128;candidates[3].block=256;candidates[4].warp="thread";
            for(const auto& k:candidates){kernel=k;unsigned s=4001;auto j=Stage(kernelCount,&s,1,target,x,false);if(ref.empty())ref=j[0].matrix;
                Matrix scale=ref;if(controls)for(size_t r=0;r<scale.size();++r)scale[r][0]=j[0].controls[r][0];
                double error=Difference(j[0].matrix,ref,scale);std::cout<<"fullauto kernel block="<<k.block<<" pipeline="<<k.pipeline<<" error="<<error<<" seconds="<<j[0].seconds<<std::endl;
                if(error<1e-7&&j[0].seconds<fastest){fastest=j[0].seconds;winner=k;}}
            Require(fastest<std::numeric_limits<double>::max(),"no validated GPU kernel");kernel=winner;}
        double best=std::numeric_limits<double>::max();int chosen=4;std::vector<std::array<double,2>> selected;std::ostringstream candidatesJson;candidatesJson<<'[';bool firstCandidate=true;
        std::vector<int> phiCandidates;for(int p:{1,2,4,8,16,32,64})if(p>=minPhi&&p<=maxPhi)phiCandidates.push_back(p);
        if(phiCandidates.empty())phiCandidates.push_back(std::max(1,minPhi));
        if(args.IsCatched("nphi"))phiCandidates={args.GetIntValue("nphi")};
        for(int p:phiCandidates)Require(p>0&&p>=minPhi&&p<=maxPhi,"phi is outside configured limits");
        for(int candidate:phiCandidates){phi=candidate;int depth=std::max(1,target-4);auto train=Stage(pilot,training,8,depth,x,false),valid=Stage(pilot,validation,8,depth,x,false);
            std::vector<std::array<double,2>> fitted(x.size(),{{1.,1.}});Samples adjusted=Data(valid);
            if(controls)for(size_t r=0;r<x.size();++r){double ym=0,cm[2]={0,0};for(int s=0;s<8;++s){ym+=train[s].controls[r][0]/8;for(int j=0;j<2;++j)cm[j]+=train[s].controls[r][j+3]/8;}
                double scale[2]={0,0};for(int s=0;s<8;++s)for(int j=0;j<2;++j)scale[j]+=std::pow(train[s].controls[r][j+3]-cm[j],2)/7;
                bool active[2];for(int j=0;j<2;++j){scale[j]=std::sqrt(scale[j]);active[j]=scale[j]>std::max(std::fabs(ym),1e-100)*1e-12;}
                double cov=0,cy[2]={0,0};for(int s=0;s<8;++s){double v[2]={0,0};for(int j=0;j<2;++j)if(active[j]){v[j]=(train[s].controls[r][j+3]-cm[j])/scale[j];cy[j]+=v[j]*(train[s].controls[r][0]-ym)/7;}cov+=v[0]*v[1]/7;}
                auto b=fitted[r];if(active[0]&&active[1]){double det=1.21-cov*cov;b[0]=(1.1*cy[0]-cov*cy[1])/det/scale[0];b[1]=(1.1*cy[1]-cov*cy[0])/det/scale[1];}
                else for(int j=0;j<2;++j)if(active[j])b[j]=cy[j]/1.1/scale[j];
                for(double& v:b)v=std::max(-4.,std::min(4.,v));
                double def[8],prop[8],dm=0,pm=0;for(int s=0;s<8;++s){const auto& c=valid[s].controls[r];def[s]=c[0]+c[1]-c[3]+c[2]-c[4];prop[s]=c[0]+b[0]*(c[1]-c[3])+b[1]*(c[2]-c[4]);dm+=def[s]/8;pm+=prop[s]/8;}
                double dv=0,pv=0;for(int s=0;s<8;++s){dv+=std::pow(def[s]-dm,2);pv+=std::pow(prop[s]-pm,2);}
                if(pv<.7*dv)fitted[r]=b;
                for(int s=0;s<8;++s)adjusted[s][r][0]=pv<.7*dv?prop[s]:def[s];}
            Statistics stat=Estimate(adjusted);double seconds=0;for(const auto& j:valid)seconds+=j.seconds;double score=seconds*stat.maxInterval*stat.maxInterval;
            candidatesJson<<(firstCandidate?"":",")<<"{\"phi\":"<<phi<<",\"native_seconds\":"<<Number(seconds)<<",\"score\":"<<Number(score)<<",\"validation95\":"<<Number(stat.maxInterval)<<"}";firstCandidate=false;
            std::cout<<"fullauto phi="<<phi<<" score="<<score<<" pilot95="<<stat.maxInterval<<std::endl;
            if(score<best){best=score;chosen=phi;selected=fitted;}}
        phi=chosen;if(controls)weights=selected;candidatesJson<<']';
        const std::string weightPath=WriteWeights(x);
        std::ostringstream c;c<<"{\"selected_phi\":"<<phi<<",\"block\":"<<kernel.block<<",\"pipeline\":"<<(kernel.pipeline?"true":"false")<<",\"warp\":"<<Quoted(kernel.warp)<<",\"pilot_count\":"<<pilot<<",\"score\":"<<Number(best)<<",\"weights\":"<<Quoted(weightPath)<<",\"phi_candidates\":"<<candidatesJson.str()<<",\"training_seeds\":[1009,1013,1019,1021,1031,1033,1039,1049],\"validation_seeds\":[2003,2011,2017,2027,2029,2039,2053,2063],\"independent_training_validation\":true,\"mirror_gamma\":"<<(mirror?"true":"false")<<",\"mirror_policy\":"<<Quoted(mirrorPolicy)<<"}\n";Atomic(root+"/calibration.json",c.str());
    }
};
std::vector<size_t> Guards(const std::vector<size_t>& knots,size_t size) {
    std::set<size_t> all(knots.begin(),knots.end());for(size_t i=1;i<knots.size();++i){size_t a=knots[i-1],b=knots[i];for(int k=1;k<=3;++k){size_t q=a+(b-a)*k/4;if(q>a&&q<b)all.insert(q);}}
    // Every unsampled interval is checked at interior positions. Once the
    // requested grid itself is exhausted the interpolation residual is zero.
    Require(!all.empty()&&*all.rbegin()<size,"invalid theta guard indices");return std::vector<size_t>(all.begin(),all.end());
}
std::vector<double> Angles(const std::vector<double>& theta,const std::vector<size_t>& ids){std::vector<double> x;for(size_t i:ids)x.push_back(theta.at(i));return x;}
Samples Support(const Samples& s,const std::vector<size_t>& evaluated,const std::vector<size_t>& knots) {
    Samples result(8);for(size_t k:knots){size_t row=std::lower_bound(evaluated.begin(),evaluated.end(),k)-evaluated.begin();Require(row<evaluated.size()&&evaluated[row]==k,"missing support node");
        for(int seed=0;seed<8;++seed)result[seed].push_back(s[seed][row]);}return result;
}
}

int Run(const ArgPP& args,const RunConfig& config,int argc,const char* argv[]) {
    Engine e(args,config);e.Init(argc,argv);Lock lock(e.root+"/controller.lock");
    const std::string identity=e.root+"/configuration.txt";
    if(Exists(identity))Require(Read(identity)==e.signature+"\n","resume inputs/binary/environment changed; choose a new output directory");else Atomic(identity,e.signature+"\n");
    WriteGrid(e.root+"/requested_theta.csv",e.theta);Directory(e.root+"/jobs");
    signal(SIGTERM,StopWorkers);signal(SIGINT,StopWorkers);
    std::vector<size_t> knots;size_t initialTheta=std::min(e.theta.size(),size_t(e.startTheta));for(size_t i=0;i<initialTheta;++i)knots.push_back(i*(e.theta.size()-1)/(initialTheta-1));
    std::vector<size_t> evaluated=Guards(knots,e.theta.size());Require(evaluated.size()<=size_t(e.maxTheta),"initial theta guards exceed max-theta-points");auto x=Angles(e.theta,evaluated);
    e.State("calibrating",0,{e.initial},evaluated.size(),0,0,0);
    auto calibrate=[&](){e.PrepareMirror(x);e.Calibrate(x);if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(x);}};
    try {try{calibrate();}catch(const std::exception& error){if(!e.DisableUnsupportedControls(error))throw;calibrate();}}
    catch(...){e.State("failed",0,{e.initial},evaluated.size(),0,0,0);throw;}
    if(e.thetaAuditCount){
        try{knots=e.AuditThetaSupport(knots);}catch(const std::exception& error){
            if(!e.DisableUnsupportedControls(error)){e.State("failed",0,{e.initial},evaluated.size(),0,0,0);throw;}
            calibrate();knots=e.AuditThetaSupport(knots);
        }
        evaluated=Guards(knots,e.theta.size());Require(evaluated.size()<=size_t(e.maxTheta),"audited theta guards exceed max-theta-points");
    }
    std::vector<int> depths;for(int d:{std::max(1,e.target-4),std::max(1,e.target-2),e.target})if(depths.empty()||depths.back()!=d)depths.push_back(d);
    std::vector<int> counts(depths.size(),e.minimum);counts[0]=e.initial;
    Matrix previous;Samples previousCoarse;int previousCoarseN=0,streak=0;bool coarseVerified=false;
    double lastConfidence=0,lastInterpolation=0;int lastRound=0;
    std::vector<int> lastCounts=counts;size_t lastEvaluated=evaluated.size();
    e.State("active",0,counts,evaluated.size(),0,0,0);
    try {
        for(int round=0;round<e.roundLimit;++round){
            try {
            x=Angles(e.theta,evaluated);
            if(!e.VerifyMirror(x)){e.SetMirror(false);e.Calibrate(x);previous.clear();previousCoarse.clear();coarseVerified=false;streak=0;}
            std::vector<Samples> levels;std::vector<double> costs;
            for(size_t l=0;l<depths.size();++l){auto upper=e.Stage(counts[l],e.seeds.data(),8,depths[l],x);Samples value=e.Data(upper);double cost=e.Cost(upper,counts[l]);
                if(l){auto lower=e.Stage(counts[l],e.seeds.data(),8,depths[l-1],x);value=Subtract(value,e.Data(lower));cost+=e.Cost(lower,counts[l]);}levels.push_back(value);costs.push_back(cost);}
            Samples data=Combine(levels);Statistics stat=Estimate(data);
            if(!previousCoarse.empty()&&counts[0]>previousCoarseN)coarseVerified=coarseVerified||Difference(Mean(levels[0]),Mean(previousCoarse),stat.mean)<=e.statistical;
            double change=previous.empty()?1:Difference(stat.mean,previous,stat.mean);
            if(stat.maxInterval<=e.statistical&&!coarseVerified&&counts[0]>=2*e.initial){auto jobs=e.Stage(counts[0]/2,e.seeds.data(),8,depths[0],x);auto predecessor=levels;predecessor[0]=e.Data(jobs);
                change=Difference(stat.mean,Mean(Combine(predecessor)),stat.mean);coarseVerified=change<=e.statistical;}
            auto sx=Angles(e.theta,knots);Samples support=Support(data,evaluated,knots),predicted(8),out(8);
            for(int s=0;s<8;++s){predicted[s]=CubicInterpolate(sx,support[s],x);out[s]=CubicInterpolate(sx,support[s],e.theta);}
            Samples delta=Subtract(predicted,data);Matrix dm=Mean(delta),dv=Variance(delta);double interpError=0;std::vector<size_t> bad;std::vector<double> guardErrors(x.size());
            for(size_t r=0;r<x.size();++r){double error=0;for(int j=0;j<16;++j)error=std::max(error,(std::fabs(dm[r][j])+t95*std::sqrt(dv[r][j]/8))/stat.mean[r][0]);
                guardErrors[r]=error;interpError=std::max(interpError,error);if(error>e.interpolation)bad.push_back(evaluated[r]);}
            bool positive=true;for(const auto& m:Mean(out))positive=positive&&m[0]>0;
            double outputCI=positive?Estimate(out).maxInterval:1.;
            bool pass=stat.maxInterval<=e.statistical&&outputCI<=e.statistical&&interpError<=e.interpolation&&change<=e.statistical&&coarseVerified;
            streak=pass?streak+1:0;e.State("active",round,counts,evaluated.size(),std::max(outputCI,stat.maxInterval),interpError,streak);
            lastConfidence=std::max(outputCI,stat.maxInterval);lastInterpolation=interpError;lastRound=round;lastCounts=counts;lastEvaluated=evaluated.size();
            std::ostringstream rec;rec<<std::setprecision(17)<<round<<' '<<counts[0]<<' '<<evaluated.size()<<' '<<knots.size()<<' '<<stat.maxInterval<<' '<<outputCI<<' '<<interpError<<' '<<change<<' '<<streak<<'\n';
            std::ofstream history((e.root+"/refinements.tsv").c_str(),std::ios::app);history<<rec.str();history.close();
            if(streak>=e.stablePasses){auto final=Estimate(out);WriteMatrix(e.root+"/mueller_fullauto.dat",e.theta,final.mean);WriteMatrix(e.root+"/mueller_evaluated.dat",x,stat.mean);WriteGrid(e.root+"/theta_support.csv",sx);WriteGrid(e.root+"/theta_evaluated.csv",x);
                std::ostringstream guards;guards<<"theta_deg,is_support,M11,estimated_guard_error_plus95_scaled_M11\n"<<std::setprecision(17);
                for(size_t i=0;i<x.size();++i)guards<<x[i]<<','<<std::binary_search(knots.begin(),knots.end(),evaluated[i])<<','<<stat.mean[i][0]<<','<<guardErrors[i]<<'\n';
                Atomic(e.root+"/theta_validation.csv",guards.str());
                std::ostringstream accuracy;accuracy<<"theta_deg";for(int j=0;j<16;++j)accuracy<<",estimated95_M"<<j/4+1<<j%4+1<<"_scaled_M11";accuracy<<'\n'<<std::setprecision(17);
                for(size_t i=0;i<e.theta.size();++i){accuracy<<e.theta[i];for(double v:final.interval[i])accuracy<<','<<v;accuracy<<'\n';}Atomic(e.root+"/all_mueller_confidence.csv",accuracy.str());
                e.State("converged",round,counts,evaluated.size(),std::max(outputCI,stat.maxInterval),interpError,streak);return 0;}
            // Resolve geometric interpolation before spending more on the
            // same unresolved sparse support. Validation uses paired samples.
            if(!bad.empty()||!positive){std::set<size_t> next(knots.begin(),knots.end());next.insert(bad.begin(),bad.end());
                if(!positive)next.insert(evaluated.begin(),evaluated.end());
                if(next.size()>knots.size()){knots.assign(next.begin(),next.end());evaluated=Guards(knots,e.theta.size());if(evaluated.size()>size_t(e.maxTheta))break;previous.clear();previousCoarse.clear();coarseVerified=false;streak=0;continue;}}
            previous=stat.mean;previousCoarse=levels[0];previousCoarseN=counts[0];std::vector<int> next=counts;
            if(std::max(outputCI,stat.maxInterval)>e.statistical){auto proposed=Allocate(levels,costs,counts,stat.mean,e.statistical*.85);
                for(size_t l=0;l<counts.size();++l)next[l]=std::max(counts[l],std::min(e.maxCount,std::min(counts[l]>e.maxCount/4?e.maxCount:4*counts[l],proposed[l])));}
            else {size_t selected=0;double largest=-1;
                if(coarseVerified&&levels.size()>1)for(size_t l=1;l<levels.size();++l){auto v=Variance(levels[l]);double score=0;for(size_t i=0;i<v.size();++i)for(int j=0;j<16;++j)score=std::max(score,v[i][j]/std::pow(stat.mean[i][0],2));if(score>largest){largest=score;selected=l;}}
                next[selected]=counts[selected]>e.maxCount/2?e.maxCount:2*counts[selected];}
            if(next==counts){bool grew=false;for(size_t l=0;l<counts.size();++l)if(counts[l]<e.maxCount){next[l]=std::min(e.maxCount,counts[l]>e.maxCount/2?e.maxCount:2*counts[l]);grew=true;}if(!grew)break;}
            counts=next;
            }catch(const std::exception& error){
                if(!e.DisableUnsupportedControls(error))throw;
                calibrate();previous.clear();previousCoarse.clear();coarseVerified=false;streak=0;
            }
        }
        e.State("requires_further_sampling",lastRound,lastCounts,lastEvaluated,lastConfidence,lastInterpolation,streak);return 3;
    } catch(...) {e.State("failed",lastRound,counts,evaluated.size(),lastConfidence,lastInterpolation,streak);throw;}
}
}
