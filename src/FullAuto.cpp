#include <memory>
#include "FullAuto.h"
#include "ArgPP.h"
#include "CliOptions.h"
#include "RunConfig.h"
#include "cuda/GpuSupport.h"
#include "RuntimeInfo.h"
#include "WorkerProtocol.h"
#include "FullAutoAuditPlan.h"
#include "IntegralCharacteristics.h"
#include <poll.h>
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
#include <dirent.h>

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
// Reuse file digests only within this controller process. Both modification
// and change times, inode and size must match; restart always rehashes files.
std::string FileHash(const std::string& path) {
    struct Entry {struct stat stamp;std::string hash;};
    static std::map<std::string,Entry> digests;
    auto same=[](const struct stat& a,const struct stat& b){return a.st_dev==b.st_dev&&a.st_ino==b.st_ino&&a.st_size==b.st_size
        &&a.st_mtim.tv_sec==b.st_mtim.tv_sec&&a.st_mtim.tv_nsec==b.st_mtim.tv_nsec
        &&a.st_ctim.tv_sec==b.st_ctim.tv_sec&&a.st_ctim.tv_nsec==b.st_ctim.tv_nsec;};
    struct stat before,after;Require(stat(path.c_str(),&before)==0,"cannot stat "+path);
    const auto found=digests.find(path);if(found!=digests.end()&&same(before,found->second.stamp))return found->second.hash;
    const std::string digest=Hash(Read(path));Require(stat(path.c_str(),&after)==0&&same(before,after),"file changed while hashing "+path);
    digests[path]=Entry{after,digest};return digest;
}
bool Exists(const std::string& p) { struct stat s;return stat(p.c_str(),&s)==0; }
// Explicit compatibility bridge for controller scheduling and CUDA scratch
// allocation changes whose native optical matrices are numerically equivalent.
// Normal resume still requires the exact executable hash. Never extend this
// release allowlist without numerical equivalence tests for the new change.
std::string ExecutableVersion(const std::string& executable){
    int fd[2];Require(pipe2(fd,O_CLOEXEC)==0,"cannot create version pipe");
    const pid_t child=fork();
    if(child<0){close(fd[0]);close(fd[1]);Require(false,"cannot start version check");}
    if(child==0){close(fd[0]);dup2(fd[1],STDOUT_FILENO);dup2(fd[1],STDERR_FILENO);close(fd[1]);
        execl(executable.c_str(),executable.c_str(),"--version",static_cast<char*>(nullptr));_exit(127);}
    close(fd[1]);std::string text;char buffer[1024];bool valid=true;
    for(;;){pollfd p={fd[0],POLLIN,0};int ready;do{ready=poll(&p,1,10000);}while(ready<0&&errno==EINTR);
        if(ready<=0){valid=false;break;}
        ssize_t n;do{n=read(fd[0],buffer,sizeof(buffer));}while(n<0&&errno==EINTR);
        if(n==0)break;if(n<0||text.size()+n>16384){valid=false;break;}text.append(buffer,n);}
    close(fd[0]);if(!valid)kill(child,SIGKILL);int status=0;while(waitpid(child,&status,0)<0&&errno==EINTR){}
    Require(valid&&WIFEXITED(status)&&WEXITSTATUS(status)==0,"cannot verify compatible executable "+executable);return text;
}
void Directory(const std::string& p) {
    if(p.empty() || p=="/" || Exists(p))return;
    size_t slash=p.find_last_of('/');if(slash!=std::string::npos)Directory(p.substr(0,slash));
    Require(mkdir(p.c_str(),0755)==0 || errno==EEXIST,"cannot create "+p);
}
void Atomic(const std::string& p,const std::string& text) {
    if(Exists(p)&&Read(p)==text)return;
    const std::string tmp=p+".tmp";std::ofstream f(tmp.c_str());Require(bool(f),"cannot write "+tmp);
    f<<text;f.close();Require(bool(f),"write failed "+tmp);
    int fd=open(tmp.c_str(),O_RDONLY);Require(fd>=0&&fsync(fd)==0,"checkpoint flush failed "+tmp);close(fd);
    Require(rename(tmp.c_str(),p.c_str())==0,"rename failed "+p);
    const size_t slash=p.find_last_of('/');const std::string parent=slash==std::string::npos?".":p.substr(0,slash);
    fd=open(parent.c_str(),O_RDONLY|O_DIRECTORY);if(fd>=0){fsync(fd);close(fd);}
}
std::string ProcessStart(pid_t pid){
    const std::string path="/proc/"+std::to_string(pid)+"/stat";if(!Exists(path))return "";
    const std::string stat=Read(path);const size_t end=stat.rfind(')');if(end==std::string::npos)return "";
    std::istringstream fields(stat.substr(end+2));std::string value;for(int i=0;i<=19;++i)if(!(fields>>value))return "";return value;
}
bool WorkerAlive(const std::string& path){
    std::istringstream file(Read(path));int pid=0;std::string start,boot;file>>pid>>start>>boot;
    if(pid<=0||kill(pid,0)!=0)return false;
    const std::string state=Read("/proc/"+std::to_string(pid)+"/stat");const auto stateEnd=state.rfind(')');
    if(stateEnd!=std::string::npos&&stateEnd+2<state.size()&&state[stateEnd+2]=='Z')return false;
    if(start.empty()||boot.empty())return true; // Legacy record: refuse ambiguous ownership.
    std::string currentBoot=Read("/proc/sys/kernel/random/boot_id");currentBoot.erase(currentBoot.find_last_not_of("\n\r ")+1);
    return boot==currentBoot&&start==ProcessStart(pid);
}
std::string Quoted(const std::string& x) {
    std::ostringstream o;o<<'"';for(unsigned char c:x){if(c=='"'||c=='\\')o<<'\\'<<c;
        else if(c=='\n')o<<"\\n";else if(c<32)o<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;else o<<c;}o<<'"';return o.str();
}
std::string CompactJson(const std::string& text){
    std::string result;result.reserve(text.size());bool quoted=false,escaped=false;
    for(char c:text){if(quoted){result+=c;if(escaped)escaped=false;else if(c=='\\')escaped=true;else if(c=='"')quoted=false;}
        else if(c=='"'){quoted=true;result+=c;}else if(c!=' '&&c!='\n'&&c!='\r'&&c!='\t')result+=c;}
    return result;
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
    Statistics r;r.mean=Mean(s);Matrix v=Variance(s);r.interval.resize(v.size());r.ratioInterval.resize(v.size());
    for(size_t i=0;i<v.size();++i){Require(r.mean[i][0]>0,"mean M11 must be positive");
        for(int j=0;j<16;++j){r.interval[i][j]=t95*std::sqrt(v[i][j]/8)/r.mean[i][0];r.maxInterval=std::max(r.maxInterval,r.interval[i][j]);
            if(!j)continue;
            double covariance=0;for(const auto& a:s)covariance+=(a[i][0]-r.mean[i][0])*(a[i][j]-r.mean[i][j])/7;
            // Fieller interval for the ratio of means, with paired covariance.
            const double q=t95*t95/8,m0=r.mean[i][0],mj=r.mean[i][j],ratio=mj/m0;
            const double aa=m0*m0-q*v[i][0];
            double error=1e100;
            if(aa>0){double transformed=0;for(const auto& a:s){const double z=(a[i][j]-mj)-ratio*(a[i][0]-m0);transformed+=z*z/7;}
                const double disc=std::max(0.,q*m0*m0*transformed-q*q*std::max(0.,v[i][0]*v[i][j]-covariance*covariance));
                const double centre=(m0*mj-q*covariance)/aa;error=std::fabs(centre-ratio)+std::sqrt(disc)/aa;}
            r.ratioInterval[i][j]=error;r.maxRatioInterval=std::max(r.maxRatioInterval,error);
        }}
    return r;
}
double RotationEquivalenceBound(const Samples& candidate,const Samples& reference,bool ratios){
    if(candidate.size()!=8||reference.size()!=8||reference[0].empty())return 1e100;
    double relative=0,largestRatio=0;
    for(size_t s=0;s<8;++s){
        if(candidate[s].size()!=reference[0].size()||reference[s].size()!=reference[0].size())return 1e100;
        for(size_t r=0;r<reference[s].size();++r){
            const auto& a=candidate[s][r];const auto& b=reference[s][r];
            for(int j=0;j<16;++j)if(!std::isfinite(a[j])||!std::isfinite(b[j]))return 1e100;
            if(!(a[0]>0&&b[0]>0)){
                // Identically dark pairs carry no intensity or ratio bias.
                bool dark=true;for(int j=0;j<16;++j)dark=dark&&a[j]==0&&b[j]==0;
                if(dark)continue;
                return 1e100;
            }
            for(int j=0;j<16;++j){
                relative=std::max(relative,std::fabs(a[j]-b[j])/b[0]);
                if(j)largestRatio=std::max(largestRatio,std::fabs(b[j])/b[0]);
            }
        }
    }
    // If |delta Mij| <= e*M11 in every paired block, the same bound holds
    // for any positive weighted sum of those blocks. Its ratio difference
    // is bounded by e*(1+max|Mij/M11|)/(1-e), even when the individual
    // orientation means have unbounded Fieller intervals. Geometry must
    // independently establish the symmetry; this is a numerical audit only.
    if(!std::isfinite(relative)||!std::isfinite(largestRatio)||relative>=1)return 1e100;
    return ratios?std::max(relative,relative*(1+largestRatio)/(1-relative)):relative;
}
Matrix LocalCubicInterpolate(const std::vector<double>& x,const Matrix& y,const std::vector<double>& query){
    const size_t n=x.size();Require(n>=2&&y.size()==n,"invalid local spline dimensions");for(size_t i=1;i<n;++i)Require(x[i]>x[i-1],"nonincreasing local spline coordinates");
    Matrix derivative(n);const size_t width=std::min(size_t(4),n);
    for(size_t i=0;i<n;++i){const size_t first=std::min(i?i-1:0,n-width);
        for(size_t j=first;j<first+width;++j){double coefficient;
            if(j==i){coefficient=0;for(size_t k=first;k<first+width;++k)if(k!=i)coefficient+=1/(x[i]-x[k]);}
            else{coefficient=1/(x[j]-x[i]);for(size_t k=first;k<first+width;++k)if(k!=i&&k!=j)coefficient*=(x[i]-x[k])/(x[j]-x[k]);}
            for(int c=0;c<16;++c)derivative[i][c]+=coefficient*y[j][c];}}
    Matrix output(query.size());for(size_t r=0;r<query.size();++r){const double q=query[r];Require(q>=x.front()-1e-12&&q<=x.back()+1e-12,"local spline extrapolation forbidden");
        size_t i=std::upper_bound(x.begin(),x.end(),q)-x.begin();i=i?i-1:0;i=std::min(i,n-2);const double h=x[i+1]-x[i],u=(q-x[i])/h;
        for(int c=0;c<16;++c)output[r][c]=(2*u*u*u-3*u*u+1)*y[i][c]+(u*u*u-2*u*u+u)*h*derivative[i][c]+(-2*u*u*u+3*u*u)*y[i+1][c]+(u*u*u-u*u)*h*derivative[i+1][c];}
    return output;
}
std::vector<std::pair<size_t,double>> LocalCubicWeights(const std::vector<double>& x,double q){
    const size_t n=x.size();Require(n>=2&&q>=x.front()-1e-12&&q<=x.back()+1e-12,"invalid local spline query");
    size_t i=std::upper_bound(x.begin(),x.end(),q)-x.begin();i=std::min(i?i-1:0,n-2);
    const double h=x[i+1]-x[i],u=(q-x[i])/h;Require(h>0,"invalid local spline interval");
    std::map<size_t,double> weights;weights[i]=2*u*u*u-3*u*u+1;weights[i+1]=-2*u*u*u+3*u*u;
    const size_t width=std::min(size_t(4),n);
    for(size_t at:{i,i+1}){const size_t first=std::min(at?at-1:0,n-width);const double factor=h*(at==i?u*u*u-2*u*u+u:u*u*u-u*u);
        for(size_t j=first;j<first+width;++j){double coefficient=0;
            if(j==at){for(size_t k=first;k<first+width;++k)if(k!=at)coefficient+=1/(x[at]-x[k]);}
            else{coefficient=1/(x[j]-x[at]);for(size_t k=first;k<first+width;++k)if(k!=at&&k!=j)coefficient*=(x[at]-x[k])/(x[j]-x[k]);}
            weights[j]+=factor*coefficient;}}
    std::vector<std::pair<size_t,double>> result;for(const auto& w:weights)if(w.second!=0)result.push_back(w);return result;
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

#include "FullAutoAllocation.inc"

// Schedule, never accept, an estimator. First learn from a cheap dominant
// correction before multiplying a much larger prefix whose noise is small.
// The 1/N work model affects only the order of refinements; all observed
// confidence, prefix, interpolation and independent-final gates remain live.
std::vector<int> PrioritizeAllocation(const std::vector<Samples>& levels,
        const std::vector<double>& costs,const std::vector<int>& counts,
        const std::vector<int>& proposed,double tolerance,bool ratios,bool useCovariance){
    Require(levels.size()==costs.size()&&levels.size()==counts.size()&&levels.size()==proposed.size()
            &&tolerance>0,"invalid priority allocation dimensions");
    if(levels.size()<2)return proposed;
    const auto stats=Estimate(Combine(levels));Require(stats.mean.size()==1,"priority allocation needs one angular row");
    int element=0;bool ratio=false;double worst=0;
    for(int j=0;j<16;++j){if(stats.interval[0][j]>worst){worst=stats.interval[0][j];element=j;ratio=false;}
        if(ratios&&stats.ratioInterval[0][j]>worst){worst=stats.ratioInterval[0][j];element=j;ratio=true;}}
    if(worst<=tolerance)return proposed; // Retain stability/interpolation work.
    const double r=ratio?stats.mean[0][element]/stats.mean[0][0]:0;
    std::vector<double> variance(levels.size()),score(levels.size());double sum=0;
    size_t best=0;
    for(size_t l=0;l<levels.size();++l){
        Require(counts[l]>0&&proposed[l]>=counts[l]&&costs[l]>0&&std::isfinite(costs[l]),"invalid priority work estimate");
        double mean=0;for(const auto& s:levels[l])mean+=(s[0][element]-r*s[0][0])/8;
        for(const auto& s:levels[l])variance[l]+=std::pow(s[0][element]-r*s[0][0]-mean,2)/7;
        sum+=variance[l];
        if(proposed[l]>counts[l]){
            // The gain includes covariance with the other levels. Refining
            // one member of a cancelling pair can otherwise increase noise.
            double covariance=0;
            for(int s=0;s<8;++s){double total=0;for(size_t k=0;k<levels.size();++k){double km=0;for(int t=0;t<8;++t)km+=(levels[k][t][0][element]-r*levels[k][t][0][0])/8;
                    total+=levels[k][s][0][element]-r*levels[k][s][0][0]-km;}
                covariance+=(levels[l][s][0][element]-r*levels[l][s][0][0]-mean)*total/7;}
            const double shrink=std::sqrt(double(counts[l])/proposed[l]);
            const double gain=useCovariance?std::max(0.,2*(1-shrink)*covariance-std::pow(1-shrink,2)*variance[l]):variance[l]*(1-double(counts[l])/proposed[l]);
            score[l]=gain/(costs[l]*(proposed[l]-counts[l]));}
        if(score[l]>score[best])best=l;
    }
    auto result=proposed;
    if(score[best]>0&&variance[best]>.5*sum)
        for(size_t l=0;l<levels.size();++l)
            if(l!=best&&score[l]*8<score[best]&&variance[l]<.25*sum)result[l]=counts[l];
    return result;
}

namespace {
struct Kernel {bool pipeline=false;int block=64;std::string warp="auto";};
struct Job {
    int offset=0,pairedDepth=0;
    std::string folder,data,summary,identity,weightHash,executionHash;
    Matrix matrix;std::vector<double> angles,inputAngles;std::vector<std::array<double,5>> controls;
    int preparationThreads=0;double assemblySeconds=0;
    double seconds=0,controlRefinement=0,controlSetup=0,traceSeconds=0,diffractionSeconds=0,processingSeconds=-1;pid_t pid=-1;
    int count=0,depth=0,phi=0,block=0;unsigned seed=0;
    bool pipeline=false,mirrored=false,paired=false,assembled=false,hasControls=false;std::string warp;
    // Phase sums overlap in the producer/consumer pipeline. Only the measured
    // processing wall clock can separate worker launch/output costs from work.
    double NativeWall()const{return processingSeconds>=0?processingSeconds:traceSeconds+diffractionSeconds;}
    double OutsideProcessing()const{return std::max(0.,seconds-NativeWall());}
    double OverlapFraction()const{
        const double smaller=std::min(traceSeconds,diffractionSeconds);
        return pipeline&&processingSeconds>=0&&smaller>0
            ?std::max(0.,std::min(1.,(traceSeconds+diffractionSeconds-processingSeconds)/smaller)):0.;
    }
};
void ReadPhaseTimings(Job& j,const std::string& text){
    std::istringstream profile(text);std::string line;std::getline(profile,line);
    while(std::getline(profile,line)){
        std::istringstream row(line);std::string phase,clock;double seconds;
        if(!(row>>phase>>seconds)||!std::isfinite(seconds)||seconds<0)continue;
        row>>clock;
        if(phase=="trace_prepare")j.traceSeconds=seconds;
        if(phase=="diffraction_and_controls")j.diffractionSeconds=seconds;
        if(phase=="processing"&&clock=="wall")j.processingSeconds=seconds;
    }
}
void ValidateStoredJob(const std::string& folder,std::set<std::string>& visited){
    if(!visited.insert(folder).second)return;
    std::istringstream metadata(Read(folder+"/complete.txt"));std::string id,dataHash,summaryHash,sourceHash,integralHash;double seconds=0;
    Require(bool(metadata>>id>>seconds>>dataHash>>summaryHash)&&seconds>0,"invalid cached metadata");metadata>>sourceHash>>integralHash;
    Require(FileHash(folder+"/result/result.dat")==dataHash&&(summaryHash=="none"||FileHash(folder+"/result/result_analytic_facets.tsv")==summaryHash),"cached job changed; do not reuse corrupted files");
    if(!integralHash.empty()){
        const auto manifest=Read(folder+"/integral_outputs.fnv");Require(Hash(manifest)==integralHash,"worker integral manifest changed");
        std::istringstream rows(manifest);std::string path,hash;while(rows>>path>>hash)Require(FileHash(path)==hash,"worker integral output changed");
    }
    if(Exists(folder+"/paired_depth.tsv")){const std::string manifest=Read(folder+"/paired_depth.tsv");Require(Hash(manifest)==sourceHash,"paired depth manifest changed");
        std::istringstream rows(manifest);std::string path,hash;while(rows>>path>>hash)Require(FileHash(path)==hash,"paired depth matrix changed");}
    if(Exists(folder+"/angular_sources.tsv")){
        const std::string sources=Read(folder+"/angular_sources.tsv");Require(Hash(sources)==sourceHash,"angular cache provenance changed");
        std::istringstream rows(sources);std::string row;while(std::getline(rows,row)){
            size_t tab=row.find('\t');Require(tab!=std::string::npos,"invalid angular cache source");
            const std::string source=row.substr(0,tab);Require(FileHash(source+"/complete.txt")==row.substr(tab+1),"angular cache source metadata changed");
            ValidateStoredJob(source,visited);
        }
    }
}
struct WorkerService {
    pid_t pid;int fd;std::string identity;unsigned long long requests=0;
    WorkerService(pid_t p,int f,const std::string& key):pid(p),fd(f),identity(key){}
    ~WorkerService(){
        close(fd);kill(pid,SIGTERM);int status;while(waitpid(pid,&status,0)<0&&errno==EINTR){}
        for(auto& active:activeWorkers)if(active==pid)active=0;
    }
};
struct Engine {
    const ArgPP& args;const RunConfig& config;std::string exe,root,signature,binaryHash,weightFile;
    std::string compatibleBinary,compatibleHash;
    std::set<bool> warmedCacheControls;
    long pendingLegacyAuditSerial=-1;int pendingLegacyAuditRound=-1;std::string pendingLegacyAuditSources;
    std::ofstream journal;std::string requestedOptions;size_t journalSequence=0;double runStarted=Seconds();
    std::vector<std::string> base,gpus;std::vector<double> theta,weightTheta;
    std::vector<std::array<double,2>> weights;std::vector<unsigned> seeds;
    std::map<std::string,Job> cache;Kernel kernel;int phi=4,target=12,pilot=512,kernelCount=256;
    // This is only a candidate index. SameSampling still verifies all optical
    // and execution options. Keep identity order identical to the old map scan:
    // overlapping angular rows must select the very same saved estimator.
    using SamplingKey=std::array<uint64_t,4>;
    using JobBucket=std::map<std::string,const Job*>;
    std::map<SamplingKey,JobBucket> samplingIndex;
    std::map<std::string,const Job*> folderIndex;
    static SamplingKey SamplingSlot(int n,unsigned seed,int depth,int offset){
        return {{uint64_t(n),uint64_t(seed),uint64_t(depth),uint64_t(offset)}};
    }
    void CacheJob(const Job& job){
        const auto key=SamplingSlot(job.count,job.seed,job.depth,job.offset);
        auto found=cache.find(job.identity);
        if(found!=cache.end()){
            const auto& old=found->second;
            Require(SamplingSlot(old.count,old.seed,old.depth,old.offset)==key&&old.folder==job.folder,
                    "cached identity changed sampling coordinates");
        }
        Job& stored=cache[job.identity];stored=job;
        samplingIndex[key][job.identity]=&stored;folderIndex[job.folder]=&stored;
    }
    const JobBucket& SamplingCandidates(int n,unsigned seed,int depth,int offset=0)const{
        static const JobBucket empty;
        const auto found=samplingIndex.find(SamplingSlot(n,seed,depth,offset));
        return found==samplingIndex.end()?empty:found->second;
    }
    int initial=1024,minimum=256,maxCount=67108864,roundLimit=32,startTheta=33;
    int maxTheta=2147483647,stablePasses=2,minPhi=1,maxPhi=2147483647;
    int prepareThreadCap=1;std::map<std::string,std::array<double,3>> preparationProfiles;
    int workerThreads=4,workerChunk=0,analyticMeanThreads=1;long long hostBudgetMb=0;
    std::string resourcePolicy="auto",rotationPolicy="auto",cutoffPolicy="auto",selectedCutoff="requested";
    bool tuneCutoffs=false,autoGpus=false;int rotationBeta=1,rotationGamma=1;double rotationError=0,cutoffError=0;
    std::vector<std::string> gpuCandidates;
    double eps=.01,interpolation=.0025,statistical=.0075;bool controls=true;
    bool mirror=false;std::string mirrorPolicy="auto";double mirrorError=0;
    std::set<std::string> verifiedMirrorGrids;
    bool angularCache=true,sharedMeans=false;size_t assembledJobs=0,reusedRows=0,newRows=0;
    bool incremental=true;int orientationBlock=65536;size_t savedOrientationBlocks=0;
    std::string controlPolicy="auto";int thetaAuditCount=0;
    bool local=true,ratios=true,phiRecheck=true,polarControls=true,finalAudit=true,physicsAudit=false;
    std::vector<std::array<Mueller,2>> polarWeights;
    double lastRatio=0,lastAnalyticError=0;bool independentPassed=false,independentNeedsSampling=true;
    int fixedPhi=0;std::string requestedMeanCache;bool coherent=true,fixedPipeline=false,fixedBlock=false,fixedWarp=false;
    bool depthSharing=true,depthVerified=false;
    bool completedCheckpoint=false;
    bool autoDepth=false;int depthStart=6,depthMaximum=18;double depthBudget=0,maximumEdge=0,maximumDiameter=0,particleVolume=0;
    double estimatedOrientations=0,estimatedBeta=0,estimatedGamma=0,empiricalFinal=0;
    int depthPilot=512;std::string initialModel="empirical",depthOrder="auto";bool startupKernel=false;std::string initialSource="explicit";
    bool collecting=false;bool go=false;std::string workQueuePolicy="auto";
    std::vector<Job> pendingJobs;
    std::vector<std::vector<std::string>> pendingCommands;
    std::string workerPolicy="auto",tracePolicy="auto";
    std::array<std::unique_ptr<WorkerService>,8> services;
    std::string localPolicy="auto";
    Engine(const ArgPP& a,const RunConfig& c):args(a),config(c){}
    void Log(const std::string& event,const std::string& payload){
        if(!journal.is_open()){
            const std::string path=root+"/fullauto_history.jsonl";
            if(Exists(path)){std::ifstream previous(path.c_str(),std::ios::binary);previous.seekg(0,std::ios::end);
                if(previous.tellg()>0){previous.seekg(-1,std::ios::end);char last;previous.get(last);if(last!='\n'){std::ofstream repair(path.c_str(),std::ios::app);repair<<'\n';}}}
            journal.open(path.c_str(),std::ios::app);Require(bool(journal),"cannot open adaptive history");
        }
        const double utc=std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
        journal<<"{\"schema\":1,\"unix_seconds\":"<<Number(utc)<<",\"elapsed_seconds\":"<<Number(Seconds()-runStarted)
            <<",\"method\":"<<Quoted(go?"go":"po")<<",\"cpu_workers\":"<<(config.useGpu?0:gpus.size())<<",\"controller_pid\":"<<getpid()<<",\"sequence\":"<<++journalSequence<<",\"event\":"<<Quoted(event)<<",\"data\":"<<CompactJson(payload)<<"}\n";
        journal.flush();Require(bool(journal),"adaptive history write failed");
    }
    void Record(const std::string& path,const std::string& text){
        Atomic(path,text);const auto slash=path.find_last_of('/');const std::string name=path.substr(slash==std::string::npos?0:slash+1);
        const bool json=name.size()>5&&name.substr(name.size()-5)==".json";
        Log("parameter_snapshot","{\"file\":"+Quoted(name)+",\"value\":"+(json?text:Quoted(text))+"}");
    }
    std::string Parameters(const std::string& state)const{
        std::ostringstream out;out<<"{\"configuration_fnv1a64\":"<<Quoted(signature)<<",\"binary_fnv1a64\":"<<Quoted(binaryHash)
            <<",\"compatible_origin_binary_fnv1a64\":"<<Quoted(compatibleHash)
            <<",\"requested_options\":"<<requestedOptions<<",\"effective_worker_base\":[";
        for(size_t i=0;i<base.size();++i)out<<(i?",":"")<<Quoted(base[i]);
        out<<"],\"epsilon\":"<<Number(eps)<<",\"interpolation_budget\":"<<Number(interpolation)<<",\"statistical_budget\":"<<Number(statistical)
            <<",\"pilot_count\":"<<pilot<<",\"kernel_pilot_count\":"<<kernelCount<<",\"dense_theta_pilot_count\":"<<thetaAuditCount
            <<",\"initial_model\":"<<Quoted(initialModel)<<",\"estimated_final_orientation_count\":"<<Number(empiricalFinal)<<",\"depth_pilot_count\":"<<depthPilot<<",\"depth_search_order\":"<<Quoted(depthOrder)<<",\"initial_source\":"<<Quoted(initialSource)<<",\"maximum_edge_um\":"<<Number(maximumEdge)<<",\"maximum_diameter_um\":"<<Number(maximumDiameter)<<",\"oldauto_N_beta_estimate\":"<<Number(estimatedBeta)<<",\"oldauto_N_gamma_estimate\":"<<Number(estimatedGamma)<<",\"oldauto_orientation_estimate\":"<<Number(estimatedOrientations)<<",\"initial_count\":"<<initial<<",\"minimum_correction_count\":"<<minimum<<",\"max_count\":"<<maxCount<<",\"max_rounds\":"<<roundLimit
            <<",\"initial_theta_support\":"<<startTheta<<",\"max_theta_evaluations\":"<<maxTheta<<",\"min_phi\":"<<minPhi<<",\"max_phi\":"<<maxPhi
            <<",\"stable_passes\":"<<stablePasses<<",\"automatic_depth\":"<<(autoDepth?"true":"false")<<",\"depth_start\":"<<depthStart<<",\"depth_maximum\":"<<depthMaximum<<",\"depth_budget\":"<<Number(depthBudget)<<",\"physical_target_depth\":"<<target<<",\"phi\":"<<phi
            <<",\"worker_lifetime\":"<<Quoted(workerPolicy)<<",\"trace_policy\":"<<Quoted(tracePolicy)<<",\"analytic_mean_threads\":"<<analyticMeanThreads<<",\"worker_threads\":"<<workerThreads<<",\"worker_chunk\":"<<workerChunk<<",\"orientation_block\":"<<orientationBlock
            <<",\"host_memory_budget_mib\":"<<hostBudgetMb<<",\"selected_gpus\":[";
        for(size_t i=0;i<gpus.size();++i)out<<(i?",":"")<<Quoted(gpus[i]);
        out<<"],\"resource_policy\":"<<Quoted(resourcePolicy)<<",\"rotation_policy\":"<<Quoted(rotationPolicy)<<",\"cutoff_policy\":"<<Quoted(cutoffPolicy)
            <<",\"mirror_policy\":"<<Quoted(mirrorPolicy)<<",\"controls_policy\":"<<Quoted(controlPolicy)<<",\"local_policy\":"<<Quoted(localPolicy)
            <<",\"work_queue_policy\":"<<Quoted(workQueuePolicy)<<",\"angular_cache\":"<<(angularCache?"true":"false")<<",\"incremental\":"<<(incremental?"true":"false")
            <<",\"coherent\":"<<(coherent?"true":"false")<<",\"depth_sharing\":"<<(depthSharing?"true":"false")<<",\"ratio_accuracy\":"<<(ratios?"true":"false")
            <<",\"phi_recheck\":"<<(phiRecheck?"true":"false")<<",\"polarization_controls\":"<<(polarControls?"true":"false")
            <<",\"independent_final_audit\":"<<(finalAudit?"true":"false")<<",\"physics_audit\":"<<(physicsAudit?"true":"false")
            <<",\"requested_theta_rows\":"<<theta.size()<<",\"requested_theta_file\":"<<Quoted(root+"/requested_theta.csv")
            <<",\"final_orientation_budget_file\":"<<Quoted(root+"/orientation_budget.csv")<<",\"selected_state\":"<<state<<"}";return out.str();
    }
    int Int(const std::string& key,int fallback){return args.IsCatched(key)?args.GetIntValue(key):fallback;}
    bool On(const char* key,bool fallback=true){return args.IsCatched(key)?args.GetStringValue(key)=="on":fallback;}
    void Init(int argc,const char* argv[],double edge=0,double diameter=0,double volume=0) {
        maximumEdge=edge;maximumDiameter=diameter;particleVolume=volume;
        go=config.method==RunMethod::GeometricalOptics;Require(!config.useFft,"requires a direct solver");
        for(const char* key:{"analytic_mean_reference","haar_mirror_audit","multigrid","multikeq","multikeq_list","orientfile_dedup","analytic_facet_samples","analytic_control_weights","analytic_backscatter","analytic_azimuth_gaussian","owen_avg"})
            Require(!args.IsCatched(key),std::string("unsupported modifier --")+key+"; run one coherent size without external controls");
        const RuntimeEnvironmentSnapshot environment=QueryRuntimeEnvironmentSnapshot();
        Require(environment.experimental.empty() || args.IsCatched("allow_experimental_environment"),"unset experimental MBS_* variables or acknowledge --allow-experimental-environment");
        for(const char* name:{"OMPI_COMM_WORLD_SIZE","PMI_SIZE","MV2_COMM_WORLD_SIZE"}){const char* v=getenv(name);Require(!v||std::atoi(v)<=1,"use one MPI rank; controller launches its own GPU workers");}
        char path[4096];ssize_t n=readlink("/proc/self/exe",path,sizeof(path)-1);Require(n>0,"cannot locate executable");path[n]=0;exe=path;binaryHash=FileHash(exe);
        root=args.IsCatched("o")?args.GetStringValue("o"):"fullauto_result";
        Require(root.find('%')==std::string::npos,"output placeholders are not supported; use an explicit directory");Directory(root);
        char* full=realpath(root.c_str(),nullptr);Require(full!=nullptr,"cannot resolve output directory");root=full;free(full);
        eps=args.GetDoubleValue("fullauto");interpolation=eps*.25;statistical=eps-interpolation;
        localPolicy=args.IsCatched("fullauto_local")?args.GetStringValue("fullauto_local"):"auto";local=localPolicy!="off";ratios=On("fullauto_ratios");phiRecheck=On("fullauto_phi_recheck");
        polarControls=On("fullauto_polar_controls");finalAudit=On("fullauto_final_audit");physicsAudit=On("fullauto_physics_audit",false);
        coherent=!go&&!args.IsCatched("incoh");depthSharing=!go&&coherent&&On("fullauto_depth_share");
        fixedPhi=args.IsCatched("nphi")?args.GetIntValue("nphi"):(args.IsCatched("grid")?args.GetIntValue("grid",args.GetArgNumber("grid")-2):0);
        fixedPipeline=args.IsCatched("orientation_pipeline");kernel.pipeline=fixedPipeline;
        if(const char* v=getenv("MBS_GPU_BLOCK")){int n=std::atoi(v);Require(std::string(v)==std::to_string(n)&&(n==64||n==128||n==256),"MBS_GPU_BLOCK must be 64, 128 or 256");kernel.block=n;fixedBlock=true;}
        if(const char* v=getenv("MBS_GPU_WARP_BEAMS")){Require(std::string(v)=="0"||std::string(v)=="1","MBS_GPU_WARP_BEAMS must be 0 or 1");kernel.warp=std::string(v)=="1"?"warp":"thread";fixedWarp=true;}
        if(args.IsCatched("analytic_mean_cache"))requestedMeanCache=args.GetStringValue("analytic_mean_cache");
        angularCache=!args.IsCatched("fullauto_angular_cache")||args.GetStringValue("fullauto_angular_cache")!="off";
        workerPolicy=args.IsCatched("fullauto_workers")?args.GetStringValue("fullauto_workers"):"auto";
        tracePolicy=args.IsCatched("po_trace_backend")?args.GetStringValue("po_trace_backend"):(args.IsCatched("fullauto_trace")?args.GetStringValue("fullauto_trace"):"auto");
        if(tracePolicy=="hybrid")kernel.pipeline=fixedPipeline=true;
        workQueuePolicy=args.IsCatched("fullauto_work_queue")?args.GetStringValue("fullauto_work_queue"):"auto";
        incremental=!args.IsCatched("fullauto_incremental")||args.GetStringValue("fullauto_incremental")!="off";
        orientationBlock=Int("fullauto_orientation_block",65536);
        Require(orientationBlock>0&&(orientationBlock&(orientationBlock-1))==0,"orientation block must be a positive power of two");
        resourcePolicy=args.IsCatched("fullauto_resources")?args.GetStringValue("fullauto_resources"):"auto";
        rotationPolicy=args.IsCatched("fullauto_rotation")?args.GetStringValue("fullauto_rotation"):"auto";
        cutoffPolicy=args.IsCatched("fullauto_cutoffs")?args.GetStringValue("fullauto_cutoffs"):"auto";
        workerChunk=Int("chunk",0);
        mirrorPolicy=args.IsCatched("mirror_gamma")?"on":(args.IsCatched("fullauto_mirror")?args.GetStringValue("fullauto_mirror"):"auto");
        pilot=Int("fullauto_pilot",512);kernelCount=Int("fullauto_kernel",256);initial=Int("fullauto_initial",1024);
        thetaAuditCount=Int("fullauto_theta_audit",(!args.IsCatched("grid")&&!args.IsCatched("tgrid"))?std::min(pilot,256):0);
        Require(thetaAuditCount>=0&&thetaAuditCount<=2147483646,"invalid theta audit count");
        minimum=Int("fullauto_correction",256);maxCount=Int("maxorient",67108864);roundLimit=Int("fullauto_rounds",32);startTheta=Int("fullauto_theta_start",33);
        for(int v:{pilot,kernelCount,initial,minimum,maxCount})Require(v>0&&v<2147483647,"counts must fit positive signed32bit");
        Require(roundLimit>0&&startTheta>=2,"invalid round/theta limits");
        if(config.adaptive.loadedFromFile){
            pilot=Int("fullauto_pilot",std::max(config.adaptive.minPilotOrientations,std::min(512,config.adaptive.maxPilotOrientations)));
            initial=Int("fullauto_initial",config.adaptive.minOrientations);
            if(!args.IsCatched("maxorient")&&config.adaptive.maxOrientations>0)maxCount=config.adaptive.maxOrientations;
            roundLimit=Int("fullauto_rounds",config.adaptive.maxJointSweeps);
            startTheta=Int("fullauto_theta_start",config.adaptive.minThetaPoints);maxTheta=config.adaptive.maxThetaPoints;
            minPhi=config.adaptive.minPhiPoints;maxPhi=config.adaptive.maxPhiPoints;
        }
        if(args.IsCatched("max_theta_points"))maxTheta=config.adaptive.maxThetaPoints;
        if(args.IsCatched("max_phi_points"))maxPhi=config.adaptive.maxPhiPoints;
        if(fixedPhi){phi=fixedPhi;minPhi=fixedPhi;if(!args.IsCatched("max_phi_points"))maxPhi=fixedPhi;}
        if(!args.IsCatched("fullauto_theta_start"))startTheta=std::min(startTheta,std::max(2,(maxTheta<1000?(maxTheta+3)/4:33)));
        Require(startTheta<=maxTheta,"initial theta support exceeds explicit theta limit");
        stablePasses=config.adaptive.stablePasses;
        if(!args.IsCatched("fullauto_initial"))initial=std::min(initial,maxCount);
        if(!args.IsCatched("fullauto_correction"))minimum=std::min(minimum,maxCount);
        if(!args.IsCatched("fullauto_pilot"))pilot=std::min(pilot,maxCount);
        if(!args.IsCatched("fullauto_kernel"))kernelCount=std::min(kernelCount,maxCount);
        initial=Power2(initial);minimum=Power2(minimum);Require(initial<=maxCount&&minimum<=maxCount,"rounded counts exceed max-orientations");
        for(int v:{pilot,kernelCount,initial,minimum,maxCount})Require(v>0&&v<2147483647,"invalid configured counts");
        Require(roundLimit>0&&startTheta>=2&&maxTheta>=2&&minPhi>0&&maxPhi>=minPhi,"invalid configured limits");
        depthStart=Int("fullauto_depth_start",6);depthMaximum=Int("fullauto_depth_max",18);
        autoDepth=!args.IsCatched("n")&&On("fullauto_depth_search");
        target=args.IsCatched("n")?config.maxReflections:depthStart;depthBudget=autoDepth?eps/4:0;
        Require(target>0&&depthStart>0&&depthMaximum>=depthStart&&depthMaximum<=28,"invalid depth search range");
        if(maximumEdge>0||(args.IsCatched("p")&&args.GetIntValue("p")==1)){
            if(!(maximumEdge>0)&&args.IsCatched("p")&&args.GetIntValue("p")==1)maximumEdge=std::max(args.GetDoubleValue("p",1),args.GetDoubleValue("p",2)/2);
            const double xi=.69*config.wavelengthUm/maximumEdge,step=xi/16;
            const double betaFactor=args.IsCatched("sym")?args.GetDoubleValue("sym",0):1,gammaFactor=args.IsCatched("sym")?args.GetDoubleValue("sym",1):1;
            estimatedBeta=std::max(3.,std::ceil(3.14159265358979323846/betaFactor/step));
            estimatedGamma=std::max(3.,std::ceil(2*3.14159265358979323846/gammaFactor/step));
            estimatedOrientations=estimatedBeta*estimatedGamma;
        }
        initialModel=args.IsCatched("fullauto_initial_model")?args.GetStringValue("fullauto_initial_model"):"empirical";
        depthOrder=args.IsCatched("fullauto_depth_order")?args.GetStringValue("fullauto_depth_order"):"auto";
        initialSource=args.IsCatched("fullauto_initial")?"explicit CLI":(config.adaptive.loadedFromFile?"adaptive-config":"pending empirical geometry model");

        controlPolicy=args.IsCatched("fullauto_controls")?args.GetStringValue("fullauto_controls"):((args.IsCatched("analytic_facet_average")||args.IsCatched("analytic_shadow_control"))?"on":"auto");
        if(controlPolicy=="auto"&&(args.IsCatched("analytic_facet_average")||args.IsCatched("analytic_shadow_control")))controlPolicy="on";
        controls=!go&&coherent&&controlPolicy!="off"&&config.refractiveReal>=1.;
        Require(controlPolicy!="on"||controls,"analytic facet controls require refractive REAL>=1");
        if(go){angularCache=false;pilot=0;kernelCount=0;thetaAuditCount=0;minimum=0;startTheta=0;controlPolicy="off";resourcePolicy="off";rotationPolicy="off";cutoffPolicy="off";mirrorPolicy="off";
            interpolation=0;statistical=eps;local=false;localPolicy="off";phi=1;fixedPhi=1;phiRecheck=false;polarControls=false;depthSharing=false;}
        SetMirror(false);
        seeds.assign(production,production+8);if(args.IsCatched("owen_seeds")){Require(args.GetArgNumber("owen_seeds")==8,"need exactly eight production seeds");
            for(int i=0;i<8;++i){int s=args.GetIntValue("owen_seeds",i);Require(s>=0,"negative seed");seeds[i]=s;}}
        std::set<unsigned> unique(seeds.begin(),seeds.end());Require(unique.size()==8,"production seeds must be distinct");
        for(unsigned s:seeds)for(int i=0;i<8;++i)Require(s!=training[i]&&s!=validation[i]&&s!=4001,"production/calibration seed overlap");
        for(unsigned s:seeds){for(unsigned reserved:{8001u,8023u,8101u,8111u,8117u,8123u,8147u,8161u,8167u,8171u,9001u,9007u,9011u,9013u,9029u,9041u,9043u,9049u})Require(s!=reserved,"production/tuning seed overlap");Require(s!=3001&&s!=5003,"production/geometry/mean-provider seed overlap");for(unsigned v:mirrorValidation)Require(s!=v,"production/mirror verification seed overlap");for(unsigned v:angularValidation)Require(s!=v,"production/theta pilot seed overlap");}
        if(go){const int intervals=args.IsCatched("grid")?args.GetIntValue("grid",3):36;
            Require(intervals>0&&intervals<1000000,"invalid GO bin count");Require(intervals+1<=maxTheta,"GO angular bins exceed explicit theta limit");
            for(int i=0;i<=intervals;++i)theta.push_back(180.*i/intervals);
        }else if(args.IsCatched("tgrid"))theta=ReadGrid(args.GetStringValue("tgrid"));
        else {double lo=0,hi=180;int intervals=720;
            if(args.IsCatched("fullauto_theta_range")){lo=args.GetDoubleValue("fullauto_theta_range",0);hi=args.GetDoubleValue("fullauto_theta_range",1);}
            if(!args.IsCatched("grid")){
                // A default 0.25-degree reference grid would miss narrow
                // forward lobes of large particles. Bound phase variation by
                // a geometric diameter before the interpolation search.
                double extent=0,volume=0;
                if(args.IsCatched("p")){
                    if(args.GetIntValue("p")==1){double height=args.GetDoubleValue("p",1),diameter=args.GetDoubleValue("p",2);
                        extent=std::hypot(height,diameter);volume=3*std::sqrt(3.)*diameter*diameter*height/8.;
                    }else{Require(maximumDiameter>0,"supply a scattering grid for this built-in geometry");extent=maximumDiameter;volume=particleVolume;}
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
        if(!args.IsCatched("fullauto_initial")&&!config.adaptive.loadedFromFile){
            Require(maximumEdge>0,"cannot determine geometric sampling estimate; supply fullauto-initial");
            Require(maxCount>=512,"automatic initial count needs max-orientations >=512 or an explicit fullauto-initial");
            double proposed=512;
            if(initialModel=="xi"){proposed=std::max(512.,estimatedOrientations);initialSource="legacy max(512, beta*gamma at xi/16)";}
            else if(go||!coherent){
                // With real index GO is invariant to uniform body scaling.
                // Absorption introduces an optical path / wavelength scale.
                const double tau=4*3.14159265358979323846*std::fabs(config.refractiveImag)*maximumDiameter/config.wavelengthUm;
                proposed=512*(1+std::sqrt(std::max(0.,tau)));empiricalFinal=proposed*16;initialSource=go?"GO scale invariance; optical-depth pilot":"incoherent fallback: bounded optical-depth pilot";
            }else{
                const bool prism=args.IsCatched("p")&&args.GetIntValue("p")==1;
                const double relative=maximumEdge/config.wavelengthUm;
                // empirical-v1: 234 completed records, 73 unique coherent
                // settings, measured variance removes excessive old floors.
                const double score=8.512893398898298+2.026426731719961*std::log2(1+relative/8)
                    +.6416260243693346*(theta.back()>90)+4.040351838305495*(theta.back()-theta.front()>90)
                    +.8222766179222369*(!controls)-2.703521520753631*(!prism);
                empiricalFinal=std::pow(2.,score)*std::pow(.01/eps,2);
                proposed=std::max(512.,empiricalFinal/16);initialSource="empirical-v1: completed-data variance model, predicted final/16";
            }
            const int cap=std::max(512,maxCount/4);initial=Power2(int(std::min(double(cap),proposed)));
        }
        depthPilot=Int("fullauto_depth_pilot",std::min(512,initial));
        Require(depthPilot>0&&depthPilot<=maxCount&&(depthPilot&(depthPilot-1))==0,"depth pilot must be a power of two within max-orientations");
        // The generated reference grid already bounds the angular search.
        // Share authoritative, lazily computed means. Independent per-row
        // refinement gives the same values for every requested subset; a
        // large reference grid no longer forces a cold calculation of all rows.
        sharedMeans=controls;
        if(config.useGpu){int count=VisibleGpuDeviceCount();Require(count>0,"no visible CUDA GPU");
            std::vector<std::string> visible;const char* env=getenv("CUDA_VISIBLE_DEVICES");if(env){std::istringstream f(env);std::string v;while(std::getline(f,v,','))visible.push_back(v);}
            std::string list=args.IsCatched("fullauto_gpus")?args.GetStringValue("fullauto_gpus"):"0";autoGpus=list=="auto";if(autoGpus){list="";for(int i=0;i<std::min(count,8);++i)list+=(i?",":"")+std::to_string(i);}
            std::istringstream f(list);std::string v;std::set<int> used;
            while(std::getline(f,v,',')){size_t pos=0;int id=std::stoi(v,&pos);Require(pos==v.size()&&id>=0&&id<count&&used.insert(id).second,"invalid/duplicate visible GPU index");gpus.push_back(visible.empty()?v:visible.at(id));}
        Require(!gpus.empty(),"empty GPU list");}else{Require(!args.IsCatched("fullauto_gpus"),"GPU list requires CUDA backend");gpus.assign(size_t(Int("fullauto_cpu_workers",1)),"");}
        workerThreads=config.threads>0?config.threads:std::max(1,std::min(16,QueryRuntimeResourceSnapshot().physicalCores/int(gpus.size())));
        // A shared mean cache has one CPU owner per controller, not per GPU.
        // Reserve the selected GPUs' share of host cores for its independent
        // theta rows. Small GPU-pilot thread tuning must not shrink this pool.
        const int physical=std::max(1,QueryRuntimeResourceSnapshot().physicalCores);
        const int perGpu=config.useGpu?physical/std::max(1,VisibleGpuDeviceCount()):1;
        while(prepareThreadCap<8&&prepareThreadCap*2<=std::max(1,perGpu/2))prepareThreadCap*=2;
        const int share=config.useGpu?std::max(1,physical*int(gpus.size())/VisibleGpuDeviceCount())
                                     :std::max(1,physical/int(gpus.size()));
        analyticMeanThreads=Int("analytic_mean_threads",config.threads>0?config.threads:std::min(32,share));
        Require(analyticMeanThreads>0,"analytic mean threads must be positive");
        std::set<std::string> drop={"fullauto","autofull","o","n","grid","tgrid","nphi","haar_alpha","mirror_gamma","analytic_facet_average","analytic_shadow_control","analytic_mean_cache","analytic_control_weights","close","profile_phases","orientation_pipeline","checkpoint","maxorient","max_theta_points","max_phi_points","adaptive_config","owen_seeds","convergence_passes","threads","chunk","analytic_mean_threads"};
        if(args.IsCatched("fullauto_resume_compatible")){
            compatibleBinary=args.GetStringValue("fullauto_resume_compatible");
            Require(Exists(root+"/configuration.txt")&&Exists(root+"/controller_command.json")&&Exists(root+"/fullauto_parameters.json"),
                "compatible resume requires an existing fullauto result");
            compatibleHash=FileHash(compatibleBinary);
            const std::string oldVersion=ExecutableVersion(compatibleBinary),newVersion=ExecutableVersion(exe);
            const size_t oldLine=oldVersion.find('\n'),newLine=newVersion.find('\n');
            Require(oldLine!=std::string::npos&&(oldVersion.substr(0,oldLine)=="MBS-fast 8ae9e85"||oldVersion.substr(0,oldLine)=="MBS-fast 89f83a3"||oldVersion.substr(0,oldLine)=="MBS-fast 336ef78"),
                "compatible resume is verified only from releases 8ae9e85, 89f83a3 and 336ef78");
            Require(newLine!=std::string::npos&&oldVersion.substr(oldLine)==newVersion.substr(newLine),
                "compatible resume requires identical CUDA architecture, precision, phase math, compiler, OpenMP and MPI profiles");
        }
        std::ostringstream fingerprint;fingerprint<<(compatibleHash.empty()?binaryHash:compatibleHash)<<'\n';
        for(int i=1;i<argc;++i){
            std::string option=argv[i];while(!option.empty()&&option[0]=='-')option.erase(0,1);
            if(argv[i][0]=='-' && (option=="fullauto-resume-compatible"||option=="fullauto_resume_compatible")){++i;continue;}
            fingerprint<<argv[i]<<'\n';
        }
        // Canonical reconstruction avoids shell quoting and option aliases.
        std::set<std::string> seen;for(const auto& spec:GetCliOptionSpecs())if(args.IsCatched(spec.key)&&seen.insert(spec.key).second){
            if(spec.key=="pf"||spec.key=="tgrid"||spec.key=="adaptive_config")fingerprint<<FileHash(args.GetStringValue(spec.key))<<'\n';
            if(drop.count(spec.key)||spec.key.compare(0,9,"fullauto_")==0)continue;
            base.push_back("--"+spec.canonical);for(unsigned i=0;i<args.GetArgNumber(spec.key);++i)base.push_back(args.GetStringValue(spec.key,i));}
        // File headers can declare an unverified fundamental domain. Use
        // the full Haar measure for arbitrary files unless explicitly set.
        if(!args.IsCatched("sym")){base.push_back("--symmetry");base.push_back("1");base.push_back("1");fingerprint<<"verified_orientation_domain=full_SO3_initial\n";}
        if(tracePolicy!="auto"&&!args.IsCatched("po_trace_backend")){base.push_back("--po-trace-backend");base.push_back(tracePolicy);}
        gpuCandidates=gpus;
        hostBudgetMb=std::max(64LL,QueryRuntimeResourceSnapshot().ramAvailableKb*65/100/1024/static_cast<long long>(std::min(gpus.size(),size_t(8))));
        // Result-affecting environment forms part of resume/cache identity.
        fingerprint<<FormatRuntimeEnvironmentReport(environment);
        char host[256]={};Require(gethostname(host,sizeof(host)-1)==0,"cannot determine host identity");
        fingerprint<<"host="<<host<<"\nthreads="<<workerThreads<<"\nanalytic_mean_threads="<<analyticMeanThreads<<'\n';
        const char* visible=getenv("CUDA_VISIBLE_DEVICES");fingerprint<<"visible="<<(visible?visible:"")<<'\n';
        for(const auto& gpu:gpus)fingerprint<<"gpu="<<gpu<<'\n';
        signature=Hash(fingerprint.str());
        std::ostringstream supplied;supplied<<"{\"priority\":\"explicit CLI > config file > automatic selection\",\"options\":[";bool first=true;
        for(const auto& spec:GetCliOptionSpecs())if(args.IsCatched(spec.key)){
            supplied<<(first?"":",")<<"{\"flag\":"<<Quoted(CliCanonicalFlag(spec.key))<<",\"policy\":"<<Quoted(FullAutoOptionPolicy(spec.key))<<",\"values\":[";
            for(unsigned i=0;i<args.GetArgNumber(spec.key);++i)supplied<<(i?",":"")<<Quoted(args.GetStringValue(spec.key,i));supplied<<"]}";first=false;
        }supplied<<"]}";requestedOptions=supplied.str();
    }
    void State(const std::string& status,int round,const std::vector<int>& counts,size_t nodes,double ci,double ie,int streak) {
        std::ostringstream s;s<<std::setprecision(17)<<"{\n\"status\":"<<Quoted(status)<<",\"configuration_fnv1a64\":"<<Quoted(signature)
            <<",\"method\":"<<Quoted(go?"go":"po")<<",\"cpu_workers\":"<<(config.useGpu?0:gpus.size())<<",\"controller_pid\":"<<getpid()<<",\"round\":"<<round<<",\"counts_per_seed\":[";
        for(size_t i=0;i<counts.size();++i)s<<(i?",":"")<<counts[i];
        s<<"],\"seeds\":[";
        for(size_t i=0;i<seeds.size();++i)s<<(i?",":"")<<seeds[i];
        s<<"],\"phi\":"<<phi<<",\"target_depth\":"<<target<<",\"requested_theta\":"<<theta.size()<<",\"evaluated_theta\":"<<nodes
            <<",\"automatic_depth\":"<<(autoDepth?"true":"false")<<",\"depth_budget\":"<<depthBudget<<",\"statistical_budget\":"<<statistical<<",\"interpolation_budget\":"<<interpolation<<",\"max_pointwise95_scaled_M11\":"<<ci
            <<",\"max_guard_residual95_scaled_M11\":"<<ie<<",\"pass_streak\":"<<streak<<",\"pipeline\":"<<(kernel.pipeline?"true":"false")
            <<",\"block\":"<<kernel.block<<",\"warp\":"<<Quoted(kernel.warp)
            <<",\"mirror_gamma\":"<<(mirror?"true":"false")<<",\"mirror_policy\":"<<Quoted(mirrorPolicy)
            <<",\"analytic_controls\":"<<(controls?"true":"false")<<",\"analytic_control_policy\":"<<Quoted(controlPolicy)<<",\"shared_analytic_reference\":"<<(sharedMeans?"true":"false")<<",\"angular_cache\":"<<(angularCache?"true":"false")<<",\"work_queue\":"<<(QueueEnabled()?"true":"false")<<",\"work_queue_policy\":"<<Quoted(workQueuePolicy)<<",\"assembled_jobs_this_invocation\":"<<assembledJobs<<",\"reused_angular_rows_this_invocation\":"<<reusedRows<<",\"computed_new_angular_rows_this_invocation\":"<<newRows
            <<",\"mirror_budget\":"<<(mirror?eps/256:0)<<",\"max_mirror_residual95_scaled_M11\":"<<mirrorError
            <<",\"prepare_thread_cap\":"<<prepareThreadCap<<",\"accuracy_scope\":\"production\",\"adaptive_schema\":2,\"local_orientations\":"<<(local?"true":"false")<<",\"ratio_accuracy\":"<<(ratios?"true":"false")<<",\"max_pointwise_ratio95\":"<<lastRatio<<",\"independent_final_validation\":"<<(independentPassed?"true":"false")<<",\"saved_orientation_blocks_this_invocation\":"<<savedOrientationBlocks
            <<",\"analytic_mean_budget\":"<<(controls?eps/64:0)<<",\"max_estimated_analytic_mean_error_scaled_M11\":"<<lastAnalyticError
            <<",\"execution_budget\":"<<((resourcePolicy=="auto"||config.useGpu)?eps/1024:0)<<",\"host_memory_budget_mib\":"<<hostBudgetMb<<",\"worker_lifetime\":"<<Quoted(workerPolicy)<<",\"trace_policy\":"<<Quoted(tracePolicy)<<",\"analytic_mean_threads\":"<<analyticMeanThreads<<",\"worker_threads\":"<<workerThreads<<",\"worker_chunk\":"<<workerChunk<<",\"orientation_block\":"<<orientationBlock
            <<",\"rotation_beta_factor\":"<<rotationBeta<<",\"rotation_gamma_factor\":"<<rotationGamma<<",\"rotation_budget\":"<<((rotationBeta>1||rotationGamma>1)?eps/256:0)
            <<",\"selected_cutoff\":"<<Quoted(selectedCutoff)<<",\"cutoff_budget\":"<<(selectedCutoff!="requested"&&selectedCutoff!="off"?eps/8:0)<<",\"max_cutoff_residual95\":"<<cutoffError
            <<",\"coherent\":"<<(coherent?"true":"false")<<",\"fixed_phi\":"<<fixedPhi<<",\"scope\":\"eight-seed pointwise intervals; interpolation guards on requested grid; mirror paired pilot when enabled; current target physics; automatic depth has fresh paired validation; no global simultaneous confidence claim\"}\n";
        Record(root+"/fullauto_status.json",s.str());
        if(status=="converged")Record(root+"/final_parameters.json",Parameters(s.str()));
        std::cout<<"fullauto "<<status<<" round="<<round<<" theta="<<nodes;
        if(status=="evaluating")std::cout<<" accuracy=pending";
        else if(status=="validating")std::cout<<" production_ci="<<ci<<" production_interpolation="<<ie<<" production_passes="<<streak<<'/'<<stablePasses<<" independent_audit=pending";
        else std::cout<<" ci="<<ci<<" interpolation="<<ie<<" streak="<<streak;
        std::cout<<std::endl;
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
        const std::string log=Read(j.folder+"/stdout.log");const auto preparation=log.find("PO preparation threads: ");
        if(preparation!=std::string::npos){std::istringstream count(log.substr(preparation+24));count>>j.preparationThreads;}
        for(const char* label:{"Hard tree-limit hits:","Orientations still incomplete:"}){
            size_t pos=log.find(label);Require(pos!=std::string::npos,"missing tree diagnostics");std::istringstream r(log.substr(pos+std::string(label).size()));int n=-1;r>>n;Require(n==0,"worker tree limit or incomplete orientation");}
        if(Exists(j.folder+"/result/result_phase_timings.tsv"))ReadPhaseTimings(j,Read(j.folder+"/result/result_phase_timings.tsv"));
        if(Exists(j.folder+"/angular_reuse.json")){const auto text=Read(j.folder+"/angular_reuse.json");const std::string key="\"controller_assembly_seconds\":";const auto pos=text.find(key);
            if(pos!=std::string::npos){std::istringstream timing(text.substr(pos+key.size()));double value=0;if(timing>>value&&std::isfinite(value)&&value>=0)j.assemblySeconds=value;}}

        std::ostringstream index;index<<signature<<' '<<j.identity<<' '<<j.count<<' '<<j.offset<<' '<<j.seed<<' '<<j.depth<<' '<<j.phi<<' '<<j.block<<' '<<j.pipeline<<' '<<j.mirrored<<' '<<j.paired<<' '<<j.assembled<<' '<<j.hasControls<<' '<<j.pairedDepth<<' '<<j.warp<<' '<<j.weightHash<<' '<<j.executionHash<<'\n';
        for(double angle:j.inputAngles)index<<Number(angle)<<'\n';
        Atomic(j.folder+"/sampling.index",index.str());Atomic(j.folder+"/sampling.index.fnv",Hash(index.str()));
        return j;
    }
    void WarmCache(){
        // The controller lock excludes outside writers. Completed workers are
        // inserted by CompleteNative/CacheJob, so each controls domain needs
        // one verified disk pass per process, independent of tuning/depth.
        if(warmedCacheControls.count(controls)){
            Log("cache_restore_reuse","{\"verified_jobs\":"+Number(cache.size())+",\"controls\":"+(controls?"true":"false")+"}");return;
        }
        DIR* directory=opendir((root+"/jobs").c_str());if(!directory)return;std::vector<std::string> folders;
        while(auto* entry=readdir(directory)){if(entry->d_name[0]!='.')folders.push_back(root+"/jobs/"+entry->d_name);}closedir(directory);
        std::set<std::string> visited;for(const auto& folder:folders){if(!Exists(folder+"/complete.txt")||!Exists(folder+"/sampling.index"))continue;
            const std::string text=Read(folder+"/sampling.index");Require(Hash(text)==Read(folder+"/sampling.index.fnv"),"sampling index changed");
            std::istringstream input(text);std::string configuration;Job j;
            Require(bool(input>>configuration>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block>>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash),"invalid sampling index");
            if(configuration!=signature||j.hasControls!=controls)continue;
            double angle;while(input>>angle)j.inputAngles.push_back(angle);Require(input.eof()&&j.inputAngles.size()>=2,"invalid indexed grid");
            j.folder=folder;j.data=folder+"/result/result.dat";j.summary=folder+"/result/result_analytic_facets.tsv";
            ValidateStoredJob(folder,visited);std::istringstream meta(Read(folder+"/complete.txt"));std::string identity;meta>>identity>>j.seconds;Require(identity==j.identity,"sampling identity changed");
            j=Load(j);Require(j.angles.size()==j.inputAngles.size(),"sampling coordinates changed");for(size_t r=0;r<j.angles.size();++r)Require(std::fabs(j.angles[r]-j.inputAngles[r])<1e-10,"sampling coordinates changed");CacheJob(j);
        }
        warmedCacheControls.insert(controls);
        Log("cache_restored","{\"verified_jobs\":"+Number(cache.size())+"}");
        std::cout<<"fullauto restored verified cached jobs: "<<cache.size()<<std::endl;
    }
    std::string WriteWeights(const std::vector<double>& x) {
        if(!controls||weights.empty())return "";
        std::ostringstream f;f<<"theta_deg reflection_weight shadow_weight\n"<<std::setprecision(17);
        for(double q:x){size_t k=std::upper_bound(weightTheta.begin(),weightTheta.end(),q)-weightTheta.begin();k=k?k-1:0;k=std::min(k,weightTheta.size()-2);
            double t=(q-weightTheta[k])/(weightTheta[k+1]-weightTheta[k]);f<<q;for(int j=0;j<2;++j)f<<' '<<((1-t)*weights[k][j]+t*weights[k+1][j]);f<<'\n';}
        std::string path=root+"/weights_"+Hash(f.str())+".tsv";Atomic(path,f.str());return path;
    }
    std::string MeanCache(const std::string& gridHash)const{return requestedMeanCache.empty()?(sharedMeans?root+"/means_shared.cache":root+"/mean_"+gridHash+".cache"):requestedMeanCache;}
    void SetSampling(Job& j,int n,unsigned seed,int depth,const std::vector<double>& x,bool paired){
        j.executionHash=ExecutionHash();j.count=n;j.seed=seed;j.depth=depth;j.phi=phi;j.block=kernel.block;j.pipeline=kernel.pipeline;j.warp=kernel.warp;
        j.mirrored=mirror&&!paired;j.paired=paired;j.inputAngles=x;j.hasControls=controls;
    }
    bool SameSampling(const Job& j,int n,unsigned seed,int depth,int offset=0)const{
        return j.offset==offset&&j.count==n&&j.seed==seed&&j.depth==depth&&j.phi==phi&&j.block==kernel.block
            &&j.pipeline==kernel.pipeline&&j.warp==kernel.warp&&j.mirrored==mirror&&!j.paired&&j.hasControls==controls&&j.executionHash==ExecutionHash();
    }
    int AutoPrepareThreads(int n,int depth,int paired,size_t rows){
        if(!config.useGpu||go||resourcePolicy!="auto"||args.IsCatched("threads")||args.IsCatched("po_prepare_threads")||n<1024||prepareThreadCap<=workerThreads)return 0;
        const auto recipe=ExecutionHash();const std::string key=recipe+" "+Number(depth)+" "+Number(paired)+" "+Number(phi)+" "+Number(mirror)+" "+Number(controls);
        auto saved=preparationProfiles.find(key);std::array<double,3> profile{};
        if(saved!=preparationProfiles.end())profile=saved->second;
        else{for(const auto& item:cache){const auto& j=item.second;
            if(j.assembled||j.paired||j.executionHash!=recipe||j.depth!=depth||j.pairedDepth!=paired||j.phi!=phi||j.mirrored!=mirror||j.hasControls!=controls||j.count<=0||j.inputAngles.empty()||j.traceSeconds<=0)continue;
            profile[0]+=j.traceSeconds*std::max(workerThreads,j.preparationThreads)/j.count;
            profile[1]+=j.diffractionSeconds/j.count/j.inputAngles.size();++profile[2];}
            if(profile[2]){profile[0]/=profile[2];profile[1]/=profile[2];preparationProfiles[key]=profile;}}
        if(!profile[2]||profile[0]/workerThreads<=1.25*profile[1]*rows)return 0;
        Log("worker_preparation_selected","{\"count\":"+Number(n)+",\"depth\":"+Number(depth)+",\"theta_rows\":"+Number(rows)+",\"threads\":"+Number(prepareThreadCap)+",\"one_thread_trace_seconds_forecast\":"+Number(profile[0]*n)+",\"diffraction_seconds_forecast\":"+Number(profile[1]*n*rows)+",\"ordered_accumulation_unchanged\":true}");
        return prepareThreadCap;
    }
    std::vector<Job> DirectStage(int n,unsigned const* seedSet,size_t ns,int depth,const std::vector<double>& x,bool weighted=true,bool pairedMirror=false,int offset=0,int pairedDepth=0) {
        Require(x.size()<=size_t(maxTheta)+(go?2:0),"explicit theta evaluation limit exceeded");Require(FileHash(exe)==binaryHash,"binary changed during run");std::ostringstream angles;for(double q:x)angles<<Number(q)<<'\n';std::string gridHash=Hash(angles.str());
        std::string grid=root+"/theta_"+gridHash+".csv";Atomic(grid,angles.str());std::string wf=weighted?WriteWeights(x):"";
        const int preparationThreads=AutoPrepareThreads(n,depth,pairedDepth,x.size());
        std::vector<Job> result(ns);std::vector<std::vector<std::string>> commands(ns);std::vector<size_t> todo;
        for(size_t s=0;s<ns;++s){std::vector<std::string> command=base;
            for(const std::string& v:std::vector<std::string>{"--max-reflections",std::to_string(depth),"--sobol-seed",std::to_string(n),std::to_string(seedSet[s]),"--threads",std::to_string(workerThreads),"--close"})command.push_back(v);
            if(go){for(const auto& v:std::vector<std::string>{"--scattering-grid","0","180","1",Number(theta.size()-1),"--go-fullauto-worker"})command.push_back(v);}
            else{for(const auto& v:std::vector<std::string>{"--haar-alpha","--theta-grid-file",grid,"--phi-points",std::to_string(phi)})command.push_back(v);}
            if(pairedMirror)command.push_back("--haar-mirror-audit");else if(mirror)command.push_back("--mirror-gamma");
            if(offset){command.push_back("--sobol-offset");command.push_back(std::to_string(offset));}
            if(pairedDepth){command.push_back("--paired-depth");command.push_back(std::to_string(pairedDepth));}
            if(!go)command.push_back("--profile-phases");
            if(controls){command.push_back("--analytic-mean-threads");command.push_back(Number(analyticMeanThreads));command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(MeanCache(gridHash));if(sharedMeans){command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");}}
            if(workerChunk>0){command.push_back("--orientation-chunk");command.push_back(Number(workerChunk));}
            if(kernel.pipeline)command.push_back("--orientation-pipeline");
            if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}
            std::ostringstream identity;identity<<signature<<ExecutionHash()<<'\n'<<kernel.block<<'\n'<<kernel.warp<<'\n';for(const auto& v:command)identity<<v<<'\n';if(!wf.empty())identity<<Read(wf);
            Job j;SetSampling(j,n,seedSet[s],depth,x,pairedMirror);j.offset=offset;j.pairedDepth=pairedDepth;j.weightHash=Hash(wf.empty()?"":Read(wf));j.identity=Hash(identity.str());j.folder=root+"/jobs/j"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";Directory(j.folder);
            auto correctAngles=[&](const Job& job){Require(job.angles.size()==x.size(),"wrong theta count");for(size_t r=0;r<x.size();++r)Require(std::fabs(job.angles[r]-x[r])<1e-10,"worker theta grid differs");};
            if(cache.count(j.identity)){result[s]=cache.at(j.identity);correctAngles(result[s]);continue;}
            const std::string meta=j.folder+"/complete.txt";
            if(Exists(meta)){std::istringstream m(Read(meta));std::string id,dataHash,summaryHash;m>>id>>j.seconds>>dataHash>>summaryHash;
                Require(id==j.identity&&j.seconds>0&&std::isfinite(j.seconds)&&FileHash(j.data)==dataHash&&(!controls||FileHash(j.summary)==summaryHash),"cached job changed; do not reuse corrupted files");
                {std::set<std::string> visited;ValidateStoredJob(j.folder,visited);}j=Load(j);correctAngles(j);CacheJob(j);result[s]=j;Log("worker_cache_reuse","{\"folder\":"+Quoted(j.folder)+",\"verified_from_disk\":true}");continue;}
            if(Exists(j.folder+"/worker.pid"))Require(!WorkerAlive(j.folder+"/worker.pid"),"an incomplete cached worker is still running; avoid duplicate launch");
            if(Exists(j.folder+"/result")){std::string saved=j.folder+"/result_incomplete_"+std::to_string(getpid())+"_"+std::to_string(s);
                Require(rename((j.folder+"/result").c_str(),saved.c_str())==0,"cannot preserve incomplete worker output");}
            if(preparationThreads){j.preparationThreads=preparationThreads;command.push_back("--po-prepare-threads");command.push_back(Number(preparationThreads));}
            command.push_back("--output");command.push_back(j.folder+"/result");commands[s]=command;result[s]=j;todo.push_back(s);
        }
        if(collecting){for(size_t index:todo){pendingJobs.push_back(result[index]);pendingCommands.push_back(commands[index]);}return result;}
        ExecutePrepared(result,commands,todo);
        return result;
    }
    #include "FullAutoWorkers.inc"
    void ExecutePrepared(std::vector<Job>& result,std::vector<std::vector<std::string>>& commands,const std::vector<size_t>& todo){
        if(PersistentEnabled()){ExecutePersistent(result,commands,todo);return;}
        for(auto& service:services)service.reset();
        WorkerCleanup cleanup;
        const size_t slots=std::min(gpus.size(),size_t(8));
        std::vector<int> slotSeed(slots,-1);std::vector<double> begun(result.size());
        size_t next=0,completed=0;
        auto launch=[&](size_t slot,size_t s){Job& j=result[s];begun[s]=Seconds();std::string text=exe+"\n";for(const auto& v:commands[s])text+=v+"\n";Atomic(j.folder+"/command.args",text);
                std::vector<char*> argv;argv.push_back(const_cast<char*>(exe.c_str()));for(auto& v:commands[s])argv.push_back(const_cast<char*>(v.c_str()));argv.push_back(nullptr);
                j.pid=fork();Require(j.pid>=0,"cannot fork worker");
                if(j.pid==0){signal(SIGTERM,SIG_DFL);signal(SIGINT,SIG_DFL);int log=open((j.folder+"/stdout.log").c_str(),O_WRONLY|O_CREAT|O_TRUNC,0644);if(log<0)_exit(126);dup2(log,1);dup2(log,2);close(log);
                    if(config.useGpu){setenv("CUDA_VISIBLE_DEVICES",gpus[slot].c_str(),1);unsetenv("MBS_GPU_DEVICE");}
                    if(!getenv("MBS_HOST_MEM_BUDGET_MB"))setenv("MBS_HOST_MEM_BUDGET_MB",std::to_string(hostBudgetMb).c_str(),1);
                    setenv("MBS_GPU_BLOCK",std::to_string(j.block).c_str(),1);
                    if(j.warp=="auto")unsetenv("MBS_GPU_WARP_BEAMS");else setenv("MBS_GPU_WARP_BEAMS",j.warp=="warp"?"1":"0",1);
                    execv(exe.c_str(),argv.data());_exit(127);}
                activeWorkers[slot]=j.pid;slotSeed[slot]=int(s);
                Log("worker_start","{\"folder\":"+Quoted(j.folder)+",\"command_file\":"+Quoted(j.folder+"/command.args")+",\"pid\":"+Number(j.pid)
                    +",\"gpu\":"+Quoted(gpus[slot])+",\"count\":"+Number(j.count)+",\"offset\":"+Number(j.offset)+",\"seed\":"+Number(j.seed)
                    +",\"depth\":"+Number(j.depth)+",\"paired_depth\":"+Number(j.pairedDepth)+",\"phi\":"+Number(j.phi)+",\"theta_rows\":"+Number(j.inputAngles.size())
                    +",\"block\":"+Number(j.block)+",\"warp\":"+Quoted(j.warp)+",\"pipeline\":"+(j.pipeline?"true":"false")+",\"threads\":"+Number(workerThreads)+",\"chunk\":"+Number(workerChunk)+"}");std::string boot=Read("/proc/sys/kernel/random/boot_id");
                Atomic(j.folder+"/worker.pid",std::to_string(j.pid)+" "+ProcessStart(j.pid)+" "+boot);};
        for(size_t slot=0;slot<slots&&next<todo.size();++slot)launch(slot,todo[next++]);
        while(completed<todo.size()){
                int status=0;pid_t waited;do{waited=waitpid(-1,&status,0);}while(waited<0&&errno==EINTR);
                Require(waited>0,"cannot reap worker");
                size_t slot=0;while(slot<slots&&activeWorkers[slot]!=waited)++slot;Require(slot<slots,"unknown worker PID");
                size_t s=size_t(slotSeed[slot]);Job& j=result[s];j.pid=-1;activeWorkers[slot]=0;slotSeed[slot]=-1;j.seconds=Seconds()-begun[s];std::remove((j.folder+"/worker.pid").c_str());
                Require(WIFEXITED(status)&&WEXITSTATUS(status)==0,"worker failed; see "+j.folder+"/stdout.log");
                if(next<todo.size())launch(slot,todo[next++]);
                CompleteNative(j);++completed;
        }
    }
    bool QueueEnabled()const{return (go||angularCache)&&workQueuePolicy!="off"&&(workQueuePolicy=="on"||config.useGpu||gpus.size()>1);}
    struct WorkRequest {int count,offset,depth,pairedDepth;std::vector<double> angles;};
    void Prefetch(const std::vector<WorkRequest>& requests,const unsigned* seedSet){
        if(!QueueEnabled()||requests.empty())return;
        auto flush=[&](){
            std::vector<size_t> todo(pendingJobs.size());for(size_t i=0;i<todo.size();++i)todo[i]=i;
            std::stable_sort(todo.begin(),todo.end(),[&](size_t a,size_t b){const auto& x=pendingJobs[a];const auto& y=pendingJobs[b];
                return double(x.count)*x.depth*x.inputAngles.size()>double(y.count)*y.depth*y.inputAngles.size();});
            const double begin=Seconds();ExecutePrepared(pendingJobs,pendingCommands,todo);
            if(!todo.empty()){double seconds=0;for(const auto& job:pendingJobs)seconds+=job.seconds;
                std::ofstream log((root+"/work_queue.tsv").c_str(),std::ios::app);
                log<<todo.size()<<'\t'<<gpus.size()<<'\t'<<Number(Seconds()-begin)<<'\t'<<Number(seconds)<<'\n';}
            pendingJobs.clear();pendingCommands.clear();
        };
        int maximum=0;for(const auto& request:requests)maximum=std::max(maximum,request.count);
        const int block=incremental?orientationBlock:maximum;
        collecting=true;
        try{for(int offset=0;offset<maximum;offset+=block){
                // Keep only one window's descriptors and at most 128 native
                // jobs in memory; even huge orientation budgets stream.
                std::map<std::string,WorkRequest> leaves;
                for(const auto& request:requests){if(offset>=request.count)continue;const int count=std::min(block,request.count-offset);
                    const std::string key=Number(count)+" "+Number(request.depth)+" "+Number(request.pairedDepth);
                    auto found=leaves.find(key);if(found==leaves.end())found=leaves.emplace(key,WorkRequest{count,offset,request.depth,request.pairedDepth,{}}).first;
                    auto& angles=found->second.angles;angles.insert(angles.end(),request.angles.begin(),request.angles.end());}
                for(auto& entry:leaves){auto& request=entry.second;auto& x=request.angles;std::sort(x.begin(),x.end());x.erase(std::unique(x.begin(),x.end()),x.end());
                    if(x.size()>size_t(maxTheta))continue; // Existing per-group execution respects the explicit cap.
                    std::vector<std::set<double>> covered(8);
                    for(size_t s=0;s<8;++s)for(const auto& item:SamplingCandidates(request.count,seedSet[s],request.depth,request.offset)){
                        const Job& job=*item.second;if(!SameSampling(job,request.count,seedSet[s],request.depth,request.offset))continue;
                        if(request.pairedDepth&&job.pairedDepth!=request.pairedDepth)continue;
                        covered[s].insert(job.inputAngles.begin(),job.inputAngles.end());}
                    std::vector<double> missing;for(double angle:x){bool complete=true;for(const auto& rows:covered)complete=complete&&rows.count(angle);if(!complete)missing.push_back(angle);}
                    if(missing.size()==1){missing.push_back(x.front()==missing.front()?x.back():x.front());std::sort(missing.begin(),missing.end());}
                    if(!missing.empty())DirectStage(request.count,seedSet,8,request.depth,missing,true,false,request.offset,request.pairedDepth);
                    if(pendingJobs.size()>=128)flush();
                }
            }flush();
        }catch(...){collecting=false;pendingJobs.clear();pendingCommands.clear();throw;}
        collecting=false;

    }
    std::string Replay(const Job& job,const std::vector<double>& x){
        std::ostringstream coordinates;for(double v:x)coordinates<<Number(v)<<'\n';const std::string grid=root+"/theta_"+Hash(coordinates.str())+".csv";Atomic(grid,coordinates.str());
        std::vector<std::string> command=base;for(const auto& v:std::vector<std::string>{"--max-reflections",Number(job.depth),"--sobol-seed",Number(job.count),Number(job.seed),"--haar-alpha","--theta-grid-file",grid,"--phi-points",Number(phi),"--threads",Number(workerThreads),"--close"})command.push_back(v);
        if(job.offset){command.push_back("--sobol-offset");command.push_back(Number(job.offset));}if(mirror)command.push_back("--mirror-gamma");
        if(controls){command.push_back("--analytic-mean-threads");command.push_back(Number(analyticMeanThreads));command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(MeanCache("shared"));command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");const auto wf=WriteWeights(x);if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}}
        if(workerChunk>0){command.push_back("--orientation-chunk");command.push_back(Number(workerChunk));}
        command.push_back("--output");command.push_back(job.folder+"/result");std::ostringstream recipe;recipe<<exe<<'\n';for(const auto& v:command)recipe<<v<<'\n';return recipe.str();
    }
    Job LowerView(const Job& source,int depth){
        const std::string identity=Hash(signature+"\npaired view\n"+source.identity+" "+Number(depth));
        const auto resident=cache.find(identity);if(resident!=cache.end())return resident->second;
        Job job=source;job.depth=depth;job.pairedDepth=0;job.assembled=true;job.identity=Hash(signature+"\npaired view\n"+source.identity+" "+Number(depth));job.folder=root+"/jobs/d"+job.identity;
        job.data=job.folder+"/result/result.dat";job.summary=job.folder+"/result/result_analytic_facets.tsv";job.matrix.clear();job.controls.clear();job.angles.clear();job.seconds=source.seconds*.5;
        if(!Exists(job.folder+"/complete.txt")){Directory(job.folder+"/result");const auto name=source.folder+"/result/result_depth"+Number(depth);
            Atomic(job.data,Read(name+".dat"));if(controls)Atomic(job.summary,Read(name+"_analytic_facets.tsv"));Atomic(job.folder+"/command.args",Replay(job,job.inputAngles));
            Atomic(job.folder+"/stdout.log","Validated paired-depth view; standalone replay recipe.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
            const std::string sources=source.folder+"\t"+FileHash(source.folder+"/complete.txt")+"\n";Atomic(job.folder+"/angular_sources.tsv",sources);
            Atomic(job.folder+"/complete.txt",job.identity+" "+Number(job.seconds)+" "+FileHash(job.data)+" "+(controls?FileHash(job.summary):"none")+" "+Hash(sources)+"\n");}
        std::set<std::string> visited;ValidateStoredJob(job.folder,visited);job=Load(job);CacheJob(job);return job;
    }
    std::vector<Job> AssemblePairParts(const std::vector<std::vector<Job>>& parts,int n,const std::vector<double>& x,const std::string& tag){
        if(parts.size()==1)return parts[0];
        const size_t ns=parts[0].size();std::vector<Job> result(ns);const auto wf=WriteWeights(x);
        std::vector<std::array<double,2>> coefficients(x.size(),{{1,1}});if(!wf.empty()){std::istringstream input(Read(wf));std::string line;std::getline(input,line);for(auto& b:coefficients){double angle;input>>angle>>b[0]>>b[1];}}
        for(size_t s=0;s<ns;++s){Job job=parts[0][s];job.count=n;job.offset=0;job.pairedDepth=0;job.assembled=true;std::ostringstream dependencies,key;key<<signature<<tag<<n;
            job.matrix.assign(x.size(),Mueller{});job.controls.assign(controls?x.size():0,std::array<double,5>{});job.angles=x;job.seconds=0;
            for(const auto& group:parts){const auto& p=group[s];const double fraction=double(p.count)/n;key<<p.identity;dependencies<<p.folder<<'\t'<<FileHash(p.folder+"/complete.txt")<<'\n';job.seconds+=p.seconds;
                for(size_t r=0;r<x.size();++r){for(int j=0;j<16;++j)job.matrix[r][j]+=fraction*p.matrix[r][j];if(controls)for(int j=0;j<5;++j)job.controls[r][j]+=fraction*p.controls[r][j];}}
            job.identity=Hash(key.str());job.folder=root+"/jobs/q"+job.identity;job.data=job.folder+"/result/result.dat";job.summary=job.folder+"/result/result_analytic_facets.tsv";Directory(job.folder+"/result");
            if(!Exists(job.folder+"/complete.txt")){WriteMatrix(job.data,x,job.matrix);if(controls){std::ostringstream summary;summary<<std::setprecision(17)<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n";
                    for(size_t r=0;r<x.size();++r){const auto& c=job.controls[r];const auto& b=coefficients[r];summary<<x[r];for(double v:c)summary<<'\t'<<v;summary<<'\t'<<c[0]-b[0]*c[3]-b[1]*c[4]<<'\t'<<job.matrix[r][0]<<'\t'<<job.controlRefinement<<'\t'<<job.controlSetup<<'\t'<<b[0]<<'\t'<<b[1]<<'\n';}Atomic(job.summary,summary.str());}
                Atomic(job.folder+"/command.args",Replay(job,x));Atomic(job.folder+"/stdout.log","Exact paired Sobol-prefix assembly.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");Atomic(job.folder+"/angular_sources.tsv",dependencies.str());
                Atomic(job.folder+"/complete.txt",job.identity+" "+Number(job.seconds)+" "+FileHash(job.data)+" "+(controls?FileHash(job.summary):"none")+" "+Hash(dependencies.str())+"\n");}
            std::set<std::string> visited;ValidateStoredJob(job.folder,visited);job.matrix.clear();job.controls.clear();job.angles.clear();job=Load(job);CacheJob(job);result[s]=job;
        }return result;
    }
    std::pair<std::vector<Job>,std::vector<Job>> SharedStage(int n,const unsigned* seedSet,size_t ns,int upper,int lower,const std::vector<double>& x){
        std::vector<std::vector<Job>> uppers,lowers;const int block=incremental?orientationBlock:n;
        for(int offset=0;offset<n;offset+=block){const int count=std::min(block,n-offset);std::vector<std::set<double>> covered(ns);std::vector<const Job*> candidates;
            for(size_t s=0;s<ns;++s)for(const auto& entry:SamplingCandidates(count,seedSet[s],upper,offset)){
                const auto& job=*entry.second;if(job.pairedDepth!=lower||!SameSampling(job,count,seedSet[s],upper,offset))continue;
                candidates.push_back(&job);covered[s].insert(job.inputAngles.begin(),job.inputAngles.end());}
            for(const auto* job:candidates)LowerView(*job,lower);
            std::vector<double> missing;for(double angle:x){bool complete=true;for(const auto& rows:covered)complete=complete&&rows.count(angle);if(!complete)missing.push_back(angle);}
            if(missing.size()==1){missing.push_back(x.front()==missing.front()?x.back():x.front());std::sort(missing.begin(),missing.end());}
            if(!missing.empty()){auto fresh=DirectStage(count,seedSet,ns,upper,missing,true,false,offset,lower);for(const auto& job:fresh)LowerView(job,lower);}
            auto jobs=Stage(count,seedSet,ns,upper,x,true,false,offset);auto views=Stage(count,seedSet,ns,lower,x,true,false,offset);
            for(auto& job:jobs){job.seconds*=.5;}
            uppers.push_back(jobs);lowers.push_back(views);}
        return {AssemblePairParts(uppers,n,x,"upper paired prefix"),AssemblePairParts(lowers,n,x,"lower paired prefix")};
    }
    void VerifyDepthSharing(const std::vector<double>& x){
        if(!depthSharing||depthVerified||target<=1)return;
        const unsigned seed=7001;const int lower=std::max(1,target-2),n=std::min(pilot,64);double maximum=0;
        try{auto paired=DirectStage(n,&seed,1,target,x,true,false,0,lower);const Job view=LowerView(paired[0],lower);auto standalone=DirectStage(n,&seed,1,lower,x);
            Matrix reference=standalone[0].matrix;if(controls)for(size_t r=0;r<reference.size();++r)reference[r][0]=standalone[0].controls[r][0];maximum=Difference(view.matrix,standalone[0].matrix,reference);depthSharing=maximum<1e-7;
        }catch(const std::exception& error){const std::string message=error.what();const size_t pos=message.find("see ");const std::string log=pos==std::string::npos?message:Read(message.substr(pos+4));if(log.find("paired-depth cache unsupported")==std::string::npos)throw;depthSharing=false;}
        depthVerified=true;Record(root+"/depth_sharing_verification.json","{\"enabled\":"+std::string(depthSharing?"true":"false")+",\"max_all16_difference_scaled_M11\":"+Number(maximum)+",\"probe_count\":"+Number(n)+",\"scope\":\"convex trace prefix with separate coherent diffraction; unsupported shapes use independent traces\"}");
    }
    #include "FullAutoPrefixCache.inc"
    std::vector<Job> PrefixStage(int n,unsigned const* seedSet,size_t ns,int depth,const std::vector<double>& x,bool weighted=true,bool paired=false){
        bool production=weighted&&!paired;
        if(!incremental||!production||n<=orientationBlock)return DirectStage(n,seedSet,ns,depth,x,weighted,paired);
        Require(FileHash(exe)==binaryHash,"binary changed during run");
        const std::string wf=WriteWeights(x),weightHash=Hash(wf.empty()?"":Read(wf));
        std::ostringstream coordinates;for(double q:x)coordinates<<Number(q)<<'\n';
        const std::string gridHash=Hash(coordinates.str()),grid=root+"/theta_"+gridHash+".csv";Atomic(grid,coordinates.str());
        std::vector<Job> result(ns);bool ready=true;
        for(size_t s=0;s<ns;++s){Job j;SetSampling(j,n,seedSet[s],depth,x,false);j.assembled=true;j.weightHash=weightHash;
            j.identity=Hash(signature+ExecutionHash()+"\norientation prefix v1\n"+Number(n)+" "+Number(seedSet[s])+" "+Number(depth)+" "+Number(phi)+" "+Number(kernel.block)+" "+Number(kernel.pipeline)+" "+kernel.warp+" "+Number(mirror)+" "+Number(controls)+"\n"+gridHash+"\n"+weightHash);
            j.folder=root+"/jobs/p"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";
            if(cache.count(j.identity))j=cache.at(j.identity);
            else if(Exists(j.folder+"/complete.txt")){std::set<std::string> visited;ValidateStoredJob(j.folder,visited);
                std::istringstream meta(Read(j.folder+"/complete.txt"));std::string id;meta>>id>>j.seconds;Require(id==j.identity,"prefix cache identity changed");j=Load(j);CacheJob(j);
            }else{ready=false;j.matrix.assign(x.size(),Mueller{});if(controls)j.controls.assign(x.size(),std::array<double,5>{{0,0,0,0,0}});}
            result[s]=j;
        }
        if(ready)return result;
        // Every part is a mean over its own exact contiguous Sobol window.
        // Combining by part_count/total_count preserves the full estimator,
        // including its analytic corrections and every Mueller entry.
        for(auto& j:result){j.matrix.assign(x.size(),Mueller{});if(controls)j.controls.assign(x.size(),std::array<double,5>{{0,0,0,0,0}});j.seconds=0;}
        std::vector<std::ostringstream> sources(ns);std::set<std::string> verified;
        std::vector<std::array<double,2>> coefficients(x.size(),{{1.,1.}});
        if(!wf.empty()){std::istringstream f(Read(wf));std::string header;std::getline(f,header);for(auto& b:coefficients){double theta;Require(bool(f>>theta>>b[0]>>b[1]),"incomplete prefix weights");}}
        const double assemblyStarted=Seconds();double progressAt=-1;size_t cachedBlocks=0;
        for(int start=0;start<n;start+=orientationBlock){const int count=std::min(orientationBlock,n-start);
            const bool cached=AccumulateCachedPrefixBlock(count,start,n,seedSet,ns,depth,x,weightHash,coefficients,result,sources,verified);
            if(cached)++cachedBlocks;
            else{auto parts=Stage(count,seedSet,ns,depth,x,true,false,start);const double fraction=double(count)/n;
            for(size_t s=0;s<ns;++s){Job& j=result[s];const Job& p=parts[s];
                for(size_t r=0;r<x.size();++r){for(int k=0;k<16;++k)j.matrix[r][k]+=fraction*p.matrix[r][k];
                    if(controls)for(int k=0;k<5;++k)j.controls[r][k]+=fraction*p.controls[r][k];}
                j.seconds+=p.seconds;j.controlRefinement=std::max(j.controlRefinement,p.controlRefinement);j.controlSetup+=p.controlSetup;
                sources[s]<<p.folder<<'\t'<<FileHash(p.folder+"/complete.txt")<<'\n';
            }}
            savedOrientationBlocks+=ns;
            // Native completion markers remain durable for every block. This
            // UI snapshot can be reconstructed and must not fsync each cache hit.
            const double now=Seconds();if(progressAt<0||now-progressAt>=2||start+count==n){
                Record(root+"/orientation_progress.json","{\"completed_prefix\":"+Number(start+count)+",\"requested_prefix\":"+Number(n)+",\"depth\":"+Number(depth)+",\"seeds\":"+Number(ns)+",\"block_size\":"+Number(orientationBlock)+",\"theta_rows\":"+Number(x.size())+",\"theta_grid_fnv1a64\":"+Quoted(gridHash)+",\"phase\":\"prefix_assembly\",\"scope\":\"verified saved blocks, per seed; intermediate status throttled to two seconds\"}\n");progressAt=now;}
        }
        for(size_t s=0;s<ns;++s){Job& j=result[s];Directory(j.folder+"/result");WriteMatrix(j.data,x,j.matrix);
            if(controls){std::ostringstream f;f<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n"<<std::setprecision(17);
                for(size_t r=0;r<x.size();++r){const auto& c=j.controls[r];const auto& b=coefficients[r];f<<x[r];for(double v:c)f<<'\t'<<v;
                    f<<'\t'<<c[0]-b[0]*c[3]-b[1]*c[4]<<'\t'<<j.matrix[r][0]<<'\t'<<j.controlRefinement<<'\t'<<j.controlSetup<<'\t'<<b[0]<<'\t'<<b[1]<<'\n';}Atomic(j.summary,f.str());}
            std::vector<std::string> command=base;
            for(const auto& v:std::vector<std::string>{"--max-reflections",std::to_string(depth),"--sobol-seed",std::to_string(n),std::to_string(seedSet[s]),"--haar-alpha","--theta-grid-file",grid,"--phi-points",std::to_string(phi),"--threads",std::to_string(workerThreads),"--close"})command.push_back(v);
            if(mirror)command.push_back("--mirror-gamma");
            if(workerChunk>0){command.push_back("--orientation-chunk");command.push_back(Number(workerChunk));}
            if(kernel.pipeline)command.push_back("--orientation-pipeline");
            if(controls){command.push_back("--analytic-mean-threads");command.push_back(Number(analyticMeanThreads));command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(MeanCache("shared"));command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");}
            if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}command.push_back("--output");command.push_back(j.folder+"/result");
            std::string recipe=exe+"\n";for(const auto& v:command)recipe+=v+"\n";Atomic(j.folder+"/command.args",recipe);
            Atomic(j.folder+"/stdout.log","Controller exact Sobol-prefix assembly; command.args is a full-prefix replay recipe.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
            Atomic(j.folder+"/angular_sources.tsv",sources[s].str());
            Atomic(j.folder+"/orientation_reuse.json","{\"kind\":\"exact_sobol_prefix\",\"count\":"+Number(n)+",\"block\":"+Number(orientationBlock)+",\"seconds_scope\":\"represented full native work, including reused blocks\"}\n");
            Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(controls?FileHash(j.summary):"none")+" "+Hash(sources[s].str())+"\n");
            j.matrix.clear();j.controls.clear();j.angles.clear();j=Load(j);CacheJob(j);
        }
        Log("prefix_assembly_profile","{\"count\":"+Number(n)+",\"depth\":"+Number(depth)+",\"theta_rows\":"+Number(x.size())+",\"seeds\":"+Number(ns)+",\"cached_blocks_without_intermediate_views\":"+Number(cachedBlocks)+",\"elapsed_seconds\":"+Number(Seconds()-assemblyStarted)+",\"native_block_checkpoints_unchanged\":true}");
        return result;
    }
    double ForecastFullGridCost(int n,unsigned seed,int depth,size_t points,double setup,double conservative)const{
        double sx=0,sy=0,sxx=0,sxy=0;size_t count=0;
        for(const auto& entry:SamplingCandidates(n,seed,depth)){const Job& j=*entry.second;
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
    std::vector<Job> Stage(int n,unsigned const* seedSet,size_t ns,int depth,const std::vector<double>& x,bool weighted=true,bool pairedMirror=false,int offset=0){
        if(!angularCache||pairedMirror)return offset?DirectStage(n,seedSet,ns,depth,x,weighted,pairedMirror,offset):PrefixStage(n,seedSet,ns,depth,x,weighted,pairedMirror);
        Require(x.size()<=size_t(maxTheta),"explicit theta evaluation limit exceeded");Require(FileHash(exe)==binaryHash,"binary changed during run");
        const std::string wf=weighted?WriteWeights(x):"",weightHash=Hash(wf.empty()?"":Read(wf));
        typedef std::pair<const Job*,size_t> SourceRow;std::vector<std::map<double,SourceRow>> rows(ns);
        std::vector<Job> exactJobs(ns);
        std::set<double> requested(x.begin(),x.end());bool everySeedHasRows=true,exact=true;
        for(size_t s=0;s<ns;++s){bool exactSeed=false;
            for(const auto& cached:SamplingCandidates(n,seedSet[s],depth,offset)){const Job& j=*cached.second;if(!SameSampling(j,n,seedSet[s],depth,offset))continue;
                if(j.inputAngles==x&&j.weightHash==weightHash){exactSeed=true;exactJobs[s]=j;}
                for(size_t r=0;r<j.inputAngles.size();++r)if(requested.count(j.inputAngles[r]))rows[s].insert({j.inputAngles[r],{&j,r}});
            }
            exact=exact&&exactSeed;everySeedHasRows=everySeedHasRows&&!rows[s].empty();
        }
        // Keep all eight GPU workers available when a complete new sample
        // budget is needed. Only the theta-refinement path assembles rows.
        if(exact)return exactJobs;
        if(!everySeedHasRows)return offset?DirectStage(n,seedSet,ns,depth,x,weighted,pairedMirror,offset):PrefixStage(n,seedSet,ns,depth,x,weighted,pairedMirror);
        std::ostringstream coordinates;for(double q:x)coordinates<<Number(q)<<'\n';const std::string gridHash=Hash(coordinates.str());
        const std::string grid=root+"/theta_"+gridHash+".csv";Atomic(grid,coordinates.str());
        std::vector<Job> result(ns);bool allCached=true;
        for(size_t s=0;s<ns;++s){Job j;SetSampling(j,n,seedSet[s],depth,x,false);j.offset=offset;j.weightHash=weightHash;j.assembled=true;
            j.identity=Hash(signature+ExecutionHash()+"\nangular assembly v2\n"+Number(offset)+" "+Number(n)+" "+Number(seedSet[s])+" "+Number(depth)+" "+Number(phi)+" "+Number(kernel.block)+" "+Number(kernel.pipeline)+" "+kernel.warp+" "+Number(mirror)+" "+Number(controls)+" "+Number(sharedMeans)+"\n"+gridHash+"\n"+weightHash);
            j.folder=root+"/jobs/a"+j.identity;j.data=j.folder+"/result/result.dat";j.summary=j.folder+"/result/result_analytic_facets.tsv";Directory(j.folder);
            if(cache.count(j.identity)){result[s]=cache.at(j.identity);continue;}
            if(Exists(j.folder+"/complete.txt")){
                std::set<std::string> visited;ValidateStoredJob(j.folder,visited);
                std::istringstream meta(Read(j.folder+"/complete.txt"));std::string id;meta>>id>>j.seconds;Require(id==j.identity,"angular cache identity changed");
                j=Load(j);Require(j.matrix.size()==x.size(),"wrong assembled theta count");CacheJob(j);result[s]=j;
            }else{allCached=false;result[s]=j;}
        }
        if(allCached)return result;
        std::vector<double> missing;
        for(double q:x){bool found=true;for(const auto& r:rows)found=found&&r.count(q);if(!found)missing.push_back(q);}
        // Native theta files need two endpoints. A single missing point uses
        // one existing row as an anchor; its cached estimate remains valid.
        if(missing.size()==1){for(double q:x)if(q!=missing[0]){missing.push_back(q);break;}std::sort(missing.begin(),missing.end());}
        std::vector<Job> fresh;if(!missing.empty())fresh=offset?DirectStage(n,seedSet,ns,depth,missing,weighted,false,offset):PrefixStage(n,seedSet,ns,depth,missing,weighted);
        Job means;if(controls&&!sharedMeans){unsigned seed=5003;means=DirectStage(1,&seed,1,1,x,false)[0];}
        std::vector<std::array<double,2>> coefficients(x.size(),{{1.,1.}});
        if(!wf.empty()){std::istringstream file(Read(wf));std::string header;std::getline(file,header);
            for(size_t r=0;r<x.size();++r){double angle;Require(bool(file>>angle>>coefficients[r][0]>>coefficients[r][1]),"invalid assembled control weights");}}
        for(size_t s=0;s<ns;++s){Job& j=result[s];if(cache.count(j.identity))continue;
            const double assemblyStarted=Seconds();
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
            for(const std::string& source:sources){ValidateStoredJob(source,visited);dependencies<<source<<'\t'<<FileHash(source+"/complete.txt")<<'\n';
                const auto cached=folderIndex.find(source);if(cached!=folderIndex.end())j.seconds+=cached->second->seconds;}
            Require(j.seconds>0,"missing angular source cost");const double representedSeconds=j.seconds;
            j.seconds=ForecastFullGridCost(n,seedSet[s],depth,x.size(),j.controlSetup,j.seconds);Directory(j.folder+"/result");WriteMatrix(j.data,x,j.matrix);
            if(controls){std::ostringstream table;table<<"theta_deg\traw_M11\treflection_mean\tshadow_mean\tsampled_reflection\tsampled_shadow\tresidual_mean\thybrid_M11\trefinement\tsetup_seconds\treflection_weight\tshadow_weight\n"<<std::setprecision(17);
                for(size_t r=0;r<x.size();++r){const auto& c=j.controls[r];const auto& b=coefficients[r];table<<x[r];for(double v:c)table<<'\t'<<v;
                    table<<'\t'<<c[0]-b[0]*c[3]-b[1]*c[4]<<'\t'<<j.matrix[r][0]<<'\t'<<j.controlRefinement<<'\t'<<j.controlSetup<<'\t'<<b[0]<<'\t'<<b[1]<<'\n';}Atomic(j.summary,table.str());}
            std::vector<std::string> command=base;
            for(const std::string& v:std::vector<std::string>{"--max-reflections",std::to_string(depth),"--sobol-seed",std::to_string(n),std::to_string(seedSet[s]),"--haar-alpha","--theta-grid-file",grid,"--phi-points",std::to_string(phi),"--threads",std::to_string(workerThreads),"--close"})command.push_back(v);
            if(mirror)command.push_back("--mirror-gamma");
            if(offset){command.push_back("--sobol-offset");command.push_back(std::to_string(offset));}
            if(workerChunk>0){command.push_back("--orientation-chunk");command.push_back(Number(workerChunk));}
            if(kernel.pipeline)command.push_back("--orientation-pipeline");
            if(controls){command.push_back("--analytic-mean-threads");command.push_back(Number(analyticMeanThreads));command.push_back("--analytic-facet-average");command.push_back("--analytic-shadow-control");command.push_back(args.IsCatched("analytic_shadow_control")?args.GetStringValue("analytic_shadow_control"):"circular");command.push_back("--analytic-mean-cache");command.push_back(MeanCache(gridHash));if(sharedMeans){command.push_back("--analytic-mean-reference-grid");command.push_back(root+"/requested_theta.csv");}}
            if(!wf.empty()){command.push_back("--analytic-control-weights");command.push_back(wf);}command.push_back("--output");command.push_back(j.folder+"/result");
            std::string recipe=exe+"\n";for(const auto& v:command)recipe+=v+"\n";Atomic(j.folder+"/command.args",recipe);
            Atomic(j.folder+"/stdout.log","Controller angular assembly; command.args is a reproducible full-grid recipe, not a launched worker.\nHard tree-limit hits: 0\nOrientations still incomplete: 0\n");
            Atomic(j.folder+"/angular_sources.tsv",dependencies.str());
            j.assemblySeconds=Seconds()-assemblyStarted;
            std::ostringstream provenance;provenance<<"{\"controller_assembly_seconds\":"<<Number(j.assemblySeconds)<<",\"kind\":\"angular_assembly\",\"requested_rows\":"<<x.size()<<",\"reused_rows\":"<<reused<<",\"new_worker_rows\":"<<missing.size()<<",\"cost_scope\":\"full-grid forecast from matching native jobs; not assembly wall time\",\"represented_native_seconds\":"<<Number(representedSeconds)<<",\"forecast_full_grid_seconds\":"<<Number(j.seconds)<<",\"M11_rebased_to_authoritative_full_grid_means\":"<<(controls?"true":"false")<<"}\n";
            Atomic(j.folder+"/angular_reuse.json",provenance.str());
            Atomic(j.folder+"/complete.txt",j.identity+" "+Number(j.seconds)+" "+FileHash(j.data)+" "+(controls?FileHash(j.summary):"none")+" "+Hash(dependencies.str())+"\n");
            j.matrix.clear();j.controls.clear();j.angles.clear();j=Load(j);CacheJob(j);++assembledJobs;reusedRows+=reused;newRows+=missing.size();
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
           &&log.find("mean M11 must be positive")==std::string::npos
           &&log.find("analytic mean numerical budget")==std::string::npos
           &&log.find("analytic combined radial kernel exceeds memory budget")==std::string::npos)return false;
        controls=false;sharedMeans=false;weights.clear();polarWeights.clear();SetMirror(false);verifiedMirrorGrids.clear();
        Record(root+"/analytic_control_fallback.json","{\"enabled\":false,\"reason\":\"analytic mean budget/convergence limit; full sampled PO retained\",\"exception\":"+Quoted(message)+",\"diagnostic\":"+Quoted(log)+"}\n");
        std::cout<<"fullauto: analytic mean unavailable; retaining complete sampled PO without facet controls"<<std::endl;return true;
    }
    std::array<Mueller,2> Polar(double angle)const{
        std::array<Mueller,2> result{};if(polarWeights.empty())return result;
        size_t k=std::upper_bound(weightTheta.begin(),weightTheta.end(),angle)-weightTheta.begin();k=k?k-1:0;k=std::min(k,weightTheta.size()-2);
        const double f=(angle-weightTheta[k])/(weightTheta[k+1]-weightTheta[k]);
        for(int c=0;c<2;++c)for(int j=0;j<16;++j)result[c][j]=(1-f)*polarWeights[k][c][j]+f*polarWeights[k+1][c][j];
        return result;
    }
    std::array<double,2> ScalarWeights(double angle)const{std::array<double,2> b{{1,1}};if(weights.empty())return b;
        size_t k=std::upper_bound(weightTheta.begin(),weightTheta.end(),angle)-weightTheta.begin();k=k?k-1:0;k=std::min(k,weightTheta.size()-2);double t=(angle-weightTheta[k])/(weightTheta[k+1]-weightTheta[k]);for(int c=0;c<2;++c)b[c]=(1-t)*weights[k][c]+t*weights[k+1][c];return b;}
    Samples Data(const std::vector<Job>& jobs,bool correctPolar=true){Samples s;for(const auto& job:jobs){Matrix m=job.matrix;
        if(correctPolar&&controls&&!polarWeights.empty())for(size_t r=0;r<m.size();++r){const auto b=Polar(job.angles[r]);const auto& c=job.controls[r];
            for(int j=1;j<16;++j)m[r][j]+=b[0][j]*(c[1]-c[3])+b[1][j]*(c[2]-c[4]);}s.push_back(m);}return s;}
    double Cost(const std::vector<Job>& j,int n){double s=0;for(const auto& v:j)s+=v.seconds;return s/j.size()/n;}
    void SetMirror(bool value){if(mirror!=value)mirrorError=0;mirror=value;statistical=eps-interpolation-depthBudget-(mirror?eps/256:0)-(controls?eps/64:0)-((rotationBeta>1||rotationGamma>1)?eps/256:0)-(selectedCutoff!="requested"&&selectedCutoff!="off"?eps/8:0)-((resourcePolicy=="auto"||config.useGpu)?eps/1024:0);}
    void PrepareMirror(const std::vector<double>& x){
        if(mirrorPolicy=="off")return;
        SetMirror(true);unsigned seed=3001;
        try{Stage(16,&seed,1,std::max(1,target-4),x,false);}
        catch(const std::exception& error){
            const std::string message=error.what();size_t start=message.find("see ");
            bool geometryFailure=start!=std::string::npos&&Read(message.substr(start+4)).find("geometry has no verified body xz reflection")!=std::string::npos;
            if(mirrorPolicy=="on"||!geometryFailure)throw;
            SetMirror(false);Record(root+"/mirror_verification.json","{\"enabled\":false,\"geometry_verified\":false,\"reason\":\"no body xz reflection; full gamma selected\"}\n");
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
        if(!positive){Record(root+"/mirror_verification.json","{\"enabled\":false,\"geometry_verified\":true,\"paired_verified\":false,\"reason\":\"nonpositive corrected mirror pilot M11; full gamma selected\"}\n");
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
        Record(root+"/mirror_verification_"+key+".json",report.str());Record(root+"/mirror_verification.json",report.str());
        std::cout<<"fullauto mirror: bound="<<bound<<" budget="<<eps/256<<" verified="<<pass<<std::endl;
        if(pass)verifiedMirrorGrids.insert(key);
        else if(mirrorPolicy=="on")throw std::runtime_error("fullauto: paired Haar mirror verification failed; use --fullauto-mirror off");
        return pass;
    }
    std::vector<size_t> AuditThetaSupport(std::vector<size_t> knots){
        if(!thetaAuditCount)return knots;
        if(theta.size()>size_t(maxTheta)){
            Require(!args.IsCatched("fullauto_theta_audit"),"explicit dense theta pilot exceeds explicit theta limit");
            Record(root+"/dense_theta_pilot.json","{\"enabled\":false,\"reason\":\"explicit theta limit\"}\n");return knots;
        }
        std::cout<<"fullauto: independent dense theta pilot N="<<thetaAuditCount<<" rows="<<theta.size()<<std::endl;
        auto jobs=DirectStage(thetaAuditCount,angularValidation,8,target,theta);Samples data=Data(jobs);Matrix reference=Mean(data);
        // A finite control-variate pilot is a signed estimator. A negative
        // corrected mean does not make the underlying physical M11 negative.
        // Scale interpolation checks by native positive M11 from the SAME
        // independent jobs; retain the corrected matrices in the residual.
        if(controls)for(size_t row=0;row<reference.size();++row){
            double raw=0;for(const auto& job:jobs)raw+=job.controls.at(row)[0]/8;
            reference[row][0]=raw;
        }
        for(const auto& r:reference)Require(r[0]>0,"mean M11 must be positive in theta pilot");
        double maximum=0;int refinements=0;
        for(;;){Samples support(8),predicted(8);std::vector<double> kx;
            for(size_t k:knots){kx.push_back(theta[k]);for(int seed=0;seed<8;++seed)support[seed].push_back(data[seed][k]);}
            for(int seed=0;seed<8;++seed)predicted[seed]=CubicInterpolate(kx,support[seed],theta);
            Samples delta=Subtract(predicted,data);Matrix mean=Mean(delta),variance=Variance(delta);std::vector<double> errors(theta.size());maximum=0;
            for(size_t r=0;r<theta.size();++r)for(int j=0;j<16;++j){double v=(std::fabs(mean[r][j])+t95*std::sqrt(variance[r][j]/8))/reference[r][0];errors[r]=std::max(errors[r],v);maximum=std::max(maximum,v);}
            std::ostringstream step;step<<"{\"iteration\":"<<refinements<<",\"max_residual95_scaled_M11\":"<<Number(maximum)<<",\"budget\":"<<Number(interpolation*.5)<<",\"support_indices\":[";
            for(size_t i=0;i<knots.size();++i)step<<(i?",":"")<<knots[i];step<<"]}";Log("theta_pilot_refinement",step.str());
            if(maximum<=interpolation*.5||knots.size()==theta.size())break;
            std::set<size_t> next(knots.begin(),knots.end());for(size_t i=1;i<knots.size();++i){size_t worst=knots[i-1];for(size_t r=knots[i-1]+1;r<knots[i];++r)if(errors[r]>errors[worst])worst=r;
                if(errors[worst]>interpolation*.5)next.insert(worst);}
            Require(next.size()>knots.size(),"theta audit cannot refine support");knots.assign(next.begin(),next.end());++refinements;
        }
        std::ostringstream report;report<<"{\"count_per_seed\":"<<thetaAuditCount<<",\"reference_rows\":"<<theta.size()<<",\"support_rows\":"<<knots.size()<<",\"refinements\":"<<refinements<<",\"max_paired_residual95_scaled_M11\":"<<Number(maximum)<<",\"pilot_budget\":"<<Number(interpolation*.5)<<",\"independent_production_seeds\":true}\n";
        Record(root+"/dense_theta_pilot.json",report.str());
        Log("theta_pilot_normalization","{\"native_M11\":"+std::string(controls?"true":"false")+",\"signed_control_residuals_retained\":"+(controls?"true":"false")+"}");return knots;
    }
    #include "FullAutoTuning.inc"
    void Calibrate(const std::vector<double>& x,bool working=false) {
        const int previousPhi=phi;weightTheta=x;weights.clear();polarWeights.clear();phi=fixedPhi?fixedPhi:4;Matrix ref;double fastest=std::numeric_limits<double>::max();Kernel winner=kernel;
        if(config.useGpu&&!working&&!startupKernel){std::vector<Kernel> candidates(5);candidates[1].pipeline=candidates[2].pipeline=candidates[3].pipeline=candidates[4].pipeline=true;candidates[2].block=128;candidates[3].block=256;candidates[4].warp="thread";
            for(auto k:candidates){
                if(fixedPipeline&&!k.pipeline)continue;if(fixedBlock)k.block=kernel.block;if(fixedWarp)k.warp=kernel.warp;
                kernel=k;unsigned seed=4001;auto jobs=Stage(kernelCount,&seed,1,target,x,false);if(ref.empty())ref=jobs[0].matrix;
                Matrix scale=ref;if(controls)for(size_t r=0;r<scale.size();++r)scale[r][0]=jobs[0].controls[r][0];
                double error=Difference(jobs[0].matrix,ref,scale);
                Log("kernel_candidate","{\"pipeline\":"+std::string(k.pipeline?"true":"false")+",\"block\":"+Number(k.block)+",\"warp\":"+Quoted(k.warp)
                    +",\"seconds\":"+Number(jobs[0].seconds)+",\"all16_difference_scaled_M11\":"+Number(error)+",\"budget\":"+Number(eps/1024)
                    +",\"numerically_accepted\":"+(error<=eps/1024?"true":"false")+",\"fixed_pipeline\":"+(fixedPipeline?"true":"false")+",\"fixed_block\":"+(fixedBlock?"true":"false")+",\"fixed_warp\":"+(fixedWarp?"true":"false")+"}");
                const double native=std::max(1e-6,jobs[0].seconds-jobs[0].controlSetup);
                if(error<=eps/1024&&native<fastest){fastest=native;winner=k;}}
            Require(fastest<std::numeric_limits<double>::max(),"no validated GPU kernel");kernel=winner;}
        double best=std::numeric_limits<double>::max();int chosen=4;std::vector<std::array<double,2>> selected;
        std::vector<std::array<Mueller,2>> selectedPolar;std::vector<int> phiCandidates;
        for(int p:{1,2,4,8,16,32,64})if(p>=minPhi&&p<=maxPhi)phiCandidates.push_back(p);
        if(working&&!fixedPhi){
            std::vector<int> nearby;for(int p:phiCandidates)if(p==previousPhi||p==previousPhi/2||p==previousPhi*2||p==minPhi)nearby.push_back(p);
            if(!nearby.empty())phiCandidates=nearby;
        }
        if(phiCandidates.empty())phiCandidates.push_back(std::max(1,minPhi));
        if(fixedPhi)phiCandidates={fixedPhi};
        std::ostringstream candidatesJson;candidatesJson<<'[';bool first=true;
        const int depth=working?target:std::max(1,target-4);
        for(int candidate:phiCandidates){Require(candidate>=minPhi&&candidate<=maxPhi,"phi outside configured limits");phi=candidate;
            // Training fits control coefficients only. With controls disabled
            // its eight native solves were discarded; keep the independent
            // validation scrambles and the complete accuracy contract.
            std::vector<Job> train;if(controls)train=Stage(pilot,training,8,depth,x,false);
            auto valid=Stage(pilot,validation,8,depth,x,false);
            Samples adjusted=Data(valid,false);std::vector<std::array<double,2>> fitted(x.size(),{{0.,0.}});
            std::vector<std::array<Mueller,2>> fittedPolar(x.size());std::array<int,16> accepted{};
            if(controls)for(size_t r=0;r<x.size();++r){
                double cm[2]={};for(int seed=0;seed<8;++seed)for(int c=0;c<2;++c)cm[c]+=train[seed].controls[r][c+3]/8;
                double scale[2]={};for(int seed=0;seed<8;++seed)for(int c=0;c<2;++c)scale[c]+=std::pow(train[seed].controls[r][c+3]-cm[c],2)/7;
                for(double& v:scale)v=std::sqrt(v);
                for(int element=0;element<(polarControls?16:1);++element){
                    double ym=0;for(int seed=0;seed<8;++seed)ym+=(element?train[seed].matrix[r][element]:train[seed].controls[r][0])/8;
                    bool active[2];for(int c=0;c<2;++c)active[c]=scale[c]>std::max(std::fabs(ym),1e-100)*1e-12;
                    double cov=0,cy[2]={};for(int seed=0;seed<8;++seed){double v[2]={};const double y=element?train[seed].matrix[r][element]:train[seed].controls[r][0];
                        for(int c=0;c<2;++c)if(active[c]){v[c]=(train[seed].controls[r][c+3]-cm[c])/scale[c];cy[c]+=v[c]*(y-ym)/7;}
                        cov+=v[0]*v[1]/7;}
                    std::array<double,2> b=element?std::array<double,2>{{0,0}}:std::array<double,2>{{1,1}};
                    if(active[0]&&active[1]){const double det=1.21-cov*cov;b[0]=(1.1*cy[0]-cov*cy[1])/det/scale[0];b[1]=(1.1*cy[1]-cov*cy[0])/det/scale[1];}
                    else {for(int c=0;c<2;++c)if(active[c])b[c]=cy[c]/1.1/scale[c];}
                    for(double& v:b)v=std::max(-4.,std::min(4.,v));
                    double def[8],prop[8],rawValues[8],dm=0,pm=0,rm=0;
                    for(int seed=0;seed<8;++seed){const auto& c=valid[seed].controls[r];const double raw=element?valid[seed].matrix[r][element]:c[0];
                        rawValues[seed]=raw;rm+=raw/8;
                        def[seed]=element?raw:raw+c[1]-c[3]+c[2]-c[4];prop[seed]=raw+b[0]*(c[1]-c[3])+b[1]*(c[2]-c[4]);dm+=def[seed]/8;pm+=prop[seed]/8;}
                    double dv=0,pv=0,rv=0;for(int seed=0;seed<8;++seed){dv+=std::pow(def[seed]-dm,2);pv+=std::pow(prop[seed]-pm,2);rv+=std::pow(rawValues[seed]-rm,2);}
                    const bool useDefault=!element&&dm>0&&dv<.7*rv;
                    const bool accept=(element||pm>0)&&pv<.7*rv&&(!useDefault||pv<dv);
                    if(accept||useDefault){++accepted[element];const auto chosen=accept?b:std::array<double,2>{{1.,1.}};
                        if(!element)fitted[r]=chosen;else for(int c=0;c<2;++c)fittedPolar[r][c][element]=chosen[c];}
                    for(int seed=0;seed<8;++seed)adjusted[seed][r][element]=accept?prop[seed]:(useDefault?def[seed]:rawValues[seed]);
                }
            }
            Statistics stats=Estimate(adjusted);double seconds=0,setup=0;
            for(const auto& job:valid){seconds+=std::max(1e-6,job.seconds-job.controlSetup);setup=std::max(setup,job.controlSetup);}
            const double uncertainty=ratios?std::max(stats.maxInterval,stats.maxRatioInterval):stats.maxInterval;
            // Predict finite work to reach the requested budget. Once the
            // minimum production sample count suffices, lower variance
            // alone must not justify a more expensive azimuth rule.
            const double multiplier=std::max(double(initial)/pilot,std::pow(uncertainty/statistical,2));
            // Cache construction/lock waiting happens once, not N/pilot
            // times. Do not make the first cold phi/kernel candidate lose
            // solely because it prepared authoritative analytic rows.
            const double score=seconds*multiplier+setup;
            Log("phi_candidate","{\"phi\":"+Number(phi)+",\"validation_seconds\":"+Number(seconds)+",\"projected_seconds\":"+Number(score)
                +",\"analytic_setup_once_seconds\":"+Number(setup)+",\"statistical95_scaled_M11\":"+Number(stats.maxInterval)+",\"ratio95\":"+Number(stats.maxRatioInterval)+",\"pilot_count\":"+Number(pilot)
                +",\"depth\":"+Number(depth)+",\"theta_rows\":"+Number(x.size())+",\"fixed_phi\":"+Number(fixedPhi)+",\"working_recheck\":"+(working?"true":"false")+"}");
            candidatesJson<<(first?"":",")<<"{\"phi\":"<<phi<<",\"native_seconds\":"<<Number(seconds)<<",\"analytic_setup_once_seconds\":"<<Number(setup)<<",\"score\":"<<Number(score)<<",\"validation95\":"<<Number(stats.maxInterval)<<",\"ratio95\":"<<Number(stats.maxRatioInterval)<<",\"accepted_control_rows\":[";
            for(int j=0;j<16;++j)candidatesJson<<(j?",":"")<<accepted[j];
            candidatesJson<<"]}";first=false;
            if(score<best){best=score;chosen=phi;selected=fitted;selectedPolar=fittedPolar;}
        }
        phi=chosen;if(controls){weights=selected;polarWeights=selectedPolar;}candidatesJson<<']';
        std::ostringstream polar;polar<<"theta_deg";for(int j=0;j<16;++j)polar<<" R_M"<<j/4+1<<j%4+1<<" S_M"<<j/4+1<<j%4+1;polar<<'\n'<<std::setprecision(17);
        for(size_t r=0;r<x.size();++r){polar<<x[r];for(int j=0;j<16;++j)polar<<' '<<(j?selectedPolar[r][0][j]:selected[r][0])<<' '<<(j?selectedPolar[r][1][j]:selected[r][1]);polar<<'\n';}
        Record(root+"/polarization_weights.tsv",polar.str());const std::string weightPath=WriteWeights(x);
        std::ostringstream c;c<<"{\"selected_phi\":"<<phi<<",\"block\":"<<kernel.block<<",\"pipeline\":"<<(kernel.pipeline?"true":"false")<<",\"warp\":"<<Quoted(kernel.warp)<<",\"pilot_count\":"<<pilot<<",\"score\":"<<Number(best)<<",\"weights\":"<<Quoted(weightPath)<<",\"phi_candidates\":"<<candidatesJson.str()
          <<",\"score_model\":\"finite work to statistical budget, with initial-count floor\",\"calibration_depth\":"<<depth<<",\"calibration_theta\":"<<x.size()<<",\"working_grid_recheck\":"<<(working?"true":"false")<<",\"score_includes_ratios\":"<<(ratios?"true":"false")<<",\"polarization_controls\":"<<(polarControls?"true":"false")
          <<",\"control_training_required\":"<<(controls?"true":"false")<<",\"training_scrambles_per_phi\":"<<(controls?8:0)
          <<",\"validation_scrambles_per_phi\":8,\"independent_training_validation\":true,\"mirror_gamma\":"<<(mirror?"true":"false")<<",\"mirror_policy\":"<<Quoted(mirrorPolicy)<<"}";
        Record(root+"/calibration.json",c.str());
        if(working)TuneTracing(x);
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

namespace { double PhysicalDepthBound(const Samples&,const Samples&,const Matrix&,bool,size_t,bool); }
#include "FullAutoIntegrals.inc"
#include "FullAutoAdaptive.inc"
#include "FullAutoGO.inc"
#include "FullAutoPhysics.inc"

int Run(const ArgPP& args,const RunConfig& config,int argc,const char* argv[],double maximumEdge,double maximumDiameter,double volume) {
    std::unique_ptr<Engine> controller(new Engine(args,config));controller->Init(argc,argv,maximumEdge,maximumDiameter,volume);
    Lock lock(controller->root+"/controller.lock");
    const std::string identity=controller->root+"/configuration.txt";
    if(Exists(identity))Require(Read(identity)==controller->signature+"\n","resume inputs/binary/environment changed; choose a new output directory");else Atomic(identity,controller->signature+"\n");
    // A legacy audit may have committed native blocks but never assembled or
    // inspected its holdout. Preserve that reserved namespace during the
    // verified controller upgrade, only while the production sources agree.
    if(!controller->compatibleHash.empty()&&!Exists(controller->root+"/independent_audit_reservation.tsv")
        &&Exists(controller->root+"/adaptive_checkpoint.tsv")&&Exists(controller->root+"/fullauto_status.json")
        &&Exists(controller->root+"/final_stage_sources.tsv")){
        const auto status=CompactJson(Read(controller->root+"/fullauto_status.json"));
        if(status.find("\"status\":\"validating\"")!=std::string::npos){
            std::istringstream saved(Read(controller->root+"/adaptive_checkpoint.tsv"));int version=0;std::string hash,signature,phase;int round=0,streak=0;unsigned serial=0;
            saved>>version>>hash;std::string header;std::getline(saved,header);std::getline(saved,header);std::istringstream fields(header);fields>>signature>>phase>>round>>streak;
            if(version>=3&&signature==controller->signature&&streak>=controller->stablePasses&&(saved>>serial)&&serial>0
                &&!Exists(controller->root+"/independent_final_validation_serial"+Number(serial-1)+".json")){
                controller->pendingLegacyAuditSerial=serial-1;controller->pendingLegacyAuditRound=round;
                controller->pendingLegacyAuditSources=Hash(Read(controller->root+"/final_stage_sources.tsv"));
            }
        }
    }
    if(!controller->compatibleHash.empty()){
        const std::string previous=controller->root+"/fullauto_parameters.before_compatible_resume.json";
        if(!Exists(previous))Atomic(previous,Read(controller->root+"/fullauto_parameters.json"));
        controller->Record(controller->root+"/binary_compatibility_resume.json",
            "{\"reason\":\"verified native-kernel compatibility; local independent audit and bounded CUDA workspace\",\"origin_executable\":"+Quoted(controller->compatibleBinary)
            +",\"origin_binary_fnv1a64\":"+Quoted(controller->compatibleHash)+",\"executable\":"+Quoted(controller->exe)
            +",\"binary_fnv1a64\":"+Quoted(controller->binaryHash)+"}");
    }
    char originalDirectory[4096];Require(getcwd(originalDirectory,sizeof(originalDirectory))!=nullptr,"cannot capture original working directory");
    std::ostringstream replay;replay<<"{\"schema\":1,\"executable\":"<<Quoted(controller->exe)<<",\"cwd\":"<<Quoted(originalDirectory)<<",\"argv\":[";
    for(int i=1;i<argc;++i)replay<<(i>1?",":"")<<Quoted(argv[i]);replay<<"],\"environment\":{";
    bool firstVariable=true;
    for(char** variable=environ;*variable;++variable){std::string value=*variable;const auto equal=value.find('=');const auto name=value.substr(0,equal);
        if(name.compare(0,4,"MBS_")!=0&&name!="CUDA_VISIBLE_DEVICES"&&name!="LD_LIBRARY_PATH"&&name!="OMP_NUM_THREADS")continue;
        replay<<(firstVariable?"":",")<<Quoted(name)<<':'<<Quoted(value.substr(equal+1));firstVariable=false;}
    replay<<"}}";
    const std::string commandFile=controller->root+"/controller_command.json";
    // Environment iteration order is not stable across exec/reboots. Keep
    // the original JSON; signature validation checks result-affecting values.
    if(!Exists(commandFile))Atomic(commandFile,replay.str());
    if(!controller->compatibleHash.empty())Atomic(controller->root+"/controller_resume_command.json",replay.str());
    signal(SIGTERM,StopWorkers);signal(SIGINT,StopWorkers);
    for(;;){Engine& e=*controller;
        e.Record(e.root+"/explicit_options.json",e.requestedOptions);e.Record(e.root+"/fullauto_parameters.json",e.Parameters("{}"));
        WriteGrid(e.root+"/requested_theta.csv",e.theta);Directory(e.root+"/jobs");
        try {
            PrepareDepthStartup(e);SelectAutomaticDepth(e);e.Record(e.root+"/fullauto_parameters.json",e.Parameters("{}"));
            const int result=e.go?RunGO(e):RunAdaptive(e);
            if(!PromoteAutomaticDepth(e,result))return result;
        }catch(...) {
            if(Exists(e.root+"/fullauto_status.json")){std::string text=Read(e.root+"/fullauto_status.json");const size_t key=text.find("\"status\":\"");
                if(key!=std::string::npos){const size_t begin=key+10,end=text.find('"',begin);text.replace(begin,end-begin,"failed");Atomic(e.root+"/fullauto_status.json",text);}}
            else e.State("failed",0,{},0,0,0,0);
            throw;
        }
        controller.reset(new Engine(args,config));controller->Init(argc,argv,maximumEdge,maximumDiameter,volume);
    }
}
}
