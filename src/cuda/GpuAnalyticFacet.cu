#include "GpuAnalyticFacet.h"
#include <cuda_runtime.h>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <map>
#include <memory>
#include <algorithm>
#include <limits>
#include <sstream>
#include <chrono>
#include <cstdlib>
#include <thread>

namespace
{
const double pi=3.1415926535897932384626433832795;
__device__ double sinc_value(double x)
{
    if(fabs(x)<1e-3){double q=x*x;return 1-q/6+q*q/120-q*q*q/5040;}
    return sin(x)/x;
}
__device__ double reflectance(double mu,double re,double im)
{
    mu=fmin(1.,fmax(0.,mu));double nr=re*re-im*im,ni=2*re*im;
    double a=nr-1+mu*mu,b=ni,length=hypot(a,b);
    double gr,gi;
    if(a>=0){gr=sqrt(.5*(length+a));gi=gr>0 ? b/(2*gr):0;}
    else{gi=copysign(sqrt(.5*(length-a)),b);gr=gi!=0 ? b/(2*gi):0;}
    double rs=((mu-gr)*(mu-gr)+gi*gi)/((mu+gr)*(mu+gr)+gi*gi);
    double pr=nr*mu-gr,ps=ni*mu-gi,qr=nr*mu+gr,qs=ni*mu+gi;
    return .5*(rs+(pr*pr+ps*ps)/(qr*qr+qs*qs));
}
__device__ double dot3(const double*a,const double*b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
__device__ double fourier_norm(const AnalyticGpuPatch &p,const double*xy,double qx,double qy)
{
    xy+=2*p.firstVertex;double maximum=0;
    for(int i=0;i<p.vertices;++i)maximum=fmax(maximum,fabs(qx*xy[2*i]+qy*xy[2*i+1]));
    double sr=0,si=0;
    if(maximum<1e-3)
    {
        double p0=qx*xy[0]+qy*xy[1];
        for(int i=1;i+1<p.vertices;++i)
        {
            double area=.5*((xy[2*i]-xy[0])*(xy[2*(i+1)+1]-xy[1])-(xy[2*i+1]-xy[1])*(xy[2*(i+1)]-xy[0]));
            double p1=qx*xy[2*i]+qy*xy[2*i+1],p2=qx*xy[2*(i+1)]+qy*xy[2*(i+1)+1];
            sr+=area*(1-(p0*p0+p1*p1+p2*p2+p0*p1+p0*p2+p1*p2)/12);
            si+=area*(p0+p1+p2)/3;
        }
        return sr*sr+si*si;
    }
    for(int i=0;i<p.vertices;++i)
    {
        int j=(i+1)%p.vertices;
        double dx=xy[2*j]-xy[2*i],dy=xy[2*j+1]-xy[2*i+1];
        double phase=.5*(qx*(xy[2*i]+xy[2*j])+qy*(xy[2*i+1]+xy[2*j+1]));
        double value=(qx*dy-qy*dx)*sinc_value(.5*(qx*dx+qy*dy));
        double s,c;sincos(phase,&s,&c);sr+=value*c;si+=value*s;
    }
    double q2=qx*qx+qy*qy;return (sr*sr+si*si)/(q2*q2);
}
__global__ void controls_kernel(const AnalyticGpuPatch*patches,int count,const double*xy,
    const double*frames,const double*theta,int rows,int phi,int orientations,int firstCell,int tileCells,
    double nr,double ni,double wave,double radius,int shadowMode,double*reflection,double*shadow)
{
    size_t index=static_cast<size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    size_t total=static_cast<size_t>(tileCells)*orientations;if(index>=total)return;
    int pose=index/tileCells,cell=firstCell+index%tileCells,row=cell%rows,p=cell/rows;
    const double*f=frames+pose*9;const double*u=f+6;double angle=theta[row];
    if(nr==1 && ni==0){reflection[index]=0;shadow[index]=0;return;}
    double st,ct,sp,cp;sincos(angle,&st,&ct);sincos(2*pi*p/phi,&sp,&cp);
    if(angle==0 || angle==pi)st=0;
    double v[3],rq[3],sq[3],k=2*pi/wave,h2=0,c2=0;
    for(int i=0;i<3;++i)
    {
        v[i]=st*(cp*f[i]+sp*f[3+i])-ct*u[i];
        double sum=u[i]+v[i],difference=u[i]-v[i];
        h2+=.25*sum*sum;c2+=.25*difference*difference;rq[i]=k*sum;
    }
    h2=fmin(1.,h2);c2=fmin(1.,c2);
    const double t=dot3(u,v);
    for(int i=0;i<3;++i)sq[i]=k*(v[i]-t*u[i]);
    double r=0,s=0,area=0;
    for(int i=0;i<count;++i)
    {
        const AnalyticGpuPatch &patch=patches[i];double mu=dot3(patch.normal,u);if(mu<=AnalyticFacetAverage::GrazingTolerance)continue;
        area+=mu*patch.area;
        if(h2>0)
        {
            double nu=dot3(patch.normal,v),weight=h2*fmax(0.,c2+mu*nu)*reflectance(mu,nr,ni);
            r+=weight*fourier_norm(patch,xy,dot3(rq,patch.horizontal),dot3(rq,patch.vertical))/(wave*wave);
        }
        if(shadowMode==0 && c2>0 && h2>0)
            s+=c2*c2*mu*mu*fourier_norm(patch,xy,dot3(sq,patch.horizontal),dot3(sq,patch.vertical))/(wave*wave);
    }
    if(shadowMode!=2 && (h2==0 || shadowMode==1))
    {
        double x=k*radius*st,env;
        if(fabs(x)<1e-3){double q=x*x;env=1-q/8+q*q/192-q*q*q/9216;}
        else env=2*j1(x)/x;
        s=c2*c2*env*env*area*area/(wave*wave);
    }
    reflection[index]=r;shadow[index]=s;
}
__global__ void sum_kernel(const double*r,const double*s,const double*weights,int orientations,int cells,double*sumR,double*sumS)
{
    int cell=blockIdx.x*blockDim.x+threadIdx.x;if(cell>=cells)return;double a=0,b=0;
    for(int i=0;i<orientations;++i){size_t j=static_cast<size_t>(i)*cells+cell;a+=weights[i]*r[j];b+=weights[i]*s[j];}
    sumR[cell]=a;sumS[cell]=b;
}
// Preserve the former CPU order: increasing angular cell, i.e. increasing
// azimuth for each theta. One thread owns each orientation/theta accumulator,
// including across tiles, so neither atomics nor a parallel reduction tree
// changes cancellation or rounding. Explicit operations prohibit FMA changes.
__global__ void sample_rows_kernel(const double*r,const double*s,int orientations,
    int rows,int phi,int firstCell,int tileCells,double*sampleR,double*sampleS)
{
    const size_t index=static_cast<size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    if(index>=static_cast<size_t>(orientations)*rows)return;
    const int pose=index/rows,row=index%rows;
    int delta=row-firstCell%rows;if(delta<0)delta+=rows;
    size_t cell=static_cast<size_t>(firstCell)+delta;
    double a=sampleR[index],b=sampleS[index];
    for(;cell<static_cast<size_t>(firstCell)+tileCells;cell+=rows){
        const size_t j=static_cast<size_t>(pose)*tileCells+cell-firstCell;
        a=__dadd_rn(a,__ddiv_rn(r[j],static_cast<double>(phi)));
        b=__dadd_rn(b,__ddiv_rn(s[j],static_cast<double>(phi)));
    }
    sampleR[index]=a;sampleS[index]=b;
}
__device__ double scaled_i0(double x)
{
    x=fabs(x);double term=1,sum=1;
    if(x<50)
    {
        for(int j=1;j<140;++j){term*=x*x/(4*j*j);sum+=term;if(term<sum*1e-16)break;}
        return exp(-x)*sum;
    }
    for(int j=1;j<=12;++j){term*=(2*j-1.)*(2*j-1.)/(8*j*x);sum+=term;}
    return sum/sqrt(2*pi*x);
}
__global__ void gaussian_kernel(const AnalyticAzimuthGaussian::Beam*beams,const int*offset,
    const double*theta,int rows,int phi,int orientations,double*point,double*mean)
{
    const size_t index=static_cast<size_t>(blockIdx.x)*blockDim.x+threadIdx.x;
    const size_t cells=static_cast<size_t>(rows)*phi;if(index>=cells*orientations)return;
    const int pose=index/cells,cell=index%cells,t=cell%rows,p=cell/rows;
    double st,ct,sp,cp;sincos(theta[t],&st,&ct);sincos(2*pi*p/phi,&sp,&cp);
    if(theta[t]==0 || theta[t]==pi)st=0;
    double a=0,m=0;
    for(int i=offset[pose];i<offset[pose+1];++i)
    {
        const auto &b=beams[i];const double x=b.dx-st*cp,y=b.dy-st*sp,z=b.dz+ct;
        const double half=sin((theta[t]-b.theta)/2);
        if(2*b.kappa*half*half>AnalyticAzimuthGaussian::PolarExponentLimit)continue;
        a+=b.peak*exp(-.5*b.kappa*(x*x+y*y+z*z));
        m+=b.peak*exp(-2*b.kappa*half*half)*scaled_i0(b.kappa*st*b.sine);
    }
    point[index]=a;mean[index]=m;
}
template<class T>struct Memory
{
    T*pointer=nullptr;
    size_t capacity=0;
    size_t allocations=0,releases=0;
    cudaError_t status=cudaSuccess;
    ~Memory(){if(pointer)cudaFree(pointer);}
    void Release(){if(pointer){cudaFree(pointer);++releases;}pointer=nullptr;capacity=0;}
    bool Allocate(size_t count)
    {
        status=cudaSuccess;
        if(count<=capacity)return true;
        Release();
        status=cudaMalloc(reinterpret_cast<void**>(&pointer),count*sizeof(T));
        if(status!=cudaSuccess)return false;
        ++allocations;capacity=count;return true;
    }
    bool Load(const T*data,size_t count)
    {return Allocate(count) && (!count || (status=cudaMemcpy(pointer,data,count*sizeof(T),cudaMemcpyHostToDevice))==cudaSuccess);}
};
struct FacetWorkspace
{
    bool pending=false;
    Memory<AnalyticGpuPatch>patch;
    Memory<double>xy,angles,frames,weights,reflection,shadow,sumReflection,sumShadow;
    Memory<double>sampleReflection,sampleShadow;
    std::vector<AnalyticGpuPatch>hostPatches;
    std::vector<double>hostXY,hostTheta;
    size_t reportedCount=0,reportedCells=0,reportedTile=0;
    bool reportedGpuSamples=false;
    size_t Allocations()const{return patch.allocations+xy.allocations+angles.allocations+frames.allocations+weights.allocations
        +reflection.allocations+shadow.allocations+sumReflection.allocations+sumShadow.allocations+sampleReflection.allocations+sampleShadow.allocations;}
    size_t Releases()const{return patch.releases+xy.releases+angles.releases+frames.releases+weights.releases
        +reflection.releases+shadow.releases+sumReflection.releases+sumShadow.releases+sampleReflection.releases+sampleShadow.releases;}
    bool Bind(const AnalyticGpuModel &model,const std::vector<double> &theta)
    {
        const bool same=hostPatches.size()==model.patches.size()
            && (hostPatches.empty() || std::memcmp(hostPatches.data(),model.patches.data(),hostPatches.size()*sizeof(AnalyticGpuPatch))==0)
            && hostXY==model.verticesXY;
        if(!same)
        {
            if(!patch.Load(model.patches.data(),model.patches.size()) || !xy.Load(model.verticesXY.data(),model.verticesXY.size()))return false;
            hostPatches=model.patches;hostXY=model.verticesXY;
        }
        if(hostTheta!=theta)
        {
            if(!angles.Load(theta.data(),theta.size()))return false;
            hostTheta=theta;
        }
        return true;
    }
};
bool enabled(const char*name){const char*value=std::getenv(name);return value && std::strcmp(value,"1")==0;}
struct FacetProfile
{
    typedef std::chrono::steady_clock Clock;
    bool active=enabled("MBS_CUDA_ANALYTIC_PROFILE");
    cudaEvent_t events[5]={};
    double controlsMs=0,sumMs=0,samplesMs=0,copyMs=0,hostSamplesMs=0;
    size_t downloadedBytes=0,tiles=0;
    Clock::time_point started=Clock::now();
    FacetProfile(){if(active)for(auto &event:events)if(cudaEventCreate(&event)!=cudaSuccess){active=false;break;}}
    ~FacetProfile(){for(auto event:events)if(event)cudaEventDestroy(event);}
    void Mark(int i){if(active && cudaEventRecord(events[i])!=cudaSuccess)active=false;}
    void EndTile(){if(!active)return;
        if(cudaEventSynchronize(events[4])!=cudaSuccess){active=false;return;}
        double*values[]={&controlsMs,&sumMs,&samplesMs,&copyMs};
        for(int i=0;i<4;++i){float ms=0;if(cudaEventElapsedTime(&ms,events[i],events[i+1])!=cudaSuccess){active=false;return;}*values[i]+=ms;}
        ++tiles;
    }
    void EndCopy(){if(!active)return;
        if(cudaEventSynchronize(events[4])!=cudaSuccess){active=false;return;}
        float ms=0;if(cudaEventElapsedTime(&ms,events[3],events[4])!=cudaSuccess){active=false;return;}copyMs+=ms;
    }
    void Print(size_t count,size_t rows,int phi,bool gpuSamples,size_t allocations,size_t releases,size_t bytes)const{
        if(!active)return;
        std::printf("CUDA analytic profile: orientations=%zu theta=%zu phi=%d sample_reducer=%s tiles=%zu "
            "controls_ms=%.6f weighted_sum_ms=%.6f sample_sum_ms=%.6f copy_ms=%.6f host_samples_ms=%.6f "
            "wall_ms=%.6f d2h_bytes=%zu allocations=%zu releases=%zu scratch_bytes=%zu\n",count,rows,phi,
            gpuSamples ? "gpu":"cpu",tiles,controlsMs,sumMs,samplesMs,copyMs,hostSamplesMs,
            std::chrono::duration<double,std::milli>(Clock::now()-started).count(),downloadedBytes,allocations,releases,bytes);
    }
};

// Shared launch sequence: synchronous evaluation and the overlap experiment
// execute exactly the same kernels, arguments and orientation reduction.
cudaError_t LaunchFacetTile(FacetWorkspace&cache,const AnalyticGpuModel&model,
    size_t rows,int phi,size_t count,size_t first,size_t length,
    FacetProfile*profile,const char*&operation)
{
    const size_t elements=count*length;
    if(profile)profile->Mark(0);
    controls_kernel<<<(elements+127)/128,128>>>(cache.patch.pointer,model.patches.size(),cache.xy.pointer,
        cache.frames.pointer,cache.angles.pointer,rows,phi,count,first,length,model.indexReal,model.indexImag,
        model.wave,model.circleRadius,model.shadowMode,cache.reflection.pointer,cache.shadow.pointer);
    operation="launch analytic controls";cudaError_t status=cudaGetLastError();
    if(status!=cudaSuccess)return status;
    if(profile)profile->Mark(1);
    sum_kernel<<<(length+127)/128,128>>>(cache.reflection.pointer,cache.shadow.pointer,cache.weights.pointer,
        count,length,cache.sumReflection.pointer,cache.sumShadow.pointer);
    operation="launch analytic sum";status=cudaGetLastError();
    if(status==cudaSuccess&&profile)profile->Mark(2);
    return status;
}
}

struct AnalyticFacetGpuPending::State
{
    FacetWorkspace*workspace=nullptr;
    AnalyticGpuModel model;
    size_t rows=0,count=0,cells=0;
    int phi=0,device=0;
    std::thread::id owner=std::this_thread::get_id();
    cudaEvent_t done=nullptr;
    bool held=false,submitted=false,recorded=false,collected=false;
    std::vector<double>reflection,shadow;
    bool Check(std::string*error)const{
        int current=-1;
        if(owner!=std::this_thread::get_id()||cudaGetDevice(&current)!=cudaSuccess||current!=device){
            if(error)*error="analytic pending operation changed host thread/device";return false;}
        return true;
    }
    void Release(){if(held){workspace->pending=false;held=false;}}
    ~State(){
        // Destruction cannot report a failed ownership check to its caller.
        // Never release a live workspace from another thread/device.
        if(held&&!Check(nullptr))std::terminate();
        if(submitted&&!collected){
            // A launch error may precede recording the completion event.
            if(recorded)cudaEventSynchronize(done);else cudaStreamSynchronize(0);
        }
        Release();if(done)cudaEventDestroy(done);
    }
};

static bool EvaluateAnalyticFacetGpuImpl(const AnalyticGpuModel &model,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&reflected,std::vector<double>&shadow,
    std::vector<AnalyticFacetAverage::Components>*samples,std::string*error,
    AnalyticFacetGpuPending::State*pending=nullptr)
{
    if(error)error->clear();
    FacetProfile profile;
    const auto invalid=[&](const char*message){if(error)*error=message;return false;};
    if(poses.size()!=weights.size() || poses.empty() || theta.empty() || phi<1)
        return invalid("invalid orientation, weight or angular array size");
    const size_t intMax=static_cast<size_t>(std::numeric_limits<int>::max());
    if(poses.size()>intMax || theta.size()>intMax/static_cast<size_t>(phi) || model.patches.size()>intMax
        || poses.size()>std::numeric_limits<size_t>::max()/theta.size()/sizeof(AnalyticFacetAverage::Components))
        return invalid("analytic grid exceeds addressable array dimensions");
    for(double angle:theta)if(!std::isfinite(angle) || angle<0 || angle>pi)return invalid("invalid theta angle");
    for(const auto &pose:poses)
    {
        for(int a=0;a<3;++a)for(int b=a;b<3;++b)
        {
            const auto &u=pose[a],&v=pose[b];const double dot=u.x*v.x+u.y*v.y+u.z*v.z;
            if(!std::isfinite(dot) || std::fabs(dot-(a==b ? 1.:0.))>1e-10)return invalid("non-orthonormal source frame");
        }
        const auto &a=pose[0],&b=pose[1],&u=pose[2];
        const double determinant=(a.y*b.z-a.z*b.y)*u.x+(a.z*b.x-a.x*b.z)*u.y+(a.x*b.y-a.y*b.x)*u.z;
        if(!(determinant>0) || !std::isfinite(determinant))return invalid("source frame is not right handed");
    }
    const size_t count=poses.size(),rows=theta.size(),cells=rows*phi;
    const auto failure=[&](const char*operation,cudaError_t status,size_t requested=0){
        if(error){
            std::ostringstream out;out<<operation<<": "<<cudaGetErrorName(status)<<" ("<<cudaGetErrorString(status)<<")"
                <<"; orientations="<<count<<" theta="<<rows<<" phi="<<phi;
            if(requested)out<<" requested_bytes="<<requested;
            size_t free=0,total=0;
            if(cudaMemGetInfo(&free,&total)==cudaSuccess)out<<" free_bytes="<<free<<" total_bytes="<<total;
            *error=out.str();
        }
        return false;
    };
    std::vector<double>packed(count*9);
    for(size_t i=0;i<count;++i)for(int a=0;a<3;++a)
    {packed[i*9+a*3]=poses[i][a].x;packed[i*9+a*3+1]=poses[i][a].y;packed[i*9+a*3+2]=poses[i][a].z;}
    // Geometry/theta are constant across streamed chunks; retain buffers per
    // thread/device, also avoiding malloc/free synchronization on every chunk.
    int device=0;cudaError_t status=cudaGetDevice(&device);
    if(status!=cudaSuccess)return failure("cudaGetDevice",status);
    static thread_local std::map<int,std::unique_ptr<FacetWorkspace>> workspaces;
    auto &entry=workspaces[device];if(!entry)entry.reset(new FacetWorkspace);
    auto &cache=*entry;
    if(cache.pending)return invalid("analytic workspace already has a pending operation");
    const size_t allocationsBefore=cache.Allocations(),releasesBefore=cache.Releases();
    if(!cache.Bind(model,theta)){
        if(cache.patch.status!=cudaSuccess)return failure("upload facet patches",cache.patch.status);
        if(cache.xy.status!=cudaSuccess)return failure("upload facet vertices",cache.xy.status);
        return failure("upload theta",cache.angles.status);
    }
    auto &frames=cache.frames;auto &w=cache.weights;auto &r=cache.reflection;auto &s=cache.shadow;
    auto &sr=cache.sumReflection;auto &ss=cache.sumShadow;
    auto &samplesR=cache.sampleReflection;auto &samplesS=cache.sampleShadow;
    if(!frames.Load(packed.data(),packed.size()))return failure("upload source frames",frames.status,packed.size()*sizeof(double));
    if(!w.Load(weights.data(),count))return failure("upload orientation weights",w.status,count*sizeof(double));

    // Tile angular cells, never the reduction over orientations. Each cell is
    // evaluated and summed in exactly the same order as the unbounded kernel.
    // The two old count*rows*phi arrays could exceed the entire device memory
    // on the denser independent fullauto audit grid. Keep at most 128 MiB of
    // scratch, and leave room for the resident diffraction/trace workspaces.
    size_t free=0,totalMemory=0;
    status=cudaMemGetInfo(&free,&totalMemory);
    if(status!=cudaSuccess)return failure("cudaMemGetInfo",status);
    const auto scratchBytes=[&](){return (r.capacity+s.capacity+sr.capacity+ss.capacity+samplesR.capacity+samplesS.capacity)*sizeof(double);};
    const size_t retained=scratchBytes();
    const size_t bytesPerCell=2*(count+1)*sizeof(double);
    const size_t budget=std::min<size_t>(128*1024*1024,(free+retained)/2);
    // Keep the sample accumulator inside the same 128 MiB scratch budget.
    // Very large N*theta grids retain the bounded reference CPU path instead
    // of allocating an unbounded GPU accumulator. Its arithmetic is identical.
    const size_t sampleBytes=count*rows*2*sizeof(double);
    const bool gpuSamples=samples && !enabled("MBS_CUDA_ANALYTIC_CPU_SAMPLES") && sampleBytes<=budget/2;
    const size_t sampleElements=gpuSamples ? count*rows:0;
    size_t tile=std::min(cells,std::max<size_t>(1,(budget-(gpuSamples ? sampleBytes:0))/bytesPerCell));
#ifdef MBS_ANALYTIC_GPU_TEST_TILE_CELLS
    tile=std::min<size_t>(tile,MBS_ANALYTIC_GPU_TEST_TILE_CELLS);
#endif
    const auto releaseScratch=[&](){r.Release();s.Release();sr.Release();ss.Release();samplesR.Release();samplesS.Release();};
    const std::array<size_t,6>capacities={{r.capacity,s.capacity,sr.capacity,ss.capacity,samplesR.capacity,samplesS.capacity}};
    tile=AnalyticGpuWorkspace::ReuseTile(budget,capacities,count,tile,sampleElements);
    const std::array<size_t,6>requested={{count*tile,count*tile,tile,tile,sampleElements,sampleElements}};
    if(!AnalyticGpuWorkspace::CanRetain(budget,capacities,requested))releaseScratch();
    for(;;){
        if(r.Allocate(count*tile) && s.Allocate(count*tile) && sr.Allocate(tile) && ss.Allocate(tile)
            && samplesR.Allocate(sampleElements) && samplesS.Allocate(sampleElements))break;
        status=r.status!=cudaSuccess ? r.status:s.status!=cudaSuccess ? s.status:sr.status!=cudaSuccess ? sr.status:
            ss.status!=cudaSuccess ? ss.status:samplesR.status!=cudaSuccess ? samplesR.status:samplesS.status;
        if(status!=cudaErrorMemoryAllocation || tile==1)return failure("allocate analytic tile",status,bytesPerCell*tile);
        // Only allocation failure permits retry. Never hide a kernel error.
        releaseScratch();cudaGetLastError();tile=std::max<size_t>(1,tile/2);
    }
    if(cache.reportedCount!=count || cache.reportedCells!=cells || cache.reportedTile!=tile || cache.reportedGpuSamples!=gpuSamples){
        std::printf("CUDA analytic workspace: orientations=%zu theta=%zu phi=%d tile_cells=%zu/%zu scratch_bytes=%zu sample_reducer=%s\n",
            count,rows,phi,tile,cells,scratchBytes(),gpuSamples ? "gpu":"cpu");
        cache.reportedCount=count;cache.reportedCells=cells;cache.reportedTile=tile;cache.reportedGpuSamples=gpuSamples;
    }
    reflected.resize(cells);shadow.resize(cells);
    if(pending){
        // Preserve the existing tile/capacity choice. A non-single-tile case
        // remains on the established synchronous path, before any kernel.
        if(samples||tile!=cells)return true;
        status=cudaEventCreateWithFlags(&pending->done,cudaEventDisableTiming);
        if(status!=cudaSuccess)return failure("create analytic completion event",status);
        pending->workspace=&cache;pending->model=model;
        pending->rows=rows;pending->count=count;pending->cells=cells;pending->phi=phi;pending->device=device;
        pending->reflection.swap(reflected);pending->shadow.swap(shadow);
        cache.pending=true;pending->held=true;
        return true;
    }
    std::vector<double>hr,hs;
    if(samples){
        samples->assign(count*rows,AnalyticFacetAverage::Components{0,0});
        hr.resize(gpuSamples ? count*rows:count*tile);hs.resize(hr.size());
        if(gpuSamples){
            status=cudaMemset(samplesR.pointer,0,sampleElements*sizeof(double));
            if(status!=cudaSuccess)return failure("clear reflected sample rows",status);
            status=cudaMemset(samplesS.pointer,0,sampleElements*sizeof(double));
            if(status!=cudaSuccess)return failure("clear shadow sample rows",status);
        }
    }
    for(size_t first=0;first<cells;first+=tile){
        const size_t length=std::min(tile,cells-first),elements=count*length;
        const char*operation=nullptr;
        status=LaunchFacetTile(cache,model,rows,phi,count,first,length,&profile,operation);
        if(status!=cudaSuccess)return failure(operation,status);
        if(gpuSamples){
            sample_rows_kernel<<<(sampleElements+127)/128,128>>>(r.pointer,s.pointer,count,rows,phi,first,length,samplesR.pointer,samplesS.pointer);
            status=cudaGetLastError();if(status!=cudaSuccess)return failure("launch analytic sample reduction",status);
        }
        profile.Mark(3);
        status=cudaMemcpy(reflected.data()+first,sr.pointer,length*sizeof(double),cudaMemcpyDeviceToHost);
        if(status!=cudaSuccess)return failure("copy reflected tile",status);
        status=cudaMemcpy(shadow.data()+first,ss.pointer,length*sizeof(double),cudaMemcpyDeviceToHost);
        if(status!=cudaSuccess)return failure("copy shadow tile",status);
        profile.downloadedBytes+=2*length*sizeof(double);
        if(samples && !gpuSamples){
            status=cudaMemcpy(hr.data(),r.pointer,elements*sizeof(double),cudaMemcpyDeviceToHost);
            if(status!=cudaSuccess)return failure("copy reflected orientation samples",status);
            status=cudaMemcpy(hs.data(),s.pointer,elements*sizeof(double),cudaMemcpyDeviceToHost);
            if(status!=cudaSuccess)return failure("copy shadow orientation samples",status);
            profile.downloadedBytes+=2*elements*sizeof(double);
        }
        profile.Mark(4);profile.EndTile();
        if(samples && !gpuSamples){
            const auto start=FacetProfile::Clock::now();
            for(size_t i=0;i<count;++i)for(size_t cell=0;cell<length;++cell){
                const size_t j=i*length+cell,t=(first+cell)%rows;
                (*samples)[i*rows+t].reflection+=hr[j]/phi;(*samples)[i*rows+t].shadow+=hs[j]/phi;
            }
            profile.hostSamplesMs+=std::chrono::duration<double,std::milli>(FacetProfile::Clock::now()-start).count();
        }
    }
    if(gpuSamples){
        profile.Mark(3);
        status=cudaMemcpy(hr.data(),samplesR.pointer,sampleElements*sizeof(double),cudaMemcpyDeviceToHost);
        if(status!=cudaSuccess)return failure("copy reflected sample rows",status);
        status=cudaMemcpy(hs.data(),samplesS.pointer,sampleElements*sizeof(double),cudaMemcpyDeviceToHost);
        if(status!=cudaSuccess)return failure("copy shadow sample rows",status);
        profile.downloadedBytes+=sampleBytes;
        profile.Mark(4);profile.EndCopy();
        const auto start=FacetProfile::Clock::now();
        for(size_t i=0;i<sampleElements;++i){(*samples)[i].reflection=hr[i];(*samples)[i].shadow=hs[i];}
        profile.hostSamplesMs+=std::chrono::duration<double,std::milli>(FacetProfile::Clock::now()-start).count();
    }
    profile.Print(count,rows,phi,gpuSamples,cache.Allocations()-allocationsBefore,cache.Releases()-releasesBefore,scratchBytes());
    return true;
}

bool EvaluateAnalyticFacetGpu(const AnalyticGpuModel&model,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&frames,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&reflection,std::vector<double>&shadow,
    std::vector<AnalyticFacetAverage::Components>*samples,std::string*error)
{
    return EvaluateAnalyticFacetGpuImpl(model,frames,weights,theta,phi,reflection,shadow,samples,error);
}

AnalyticFacetGpuPending::AnalyticFacetGpuPending()=default;
AnalyticFacetGpuPending::~AnalyticFacetGpuPending()=default;
bool AnalyticFacetGpuPending::Prepare(const AnalyticGpuModel&model,
    const std::vector<std::array<AnalyticBackscatter::Vec,3>>&frames,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    bool&supported,std::string*error)
{
    supported=false;
    if(state){if(error)*error="analytic pending operation was already prepared";return false;}
    state.reset(new State);std::vector<double>reflection,shadow;
    if(!EvaluateAnalyticFacetGpuImpl(model,frames,weights,theta,phi,reflection,shadow,nullptr,error,state.get()))return false;
    supported=state->held;return true;
}
bool AnalyticFacetGpuPending::Submit(std::string*error)
{
    if(error)error->clear();
    if(!state||!state->held||state->submitted){if(error)*error="invalid analytic submit state";return false;}
    if(!state->Check(error))return false;
    state->submitted=true;const char*operation=nullptr;
    cudaError_t status=LaunchFacetTile(*state->workspace,state->model,state->rows,state->phi,
        state->count,0,state->cells,nullptr,operation);
    if(status==cudaSuccess){operation="record analytic completion event";status=cudaEventRecord(state->done);state->recorded=status==cudaSuccess;}
    if(status!=cudaSuccess){if(error)*error=std::string(operation)+": "+cudaGetErrorString(status);return false;}
    return true;
}
bool AnalyticFacetGpuPending::Collect(std::vector<double>&reflection,std::vector<double>&shadow,std::string*error)
{
    if(error)error->clear();
    if(!state||!state->held||!state->recorded||state->collected){if(error)*error="invalid analytic collect state";return false;}
    if(!state->Check(error))return false;
    cudaError_t status=cudaEventSynchronize(state->done);
    if(status==cudaSuccess)status=cudaMemcpy(state->reflection.data(),state->workspace->sumReflection.pointer,
        state->cells*sizeof(double),cudaMemcpyDeviceToHost);
    if(status==cudaSuccess)status=cudaMemcpy(state->shadow.data(),state->workspace->sumShadow.pointer,
        state->cells*sizeof(double),cudaMemcpyDeviceToHost);
    if(status!=cudaSuccess){if(error)*error=std::string("collect analytic controls: ")+cudaGetErrorString(status);return false;}
    reflection.swap(state->reflection);shadow.swap(state->shadow);
    state->collected=true;state->Release();return true;
}

bool EvaluateAnalyticAzimuthGaussianGpu(
    const std::vector<std::vector<AnalyticAzimuthGaussian::Beam>>&poses,
    const std::vector<double>&weights,const std::vector<double>&theta,int phi,
    std::vector<double>&point,std::vector<double>&mean,
    std::vector<AnalyticAzimuthGaussian::Values>*samples)
{
    if(poses.empty() || poses.size()!=weights.size() || theta.empty() || phi<1)return false;
    const size_t count=poses.size(),rows=theta.size(),cells=rows*phi,total=count*cells;
    std::vector<AnalyticAzimuthGaussian::Beam> beams;std::vector<int>offset(1,0);
    for(const auto &pose:poses){beams.insert(beams.end(),pose.begin(),pose.end());offset.push_back(beams.size());}
    if(beams.empty())
    {point.assign(cells,0);mean.assign(cells,0);if(samples)samples->assign(count*rows,AnalyticAzimuthGaussian::Values{0,0});return true;}
    Memory<AnalyticAzimuthGaussian::Beam>packed;Memory<int>indices;Memory<double>angles,w,r,s,sr,ss;
    if(!packed.Load(beams.data(),beams.size()) || !indices.Load(offset.data(),offset.size())
        || !angles.Load(theta.data(),rows) || !w.Load(weights.data(),count)
        || !r.Allocate(total) || !s.Allocate(total) || !sr.Allocate(cells) || !ss.Allocate(cells))return false;
    gaussian_kernel<<<(total+127)/128,128>>>(packed.pointer,indices.pointer,angles.pointer,rows,phi,count,r.pointer,s.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    sum_kernel<<<(cells+127)/128,128>>>(r.pointer,s.pointer,w.pointer,count,cells,sr.pointer,ss.pointer);
    if(cudaGetLastError()!=cudaSuccess)return false;
    point.resize(cells);mean.resize(cells);
    if(cudaMemcpy(point.data(),sr.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
        || cudaMemcpy(mean.data(),ss.pointer,cells*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
    if(samples)
    {
        std::vector<double>hr(total),hs(total);
        if(cudaMemcpy(hr.data(),r.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess
            || cudaMemcpy(hs.data(),s.pointer,total*sizeof(double),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
        samples->assign(count*rows,AnalyticAzimuthGaussian::Values{0,0});
        for(size_t i=0;i<count;++i)for(int p=0;p<phi;++p)for(size_t t=0;t<rows;++t)
        {const size_t j=i*cells+p*rows+t;(*samples)[i*rows+t].point+=hr[j]/phi;(*samples)[i*rows+t].mean+=hs[j]/phi;}
    }
    return true;
}
