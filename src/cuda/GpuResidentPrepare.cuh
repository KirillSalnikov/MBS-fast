// Included after the diffraction kernels: the resident path deliberately reuses
// their polygon integral, Jones transform and Mueller conventions.
#include "GpuConvexQueueMath.cuh"
#include <climits>
#include <map>
#include "geometry/Intersection.h"
namespace {
struct ResidentSummary {
    double energy, extinction;
    unsigned long long candidates, rejected, rejectedJones, rejectedArea, rejectedImportance;
};
struct ResidentParameters {
    double k, wavelength, edgeTolerance, cutoff[3];
    GpuConvexQueue::Vec incidentDirection, incidentBasis;
    GpuReal eps1,eps2,cwr,cwi,icwr,icwi;
    int legacy;
};
__device__ void residentMax(double* address,double value) {
    auto* bits=reinterpret_cast<unsigned long long*>(address);
    unsigned long long old=*bits,assumed;
    do {assumed=old;if(__longlong_as_double(assumed)>=value)return;
        old=atomicCAS(bits,assumed,__double_as_longlong(value));}while(old!=assumed);
}
__global__ void residentReferences(const GpuConvexQueue::Output* raw,int total,
                                    int* counts,double* maxima) {
    const int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=total)return;
    const auto& item=raw[i];atomicAdd(counts+item.orientation,1);
    if(item.beam.facet==INT_MAX)return;
    const double j=GpuConvexQueue::jnorm(item.beam),a=GpuConvexQueue::area(item.beam);
    residentMax(maxima+3*item.orientation,j);residentMax(maxima+3*item.orientation+1,a);
    residentMax(maxima+3*item.orientation+2,j*a);
}
__global__ void residentOffsets(int* counts,int* offsets,int count) {
    if(threadIdx.x||blockIdx.x)return;
    int sum=0;for(int i=0;i<count;++i){offsets[i]=sum;sum+=counts[i];counts[i]=0;}offsets[count]=sum;
}
__global__ void residentPrepare(const GpuConvexQueue::Input* inputs,
    const GpuConvexQueue::Output* raw,int total,ResidentParameters par,
    const double* maxima,const int* offsets,int* kept,GpuBeam* prepared,
    const GpuReal* weights,ResidentSummary* summaries,int* invalid) {
    using namespace GpuConvexQueue;
    const int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=total)return;
    const int orientation=raw[i].orientation;
    GpuConvexQueue::Beam beam=raw[i].beam;const bool shadow=beam.facet==INT_MAX;
    const double cs=maximum(-1.,minimum(1.,dot(mul(par.incidentDirection,-1),beam.direction)));
    Vec basis;
    if(fabs(1-cs)<=DBL_EPSILON)basis=mul(par.incidentBasis,-1);
    else if(fabs(1+cs)<=DBL_EPSILON)basis=par.incidentBasis;
    else {double phi=hypot(beam.direction.x,beam.direction.y)<=64*DBL_EPSILON
            ?0:atan2(beam.direction.y,beam.direction.x);
        if(phi<0)phi+=2*3.14159265358979323846;
        basis={-sin(phi),cos(phi),0};}
    rotateBasis(beam,basis);
    Vec normal={0,0,0};
    for(int v=0;v<beam.count;++v){const Vec a=beam.vertices[v],b=beam.vertices[(v+1)%beam.count];
        normal.x+=(a.y-b.y)*(a.z+b.z);normal.y+=(a.z-b.z)*(a.x+b.x);normal.z+=(a.x-b.x)*(a.y+b.y);}
    normal=unit(normal);
    Vec beamBasis=cross(beam.basis,beam.direction);const double length=sqrt(norm(beamBasis));
    if(!(length>Epsilon)||!isfinite(length)){atomicExch(invalid,1);return;}
    beamBasis={beamBasis.x/length,beamBasis.y/length,beamBasis.z/length};
    const double ar=area(beam),jn=jnorm(beam);
    ResidentSummary& summary=summaries[orientation];
    if(!shadow){
        atomicAdd(&summary.candidates,1ULL);
        const Vec n=inputs[orientation].faces[beam.facet].outward;
        const double projection=fabs(dot(n,beam.direction));
        const double crossSection=projection<1e7*DBL_EPSILON?0:projection*ar/sqrt(norm(n));
        atomicAdd(&summary.energy,crossSection*.5*jn*weights[orientation]);
        const bool smallJ=par.cutoff[0]>0&&jn<par.cutoff[0]*maxima[3*orientation];
        const bool smallA=par.cutoff[1]>0&&ar<par.cutoff[1]*maxima[3*orientation+1];
        const bool smallI=par.cutoff[2]>0&&ar*jn<par.cutoff[2]*maxima[3*orientation+2];
        if(smallJ||smallA||smallI){atomicAdd(&summary.rejected,1ULL);
            if(smallJ)atomicAdd(&summary.rejectedJones,1ULL);
            if(smallA)atomicAdd(&summary.rejectedArea,1ULL);
            if(smallI)atomicAdd(&summary.rejectedImportance,1ULL);return;}
    }
    Vec horizontal,vertical;const double xy=hypot(normal.x,normal.y);
    if(xy<=64*DBL_EPSILON){horizontal={0,-normal.z,0};vertical={1,0,0};}
    else {horizontal={normal.y/xy,-normal.x/xy,0};vertical=cross(normal,horizontal);}
    const Vec cen=center(beam);
    // ChangeCoordinateSystem first projects onto the plane, including for the
    // center; retain these roundoff operations to match CPU preparation.
    const Vec projectedCenter=sub(cen,mul(normal,dot(normal,cen)));
    double x[MaxVertices],y[MaxVertices];
    for(int v=0;v<beam.count;++v){const Vec p=beam.vertices[v==0?beam.count-1:v-1];
        const Vec projected=sub(p,mul(normal,dot(normal,p)));
        x[v]=dot(projected,horizontal)-dot(projectedCenter,horizontal);
        y[v]=dot(projected,vertical)-dot(projectedCenter,vertical);}
    const int slot=offsets[orientation]+atomicAdd(kept+orientation,1);
    GpuBeam& out=prepared[slot];out.nVertices=beam.count;out.nEdgeX=out.nEdgeY=0;
    for(int v=0;v<beam.count;++v){const int next=(v+1)%beam.count;
        const double dx=x[next]-x[v],dy=y[next]-y[v];
        out.x[v]=x[v]*par.k;out.y[v]=y[v]*par.k;
        if(fabs(dx)>=par.edgeTolerance){out.slope_yx[v]=dy/dx;out.edge_valid_x[out.nEdgeX++]=v;}
        if(fabs(dy)>=par.edgeTolerance){out.slope_xy[v]=dx/dy;out.edge_valid_y[out.nEdgeY++]=v;}}
    out.isExternal=shadow;out.orientation=orientation;
    out.bdx=dot(beam.direction,horizontal);out.bdy=dot(beam.direction,vertical);out.bdz=0;
    out.horAx=horizontal.x;out.horAy=horizontal.y;out.horAz=horizontal.z;
    out.verAx=vertical.x;out.verAy=vertical.y;out.verAz=vertical.z;
    out.cenx=cen.x*par.k;out.ceny=cen.y*par.k;out.cenz=cen.z*par.k;out.beam_area=ar;
    Vec nt=cross(normal,beamBasis),np=cross(normal,beam.basis);
    Vec ndt=cross(normal,cross(beam.direction,beamBasis)),ndp=cross(normal,cross(beam.direction,beam.basis));
    out.pNTx=nt.x;out.pNTy=nt.y;out.pNTz=nt.z;out.pNPx=np.x;out.pNPy=np.y;out.pNPz=np.z;
    out.pnxDTx=ndt.x;out.pnxDTy=ndt.y;out.pnxDTz=ndt.z;out.pnxDPx=ndp.x;out.pnxDPy=ndp.y;out.pnxDPz=ndp.z;
    out.raw00r=beam.j[0].re;out.raw00i=beam.j[0].im;out.raw01r=beam.j[1].re;out.raw01i=beam.j[1].im;
    out.raw10r=beam.j[2].re;out.raw10i=beam.j[2].im;out.raw11r=beam.j[3].re;out.raw11i=beam.j[3].im;
    out.phasePath=shadow?beam.path:beam.path+dot(cen,beam.direction);
    out.phaseSign=(shadow||(beam.depth&1))?-1:1;out.absOffset=out.absCount=0;
    phase_raw_jones_scaled(out,par.k,1,nullptr,0,out.jp00r,out.jp00i,out.jp01r,out.jp01i,
        out.jp10r,out.jp10i,out.jp11r,out.jp11i);
}
__global__ void residentExtinction(const GpuBeam* beams,const int* offsets,const int* kept,
    int count,ResidentParameters par,const GpuReal* weights,ResidentSummary* summary) {
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=count)return;
    double forward=0;
    for(int j=offsets[i];j<offsets[i]+kept[i];++j){
        GpuBeam b=beams[j];
        if(b.isExternal){b.jp00r=-b.raw00r;b.jp00i=-b.raw00i;b.jp01r=-b.raw01r;b.jp01i=-b.raw01i;
            b.jp10r=-b.raw10r;b.jp10i=-b.raw10i;b.jp11r=-b.raw11r;b.jp11i=-b.raw11i;}
        GpuReal a,b0,c,d,e,f,g,h;
        if(compute_beam_jones_context_gpu<true>(b,1,0,0,1,0,0,-1,0,1,0,-1,0,0,
            1,1,par.eps1,par.eps2,par.cwr,par.cwi,par.icwr,par.icwi,par.legacy,1,a,b0,c,d,e,f,g,h))
            forward+=b0+h;
    }
    // GPU Jones amplitudes omit the CPU factor 1/2 (Mueller applies 1/4).
    summary[i].extinction=.5*par.wavelength*forward*weights[i];
}
struct ResidentPool {
    struct Slot {void* pointer=nullptr;size_t bytes=0;};
    std::vector<Slot> slots;
    int device=0;
    explicit ResidentPool(int d):device(d){}
    ~ResidentPool(){int previous=0;cudaGetDevice(&previous);cudaSetDevice(device);
        for(const auto& slot:slots)cudaFree(slot.pointer);cudaSetDevice(previous);}
};
struct ResidentAllocations {
    std::vector<void*> pointers;
    ResidentPool* pool=nullptr;size_t next=0;
    ResidentAllocations(){
        const char* enabled=std::getenv("MBS_GPU_RESIDENT_REUSE");
        const char* hybrid=std::getenv("MBS_GPU_TRACE_HYBRID");
        const bool reuse=enabled?(enabled[0]=='1'&&enabled[1]=='\0')
            :(hybrid&&hybrid[0]=='1'&&hybrid[1]=='\0');
        if(reuse){
            int device=0;if(cudaGetDevice(&device)!=cudaSuccess)return;
            static thread_local std::map<int,std::unique_ptr<ResidentPool>> pools;
            auto& item=pools[device];if(!item)item.reset(new ResidentPool(device));pool=item.get();
        }
    }
    ~ResidentAllocations(){cudaStreamSynchronize(0);for(void* p:pointers)cudaFree(p);}
    template<typename T> bool make(T*& p,size_t n,bool zero=false){
        p=nullptr;const size_t bytes=std::max(size_t(1),n)*sizeof(T);
        if(pool){
            if(next==pool->slots.size())pool->slots.emplace_back();
            auto& slot=pool->slots[next++];
            if(slot.bytes<bytes){
                // Free the superseded buffer before growth. Prior calls have
                // completed their D2H copies; no live batch references it.
                cudaFree(slot.pointer);slot.pointer=nullptr;slot.bytes=0;
                if(cudaMalloc(&slot.pointer,bytes)!=cudaSuccess)return false;
                slot.bytes=bytes;
            }
            p=static_cast<T*>(slot.pointer);
        }else{
            if(cudaMalloc(&p,bytes)!=cudaSuccess)return false;pointers.push_back(p);
        }
        return !zero||cudaMemset(p,0,n*sizeof(T))==cudaSuccess;
    }
    size_t cachedBytes() const {size_t n=0;if(pool)for(const auto& slot:pool->slots)n+=slot.bytes;return n;}
};
} // namespace

bool HandlerPO::CanHandleResidentConvexGpu() const {
#ifdef MBS_GPU_FLOAT
    // This prototype is validated in FP64 only. FP32 keeps its existing path.
    return false;
#endif
    // Explicitly decline modes whose CPU preparation has additional semantics.
    // Absorbing/nonconvex particles are also rejected by PrepareGpuInput.
    return m_gpuEnabled && isCoh && !useKarczewski && !IsFftEnabled()
        && !HasAbsorptionAccounting() && !isBackScatteringConusEnabled
        && !outputJones && m_otFarReferencePath==0 && m_otPingDistance==0 && !m_otPhaseAverage;
}
bool HandlerPO::HandleResidentConvexGpu(const GpuConvexQueue::Input* inputs,
    const GpuConvexQueue::Output* raw,size_t beamCount,const std::vector<double>& weights,
    Arr2D& full,Arr2D& noShadow,std::vector<double>& energies,std::vector<double>& extinction) {
    if(!CanHandleResidentConvexGpu()||beamCount>INT_MAX)return false;
    const int count=weights.size(),nAz=m_sphere.nAzimuth,nZen=m_sphere.nZenith;
    const size_t gridCount=size_t(nAz)*(nZen+1),mCount=gridCount*16;
    if(!count||!gridCount)return false;
    const auto start=std::chrono::steady_clock::now();
    ResidentAllocations alloc;
    GpuBeam* prepared=nullptr;int *kept=nullptr,*offsets=nullptr,*invalid=nullptr;
    double* maxima=nullptr;ResidentSummary* summary=nullptr;
    GpuReal *w=nullptr,*st=nullptr,*ct=nullptr,*sp=nullptr,*cp=nullptr,*vf=nullptr,*matrix=nullptr,*matrixNS=nullptr;
    if(!alloc.make(prepared,beamCount)||!alloc.make(kept,count,true)||!alloc.make(offsets,count+1)
        ||!alloc.make(invalid,1,true)||!alloc.make(maxima,3*count,true)||!alloc.make(summary,count,true)
        ||!alloc.make(w,count)||!alloc.make(st,nZen+1)||!alloc.make(ct,nZen+1)
        ||!alloc.make(sp,nAz)||!alloc.make(cp,nAz)||!alloc.make(vf,gridCount*3)
        ||!alloc.make(matrix,mCount,true)||(ComputeNoShadow()&&!alloc.make(matrixNS,mCount,true)))return false;
    const auto copy=[](void* dst,const void* src,size_t bytes){return cudaMemcpy(dst,src,bytes,cudaMemcpyHostToDevice)==cudaSuccess;};
    std::vector<GpuReal> hw(weights.begin(),weights.end()),hst(nZen+1),hct(nZen+1),hsp(nAz),hcp(nAz),hvf(gridCount*3);
    for(int t=0;t<=nZen;++t){const GpuReal theta=gpu_theta(m_sphere,t);hst[t]=sin(theta);hct[t]=cos(theta);}
    for(int p=0;p<nAz;++p){const GpuReal phi=p*m_sphere.azinuthStep;hsp[p]=sin(phi);hcp[p]=cos(phi);
        for(int t=0;t<=nZen;++t){const size_t g=size_t(p)*(nZen+1)+t;const auto& v=m_sphere.vf[p][t];
            hvf[3*g]=v.x;hvf[3*g+1]=v.y;hvf[3*g+2]=v.z;}}
    if(!copy(w,hw.data(),count*sizeof(GpuReal))||!copy(st,hst.data(),hst.size()*sizeof(GpuReal))
        ||!copy(ct,hct.data(),hct.size()*sizeof(GpuReal))||!copy(sp,hsp.data(),hsp.size()*sizeof(GpuReal))
        ||!copy(cp,hcp.data(),hcp.size()*sizeof(GpuReal))||!copy(vf,hvf.data(),hvf.size()*sizeof(GpuReal)))return false;
    ResidentParameters par={};par.k=m_waveIndex;par.wavelength=m_wavelength;
    par.edgeTolerance=geometry_length_tolerance(m_particle->MaximalDimention());
    par.cutoff[0]=m_beamCutoffJRel>=0?m_beamCutoffJRel:m_targetEps;
    par.cutoff[1]=m_beamCutoffAreaRel>=0?m_beamCutoffAreaRel:m_targetEps;
    par.cutoff[2]=m_beamCutoffImportanceRel;
    const auto& d=m_incidentLight->direction;const auto& b=m_incidentLight->polarizationBasis;
    par.incidentDirection={d.cx,d.cy,d.cz};par.incidentBasis={b.cx,b.cy,b.cz};
    par.eps1=gpu_effective_eps1(m_eps1);par.eps2=m_eps2;par.legacy=m_legacySign;
    par.cwr=real(m_complWave);par.cwi=imag(m_complWave);par.icwr=real(m_invComplWave);par.icwi=imag(m_invComplWave);
    if(beamCount)residentReferences<<<(beamCount+127)/128,128>>>(raw,beamCount,kept,maxima);
    residentOffsets<<<1,1>>>(kept,offsets,count);
    if(beamCount)residentPrepare<<<(beamCount+127)/128,128>>>(inputs,raw,beamCount,par,maxima,offsets,kept,prepared,w,summary,invalid);
    if(!gpu_report_cuda_error(cudaGetLastError(),"resident prepare launch")
        ||!gpu_report_cuda_error(cudaStreamSynchronize(0),"resident prepare sync"))return false;
    const double prepareSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    int bad=0;if(cudaMemcpy(&bad,invalid,sizeof(int),cudaMemcpyDeviceToHost)!=cudaSuccess||bad)return false;
    residentExtinction<<<(count+63)/64,64>>>(prepared,offsets,kept,count,par,w,summary);
    const int block=gpu_block_size();
    diffraction_grid_mueller_kernel<<<(count*gridCount+block-1)/block,block>>>(
        prepared,offsets,st,ct,sp,cp,vf,w,nAz,nZen,count,1,1,par.eps1,par.eps2,
        par.cwr,par.cwi,par.icwr,par.icwi,par.legacy,matrix,matrixNS,kept);
    if(!gpu_report_cuda_error(cudaGetLastError(),"resident diffraction launch"))return false;
    std::vector<GpuReal> hm(mCount),hns(ComputeNoShadow()?mCount:0);
    std::vector<ResidentSummary> hs(count);
    if(cudaMemcpy(hm.data(),matrix,mCount*sizeof(GpuReal),cudaMemcpyDeviceToHost)!=cudaSuccess
        ||(ComputeNoShadow()&&cudaMemcpy(hns.data(),matrixNS,mCount*sizeof(GpuReal),cudaMemcpyDeviceToHost)!=cudaSuccess)
        ||cudaMemcpy(hs.data(),summary,count*sizeof(ResidentSummary),cudaMemcpyDeviceToHost)!=cudaSuccess)return false;
    for(auto v:hm)if(!std::isfinite(v))return false;
    for(auto v:hns)if(!std::isfinite(v))return false;
    for(auto s:hs)if(!std::isfinite(s.energy)||!std::isfinite(s.extinction))return false;
    // Commit only after every kernel, copy and finite-value check succeeded.
    energies.resize(count);extinction.resize(count);unsigned long long candidates=0,rejected=0;
    for(int i=0;i<count;++i){const auto& s=hs[i];energies[i]=s.energy;extinction[i]=s.extinction;
        candidates+=s.candidates;rejected+=s.rejected;
        m_beamCutoffStatistics->candidates.fetch_add(s.candidates,std::memory_order_relaxed);
        m_beamCutoffStatistics->kept.fetch_add(s.candidates-s.rejected,std::memory_order_relaxed);
        m_beamCutoffStatistics->rejected.fetch_add(s.rejected,std::memory_order_relaxed);
        m_beamCutoffStatistics->rejectedJones.fetch_add(s.rejectedJones,std::memory_order_relaxed);
        m_beamCutoffStatistics->rejectedArea.fetch_add(s.rejectedArea,std::memory_order_relaxed);
        m_beamCutoffStatistics->rejectedImportance.fetch_add(s.rejectedImportance,std::memory_order_relaxed);}
    m_beamCutoffStatistics->prepareCalls.fetch_add(count,std::memory_order_relaxed);
    for(int p=0;p<nAz;++p)for(int t=0;t<=nZen;++t){size_t g=size_t(p)*(nZen+1)+t;
        double* cell=full.RawCell(p,t);for(int k=0;k<16;++k)cell[k]+=hm[g*16+k];
        if(ComputeNoShadow()){double* ns=noShadow.RawCell(p,t);for(int k=0;k<16;++k)ns[k]+=hns[g*16+k];}}
    std::fprintf(stderr,"GPU resident direct: orientations=%d raw_beams=%zu retained=%llu prepare_s=%.6f diffraction_download_s=%.6f raw_download_bytes=0 prepared_MiB=%.3f cached_MiB=%.3f\n",
        count,beamCount,candidates-rejected,prepareSeconds,
        std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()-prepareSeconds,
        beamCount*sizeof(GpuBeam)/1048576.,alloc.cachedBytes()/1048576.);
    return true;
}
