// Included by FullAuto.cpp. Integral diagnostics use production samples only;
// angular re-evaluations and independent audits must never add statistical weight.
namespace {
const double integralNaN=std::numeric_limits<double>::quiet_NaN();
struct IntegralScalar {double area=0,ext=0;};
struct IntegralSampling {int count=0,offset=0,depth=0;unsigned seed=0;};
struct LegacyPairedIntegral : std::runtime_error {LegacyPairedIntegral():std::runtime_error("legacy paired-depth extinction requires a direct target-depth estimator"){};};
struct IntegralReader {
    std::map<std::string,std::string> files;
    std::map<std::string,IntegralSampling> sampling;
    std::map<std::tuple<std::string,int,bool>,IntegralScalar> values;
    std::set<std::tuple<std::string,int,bool>> active;
    std::set<std::string> verified;
    std::string ReadSource(const std::string& path){const auto text=ReadValidationText(path);files[path]=Hash(text);return text;}
    IntegralSampling Sampling(const std::string& folder){
        const auto found=sampling.find(folder);if(found!=sampling.end())return found->second;
        const auto path=SamplingIndexPath(folder);
        std::istringstream in(DecodeSamplingIndex(path,ReadSource(path)));std::string signature,identity;IntegralSampling s;
        Require(bool(in>>signature>>identity>>s.count>>s.offset>>s.seed>>s.depth)&&s.count>0&&s.offset>=0,"invalid integral sampling index");
        sampling[folder]=s;return s;
    }
    IntegralScalar Prefix(const std::vector<std::string>& folders,int count,int offset,int depth,unsigned seed,bool areaOnly=false){
        // A complete prefix may have been assembled repeatedly on different
        // theta subsets. Pick one deterministic representative, never sum them.
        std::map<int,std::string> windows;
        for(const auto& folder:folders){const auto s=Sampling(folder);if(s.seed!=seed)continue;
            if(s.count==count&&s.offset==offset)return Resolve(folder,depth,areaOnly);
            Require(s.offset>=offset&&s.count<=count&&s.offset<=offset+count-s.count,"integral prefix window outside allocation");
            const auto old=windows.find(s.offset);
            if(old==windows.end()||s.count>Sampling(old->second).count)windows[s.offset]=folder;
        }
        IntegralScalar sum;int next=offset;
        for(const auto& entry:windows){const auto s=Sampling(entry.second);Require(s.offset==next,"overlapping or incomplete integral prefix");
            const auto v=Resolve(entry.second,depth,areaOnly);const double w=double(s.count)/count;sum.area+=w*v.area;sum.ext+=w*v.ext;next+=s.count;}
        Require(next==offset+count,"incomplete integral prefix");return sum;
    }
    IntegralScalar Resolve(const std::string& folder,int depth,bool areaOnly=false){
        const auto key=std::make_tuple(folder,depth,areaOnly);const auto found=values.find(key);if(found!=values.end())return found->second;
        Require(active.insert(key).second,"cycle in integral sources");ValidateStoredJob(folder,verified);
        ReadSource(folder+"/complete.txt");const auto s=Sampling(folder);
        // Projected area depends on the incident orientation and geometry, not
        // on tracing depth. Read it from the native upper result for the SAME
        // orientation window. A legacy lower result has untrusted extinction;
        // requesting area must neither read that scalar nor demand its revision.
        const auto path=folder+"/result/result"+(areaOnly||depth==s.depth?"":"_depth"+Number(depth))+"_integrals.tsv";
        IntegralScalar value;
        if(Exists(folder+"/cutoff_sources.tsv")){
            std::istringstream rows(ReadSource(folder+"/cutoff_sources.tsv"));int targetDepth,sign;std::string source;int terms=0;
            while(rows>>targetDepth>>sign>>source)if(targetDepth==depth){
                Require(sign==1||sign==-1,"invalid cutoff integral sign");
                const auto part=Resolve(source,depth,areaOnly);value.area+=sign*part.area;value.ext+=sign*part.ext;++terms;
            }
            Require(rows.eof()&&terms==3,"incomplete signed cutoff integral recipe");
        }else if(Exists(path)){
            if(!areaOnly&&depth!=s.depth){const auto revision=folder+"/result/result_depth"+Number(depth)+"_integrals_revision.txt";
                if(!Exists(revision))throw LegacyPairedIntegral();
                Require(ReadSource(revision)=="paired_extinction_v1\n","unsupported paired integral revision");}
            std::istringstream input(ReadSource(path));std::string line;Require(bool(std::getline(input,line)),"missing integral header");
            std::vector<std::string> header;std::istringstream columns(line);while(std::getline(columns,line,'\t'))header.push_back(line);
            Require(bool(std::getline(input,line)),"missing integral row");std::map<std::string,std::string> row;std::istringstream cells(line);
            for(const auto& name:header){Require(bool(std::getline(cells,line,'\t')),"incomplete integral row");row[name]=line;}
            auto number=[&](const std::string& name){Require(row.count(name),"missing integral column "+name);size_t end=0;double v=std::stod(row.at(name),&end);Require(end==row.at(name).size(),"invalid integral value");return v;};
            value.ext=areaOnly?0:number("Cext");
            // Older workers did not write projected_area explicitly.
            value.area=row.count("projected_area")?number("projected_area"):number("Cext")/number("Qext");
            Require(std::isfinite(value.area)&&value.area>0,"projected area unavailable in integral source "+path);
        }else{
            Require(Exists(folder+"/angular_sources.tsv"),"missing native integral output "+path);
            std::istringstream rows(ReadSource(folder+"/angular_sources.tsv"));std::vector<std::string> sources;std::string line;
            while(std::getline(rows,line)){const auto tab=line.find('\t');Require(tab!=std::string::npos,"invalid integral provenance");sources.push_back(line.substr(0,tab));}
            value=Prefix(sources,s.count,s.offset,depth,s.seed,areaOnly);
        }
        active.erase(key);values[key]=value;return value;
    }
};
struct IntegralProduction {
    std::array<IntegralScalar,8> samples{};std::string selectedSources;
    bool directTargetFallback=false;
    std::map<std::string,std::string> files;
};
// Keep only the largest existing prefix per level while parsing the source
// list. Angular groups often repeat the same orientation source. Index and
// deduplicate them once, including when legacy extinction needs a second pass.
// Moving to a larger prefix releases the discarded groups immediately.
struct IntegralSourcePlan {
    struct Level {
        int count=0;
        std::map<std::pair<int,unsigned>,std::vector<std::string>> groups;
    };
    std::map<size_t,Level> levels;
    explicit IntegralSourcePlan(const std::string& sources){
        std::istringstream input(sources);std::string line;
        while(std::getline(input,line)){if(line.empty())continue;
            std::istringstream row(line);size_t level;int sign,count;unsigned seed;std::string folder,extra;
            Require(bool(row>>level>>sign>>count>>seed>>folder)&&!(row>>extra)&&count>0&&(sign==1||sign==-1),"invalid integral estimator sources");
            auto& entry=levels[level];
            if(count>entry.count){entry.count=count;entry.groups.clear();}
            if(count==entry.count)entry.groups[{sign,seed}].push_back(std::move(folder));
        }
        Require(!levels.empty()&&levels.begin()->first==0,"missing integral production samples");
        for(size_t level=0;level<levels.size();++level)Require(levels.count(level),"missing integral correction level");
        for(auto& level:levels)for(auto& entry:level.second.groups){
            auto& folders=entry.second;std::sort(folders.begin(),folders.end());
            folders.erase(std::unique(folders.begin(),folders.end()),folders.end());
        }
    }
    const std::vector<std::string>& Folders(size_t level,int sign,unsigned seed)const{
        const auto& groups=levels.at(level).groups;const auto found=groups.find({sign,seed});
        Require(found!=groups.end()&&!found->second.empty(),"missing integral seed/level");return found->second;
    }
};
IntegralProduction IntegralProductionSources(const std::string& sources,const std::vector<unsigned>& seeds,int requestedThreads=1){
    const IntegralSourcePlan plan(sources);Require(seeds.size()==8,"missing integral production samples");
#ifdef _OPENMP
    const int threads=std::max(1,std::min(4,requestedThreads));
#else
    const int threads=1;(void)requestedThreads;
#endif
    auto collect=[&](bool directTarget){
        IntegralProduction result;result.directTargetFallback=directTarget;
        std::array<std::map<std::string,std::string>,8> files;
        std::vector<std::array<std::string,8>> selected(plan.levels.size());
        std::array<std::exception_ptr,8> errors{};
        // Replicas are independent. Preserve the original level/sign/window
        // accumulation order inside each replica, and the published source order.
        // At most four readers coexist; no full graph is copied to every thread.
#ifdef _OPENMP
#pragma omp parallel for num_threads(threads) schedule(dynamic,1)
#endif
        for(int i=0;i<8;++i){try{
            IntegralReader reader;
            for(size_t level=0;level<plan.levels.size();++level){
                if(directTarget&&level!=0&&level!=plan.levels.size()-1)continue;
                std::ostringstream selectedRows;
                for(int sign:{1,-1}){
                    if((level==0||directTarget)&&sign==-1)continue;
                    const auto& folders=plan.Folders(level,sign,seeds[i]);const int count=plan.levels.at(level).count;
                    const int depth=reader.Sampling(folders.front()).depth;
                    const bool areaOnly=directTarget&&level==0&&level!=plan.levels.size()-1;
                    const auto value=reader.Prefix(folders,count,0,depth,seeds[i],areaOnly);
                    if(!directTarget||level==plan.levels.size()-1)result.samples[i].ext+=sign*value.ext;
                    if(level==0)result.samples[i].area=value.area;
                    for(const auto& folder:folders)selectedRows<<level<<' '<<sign<<' '<<count<<' '<<seeds[i]<<' '<<folder<<'\n';
                }
                selected[level][i]=selectedRows.str();
            }
            files[i]=std::move(reader.files);
        }catch(...){errors[i]=std::current_exception();}}
        // A real integrity failure must never be hidden by a legacy fallback.
        bool legacy=false;
        for(const auto& error:errors)if(error){try{std::rethrow_exception(error);}catch(const LegacyPairedIntegral&){legacy=true;}}
        if(legacy)throw LegacyPairedIntegral();
        for(const auto& replica:files)for(const auto& entry:replica){
            const auto old=result.files.find(entry.first);Require(old==result.files.end()||old->second==entry.second,"integral source changed between replicas");
            result.files.insert(entry);
        }
        for(const auto& level:selected)for(const auto& replica:level)result.selectedSources+=replica;
        return result;
    };
    try{return collect(false);}catch(const LegacyPairedIntegral&){return collect(true);}
}
double FullAutoAngularIntegral(const std::vector<double>& theta,const Matrix& matrix,bool go){
    const size_t bins=theta.size();Require(bins>=2&&matrix.size()==bins+(go?2:0),"integral matrix dimensions differ");
    const double rad=std::acos(-1.)/180;double sum=0;
    for(size_t i=0;i<bins;++i){const double lower=(i?.5*(theta[i-1]+theta[i]):theta.front())*rad;
        const double upper=(i+1<bins?.5*(theta[i]+theta[i+1]):theta.back())*rad;
        sum+=matrix[i][0]*(2*std::acos(-1.))*(std::cos(lower)-std::cos(upper));}
    if(go)sum+=matrix[bins][0]+matrix[bins+1][0]; // Pole rows already have cross-section units.
    return sum;
}
struct IntegralInterval {double mean=integralNaN,lower=integralNaN,upper=integralNaN;};
IntegralInterval IntegralConfidence(const std::array<double,8>& values,const std::array<double,8>* denominator=nullptr){
    IntegralInterval result;double mean=0,den=0,var=0,denVar=0,cov=0;
    for(size_t i=0;i<8;++i){if(!std::isfinite(values[i])||(denominator&&!std::isfinite((*denominator)[i])))return result;
        mean+=values[i]/8;if(denominator)den+=(*denominator)[i]/8;}
    for(size_t i=0;i<8;++i){var+=std::pow(values[i]-mean,2)/7;
        if(denominator){denVar+=std::pow((*denominator)[i]-den,2)/7;cov+=(values[i]-mean)*((*denominator)[i]-den)/7;}}
    if(!denominator){const double half=t95*std::sqrt(var/8);result.mean=mean;result.lower=mean-half;result.upper=mean+half;return result;}
    if(!(den>0))return result;result.mean=mean/den;
    const double c=t95*t95/8,a=den*den-c*denVar,b=mean*den-c*cov,d=mean*mean-c*var;
    if(!(a>0))return result;double discriminant=b*b-a*d;
    if(discriminant<0&&discriminant>-1e-12*std::max(std::fabs(b*b),std::fabs(a*d)))discriminant=0;
    if(discriminant<0)return result;
    result.lower=(b-std::sqrt(discriminant))/a;result.upper=(b+std::sqrt(discriminant))/a;return result;
}
const char* const integralProducts[]={"mueller_fullauto_integrals.tsv","fullauto_integrals.log","integral_confidence.tsv","integral_replicas.tsv","integral_sources.fnv","integral_estimator_sources.tsv","integral_summary.json"};
void ValidateIntegralSources(const std::string& root){
    if(!Exists(root+"/integral_sources.fnv"))return;
    std::istringstream rows(Read(root+"/integral_sources.fnv"));std::string line;
    while(std::getline(rows,line)){const auto tab=line.rfind('\t');Require(tab!=std::string::npos,"invalid integral source row");
        Require(FileHash(line.substr(0,tab))==line.substr(tab+1),"final integral source changed");}
    Require(rows.eof(),"invalid final integral source manifest");
}
void WriteFullAutoIntegrals(Engine& e,const std::string& folder,const Samples& output,const std::string& sources,bool final){
    const double begun=Seconds();
    Require(output.size()==8,"integrals need eight production scrambles");Directory(folder);
    auto production=IntegralProductionSources(sources,e.seeds,e.dispatchThreads>0?e.dispatchThreads:e.workerThreads);const Matrix matrix=Mean(output);
    const double rad=std::acos(-1.)/180;std::vector<double> radians,m11;
    for(size_t i=0;i<e.theta.size();++i){radians.push_back(e.theta[i]*rad);m11.push_back(matrix[i][0]);}
    const bool full=e.go||HasFullAngularCoverage(radians),absorption=e.args.IsCatched("abs")||e.config.refractiveImag!=0;
    const double angular=FullAutoAngularIntegral(e.theta,matrix,e.go);
    const double quadrature=e.go?integralNaN:EstimateAngularIntegralRelativeError(radians,m11,angular);
    std::array<double,8> areas,ext,sca,abs,angularSamples;
    std::ostringstream replicas;replicas<<"seed\tprojected_area\tCext\tCsca\tCabs\tangular_scattering_integral\n"<<std::setprecision(17);
    double meanArea=0,meanExt=0;
    for(size_t i=0;i<8;++i){areas[i]=production.samples[i].area;ext[i]=production.samples[i].ext;
        angularSamples[i]=FullAutoAngularIntegral(e.theta,output[i],e.go);meanArea+=areas[i]/8;meanExt+=ext[i]/8;
        auto value=ComputeIntegralCharacteristics(e.go?IntegralMethod::GeometricalOptics:IntegralMethod::PhysicalOptics,
            areas[i],angularSamples[i],ext[i],std::isfinite(ext[i]),absorption,integralNaN,full);
        sca[i]=value.cSca;abs[i]=value.cAbs;
        replicas<<e.seeds[i]<<'\t'<<areas[i]<<'\t'<<ext[i]<<'\t'<<sca[i]<<'\t'<<abs[i]<<'\t'<<angularSamples[i]<<'\n';}
    const auto value=ComputeIntegralCharacteristics(e.go?IntegralMethod::GeometricalOptics:IntegralMethod::PhysicalOptics,
        meanArea,angular,meanExt,std::isfinite(meanExt),absorption,quadrature,full);
    const std::string label=final?"fullauto_final":"fullauto_intermediate";
    WriteIntegralCharacteristicsTsv(folder+"/mueller_fullauto",label,value);
    std::ostringstream confidence;confidence<<"quantity\tmean\tlower95\tupper95\tinterval_method\n"<<std::setprecision(17);
    auto interval=[&](const std::string& name,const std::array<double,8>& values,const std::array<double,8>* den=nullptr){
        const auto ci=IntegralConfidence(values,den);confidence<<name<<'\t'<<ci.mean<<'\t'<<ci.lower<<'\t'<<ci.upper<<'\t'<<(den?"paired_Fieller":"Student_t_7df")<<'\n';};
    interval("projected_area",areas);interval("Cext",ext);interval("Csca",sca);interval("Cabs",abs);
    interval("Qext",ext,&areas);interval("Qsca",sca,&areas);interval("Qabs",abs,&areas);interval("angular_scattering_integral",angularSamples);
    std::string log=FormatIntegralCharacteristicsLog(value,label);
    log+="Integral estimator: eight production scrambles; exact orientation-window weights; "+std::string(production.directTargetFallback?"direct target-depth extinction (legacy paired scalar outputs excluded)":"multilevel extinction corrections")+"; angular integral from the final production Mueller estimator.\n";
    log+="Csca source: "+std::string(absorption?(full?"integrated M11":"unavailable on a partial angular range"):"nonabsorbing balance Csca=Cext")+".\n";
    log+="Integral confidence intervals describe production-scramble variability only. Pointwise fullauto tolerance is not a certified total integral error; quadrature, interpolation, depth and cutoff errors are separate.\n";
    Atomic(folder+"/fullauto_integrals.log",log);Atomic(folder+"/integral_confidence.tsv",confidence.str());Atomic(folder+"/integral_replicas.tsv",replicas.str());
    if(final&&!e.go)for(unsigned seed:e.seeds){const auto path=e.root+"/replicas/production_interpolated_seed"+Number(seed)+".dat";
        if(Exists(path))production.files[path]=FileHash(path);}
    std::ostringstream provenance;for(const auto& entry:production.files)provenance<<entry.first<<'\t'<<entry.second<<'\n';
    Atomic(folder+"/integral_sources.fnv",provenance.str());Atomic(folder+"/integral_estimator_sources.tsv",production.selectedSources);
    e.Record(folder+"/integral_summary.json","{\"schema\":1,\"final\":"+std::string(final?"true":"false")+",\"method\":"+Quoted(e.go?"GO":"PO")
        +",\"full_angular_coverage\":"+(full?"true":"false")+",\"status\":"+Quoted(IntegralCharacteristicsStatus(value,label))
        +",\"legacy_paired_extinction_fallback\":"+(production.directTargetFallback?"true":"false")+",\"extinction_sampling\":"+Quoted(production.directTargetFallback?"direct target depth, largest existing production prefix; area from base prefix":"largest existing production prefix per level; multilevel corrections; angular duplicates excluded")+",\"angular_integral_source\":\"production Mueller estimator including interpolation; GO adds integrated poles\",\"confidence_scope\":\"production variability only, no total-integral accuracy certification\",\"fullauto_tolerance_certifies_integrals\":false}");
    e.Log("integral_summary","{\"folder\":"+Quoted(folder)+",\"final\":"+(final?"true":"false")+",\"status\":"+Quoted(IntegralCharacteristicsStatus(value,label))+",\"assembly_seconds\":"+Number(Seconds()-begun)+",\"source_files\":"+Number(production.files.size())+"}");
    if(final)std::cout<<log<<std::flush;
}
void HashIntegralProducts(const std::string& root,std::ostringstream& hashes){
    for(const char* path:integralProducts)hashes<<path<<' '<<FileHash(root+"/"+path)<<'\n';
}
Matrix ReadIntegralMatrix(const std::string& path,const std::vector<double>& theta){
    std::istringstream input(Read(path));std::string line;std::getline(input,line);Matrix result;
    while(std::getline(input,line)){if(line.empty())continue;std::istringstream row(line);double angle,area;Mueller m;
        Require(bool(row>>angle>>area)&&result.size()<theta.size()&&std::fabs(angle-theta[result.size()])<1e-10,"invalid saved integral replica grid");
        for(double& v:m)Require(bool(row>>v)&&std::isfinite(v),"invalid saved integral replica");result.push_back(m);}
    Require(result.size()==theta.size(),"incomplete saved integral replica");return result;
}
void CompleteSavedPOIntegrals(Engine& e){
    const auto previousHashes=Read(e.root+"/final_outputs.fnv");bool complete=true;
    for(const char* path:integralProducts)complete=complete&&previousHashes.find(std::string(path)+" ")!=std::string::npos;
    if(complete)return;
    Samples replicas;for(unsigned seed:e.seeds)replicas.push_back(ReadIntegralMatrix(e.root+"/replicas/production_interpolated_seed"+Number(seed)+".dat",e.theta));
    const auto saved=ReadIntegralMatrix(e.root+"/mueller_fullauto.dat",e.theta),mean=Mean(replicas);
    Require(Difference(mean,saved,saved)<1e-12,"saved production replicas differ from completed Mueller matrix");
    WriteFullAutoIntegrals(e,e.root,replicas,Read(e.root+"/final_stage_sources.tsv"),true);
    std::ostringstream hashes;std::istringstream old(previousHashes);std::string line;
    while(std::getline(old,line)){bool integral=false;for(const char* path:integralProducts)integral=integral||line.find(std::string(path)+" ")==0;if(!integral)hashes<<line<<'\n';}
    HashIntegralProducts(e.root,hashes);Atomic(e.root+"/final_outputs.fnv",hashes.str());
}
}
