// Included after AdaptiveData. Each receipt is one complete eight-replica
// group, not an accepted-accuracy decision. Frozen inputs include coefficients
// and every sampling/execution choice; source graphs remain authoritative.
std::string InterpolationGroupKey(const Engine& e,const std::vector<int>& depths,const Budgets& counts,
                                  const std::vector<size_t>& ids,const unsigned* seeds){
    std::ostringstream s;s<<std::setprecision(17)<<"interpolation-group-v1\n"<<e.signature<<'\n'<<e.binaryHash<<'\n'<<e.ExecutionHash()<<'\n'
        <<e.kernel.block<<' '<<e.kernel.pipeline<<' '<<e.kernel.warp<<' '<<e.phi<<' '<<e.controls<<' '<<e.sharedMeans<<' '
        <<e.mirror<<' '<<e.depthSharing<<' '<<e.balancedCorrection<<' '<<e.orientationBlock<<' '<<e.ratios<<'\n';
    for(int seed=0;seed<8;++seed)s<<seeds[seed]<<' ';s<<'\n';
    for(int depth:depths)s<<depth<<' ';s<<'\n';
    for(size_t id:ids){
        s<<id<<' '<<e.theta.at(id)<<' '<<e.PhiForTheta(id);
        for(const auto& level:counts)s<<' '<<level.at(id);
        if(e.controls){const auto scalar=e.ScalarWeights(e.theta[id]);const auto polar=e.Polar(e.theta[id]);
            for(double value:scalar)s<<' '<<value;
            for(int j=1;j<16;++j)for(int c=0;c<2;++c)s<<' '<<polar[c][j];}
        s<<'\n';
    }
    return Hash(s.str());
}
std::string InterpolationGroupPath(const Engine& e,const std::string& key){
    return e.root+"/interpolation_groups/"+key+".tsv";
}
std::set<std::string> InterpolationGroupRoots(const std::string& sources){
    std::set<std::string> roots;std::istringstream records(sources);std::string line;
    while(std::getline(records,line)){
        if(line.empty())continue;
        std::istringstream row(line);size_t level;int sign,n;unsigned seed;std::string path,extra;
        Require(bool(row>>level>>sign>>n>>seed>>path)&&!(row>>extra)&&(sign==1||sign==-1)&&n>0,"invalid interpolation group source");
        roots.insert(path);
    }
    Require(!roots.empty(),"interpolation group has no sources");return roots;
}
bool RestoreInterpolationGroup(Engine& e,const std::string& key,size_t levels,size_t rows,AdaptiveData& result){
    const auto path=InterpolationGroupPath(e,key);if(!Exists(path))return false;
    const double begin=Seconds();const auto text=Read(path);const auto split=text.find('\n');
    Require(split!=std::string::npos,"invalid interpolation group checkpoint");
    std::istringstream header(text.substr(0,split));int version;std::string hash,extra;
    Require(bool(header>>version>>hash)&&version==1&&!(header>>extra)&&Hash(text.substr(split+1))==hash,"interpolation group checkpoint changed");
    std::istringstream input(text.substr(split+1));std::string savedKey;size_t savedLevels,savedRows,sourceCount;
    Require(bool(input>>savedKey>>savedLevels>>savedRows)&&savedKey==key&&savedLevels==levels&&savedRows==rows,"interpolation group dimensions or inputs changed");
    size_t sourceBytes;
    Require(bool(input>>sourceBytes)&&sourceBytes>0&&sourceBytes<=text.size()&&input.get()=='\n',"invalid interpolation group sources");
    result.sources.resize(sourceBytes);input.read(&result.sources[0],sourceBytes);
    Require(bool(input)&&bool(input>>sourceCount),"truncated interpolation group sources");
    const auto roots=InterpolationGroupRoots(result.sources);Require(sourceCount==roots.size(),"interpolation group source count changed");
    std::vector<std::string> paths;std::set<std::string> unique;
    for(size_t i=0;i<sourceCount;++i){std::string source,digest;
        Require(bool(input>>source>>digest)&&roots.count(source)&&unique.insert(source).second,"invalid interpolation group source marker");
        Require(FileHash(source+"/complete.txt")==digest,"interpolation group source completion changed");paths.push_back(source);}
    // Reuse never bypasses corruption checks in a descendant. ValidationBatch
    // shares only unchanged dependency graphs across neighbouring groups.
    ValidateStoredRoots(paths,std::max(1,std::min(4,e.workerThreads)));
    result.levels.assign(levels,Samples(8,Matrix(rows)));result.costs.assign(levels,std::vector<double>(rows));result.analytic.resize(rows);
    auto scalar=[&](double& value){Require(bool(input>>value)&&std::isfinite(value),"invalid interpolation group value");};
    for(size_t row=0;row<rows;++row){scalar(result.analytic[row]);
        for(size_t level=0;level<levels;++level){scalar(result.costs[level][row]);
            for(size_t seed=0;seed<8;++seed)for(double& value:result.levels[level][seed][row])scalar(value);}}
    Require(!(input>>extra),"extra interpolation group fields");
    e.Log("interpolation_group_resume","{\"key\":"+Quoted(key)+",\"theta_rows\":"+Number(rows)+",\"levels\":"+Number(levels)
        +",\"replicas\":8,\"source_roots\":"+Number(paths.size())+",\"seconds\":"+Number(Seconds()-begin)
        +",\"source_graphs_verified\":true,\"accuracy_recomputed\":true}");return true;
}
void SaveInterpolationGroup(Engine& e,const std::string& key,const AdaptiveData& result){
    const size_t levels=result.levels.size(),rows=result.analytic.size();
    // A single snapshot remains bounded; exceptionally wide groups fall back
    // to their existing native/prefix checkpoints instead of adding RAM risk.
    const size_t limit=64*1024*1024;
    if(levels&&rows>limit/(levels*8*16*25)){
        e.Log("interpolation_group_checkpoint_skip","{\"reason\":\"snapshot exceeds bounded serialization budget\",\"theta_rows\":"+Number(rows)+"}");return;}
    const auto roots=InterpolationGroupRoots(result.sources);std::ostringstream payload;
    payload<<std::setprecision(17)<<key<<'\n'<<levels<<' '<<rows<<'\n'<<result.sources.size()<<'\n'<<result.sources<<'\n'<<roots.size()<<'\n';
    for(const auto& source:roots)payload<<source<<' '<<FileHash(source+"/complete.txt")<<'\n';
    for(size_t row=0;row<rows;++row){payload<<result.analytic[row];
        for(size_t level=0;level<levels;++level){payload<<' '<<result.costs.at(level).at(row);
            for(size_t seed=0;seed<8;++seed)for(double value:result.levels.at(level).at(seed).at(row))payload<<' '<<value;}payload<<'\n';}
    const auto text=payload.str();if(text.size()>limit){e.Log("interpolation_group_checkpoint_skip","{\"reason\":\"snapshot exceeds bounded serialization budget\"}");return;}
    Directory(e.root+"/interpolation_groups");Atomic(InterpolationGroupPath(e,key),"1 "+Hash(text)+"\n"+text);
    e.Log("interpolation_group_checkpoint","{\"key\":"+Quoted(key)+",\"theta_rows\":"+Number(rows)+",\"levels\":"+Number(levels)
        +",\"replicas\":8,\"bytes\":"+Number(text.size())+",\"source_roots\":"+Number(roots.size())+"}");
}
