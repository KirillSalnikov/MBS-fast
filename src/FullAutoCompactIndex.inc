// Included inside Engine. The journal is an acceleration index, never the
// authority for an optical result. Each referenced complete marker and source
// graph is verified before a restored payload can enter an estimator.
bool compactIndex=false,restoringIndex=false,indexJournalRestored=false;
std::string SamplingText(const Job& j)const{
    std::ostringstream out;
    out<<signature<<' '<<j.identity<<' '<<j.count<<' '<<j.offset<<' '<<j.seed<<' '<<j.depth<<' '<<j.phi<<' '<<j.block
       <<' '<<j.pipeline<<' '<<j.mirrored<<' '<<j.paired<<' '<<j.assembled<<' '<<j.hasControls<<' '<<j.pairedDepth
       <<' '<<j.warp<<' '<<j.weightHash<<' '<<j.executionHash<<'\n';
    for(double angle:j.inputAngles)out<<Number(angle)<<'\n';
    return out.str();
}
bool ParseSampling(const std::string& text,const std::string& folder,Job& j,bool currentControls=true)const{
    std::istringstream input(text);std::string configuration;
    Require(bool(input>>configuration>>j.identity>>j.count>>j.offset>>j.seed>>j.depth>>j.phi>>j.block
        >>j.pipeline>>j.mirrored>>j.paired>>j.assembled>>j.hasControls>>j.pairedDepth>>j.warp>>j.weightHash>>j.executionHash),"invalid sampling index");
    if(configuration!=signature||(currentControls&&j.hasControls!=controls))return false;
    double angle;while(input>>angle){Require(std::isfinite(angle),"nonfinite indexed angle");j.inputAngles.push_back(angle);}
    Require(input.eof()&&j.inputAngles.size()>=2&&j.count>0&&j.offset>=0,"invalid indexed grid");
    j.folder=folder;j.data=folder+"/result/result.dat";j.summary=folder+"/result/result_analytic_facets.tsv";
    return true;
}
void ConfigureIndex(bool fresh){
    const auto path=root+"/job_index_format.txt";
    if(Exists(path)){Require(Read(path)=="compact_v1\n","unsupported job index format");compactIndex=true;}
    else if(fresh){Atomic(path,"compact_v1\n");compactIndex=true;}
}
std::string DiagnosticText(const Job& j)const{
    std::ostringstream d;d<<std::setprecision(17)
        <<j.seconds<<' '<<j.controlRefinement<<' '<<j.controlSetup<<' '<<j.traceSeconds<<' '<<j.diffractionSeconds<<' '
        <<j.processingSeconds<<' '<<j.preparedWaitSeconds<<' '<<j.producerWaitSeconds<<' '<<j.assemblySeconds<<' '<<j.preparationThreads<<'\n';
    return d.str();
}
void ParseDiagnostics(const std::string& text,Job& j)const{
    std::istringstream d(text);std::string extra;
    Require(bool(d>>j.seconds>>j.controlRefinement>>j.controlSetup>>j.traceSeconds>>j.diffractionSeconds
        >>j.processingSeconds>>j.preparedWaitSeconds>>j.producerWaitSeconds>>j.assemblySeconds>>j.preparationThreads)
        &&!(d>>extra)&&j.seconds>0&&std::isfinite(j.seconds),"invalid journal diagnostics");
    for(double value:{j.controlRefinement,j.controlSetup,j.traceSeconds,j.diffractionSeconds,j.processingSeconds,
        j.preparedWaitSeconds,j.producerWaitSeconds,j.assemblySeconds})Require(std::isfinite(value),"nonfinite journal diagnostic");
}
void ReadCompactDiagnostics(const std::string& marker,const std::string& index,Job& j)const{
    const auto start=marker.find('\n',marker.find('\n')+1)+1+index.size();
    Require(start<marker.size(),"missing compact diagnostics");
    const auto end=marker.find('\n',start);Require(end!=std::string::npos,"incomplete compact diagnostics");
    std::istringstream header(marker.substr(start,end-start));std::string magic,hash,extra;
    Require(bool(header>>magic>>hash)&&!(header>>extra)&&magic=="MBS_DIAGNOSTICS_V1","invalid compact diagnostics header");
    const auto diagnostics=marker.substr(end+1);Require(Hash(diagnostics)==hash,"compact diagnostics changed");
    ParseDiagnostics(diagnostics,j);
}
void PublishComplete(const Job& j,const std::string& metadata){
    if(!compactIndex){Atomic(j.folder+"/complete.txt",metadata);return;}
    const auto index=SamplingText(j),diagnostics=DiagnosticText(j);
    Atomic(j.folder+"/complete.txt",metadata+"MBS_SAMPLING_V1 "+Hash(index)+" "+Number(index.size())+"\n"+index
        +"MBS_DIAGNOSTICS_V1 "+Hash(diagnostics)+"\n"+diagnostics);
}
void AppendIndex(const Job& j){
    if(!compactIndex||restoringIndex)return;
    std::ostringstream body;body<<j.folder.substr(root.size()+6)<<'\n'<<FileHash(j.folder+"/complete.txt")<<'\n'
        <<DiagnosticText(j)<<SamplingText(j);
    const auto data=body.str(),record=Number(data.size())+" "+Hash(data)+"\n"+data+"\n";
    const auto path=root+"/completed_jobs.journal";int fd=open(path.c_str(),O_WRONLY|O_CREAT|O_APPEND,0644);
    Require(fd>=0,"cannot append job index journal");size_t done=0;
    while(done<record.size()){ssize_t n=write(fd,record.data()+done,record.size()-done);if(n<0&&errno==EINTR)continue;
        if(n<=0){close(fd);Require(false,"job index journal write failed");}done+=size_t(n);}
    const bool okay=fsync(fd)==0;close(fd);Require(okay,"job index journal flush failed");
}
size_t RestoreIndexJournal(){
    if(!compactIndex||indexJournalRestored)return 0;
    const auto path=root+"/completed_jobs.journal";std::ifstream in(path.c_str(),std::ios::binary);
    if(!in){Require(!Exists(path),"cannot read job index journal");indexJournalRestored=true;return 0;}
    size_t restored=0;std::streamoff validEnd=0;bool torn=false;std::string header;
    while(std::getline(in,header)){
        if(in.eof()){torn=true;break;}
        size_t length=0;std::string hash,extra;std::istringstream h(header);
        Require(bool(h>>length>>hash)&&!(h>>extra)&&length>0&&length<=64*1024*1024,"invalid job index journal record");
        std::string body(length,'\0');in.read(&body[0],length);char newline=0;
        if(size_t(in.gcount())!=length||!in.get(newline)){torn=true;break;}
        Require(newline=='\n'&&Hash(body)==hash,"job index journal checksum mismatch");
        validEnd=in.tellg();std::istringstream row(body);std::string name,marker,diagnostics;
        Require(bool(std::getline(row,name)&&std::getline(row,marker)&&std::getline(row,diagnostics))
            &&!name.empty()&&name.find('/')==std::string::npos&&name!="."&&name!="..","invalid journal job path");
        std::ostringstream index;index<<row.rdbuf();Job j;
        if(ParseSampling(index.str(),root+"/jobs/"+name,j,false)){
            ParseDiagnostics(diagnostics,j);
            j.indexMarkerHash=marker;CacheJob(j);++restored;
            if(restored%1024==0)CheckHostMemory();
            if(restored%16384==0)Record(root+"/cache_restore_status.json",
                "{\"status\":\"reading_journal\",\"controller_pid\":"+Number(getpid())+",\"indexed_jobs\":"+Number(cache.size())+"}");
        }
    }
    Require(!in.bad(),"job index journal read failed");
    in.close();
    // Under the exclusive controller lock, discard only an incomplete final
    // append. Its completed directory remains discoverable in the fallback scan.
    if(torn){int fd=open(path.c_str(),O_WRONLY);Require(fd>=0,"cannot repair torn index tail");
        const bool okay=ftruncate(fd,validEnd)==0&&fsync(fd)==0;close(fd);Require(okay,"cannot truncate torn index tail");
        Log("cache_index_tail_repaired","{\"bytes_preserved\":"+Number(validEnd)+"}");}
    indexJournalRestored=true;return restored;
}
