// Included inside Engine. A presentation layer only: numerical/recovery state
// remains under root (work/). Never infer a delta mass from a PO endpoint.
void ConfigureResultLayout(){
    const bool compact=args.IsCatched("fullauto_layout")&&args.GetStringValue("fullauto_layout")=="compact";
    const bool marked=Exists(root+"/.mbs-result-layout-v1");
    Require(compact||!marked,"compact result requires --fullauto-layout compact when resuming");
    if(!compact)return;
    Require(!Exists(root+"/configuration.txt")&&!Exists(root+"/memory_storage_v1.txt"),
        "compact layout requires a new output directory; do not move existing checkpoints");
    if(marked)Require(Read(root+"/.mbs-result-layout-v1")=="mbs_compact_result_v1\n","unknown result layout marker");
    else{
        for(const char* name:{"work","results","archive","latest","config.json","status.json"}){
            struct stat st{};Require(lstat((root+"/"+name).c_str(),&st)!=0&&errno==ENOENT,"new compact output contains reserved path "+std::string(name));}
        Atomic(root+"/.mbs-result-layout-v1","mbs_compact_result_v1\n");
    }
    publicRoot=root;root+="/work";Directory(root);Directory(publicRoot+"/results");Directory(publicRoot+"/archive");
    struct stat pointer{};const auto latest=publicRoot+"/latest";
    if(lstat(latest.c_str(),&pointer)==0){
        Require(S_ISLNK(pointer.st_mode),"compact latest must be a symbolic link");
        char buffer[256];const auto count=readlink(latest.c_str(),buffer,sizeof(buffer));
        Require(count>0&&size_t(count)<sizeof(buffer),"invalid compact latest pointer");
        const std::string target(buffer,size_t(count));
        Require(target.compare(0,14,"results/round_")==0&&target.substr(8).find_first_not_of("abcdefghijklmnopqrstuvwxyz0123456789_")==std::string::npos,"unsafe compact latest pointer");
        lastCompactGeneration=target.substr(8);
    }

}
std::string CompactOptionalJson(const std::string& name)const{
    return Exists(root+"/"+name)?CompactJson(Read(root+"/"+name)):"null";
}
void ExportCompactResults(bool force=false){
    if(publicRoot.empty())return;
    const double now=Seconds();if(!force&&compactPublishedAt>=0&&now-compactPublishedAt<10)return;compactPublishedAt=now;
    const auto state=CompactOptionalJson("fullauto_status.json");
    const auto parameters=Exists(root+"/final_parameters.json")?CompactOptionalJson("final_parameters.json"):CompactOptionalJson("fullauto_parameters.json");
    const auto observation=CompactOptionalJson("execution_progress.json");
    const auto current=CompactOptionalJson("intermediate_current.json");
    const bool completed=Exists(root+"/final_outputs.fnv")&&state.find("\"status\":\"converged\"")!=std::string::npos;
    std::string source=completed?root:MemorySnapshotFolder(current);
    const bool pilot=current.find("\"phase\":\"depth_selection\"")!=std::string::npos;
    std::map<std::string,std::string> members;
    // Pilot estimates are not a converged angular result. They remain in work/.
    const std::string matrix=source+(go?"/mueller_evaluated.dat":"/mueller_fullauto.dat");
    if(force&&!source.empty()&&(completed||!pilot)&&Exists(matrix)){
        members["mueller_bins.dat"]=Read(matrix);
        if(Exists(source+"/mueller_fullauto_integrals.tsv"))members["integrals.tsv"]=Read(source+"/mueller_fullauto_integrals.tsv");
        for(const auto& entry:std::vector<std::pair<std::string,std::string>>{
            {"mueller_fullauto_annotated.csv","angles.csv"},{"parameters_accuracy.csv","accuracy_by_angle.csv"},
            {"orientation_budget.csv","orientation_budget.csv"},{"pole_status.json","pole_status.json"},
            {"matrix_output_schema.json","matrix_schema.json"}})
            if(Exists(source+"/"+entry.first))members[entry.second]=Read(source+"/"+entry.first);
        std::ostringstream poles;poles<<"# integrated delta cross sections; not densities per steradian\n# theta_deg M11 M12 M13 M14 M21 M22 M23 M24 M31 M32 M33 M34 M41 M42 M43 M44\n";
        if(go){
            for(int side=0;side<2;++side){
                const auto file=source+(side?"/mueller_backward.dat":"/mueller_forward.dat");
                Require(Exists(file),"GO result lacks separate pole contribution");
                std::istringstream in(Read(file));std::string line;std::getline(in,line);double angle,weight,value;
                Require(bool(in>>angle>>weight),"invalid GO pole header");poles<<(side?180:0);
                for(int j=0;j<16;++j){Require(bool(in>>value)&&std::isfinite(value),"invalid GO pole matrix");poles<<' '<<std::setprecision(17)<<value;}poles<<'\n';
            }
        }else poles<<"# not_applicable: PO reports angular densities; no separate delta estimator\n";
        members["pole_contributions.dat"]=poles.str();
        members["accuracy.json"]="{\"schema\":1,\"method\":"+Quoted(go?"go":"po")+",\"converged\":"+(completed?"true":"false")
            +",\"matrix_kind\":"+Quoted(go?"finite_solid_angle_bin_density":"pointwise_angular_density_with_interpolation")
            +",\"poles\":"+Quoted(go?"separate_integrated_delta_cross_sections":"not_applicable")
            +",\"integrals_available\":"+(members.count("integrals.tsv")?"true":"false")
            +",\"estimate\":"+current+",\"precision_contract\":"+CompactOptionalJson("precision_contract.json")
            +",\"final_checks\":"+CompactOptionalJson("final_checks_policy.json")
            +",\"depth_selection\":"+CompactOptionalJson("depth_selection.json")
            +",\"physics_validation\":"+(completed?CompactOptionalJson("physics_validation.json"):"null")
            +",\"independent_validation\":"+(completed?CompactOptionalJson("independent_final_validation.json"):"null")+"}\n";
        // Freeze selected settings with the matrix, not a later in-flight round.
        members["config.json"]="{\"parameters\":"+parameters+",\"selected_state\":"+state+"}\n";
        std::ostringstream manifest;manifest<<"mbs_compact_generation_v1\n";
        for(const auto& item:members)manifest<<item.first<<' '<<Hash(item.second)<<'\n';
        std::string round="unknown";
        const auto roundKey=current.find("\"round\":");
        if(roundKey!=std::string::npos){size_t at=roundKey+8;while(at<current.size()&&current[at]==' ')++at;const auto begin=at;
            while(at<current.size()&&current[at]>='0'&&current[at]<='9')++at;if(at>begin)round=current.substr(begin,at-begin);}
        const auto id="round_"+round+"_"+Hash(manifest.str()),relative="results/"+id,destination=publicRoot+"/"+relative;
        if(id!=lastCompactGeneration){
            PhysicalStorageScope disk;
            if(Exists(destination)){
                Require(Read(destination+"/manifest.fnv")==manifest.str(),"published compact manifest changed");
                for(const auto& item:members)Require(Read(destination+"/"+item.first)==item.second,"published compact result changed");
            }else{
                const auto pattern=publicRoot+"/results/.pending-XXXXXX";std::vector<char> tmp(pattern.begin(),pattern.end());tmp.push_back(0);
                char* made=mkdtemp(tmp.data());Require(made,"cannot stage compact result");const std::string stage=made;
                try{for(const auto& item:members)Atomic(stage+"/"+item.first,item.second);Atomic(stage+"/manifest.fnv",manifest.str());
                    AtomicFileBatch::FlushPath(stage,true);Require(rename(stage.c_str(),destination.c_str())==0,"cannot publish compact generation");
                    AtomicFileBatch::FlushPath(publicRoot+"/results",true);
                }catch(...){WorkerMemoryStage::RemoveTree(stage);throw;}
            }
            const auto temp=publicRoot+"/.latest-"+Number(getpid());struct stat st{};
            if(lstat(temp.c_str(),&st)==0){Require(S_ISLNK(st.st_mode),"unexpected compact pointer staging file");Require(unlink(temp.c_str())==0,"cannot remove stale pointer");}
            struct stat latestStat{};
            if(lstat((publicRoot+"/latest").c_str(),&latestStat)==0)Require(S_ISLNK(latestStat.st_mode),"refusing to replace non-symlink latest");
            Require(symlink(relative.c_str(),temp.c_str())==0,"cannot stage latest symlink");
            Require(rename(temp.c_str(),(publicRoot+"/latest").c_str())==0,"cannot publish latest symlink");AtomicFileBatch::FlushPath(publicRoot,true);
            lastCompactGeneration=id;
        }
    }
    const auto utc=std::chrono::duration<double>(std::chrono::system_clock::now().time_since_epoch()).count();
    PhysicalStorageScope disk;
    Atomic(publicRoot+"/config.json","{\"schema\":1,\"method\":"+Quoted(go?"go":"po")+",\"parameters\":"+parameters+",\"selected_state\":"+state+"}\n");
    Atomic(publicRoot+"/status.json","{\"schema\":1,\"updated_unix\":"+Number(utc)+",\"controller\":"+state+",\"execution\":"+observation
        +",\"recovery\":"+MemoryStorageReport()+",\"latest\":"+(lastCompactGeneration.empty()?"null":Quoted("results/"+lastCompactGeneration))+"}\n");
}
