// Output provenance only: never replaces an estimator or launches optical work.
std::string AnnotatedMatrixHeader(){
    std::ostringstream out;
    out<<"theta_deg,point_type,value_source,is_support,has_direct_calculation,interpolation_estimate_kind,max_statistical95_scaled_M11,max_ratio95_absolute,max_interpolation95_or_allowance,max_combined95";
    for(int j=0;j<16;++j)out<<",M"<<j/4+1<<j%4+1;
    out<<'\n';return out.str();
}
void AnnotatedMatrixRow(std::ostringstream& out,double theta,const char* kind,const char* source,
                        bool support,bool direct,const char* interpolationKind,
                        const Statistics& stats,size_t row,bool checked,bool ratios,
                        double interpolation=-1,double combined=-1){
    out<<Number(theta)<<','<<kind<<','<<source<<','<<support<<','<<direct<<','<<interpolationKind<<',';
    if(checked)out<<Number(*std::max_element(stats.interval.at(row).begin(),stats.interval.at(row).end()));
    out<<',';
    if(checked&&ratios)out<<Number(*std::max_element(stats.ratioInterval.at(row).begin()+1,stats.ratioInterval.at(row).end()));
    out<<',';if(interpolation>=0)out<<Number(interpolation);
    out<<',';if(combined>=0)out<<Number(combined);
    for(double v:stats.mean.at(row))out<<','<<Number(v);
    out<<'\n';
}
void WriteOutputGuide(Engine& e,const std::string& folder,size_t fullRows,size_t directRows,size_t supportRows,bool commonPrefix){
    Atomic(folder+"/matrix_output_schema.json","{\"schema\":1,\"method\":"+Quoted(e.go?"go":"po")
        +",\"full_matrix\":\"mueller_fullauto_annotated.csv\",\"direct_matrix\":\"mueller_evaluated.dat\",\"direct_matrix_with_errors\":\"mueller_evaluated_annotated.csv\""
        +",\"full_rows\":"+Number(fullRows)+",\"direct_rows\":"+Number(directRows)+",\"support_rows\":"+Number(supportRows)
        +",\"ratio_accuracy_enabled\":"+(e.ratios?"true":"false")
        +",\"interpolation_checks_use_common_prefix\":"+(commonPrefix?"true":"false")
        +",\"value_source\":{\"direct\":\"direct orientation average\",\"interpolated\":\"local cubic interpolation of support estimates, including guard angles\",\"direct_bin\":\"GO finite solid-angle bin density; poles are separate\"}"
        +",\"point_type\":{\"support\":\"PO interpolation anchor\",\"guard\":\"PO angle with a separate direct interpolation check\",\"interpolated\":\"no direct calculation in this production estimator; may have been evaluated in a pilot or earlier round\",\"finite_bin\":\"GO direct bin, no interpolation\"}"
        +",\"errors\":\"Dimensionless, multiply by 100 for percent of M11 or percentage points of Mij/M11. Statistical columns are pointwise 95% half-width estimates, not bounds on physical error. Empty fields mean not provided or not checked, never zero error.\""
        +",\"interpolation\":\"For the full PO matrix: componentwise guard estimates or the configured allowance at untested angles; the displayed column is their maximum. Common-prefix checks isolate interpolation on matched samples; they do not bound the difference of two DAT files with different local sample counts.\""
        +",\"combined\":\"Full PO matrix only: maximum of componentwise statistical plus interpolation estimates plus fixed reserves, as in pointwise_budget_production.csv; not the sum of the displayed maxima. Direct-matrix exports do not assert this full-grid acceptance metric.\""
        +",\"go_poles\":\"GO forward and backward integrated delta contributions remain in mueller_forward.dat and mueller_backward.dat, not in annotated finite-bin matrices.\""
        +",\"global_simultaneous_confidence\":false}\n");
}
void WritePOOutputAnnotations(Engine& e,const std::string& folder,const std::vector<size_t>& knots,
                              const std::vector<size_t>& ids,const Statistics& direct,
                              const Statistics& output,const InterpolationResult& interp,bool commonPrefix){
    Require(output.mean.size()==e.theta.size()&&direct.mean.size()==ids.size(),"annotated matrix dimensions differ");
    std::ostringstream full,evaluated;full<<AnnotatedMatrixHeader();evaluated<<AnnotatedMatrixHeader();
    size_t row=0;
    for(size_t q=0;q<e.theta.size();++q){
        const bool hasDirect=row<ids.size()&&ids[row]==q;
        const bool support=std::binary_search(knots.begin(),knots.end(),q);
        Require(!support||hasDirect,"annotated support lacks direct calculation");
        const char* kind=support?"support":(hasDirect?"guard":"interpolated");
        const auto residual=ResidualAt(e,interp,ids,q);
        double maximum=0;for(int j=0;j<16;++j){maximum=std::max(maximum,residual.matrix[j]);if(e.ratios&&j)maximum=std::max(maximum,residual.ratio[j]);}
        const char* check=support?"support_zero":(hasDirect?(commonPrefix?"common_prefix_guard95":"paired_guard95"):"unmeasured_allowance");
        AnnotatedMatrixRow(full,e.theta[q],kind,support?"direct":"interpolated",support,hasDirect,check,
            output,q,true,e.ratios,maximum,MetricTotal(e,BudgetMetric(e,output,q,residual),residual));
        if(hasDirect){
            AnnotatedMatrixRow(evaluated,e.theta[q],kind,"direct",support,true,"not_applicable",direct,row,true,e.ratios);
            ++row;
        }
    }
    Require(row==ids.size(),"annotated direct angles not fully exported");
    Atomic(folder+"/mueller_fullauto_annotated.csv",full.str());
    Atomic(folder+"/mueller_evaluated_annotated.csv",evaluated.str());
    WriteOutputGuide(e,folder,e.theta.size(),ids.size(),knots.size(),commonPrefix);
}
void WriteGOOutputAnnotations(Engine& e,const std::string& folder,const Statistics& stats,const std::vector<bool>& checked){
    std::ostringstream out;out<<AnnotatedMatrixHeader();
    for(size_t q=0;q<e.theta.size();++q)
        AnnotatedMatrixRow(out,e.theta[q],"finite_bin","direct_bin",false,true,"not_applicable",stats,q,checked.at(q),e.ratios);
    Atomic(folder+"/mueller_fullauto_annotated.csv",out.str());
    Atomic(folder+"/mueller_evaluated_annotated.csv",out.str());
    WriteOutputGuide(e,folder,e.theta.size(),e.theta.size(),0,false);
}
void HashMatrixAnnotations(const std::string& folder,std::ostringstream& hashes){
    for(const char* name:{"mueller_fullauto_annotated.csv","mueller_evaluated_annotated.csv","matrix_output_schema.json"})
        hashes<<name<<' '<<FileHash(folder+"/"+name)<<'\n';
}
