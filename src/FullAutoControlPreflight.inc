// A pilot ranks theta rows only. The actual eight-scramble production
// estimator decides whether the automatic analytic control can be dropped.
bool controlPreflight=false;
int controlPreflightPhi=0;
std::vector<double> controlPreflightAngles;
std::vector<double> RankControlPreflight(const std::vector<double>& x,const std::vector<Job>& valid,
        const std::vector<std::array<double,2>>& scalar,const std::vector<std::array<Mueller,2>>& polar)const{
    Require(valid.size()==8 && scalar.size()==x.size() && polar.size()==x.size(),"control preflight requires eight pilot samples and weights for each requested angle");
    for(const auto& job:valid)Require(job.controls.size()==x.size(),"control preflight pilot row count mismatch");
    std::vector<std::pair<double,double>> risks;
    const bool corrected=balancedCorrection&&ProductionDepths().front()>=6&&initial>=512;
    const double terms=corrected&&componentwiseCutoffMeanBound
        ?(StagedCutoff()&&(referenceCutoffMeanBound||initial>512)?5.:3.):1.;
    for(size_t row=0;row<x.size();++row){
        double mean=0,variance=0,numerator=0;
        for(const auto& job:valid)mean+=job.controls[row][0]/valid.size();
        for(const auto& job:valid){
            const auto& c=job.controls[row];const double deviation=c[0]-mean;
            variance+=deviation*deviation/(valid.size()*(valid.size()-1));
            for(size_t k=0;k<16;++k){
                const double wr=k?polar[row][0][k]:scalar[row][0],ws=k?polar[row][1][k]:scalar[row][1];
                numerator=std::max(numerator,ControlMeanRefinement(job,row,wr,ws));
            }
        }
        const double denominator=mean+t95*std::sqrt(variance);
        if(denominator>0&&terms*numerator>eps/64*denominator)risks.emplace_back(terms*numerator/denominator,x[row]);
    }
    std::stable_sort(risks.begin(),risks.end(),[](const std::pair<double,double>& a,const std::pair<double,double>& b){return a.first>b.first;});
    std::vector<double> selected;
    for(size_t i=0;i<std::min<size_t>(4,risks.size());++i)selected.push_back(risks[i].second);
    std::sort(selected.begin(),selected.end());return selected;
}
