// Pilot only. The fixed stencil applies to each independent orientation scramble.
bool CardinalThetaPilot(int count,Samples& sample,Matrix& reference,const std::vector<size_t>& initial){
    auto readRows=[&](const std::vector<size_t>& ids){
        std::vector<double> angles;angles.reserve(ids.size());
        for(size_t i:ids)angles.push_back(theta[i]);
        const auto jobs=DirectStage(count,angularValidation,8,target,angles);
        auto values=Data(jobs);std::vector<std::vector<double>> positive(8,std::vector<double>(ids.size()));
        for(size_t s=0;s<8;++s)for(size_t i=0;i<ids.size();++i)
            positive[s][i]=controls?jobs[s].controls.at(i)[0]:values[s][i][0];
        return std::make_pair(std::move(values),std::move(positive));
    };
    return RunCardinalThetaPilot(count,sample,reference,initial,readRows);
}
template<class Reader>
bool RunCardinalThetaPilot(int count,Samples& sample,Matrix& reference,const std::vector<size_t>& initial,Reader readRows){
    if(!cardinalTheta||go||!coherent||theta.size()<257||count<32||phiCutoffPilot||
       (controls&&!sharedMeans)||args.IsCatched("fullauto_theta_audit")||
       args.IsCatched("fullauto_pilot")||config.adaptive.loadedFromFile)return false;
    constexpr int stride=4,radius=12,band=24;
    const size_t rows=theta.size();
    if(cardinalAuto && double(rows)*count*phi<1024.*1024)return false;
    const double h=(theta.back()-theta.front())/(rows-1);
    if(!(h>0))return false;
    // A coarse step spanning a large scattering phase can hide an oscillation.
    // Use Dmax as a conservative phase scale when the geometry is available.
    if(!(maximumDiameter>0&&config.wavelengthUm>0)||
       (2*std::acos(-1.)*maximumDiameter/config.wavelengthUm)*(stride*h*std::acos(-1.)/180.)>2.5)return false;
    for(size_t i=1;i<rows;++i)
        if(std::fabs((theta[i]-theta[i-1])-h)>1e-9*h)return false;
    static const double coefficients[3][24]={
        {-2.6753989018446346e-05, 0.00015699691006465445, -0.00056672785020835082, 0.0015826415825194132, -0.0037347601279028478, 0.0078099858701251752, -0.014941761956588279, 0.026867588241296626, -0.046815770228987655, 0.083063621166207802, -0.17007898487366735, 0.89826934955310822, 0.29401522537286717, -0.11497872291377847, 0.061926066573610383, -0.035521797285152423, 0.020157598948294182, -0.010904698286555367, 0.0054709131831819843, -0.0024742248805638845, 0.00097114238006380587, -0.00031067543840562306, 7.0829857546183121e-05, -7.0818080567068228e-06},
        {-2.0812839851279684e-05, 0.00015169464410042573, -0.00059865099204440924, 0.001763234382245932, -0.0043163880975570527, 0.0092718716911907532, -0.018090711098651787, 0.032953649945707215, -0.05767256534215412, 0.10116430135824017, -0.19545296916331223, 0.63084734551208632, 0.63084734551208632, -0.19545296916331223, 0.10116430135824017, -0.05767256534215412, 0.032953649945707215, -0.018090711098651787, 0.0092718716911907532, -0.0043163880975570527, 0.001763234382245932, -0.00059865099204440924, 0.00015169464410042573, -2.0812839851279684e-05},
        {-7.0818080567068228e-06, 7.0829857546183121e-05, -0.00031067543840562306, 0.00097114238006380587, -0.0024742248805638845, 0.0054709131831819843, -0.010904698286555367, 0.020157598948294182, -0.035521797285152423, 0.061926066573610383, -0.11497872291377847, 0.29401522537286717, 0.89826934955310822, -0.17007898487366735, 0.083063621166207802, -0.046815770228987655, 0.026867588241296626, -0.014941761956588279, 0.0078099858701251752, -0.0037347601279028478, 0.0015826415825194132, -0.00056672785020835082, 0.00015699691006465445, -2.6753989018446346e-05}
    };
    static const double narrow[3][12]={
        {-0.0002145279719746456, 0.0024935824094670551, -0.013112366900728954, 0.046603357250601288, -0.14306749282622983, 0.89215016430442706, 0.27641140094088951, -0.081561076561389501, 0.025584051345963726, -0.0061042760627157532, 0.00084661778275358745, -2.9433711063471827e-05},
        {-0.0001267048912736189, 0.0021129192469363277, -0.012829503642883489, 0.049171630326011059, -0.15214036364354458, 0.61381202260475431, 0.61381202260475431, -0.15214036364354458, 0.049171630326011059, -0.012829503642883489, 0.0021129192469363277, -0.0001267048912736189},
        {-2.9433711063471827e-05, 0.00084661778275358745, -0.0061042760627157532, 0.025584051345963726, -0.081561076561389501, 0.27641140094088951, 0.89215016430442706, -0.14306749282622983, 0.046603357250601288, -0.013112366900728954, 0.0024935824094670551, -0.0002145279719746456}
    };
    std::vector<bool> anchor(rows),guard(rows);
    std::vector<size_t> measured;
    for(size_t r=0;r<rows;++r)
        anchor[r]=r<stride*radius||r+stride*radius>=rows||r%stride==0;
    for(size_t begin=0;begin<rows;begin+=band){
        std::vector<size_t> remaining;
        for(size_t r=begin;r<std::min(rows,begin+band);++r)if(!anchor[r])remaining.push_back(r);
        for(size_t residue=1;residue<stride;++residue){
            std::vector<size_t> group;
            for(size_t r:remaining)if(r%stride==residue)group.push_back(r);
            if(!group.empty())guard[group[(begin/band+2*(residue-1))%group.size()]]=true;
        }
    }
    for(size_t r=0;r<rows;++r)if(anchor[r]||guard[r])measured.push_back(r);
    sample.assign(8,Matrix(rows));Samples direct(8,Matrix(measured.size()));
    std::vector<std::vector<double>> intensity(8,std::vector<double>(rows)),raw(8,std::vector<double>(measured.size()));
    auto fetch=[&](const std::vector<size_t>& ids,bool first){
        const auto result=readRows(ids);const auto& values=result.first;const auto& positive=result.second;
        for(size_t s=0;s<8;++s)for(size_t i=0;i<ids.size();++i){
            const auto r=ids[i];const auto m11=positive[s][i];
            sample[s][r]=values[s][i];intensity[s][r]=m11;
            if(first){direct[s][i]=values[s][i];raw[s][i]=m11;}
        }
    };
    fetch(measured,true);
    std::vector<bool> reject((rows+band-1)/band);
    for(size_t s=0;s<8;++s)for(size_t r=0;r<rows;++r){
        if(anchor[r])continue;
        Mueller m{};double scale=0;const auto& w=coefficients[r%stride-1];
        const size_t begin=r-r%stride-stride*(radius-1);
        for(size_t j=0;j<2*radius;++j){
            const size_t i=begin+stride*j;scale+=w[j]*intensity[s][i];
            for(size_t k=0;k<16;++k)m[k]+=w[j]*sample[s][i][k];
        }
        sample[s][r]=m;intensity[s][r]=scale;
        if(!(scale>0&&std::isfinite(scale)))reject[r/band]=true;
        for(double v:m)if(!std::isfinite(v))reject[r/band]=true;
    }
    // Validate every matrix channel and the raw normalization separately.
    const double bound=std::min(2e-4,interpolation/50.);
    for(size_t s=0;s<8;++s)for(size_t i=0;i<measured.size();++i){
        const size_t r=measured[i];if(!guard[r])continue;
        const double m11=raw[s][i];double error=std::fabs(intensity[s][r]-m11);
        for(int k=0;k<16;++k)error=std::max(error,std::fabs(sample[s][r][k]-direct[s][i][k]));
        if(!(m11>0&&std::isfinite(m11)&&error<=bound*m11))reject[r/band]=true;
        sample[s][r]=direct[s][i];intensity[s][r]=m11;
    }
    std::vector<size_t> missing;
    for(size_t r=0;r<rows;++r)if(reject[r/band]&&!anchor[r]&&!guard[r])missing.push_back(r);
    if(!missing.empty())fetch(missing,false);
    std::vector<bool> exact=anchor;for(size_t i:measured)exact[i]=true;
    for(size_t i:missing)exact[i]=true;
    auto normalize=[&](const Samples& signals,const std::vector<std::vector<double>>& rawM){
        Matrix result=Mean(signals);
        if(controls)for(size_t r=0;r<rows;++r){result[r][0]=0;for(size_t s=0;s<8;++s)result[r][0]+=rawM[s][r]/8;}
        return result;
    };
    auto competing=sample;auto altRaw=intensity;
    for(size_t r=0;r<rows;++r)if(!exact[r]){
        const auto begin=r-r%stride-5*stride;
        for(size_t s=0;s<8;++s){
            Mueller m{};double rawM=0;
            for(int j=0;j<12;++j){
                const auto source=begin+stride*j;const auto w=narrow[r%stride-1][j];
                rawM+=w*intensity[s][source];
                for(int k=0;k<16;++k)m[k]+=w*sample[s][source][k];
            }
            competing[s][r]=m;altRaw[s][r]=rawM;
        }
    }
    auto support=[&](const Samples& values,const Matrix& norm){
        auto chosen=initial;const auto flag=std::getenv("MBS_FULLAUTO_MATCH_THETA_PILOT");bool match=flag&&std::string(flag)=="1";
        double maximum=0;int steps=0;
        SelectPilotKnots(chosen,values,norm,interpolation*(match?.75:.5),match,false,maximum,steps);
        return chosen;
    };
    size_t disputed=0,extraRows=0;
    for(int attempt=0;attempt<2;++attempt){
        reference=normalize(sample,intensity);const auto alt=normalize(competing,altRaw);
        bool valid=true;for(size_t r=0;r<rows;++r)if(!(reference[r][0]>0&&alt[r][0]>0))valid=false;
        std::vector<size_t> differences;
        std::set<size_t> needed;
        if(valid){
            const auto first=support(sample,reference),second=support(competing,alt);
            std::set_symmetric_difference(first.begin(),first.end(),second.begin(),second.end(),std::back_inserter(differences));
            if(differences.empty())break;
            disputed+=differences.size();
            if(attempt==0)for(size_t r:differences)for(auto set:{&first,&second}){
                const size_t at=std::lower_bound(set->begin(),set->end(),r)-set->begin();
                const size_t lo=(*set)[at>1?at-2:0],hi=(*set)[std::min(set->size()-1,at+2)];
                for(size_t j=lo;j<=hi;++j)if(!exact[j])needed.insert(j);
            }
        }
        if(needed.empty())for(size_t j=0;j<rows;++j)if(!exact[j])needed.insert(j);
        if(needed.empty())break;
        std::vector<size_t> verify(needed.begin(),needed.end());fetch(verify,false);extraRows+=verify.size();
        for(size_t j:verify){exact[j]=true;for(size_t s=0;s<8;++s){competing[s][j]=sample[s][j];altRaw[s][j]=intensity[s][j];}}
    }
    reference=normalize(sample,intensity);
    Log("cardinal_theta_pilot","{\"rows\":"+Number(rows)+",\"initial_direct\":"+Number(measured.size())+
        ",\"fallback_direct\":"+Number(missing.size())+",\"disputed_knots\":"+Number(disputed)+",\"extra_direct\":"+Number(extraRows)+",\"count_per_seed\":"+Number(count)+
        ",\"replicas\":8,\"bound\":"+Number(bound)+",\"production_accuracy_accepted\":false}");
    return true;
}
