// Plan the theta pilot's tail from its observed dense orientation prefix.
// The independent production scrambles never enter the plan or its guard.
bool PrefixAnglePlan(const Samples& head,const std::vector<double>& scale,
                     std::vector<size_t>& anchors,std::vector<size_t>& probes,
                     std::vector<ResidualStencil>& bestWeights)const{
    const size_t n=theta.size();
    if(!highOrderThetaResidual||head.size()!=8||n<256||!prefixAnglePlan)return false;
    const auto fineAnchors=anchors;
    std::vector<double> fineTheta;for(auto i:fineAnchors)fineTheta.push_back(theta[i]);
    std::vector<ResidualStencil> fine(n);
    for(size_t i=0;i<n;++i){
        fine[i]=ResidualWeights(fineTheta,theta[i]);
        for(auto& term:fine[i])term.first=fineAnchors.at(term.first);
    }
    size_t best=anchors.size()+probes.size(),original=best;
    for(size_t step:{size_t(8),size_t(12)}){
        constexpr size_t width=16;
        std::vector<ResidualStencil> wide(n);
        std::vector<double> bounds(n,1e100);
        for(size_t i=0;i<n;++i){
            if(i/step<width/2-1)continue;
            const size_t left=(i/step-width/2+1)*step;
            if(left+(width-1)*step>=n||!(scale[i]>0))continue;
            if(i%step==0){wide[i]={{i,1.}};bounds[i]=0;continue;}
            double gain=0;
            for(size_t j=0;j<width;++j){
                size_t idx=left+j*step;double weight=1;
                for(size_t k=0;k<width;++k)if(k!=j){
                    const size_t other=left+k*step;
                    weight*=(theta[i]-theta[other])/(theta[idx]-theta[other]);
                }
                wide[i].push_back({idx,weight});gain+=std::abs(weight);
            }
            if(!(gain<=3))continue;
            Mueller mean{},squares{};
            for(const auto& s:head){
                Mueller diff=s[i];
                for(auto term:wide[i])for(size_t c=0;c<16;++c)diff[c]-=term.second*s[term.first][c];
                for(size_t c=0;c<16;++c){mean[c]+=diff[c]/8;squares[c]+=diff[c]*diff[c]/8;}
            }
            double worst=0;
            for(size_t c=0;c<16;++c)worst=std::max(worst,
                (std::abs(mean[c])+t95*std::sqrt(std::max(0.,squares[c]-mean[c]*mean[c])/7))/scale[i]);
            bounds[i]=worst;
        }
        std::vector<bool> coarse(n,false),used(n,false),checked(n,false);
        for(size_t start=0;start<n;start+=4*step){
            const size_t end=std::min(n,start+4*step);
            const auto lo=bounds.begin()+(start>step?start-step:0),hi=bounds.begin()+std::min(n,end+step);
            if(*std::max_element(lo,hi)<=interpolation/4)
                std::fill(coarse.begin()+start,coarse.begin()+end,true);
        }
        if(std::find(coarse.begin(),coarse.end(),true)==coarse.end())continue;
        std::vector<ResidualStencil> maps(n);
        for(size_t i=0;i<n;++i){
            maps[i]=coarse[i]?wide[i]:fine[i];
            for(auto t:maps[i])if(t.second!=0)used[t.first]=true;
            const size_t k=coarse[i]?step:4,block=i/(4*k);
            checked[i]=i%(4*k)==k*(block%4)+1+block%(k-1);
        }
        std::vector<size_t> a,b;std::vector<size_t> index(n,n);
        for(size_t i=0;i<n;++i){
            if(used[i]){index[i]=a.size();a.push_back(i);}
            else if(checked[i])b.push_back(i);
        }
        if(b.empty()||a.size()+b.size()>=best||a.size()+b.size()>original*.95)continue;
        for(auto& row:maps){
            row.erase(std::remove_if(row.begin(),row.end(),[](const std::pair<size_t,double>& v){return v.second==0;}),row.end());
            for(auto& t:row){t.first=index.at(t.first);Require(t.first<a.size(),"missing angular interpolation source");}
        }
        best=a.size()+b.size();anchors.swap(a);probes.swap(b);bestWeights.swap(maps);
    }
    return !bestWeights.empty();
}
