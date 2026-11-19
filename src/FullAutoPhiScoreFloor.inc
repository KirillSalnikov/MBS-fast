// Bound the validation score for every admissible pair of control weights.
// Validation observations are already available. Training may choose a pair in
// [-4,4]^2, retain the default, or leave the raw element unchanged.
bool PhiValidationErrorFloor(const std::vector<Job>& valid,double& minimum) const {
    if(valid.size()!=8||valid.front().controls.empty())return false;
    using Vector=std::array<double,8>;
    auto dot=[](const Vector& a,const Vector& b){
        return std::inner_product(a.begin(),a.end(),b.begin(),0.);
    };
    auto center=[](Vector& a){
        double mean=std::accumulate(a.begin(),a.end(),0.)/8;
        for(double& v:a)v-=mean;
        return mean;
    };
    minimum=0;
    for(size_t r=0;r<valid.front().controls.size();++r){
        double mean=0;
        for(const Job& j:valid)mean+=j.controls[r][0]/8;
        if(!(mean>1e-25&&mean<1e25))return false;
        std::array<Vector,2> controls{};
        for(size_t s=0;s<8;++s){
            const auto& c=valid[s].controls[r];
            for(double value:c)if(!std::isfinite(value)||std::abs(value/mean)>1e6)return false;
            for(size_t k=0;k<2;++k)controls[k][s]=(c[1+k]-c[3+k])/mean;
        }
        const double c0=center(controls[0]),c1=center(controls[1]);
        const double meanShift=4*(std::abs(c0)+std::abs(c1))+1e-7;
        if(meanShift>=.8)return false;
        // Choose a test vector approximately orthogonal to both controls.
        // The bound below subtracts the actual remaining inner products and
        // therefore also works for collinear or nearly singular controls.
        auto basis=controls;
        for(size_t k=0;k<2;++k){
            if(k){const double projection=dot(basis[0],basis[1]);for(size_t s=0;s<8;++s)basis[1][s]-=projection*basis[0][s];}
            const double norm=std::sqrt(dot(basis[k],basis[k]));
            if(norm>1e-12){for(double& q:basis[k])q/=norm;}
            else basis[k].fill(0);
        }
        for(size_t j=0;j<16;++j){
            Vector raw;
            for(size_t s=0;s<8;++s){
                raw[s]=(j?valid[s].matrix[r][j]:valid[s].controls[r][0])/mean;
                if(!std::isfinite(raw[s])||std::abs(raw[s])>1e6)return false;
            }
            center(raw);
            // Calibrate accepts a fitted M11 only if its validation variance
            // decreases. A strictly positive raw worst-case Fieller denominator
            // also guarantees that this candidate cannot request a larger pilot.
            if(j==0&&t95*std::sqrt(dot(raw,raw)/56.)>=.9*(1-meanShift))return false;
            Vector u=raw;const bool fitted=j==0||polarControls;
            if(fitted)for(const auto& b:basis){
                const double projection=dot(u,b);for(size_t s=0;s<8;++s)u[s]-=projection*b[s];
            }
            const double norm=std::sqrt(dot(u,u));if(norm<1e-12)continue;
            for(double& v:u)v/=norm;
            const double weightPenalty=fitted?4*(std::abs(dot(u,controls[0]))+std::abs(dot(u,controls[1]))):0.;
            const double bound=std::max(0.,std::abs(dot(u,raw))-weightPenalty-1e-7);
            minimum=std::max(minimum,t95*bound/(std::sqrt(56.)*(1+meanShift)));
        }
    }
    return true;
}
