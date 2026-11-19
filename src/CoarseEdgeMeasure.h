// Cubic moments of the original Gauss nodes. Within one radial Hermite cell their
// weighted integral is identical in real arithmetic. Across a knot, bound the
// change of the second polynomial on the whole affected interval.
struct BinnedEdgeMeasure {
    struct Bin {
        double origin=0,lo=std::numeric_limits<double>::infinity(),hi=0,absolute=0;
        std::array<double,4> sums{};
    };
    double width=0;
    std::vector<Bin> bins;
    std::vector<std::pair<double,double>> diagonal;
    template<class Patches>
    void Build(const Patches& patches,const std::vector<size_t>& ids,const std::vector<size_t>& mult,
               const std::vector<double>& x,const std::vector<double>& w,double step,double maxQ,double diameter) {
        if(!(maxQ>0 && diameter>0))return;
        width=step/maxQ;
        const double intervals=std::ceil(diameter/width);
        if(!(width>0 && std::isfinite(width) && intervals>=0 && intervals<=200000))return;
        bins.resize(size_t(intervals)+2);
        for(size_t i=0;i<bins.size();++i)bins[i].origin=(double(i)+0.5)*width;
        for(size_t g=0;g<ids.size();++g){
            const auto& p=patches[ids[g]];
            const size_t m=p.x.size();
            for(size_t i=0;i<m;++i){
                size_t k=(i+1)%m;double ex=p.x[k]-p.x[i],ey=p.y[k]-p.y[i];
                diagonal.emplace_back(std::hypot(ex,ey),double(mult[g])*(ex*ex+ey*ey));
                for(size_t j=i+1;j<m;++j){
                    size_t l=(j+1)%m;double fx=p.x[l]-p.x[j],fy=p.y[l]-p.y[j];
                    double dot=2.*mult[g]*(ex*fx+ey*fy);if(dot==0)continue;
                    for(size_t a=0;a<x.size();++a)for(size_t b=0;b<x.size();++b){
                        double dx=p.x[i]+x[a]*ex-p.x[j]-x[b]*fx,dy=p.y[i]+x[a]*ey-p.y[j]-x[b]*fy;
                        double d=std::hypot(dx,dy),weight=dot*(w[a]*w[b]);
                        double index=d/width;
                        if(!(index>=0 && index<bins.size())){bins.clear();return;}
                        size_t idx=size_t(index);
                        auto& bin=bins[idx];bin.lo=std::min(bin.lo,d);bin.hi=std::max(bin.hi,d);
                        bin.absolute+=std::abs(weight);
                        const double y=(d-bin.origin)/width;
                        for(int power=0;power<4;++power){bin.sums[power]+=weight;weight*=y;}
                    }
                }
            }
        }
        bins.erase(std::remove_if(bins.begin(),bins.end(),[](const Bin& b){return b.absolute==0;}),bins.end());
    }
    bool Evaluate(const Green& green,double Q,double& reflection)const{
        if(bins.empty())return false;
        double sum=0,error=0;
        const double rate=Q*width/green.step;
        auto cell=[&](double d){return std::min(size_t(Q*d/green.step),green.g.size()-2);};
        for(const auto& b:bins){
            size_t k=cell(b.lo),last=cell(b.hi);
            double a,c1,c2,c3;green.Polynomial(k,a,c1,c2,c3);
            if(last!=k){
                if(last!=k+1)return false;
                double v,s,t,u;green.Polynomial(last,v,s,t,u);
                const double h=std::abs((Q*b.hi-last*green.step)/green.step)+1e-11;
                // Express both cubics in the right interval's local variable.
                error+=b.absolute*(std::abs((a+c1+c2+c3)-v)+h*(std::abs((c1+2*c2+3*c3)-s)
                    +h*(std::abs((c2+3*c3)-t)+h*std::abs(c3-u))));
            }
            double t=(Q*b.origin-k*green.step)/green.step;
            sum+=b.sums[0]*(a+t*(c1+t*(c2+t*c3)))
                +rate*(b.sums[1]*(c1+t*(2*c2+3*t*c3))+rate*(b.sums[2]*(c2+3*t*c3)+rate*b.sums[3]*c3));
        }
        for(const auto& d:diagonal)sum+=d.second*green.Same(Q*d.first);
        // Margin on the approximation alone. Roundoff is tested independently
        // against the same quadrature in its ungrouped accumulation order.
        if(!std::isfinite(sum) || !(error<=1.e-8*std::abs(sum)))return false;
        reflection=-sum/(Q*Q);return std::isfinite(reflection) && reflection>=0;
    }
};
