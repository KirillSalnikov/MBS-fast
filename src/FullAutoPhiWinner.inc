// Calibration-only views. Extend each stored prefix with the next Sobol
// block from the same scramble, keeping training and validation disjoint.
std::vector<Job> ExtendPhiWinner(std::vector<Job> prefix,const unsigned* scrambles,
                                 int wanted,int depth,const std::vector<double>& x){
    Require(prefix.size()==8&&sharedMeans&&controls,"invalid phi refit prefix");
    const int first=prefix.front().count;
    Require(first>0&&wanted%first==0&&((wanted/first)&(wanted/first-1))==0,"unaligned phi refit prefix");
    for(int offset=first;offset<wanted;offset*=2){
        auto extra=DirectStage(offset,scrambles,8,depth,x,false,false,offset);
        for(size_t s=0;s<8;++s){
            Job& a=prefix[s];const Job& b=extra[s];
            Require(a.count==offset&&a.offset==0&&b.offset==offset&&b.count==offset&&
                a.seed==scrambles[s]&&b.seed==a.seed&&a.phi==phi&&b.phi==phi&&
                a.depth==depth&&b.depth==depth&&a.angles==b.angles&&a.controls.size()==x.size()
                &&b.controls.size()==x.size(),"phi refit changed geometry, rule or sequence");
            for(size_t r=0;r<x.size();++r){
                for(int c:{1,2})Require(a.controls[r][c]==b.controls[r][c],"phi refit mean reference changed");
                for(int c:{0,3,4})a.controls[r][c]=.5*(a.controls[r][c]+b.controls[r][c]);
                for(size_t c=0;c<16;++c)a.matrix[r][c]=.5*(a.matrix[r][c]+b.matrix[r][c]);
            }
            a.count+=offset;a.seconds+=b.seconds;a.controlSetup+=b.controlSetup;
            a.controlRefinement=std::max(a.controlRefinement,b.controlRefinement);
        }
    }
    return prefix;
}
