#include "../src/FullAutoAuditPlan.h"
#include <fstream>
#include <iostream>
#include <cstdint>
using namespace FullAutoAuditPlan;
void require(bool x,const char*m){if(!x)throw std::runtime_error(m);}
int main(int argc,char**argv){try{
 require(argc==2,"expected measured allocation fixture");std::ifstream f(argv[1]);size_t size,levels,p,k;int initial;f>>size>>levels>>p>>k>>initial;
 std::vector<size_t> production(p),knots(k);for(auto&x:production)f>>x;for(auto&x:knots)f>>x;
 Budgets budgets(levels,std::vector<int>(size,256));budgets[0].assign(size,initial);
 for(size_t row=0;row<p;++row){size_t id;f>>id;for(auto&l:budgets)f>>l[id];}
 require(bool(f),"incomplete real allocation fixture");std::set<size_t> all(knots.begin(),knots.end());
 for(size_t i=1;i<knots.size();++i)for(size_t j=1;j<=3;++j)all.insert(knots[i-1]+(knots[i]-knots[i-1])*j/4);
 std::vector<size_t> audit(all.begin(),all.end());if(size<=4097){audit.clear();for(size_t i=0;i<size;++i)audit.push_back(i);}
 Extend(budgets,production,audit);
 uint64_t localWork=0,uniformWork=0,stencilWork=0;
 for(double fraction:{1.,.5,.0625}){
  const auto local=Scale(budgets,audit,fraction,initial,256,true),uniform=Scale(budgets,audit,fraction,initial,256,false);
  for(size_t l=0;l<levels;++l)for(size_t id:audit){require(local[l][id]<=uniform[l][id]&&local[l][id]<=budgets[l][id],"local allocation exceeds frozen budget");require(local[l][id]>0&&!(local[l][id]&(local[l][id]-1)),"not dyadic");}
  auto groups=Stencils(knots,local,audit);std::set<size_t> covered;size_t intervals=0;
  for(const auto&group:groups){intervals+=group.second.intervals.size();covered.insert(group.second.rows.begin(),group.second.rows.end());}
  require(intervals==knots.size()-1,"interpolation interval skipped");for(size_t id:audit)require(covered.count(id),"audit angle skipped");
  if(fraction==1){for(size_t l=0;l<levels;++l){for(size_t id:audit){localWork+=local[l][id];uniformWork+=uniform[l][id];}for(const auto&g:groups)stencilWork+=uint64_t(g.first[l])*g.second.rows.size();}}
 }
 require(localWork<uniformWork,"real heterogeneous fixture was globally inflated");
 bool rejected=false;try{Scale(budgets,audit,0,initial,256,true);}catch(const std::invalid_argument&){rejected=true;}require(rejected,"invalid fraction accepted");
 std::cout<<"{\"checked_theta\":"<<audit.size()<<",\"local_row_work\":"<<localWork<<",\"uniform_row_work\":"<<uniformWork<<",\"stencil_row_work\":"<<stencilWork<<",\"no_angles_skipped\":true}\n";return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
