// Pure scheduling/decision regressions. No optical calculations or files.
#include "../src/FullAutoDepthSelection.h"
#include <algorithm>
#include <climits>
#include <cmath>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>
namespace D=FullAutoDepthSelection;
namespace {
void Check(bool pass,const char* message){if(!pass)throw std::runtime_error(message);}
void ActualPilotRegressions(){
    const double budget=.01,quietGuard=.0001;
    // Recorded coherent L112.2 n6/N1024: neither pass nor statistically
    // resolved failure. These old adjacent-depth scalar diagnostics exercise
    // classification only; they are not evidence for the new anchor policy.
    Check(D::Classify(.015704045626191973,.0069374632438099692,
                     quietGuard,0,budget)==D::Unresolved,
          "L112.2 ambiguous pilot was promoted to measured failure");
    // Recorded L112.2 n10/N512 is a clear lower-bound failure.
    Check(D::Classify(.14461627685516504,.1027051440704482,
                     quietGuard,0,budget)==D::Fail,
          "resolved L112.2 residual did not fail");
    // Recorded L100 n6/N2048 was resolved below its one-percent budget.
    Check(D::Classify(.0074249833017337548,.004507890624674829,
                     quietGuard,0,budget)==D::Pass,
          "resolved residual below budget rejected");
    // L794.3 n18/N512..2048: small paired differences but unresolved M11
    // denominator. 1e100 is a sentinel, not a finite acceptance bound.
    for(double lower:{0.,.00001,.00008}){
        Check(D::Classify(1e100,lower,quietGuard,0,budget)==D::Unresolved,
              "unresolved L794.3 denominator accepted");
        Check(D::Classify(.000085533564328636009,lower,1e100,0,budget)==D::Unresolved,
              "unresolved reference denominator accepted");
    }
    Check(D::Classify(.000085533564328636009,0,quietGuard,0,budget)==D::Pass,
          "resolved tiny residual failed");
}
void AcceptanceGuards(){
    // Saved bullet L100 pilot: the old quarter-budget guard rejected this.
    Check(D::Classify(.0041059213384853866,0,.0038401379400098304,0,.01)==D::Pass,
          "L100 reference still uses obsolete quarter-budget threshold");
    // Relaxing the guard does not permit a candidate above the depth budget.
    Check(D::Classify(.024709904818480715,.02,.024692857733782478,.02,.01)==D::ReferenceFail,
          "L141 reference above full depth budget accepted");

    const double q=std::numeric_limits<double>::quiet_NaN();
    const double inf=std::numeric_limits<double>::infinity();
    for(double bad:{q,inf,-inf,-1.,1e50,1e100}){
        Check(!D::FiniteBound(bad),"invalid bound called finite");
        Check(D::Classify(bad,0,.0001,0,.01)==D::Unresolved,
              "invalid candidate upper accepted");
        Check(D::Classify(.001,0,bad,0,.01)==D::Unresolved,
              "invalid reference upper accepted");
        Check(D::Classify(.001,bad,.0001,0,.01)==D::Unresolved,
              "invalid candidate lower accepted");
        Check(D::Classify(.001,0,.0001,bad,.01)==D::Unresolved,
              "invalid reference lower accepted");
        Check(D::Classify(bad,.2,.005,.004,.01)==D::Unresolved,
              "unbounded upper with finite lower called a resolved failure");
    }
    Check(D::FiniteBound(0)&&D::FiniteBound(.01),"valid finite bound rejected");
    for(double bad:{0.,-1.,q,inf})
        Check(D::Classify(0,0,0,0,bad)==D::Unresolved,"invalid budget accepted");
    Check(D::Classify(.01,0,.01,0,.01)==D::Pass,"closed upper thresholds changed");
    Check(D::Classify(.009,0,.01000001,0,.01)==D::Unresolved,
          "unresolved reference guard accepted candidate");
    Check(D::Classify(.001,0,.014,.013,.01)==D::ReferenceFail,
          "measured reference drift not distinguished");
    Check(D::Classify(.012,.01,0,0,.01)==D::Unresolved,
          "lower bound equal to budget wrongly failed");
    Check(D::Classify(.012,.010001,0,0,.01)==D::Fail,
          "lower bound strictly above budget not failed");
    // Candidate upper is DIRECTLY n->anchor; do not add guard a second time.
    Check(D::Classify(.009,0,.002,0,.01)==D::Pass,
          "direct-to-anchor error was double-counted with guard");
}
void CandidateOrder(){
    Check(D::ProductionDepth(false,6,18)==18,"early reference failure did not use requested maximum");
    Check(D::ProductionDepth(false,12,14)==14,"explicit maximum ignored in fallback");
    Check(D::ProductionDepth(true,8,18)==8,"passing candidate promoted unnecessarily");

    Check(D::Candidates(6,18)==std::vector<int>({6,12,18,8,10,14,16}),
          "default coarse/fine order changed");
    // Nonmonotonic coherent residual: n12 fails although skipped n8 passes.
    // The candidate list must retain n8 even after a coarse n18 pass.
    const auto order=D::Candidates(6,18);int selected=0;std::vector<int> visited;
    for(int n:order){
        if(selected&&n>=selected)continue;
        visited.push_back(n);
        const double upper=(n==8||n==18)?.005:.025;
        if(D::Classify(upper,upper,0,0,.01)==D::Pass)selected=n;
    }
    Check(selected==8,"coarse bracket discarded nonmonotonic passing n8");
    Check(std::find(visited.begin(),visited.end(),8)!=visited.end(),
          "skipped lower depth never revisited");
    for(int maximum=1;maximum<=28;++maximum)for(int start=1;start<=maximum;++start){
        const auto v=D::Candidates(start,maximum);std::set<int> unique(v.begin(),v.end());
        Check(!v.empty()&&v.front()==start&&v.size()==unique.size(),"invalid/duplicate candidate order");
        Check(unique.count(maximum)==1,"maximum not inspected");
        for(int n:v)Check(n>=start&&n<=maximum,"candidate outside requested limits");
        for(int n=start;n<=maximum;n+=2)Check(unique.count(n)==1,"fine candidate omitted");
        Check(D::Reference(maximum)>=maximum&&D::Reference(maximum)<=28,
              "reference exceeds representable bound");
        Check(D::Anchor(maximum)>maximum&&D::Anchor(maximum)<=30,
              "anchor exceeds worker depth30");
        Check(D::Anchor(maximum)>D::Reference(maximum),"reference guard collapsed");
    }
    Check(D::Reference(18)==20&&D::Anchor(18)==22,"default deep guard wrong");
    Check(D::Reference(28)==28&&D::Anchor(28)==30,"explicit max28 requested unsupported32");
    Check(D::Candidates(6,6)==std::vector<int>({6}),"fixed search interval gained candidates");
}
void BoundedSampling(){
    Check(D::PilotFloor(512,1<<24)==4096,"default ambiguity floor not4096");
    Check(D::PilotLimit(512,1<<24)==131072,"default pilot cap not512x256");
    Check(D::PilotFloor(8,64)==64&&D::PilotLimit(8,64)==64,
          "explicit max64 silently raised to4096");
    Check(D::PilotFloor(8192,1<<24)==8192,"explicit initial8192 shrank");
    Check(D::PilotLimit(8,1<<24)==4096,"small initial pilot has no4096 allowance");
    Check(D::PilotLimit(INT_MAX,INT_MAX)==INT_MAX,"large counts overflowed pilot cap");
    Check(D::DoubleCount(INT_MAX/2+1,INT_MAX)==INT_MAX,"doubling overflowed");
    Check(D::DoubleCount(1024,1024)==1024,"completed pilot gained work");
    int n=512;const int limit=D::PilotLimit(n,1<<24);int steps=0;
    while(n<limit){int next=D::DoubleCount(n,limit);Check(next>n&&next<=limit,"pilot growth invalid");n=next;++steps;}
    Check(n==131072&&steps==8,"default doubled sequence not bounded");
    Check(D::DoubleCount(4096,5000)==5000,"non-power explicit cap ignored");
}
}
int main(){try{
    ActualPilotRegressions();AcceptanceGuards();CandidateOrder();BoundedSampling();
    std::cout<<"PASS depth selection: real112/794 scalar regressions; reference guard; direct anchor budget; nonmonotonic coarse/fine; all valid depth ranges; bounded count growth\n";
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}}
