#ifndef MBS_FULLAUTO_ANGULAR_COHORT_H
#define MBS_FULLAUTO_ANGULAR_COHORT_H

#include <algorithm>
#include <cstddef>
#include <utility>
#include <vector>

namespace FullAutoAngularCohort {
struct Plan {
    std::vector<size_t> nodes;
    std::vector<std::pair<size_t,int>> targets;
};

// Scheduling only. Read ORIGINAL counts, never another stencil's proposal.
// Extend a synchronization front by a small, contiguous support window. This
// neither creates a new target nor propagates a distant maximum recursively.
template<class PhiAt>
Plan Build(const std::vector<size_t>& knots,const std::vector<size_t>& evaluated,
        const std::vector<std::vector<int>>& original,const std::vector<size_t>& levels,
        const std::vector<size_t>& stencil,int maximum,PhiAt phiAt){
    Plan plan;plan.nodes=stencil;
    if(stencil.size()<2||stencil.size()>3||knots.empty()||maximum<1)return plan;
    std::vector<size_t> positions;
    for(size_t id:stencil){
        const auto at=std::lower_bound(knots.begin(),knots.end(),id);
        if(at==knots.end()||*at!=id||!std::binary_search(evaluated.begin(),evaluated.end(),id))return plan;
        if(!positions.empty()&&size_t(at-knots.begin())<=positions.back())return plan;
        positions.push_back(at-knots.begin());
    }
    const size_t first=positions.front(),last=positions.back();
    if(last-first>=8)return plan;
    for(size_t level:levels){
        if(level>=original.size())return Plan{stencil,{}};
        int low=maximum,high=0;
        for(size_t id:stencil){
            if(id>=original[level].size())return Plan{stencil,{}};
            low=std::min(low,original[level][id]);high=std::max(high,original[level][id]);
        }
        // Equal-prefix growth is deliberately excluded. Large gaps keep the
        // pre-existing direct-stencil policy without recruiting a halo.
        if(low<high&&low>0&&high<=maximum&&low>=high/2+high%2)
            plan.targets.push_back({level,high});
    }
    if(plan.targets.empty())return plan;
    const int phi=phiAt(stencil.front());
    auto eligible=[&](size_t pos){
        const size_t id=knots[pos];
        if(!std::binary_search(evaluated.begin(),evaluated.end(),id)||phiAt(id)!=phi)return false;
        for(const auto& target:plan.targets){
            if(id>=original[target.first].size())return false;
            const int n=original[target.first][id],largest=target.second;
            if(n<largest/2+largest%2||n>largest)return false;
        }
        return true;
    };
    // Span-two stencils include their intervening knots. A missing, different-
    // phi or out-of-band knot is a barrier, never a point to skip over.
    for(size_t pos=first;pos<=last;++pos)if(!eligible(pos))return Plan{stencil,{}};
    size_t left=first,right=last;
    bool leftOpen=left>0,rightOpen=right+1<knots.size();
    auto grows=[&](size_t pos){
        for(const auto& target:plan.targets)
            if(original[target.first][knots[pos]]<target.second)return true;
        return false;
    };
    while(right-left+1<8&&(leftOpen||rightOpen)){
        if(leftOpen&&!eligible(left-1))leftOpen=false;
        if(rightOpen&&!eligible(right+1))rightOpen=false;
        if(!leftOpen&&!rightOpen)break;
        const bool leftGrows=leftOpen&&grows(left-1),rightGrows=rightOpen&&grows(right+1);
        // Spend the bounded window on new work before already aligned rows.
        // Both eligibility and growth still come solely from original counts.
        const bool takeLeft=leftOpen&&(!rightOpen||(leftGrows!=rightGrows?leftGrows:first-left<=right-last));
        if(takeLeft){
            --left;leftOpen=left>0;
        }else{
            ++right;rightOpen=right+1<knots.size();
        }
    }
    plan.nodes.assign(knots.begin()+left,knots.begin()+right+1);
    return plan;
}
} // namespace FullAutoAngularCohort
#endif
