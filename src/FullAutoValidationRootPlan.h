#pragma once
#include <algorithm>
#include <cstddef>
#include <map>
#include <string>
#include <vector>

// This plan changes only validation order. Every requested root is retained;
// the actual validator must authenticate files and every incoming source hash.
struct ValidationRootPlan {
    std::vector<size_t> first,second;
    size_t manifestBytes=0,edges=0;
    bool incomplete=false;
    const char* reason="independent roots";
    bool Ordered()const{return !second.empty();}
};
// Reader(i, remaining budget, text) returns false to stop collecting hints;
// already discovered acyclic links remain usable. An absent optional angular
// manifest returns true with empty text. The reader must bound the file before
// reading; the planner additionally checks returned byte counts.
template<class Reader>
ValidationRootPlan PlanValidationRoots(const std::vector<std::string>& roots,Reader reader,
        size_t maxRoots=512,size_t byteLimit=8*1024*1024,size_t edgeLimit=8192){
    ValidationRootPlan plan;
    auto fallback=[&](const char* reason){plan.first.clear();plan.second.clear();
        for(size_t i=0;i<roots.size();++i)plan.first.push_back(i);plan.reason=reason;return plan;};
    if(roots.size()<2)return fallback("fewer than two roots");
    if(roots.size()>maxRoots){plan.incomplete=true;return fallback("root planning limit");}
    std::map<std::string,size_t> ids;
    for(size_t i=0;i<roots.size();++i)if(!ids.emplace(roots[i],i).second)return fallback("duplicate roots");
    std::vector<std::vector<size_t>> children(roots.size());std::vector<size_t> indegree(roots.size(),0);
    std::vector<bool> inspected(roots.size(),false);
    for(size_t i=0;i<roots.size();++i){std::string text;
        // A bounded ordering hint need not discard already discovered links.
        // Every root still reaches the authoritative validator. An unread
        // manifest is unknown, never evidence that its root is independent.
        if(!reader(i,byteLimit-plan.manifestBytes,text)){plan.incomplete=true;break;}
        if(text.size()>byteLimit-plan.manifestBytes){plan.incomplete=true;break;}
        plan.manifestBytes+=text.size();size_t begin=0;
        while(begin<text.size()){
            const size_t end=text.find('\n',begin),stop=end==std::string::npos?text.size():end;
            const size_t tab=text.find('\t',begin);
            if(tab==std::string::npos||tab>=stop||tab==begin||tab+1==stop)return fallback("malformed source manifest");
            const auto found=ids.find(text.substr(begin,tab-begin));
            if(found!=ids.end()){
                const size_t child=found->second;
                if(child==i)return fallback("cyclic root links");
                if(std::find(children[i].begin(),children[i].end(),child)==children[i].end()){
                    if(plan.edges==edgeLimit){plan.incomplete=true;break;}
                    children[i].push_back(child);++indegree[child];++plan.edges;
                }
            }
            if(end==std::string::npos)break;begin=end+1;
        }
        if(plan.incomplete)break;
        inspected[i]=true;
    }
    std::vector<size_t> ready;
    for(size_t i=0;i<roots.size();++i)if(!indegree[i])ready.push_back(i);
    for(size_t next=0;next<ready.size();++next)for(size_t child:children[ready[next]])if(--indegree[child]==0)ready.push_back(child);
    if(ready.size()!=roots.size())return fallback("cyclic root links");
    // Reconstruct incoming links after Kahn's cycle check. Unknown roots go
    // first: they may contain large already requested prefixes. Known child
    // roots follow, with their own direct checks still mandatory.
    std::fill(indegree.begin(),indegree.end(),0);
    for(const auto& links:children)for(size_t child:links)++indegree[child];
    if(plan.incomplete)for(size_t i=0;i<roots.size();++i)if(!inspected[i])plan.first.push_back(i);
    for(size_t i=0;i<roots.size();++i){
        if(plan.incomplete&&!inspected[i])continue;
        if(!indegree[i])plan.first.push_back(i);else plan.second.push_back(i);
    }
    if(plan.incomplete){
        if(plan.second.empty())return fallback("bounded root plan without covered roots");
        plan.reason="bounded partial root links; unread roots first";return plan;
    }
    if(plan.second.empty())return fallback("independent roots");
    plan.reason="outer roots before their requested descendants";return plan;
}
