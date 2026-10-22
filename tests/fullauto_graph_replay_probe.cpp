// Read-only replay of immutable saved graphs. The only output is stdout.
#ifndef MBS_TEST_FULLAUTO_PATH
#define MBS_TEST_FULLAUTO_PATH "../src/FullAuto.cpp"
#endif
#include MBS_TEST_FULLAUTO_PATH
using namespace FullAuto;
int main(int argc,char** argv){try{
    Require(argc==2,"GROUPS_JSONL_AS_TSV");
    std::ifstream input(argv[1]);Require(bool(input),"cannot read replay list");
    std::vector<std::vector<std::string>> groups;std::string line;
    while(std::getline(input,line)){std::istringstream row(line);std::vector<std::string> roots;std::string path;
        while(row>>path)roots.push_back(path);if(!roots.empty())groups.push_back(roots);}
    Require(!groups.empty(),"empty graph replay");ValidationBatch batch("read-only saved-graph replay");
    for(size_t k=0;k<groups.size();++k){const auto begin=Seconds();const auto checked=batch.checked,reused=batch.reused;
        ValidateStoredRoots(groups[k],4);std::string hashes;
        for(const auto& root:groups[k])hashes+=root+" "+FileHash(root+"/result/result.dat")+"\n";
        std::cout<<"{\"group\":"<<k<<",\"roots\":"<<groups[k].size()<<",\"seconds\":"<<Number(Seconds()-begin)
            <<",\"checked_nodes\":"<<batch.checked-checked<<",\"reused_graphs\":"<<batch.reused-reused
            <<",\"memo_bytes\":"<<batch.compounds.Bytes()<<",\"matrix_digests\":"<<Quoted(Hash(hashes))<<"}"<<std::endl;
    }
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
