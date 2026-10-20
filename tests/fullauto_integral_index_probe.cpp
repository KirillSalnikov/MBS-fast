// Source selection only: compare the indexed collector with the old scan in
// ordinary and legacy-fallback order. No native optical work is launched.
#include "../src/FullAuto.cpp"
#include <random>
using namespace FullAuto;
namespace {
struct ReferenceRow {size_t level;int sign,count;unsigned seed;std::string folder;};
std::string Reference(const std::string& text,const std::vector<unsigned>& seeds,bool fallback,size_t& inspected){
    std::vector<ReferenceRow> rows;std::map<size_t,int> maximum;
    std::istringstream input(text);std::string line;
    while(std::getline(input,line)){if(line.empty())continue;std::istringstream row(line);ReferenceRow r;
        Require(bool(row>>r.level>>r.sign>>r.count>>r.seed>>r.folder),"bad reference row");
        maximum[r.level]=std::max(maximum[r.level],r.count);rows.push_back(r);}
    std::ostringstream result;
    for(size_t level=0;level<maximum.size();++level){
        if(fallback&&level!=0&&level!=maximum.size()-1)continue;
        for(unsigned seed:seeds)for(int sign:{1,-1}){
            if((level==0||fallback)&&sign==-1)continue;
            std::vector<std::string> folders;
            for(const auto& r:rows){++inspected;if(r.level==level&&r.sign==sign&&r.count==maximum.at(level)&&r.seed==seed)folders.push_back(r.folder);}
            std::sort(folders.begin(),folders.end());folders.erase(std::unique(folders.begin(),folders.end()),folders.end());
            Require(!folders.empty(),"missing reference sources");
            for(const auto& folder:folders)result<<level<<' '<<sign<<' '<<maximum.at(level)<<' '<<seed<<' '<<folder<<'\n';
        }
    }
    return result.str();
}
std::string Indexed(const std::string& text,const std::vector<unsigned>& seeds,bool fallback){
    IntegralSourcePlan plan(text);std::ostringstream result;
    for(size_t level=0;level<plan.levels.size();++level){
        if(fallback&&level!=0&&level!=plan.levels.size()-1)continue;
        for(unsigned seed:seeds)for(int sign:{1,-1}){
            if((level==0||fallback)&&sign==-1)continue;
            for(const auto& folder:plan.Folders(level,sign,seed))
                result<<level<<' '<<sign<<' '<<plan.levels.at(level).count<<' '<<seed<<' '<<folder<<'\n';
        }
    }
    return result.str();
}
void Benchmark(const std::string& name,const std::string& text,const std::vector<unsigned>& seeds,int repeats){
    size_t inspected=0;const auto expected=Reference(text,seeds,false,inspected);
    Require(Indexed(text,seeds,false)==expected,"source order changed");
    size_t unused=0;Require(Indexed(text,seeds,true)==Reference(text,seeds,true,unused),"legacy fallback source order changed");
    double old=0,now=0;size_t observed=0;
    for(int k=0;k<repeats;++k){
        auto run=[&](bool indexed){const double begun=Seconds();size_t counter=0;
            const auto result=indexed?Indexed(text,seeds,false):Reference(text,seeds,false,counter);
            const auto elapsed=Seconds()-begun;Require(result==expected,"benchmark sources differ");observed+=result.size();return elapsed;};
        if(k%2){now+=run(true);old+=run(false);}else{old+=run(false);now+=run(true);}
    }
    std::cout<<std::setprecision(9)<<"{\"case\":"<<Quoted(name)<<",\"source_bytes\":"<<text.size()
             <<",\"old_row_inspections\":"<<inspected<<",\"new_row_inspections\":"<<std::count(text.begin(),text.end(),'\n')
             <<",\"repeats\":"<<repeats<<",\"old_seconds\":"<<old/repeats<<",\"new_seconds\":"<<now/repeats
             <<",\"speedup\":"<<old/now<<",\"observed_bytes\":"<<observed<<"}\n";
}
}
int main(int argc,char** argv){try{
    const std::vector<unsigned> seeds(production,production+8);
    if(argc==2){Benchmark(argv[1],Read(argv[1]),seeds,101);return 0;}
    Require(argc==1,"optional production source list expected");
    std::vector<std::string> rows;
    for(int level=0;level<12;++level)for(int sign:{1,-1})for(unsigned seed:seeds)
        for(int count:{32,64,256,128,1024,512})for(int angle=0;angle<128;++angle){
            std::ostringstream row;row<<level<<' '<<sign<<' '<<count<<' '<<seed
                <<" /immutable/optical/job/level"<<level<<"/sign"<<sign<<"/seed"<<seed
                <<"/prefix"<<count<<"/angular_group"<<angle%16<<'\n';rows.push_back(row.str());
        }
    std::mt19937 random(65537);std::shuffle(rows.begin(),rows.end(),random);
    std::ostringstream text;for(const auto& row:rows)text<<row;
    Benchmark("shuffled_12_levels_duplicate_angular_groups",text.str(),seeds,7);
    for(const auto& bad:{std::string(),std::string("1 1 1 17 /a\n"),std::string("0 1 0 17 /a\n"),
                         std::string("0 0 1 17 /a\n"),std::string("0 1 1 17 /a extra\n"),
                         std::string("0 1 1 17 /a\n2 1 1 17 /b\n")}){
        bool rejected=false;try{IntegralSourcePlan plan(bad);}catch(const std::exception&){rejected=true;}
        Require(rejected,"invalid source plan accepted");
    }
    IntegralSourcePlan missing("0 1 32 17 /a\n0 1 64 23 /b\n");
    bool rejected=false;try{missing.Folders(0,1,17);}catch(const std::exception&){rejected=true;}
    Require(rejected,"missing seed at largest prefix was replaced with a smaller prefix");
    char temporary[]="/tmp/mbs_integral_read_XXXXXX";const char* directory=mkdtemp(temporary);
    Require(directory,"cannot create integral read fixture");
    const std::string path=std::string(directory)+"/complete.txt",first(4096,'a'),second(4096,'b');
    Atomic(path,first);IntegralReader reader;const auto opens=ValidationFiles().opens;
    for(int i=0;i<100;++i)Require(reader.ReadSource(path)==first,"cached integral source differs");
    Require(ValidationFiles().opens==opens+1,"unchanged integral source reopened");
    struct stat stamp;Require(stat(path.c_str(),&stamp)==0,"cannot stat integral fixture");
    {std::ofstream file(path);file<<second;} // Same inode and size, then restore mtime.
    const timespec times[2]={stamp.st_atim,stamp.st_mtim};
    Require(utimensat(AT_FDCWD,path.c_str(),times,0)==0,"cannot restore integral fixture mtime");
    Require(reader.ReadSource(path)==second&&reader.files.at(path)==Hash(second),"same-size changed source accepted from cache");
    Require(ValidationFiles().Bytes()<=ValidationFiles().Limit(),"source text cache exceeded RAM limit");
    unlink(path.c_str());rmdir(directory);
    std::cout<<"PASS: exact source ordering; largest prefix; missing seed and malformed list rejection; 100 reads use one open; same-size/mtime mutation detected\n";
    return 0;
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}}
