#include <time.h>
// These capacity/hit tests exercise mature immutable files. A deterministic
// clock offset leaves timestamp quarantine itself to test_checked_file_cache.
static int MatureTestClock(clockid_t kind,struct timespec* result){
 const int rc=::clock_gettime(kind,result);if(rc==0)result->tv_sec+=2;return rc;
}
#define clock_gettime(...) MatureTestClock(__VA_ARGS__)
#include "../src/CheckedFileCache.h"
#undef clock_gettime
#include <cassert>
#include <fcntl.h>
#include <unistd.h>
#include <iostream>
static void Write(const std::string& p,const std::string& s){std::ofstream f(p,std::ios::binary);f<<s;assert(bool(f));}
template<class F>void Reject(F f){bool rejected=false;try{f();}catch(const std::exception&){rejected=true;}assert(rejected);}
int main(){
 char temp[]="/tmp/mbs-checked-large-XXXXXX";assert(mkdtemp(temp));const std::string r=temp,p=r+"/manifest",q=r+"/replacement";
 CheckedFileCache c(4*1024*1024);std::string data(200000,'a');Write(p,data);const auto hash=c.Hash(p);
 assert(c.VerifiedText(p,hash)==data);const auto opened=c.opens;
 assert(c.VerifiedText(p,hash)==data&&c.Hash(p)==hash&&c.opens==opened);
 Reject([&]{c.VerifiedText(p,"wrong");});
 struct stat old{};assert(stat(p.c_str(),&old)==0);usleep(2000);data[123]='b';Write(p,data);
 const timespec times[2]={old.st_atim,old.st_mtim};assert(utimensat(AT_FDCWD,p.c_str(),times,0)==0);
 Reject([&]{c.VerifiedText(p,hash);});assert(c.Text(p)==data);
 Write(q,std::string(data.size(),'c'));assert(rename(q.c_str(),p.c_str())==0);assert(c.Text(p)==std::string(data.size(),'c'));
 assert(link(p.c_str(),q.c_str())==0);assert(c.Text(q)==c.Text(p));unlink(q.c_str());
 // Large-manifest churn never exceeds the original budget and does not evict
 // the small-digest region. Oversized manifests are returned but not retained.
 const auto tiny=r+"/tiny";Write(tiny,"small");const auto tinyHash=c.Hash(tiny);
 for(int i=0;i<20;++i){const auto path=r+"/m"+std::to_string(i);Write(path,std::string(200000,char('a'+i)));c.Text(path);assert(c.Bytes()<=c.Limit());}
 const auto prior=c.opens;assert(c.Hash(tiny)==tinyHash&&c.opens==prior);
 Write(p,std::string(2*1024*1024,'z'));assert(c.Text(p).size()==2*1024*1024);const auto big=c.opens;c.Text(p);assert(c.opens==big+1);

 // Lossless common-prefix storage covers actual long directory lists and all
 // byte-level endings. An empty line must simply disable prefix factoring.
 const std::string prefix="/home/user/campaign/"+std::string(180,'p');
 for(int mode=0;mode<6;++mode){std::string manifest;
  for(int i=0;i<800;++i){manifest+=prefix+"/job"+std::to_string(i)+"\t0123456789abcdef";
   if(mode==2)manifest+=std::string("\0x",2);
   if(mode==3)manifest+='\r';manifest+='\n';}
  if(mode==1)manifest.pop_back();
  if(mode==4)manifest+=prefix; // empty encoded final suffix, no newline
  if(mode==5)manifest.insert(manifest.find('\n')+1,"\n");
  Write(p,manifest);const auto expected=c.Hash(p);assert(c.VerifiedText(p,expected)==manifest);const auto oldOpens=c.opens;
  for(int repeat=0;repeat<3;++repeat)assert(c.VerifiedText(p,expected)==manifest&&c.Hash(p)==expected&&c.opens==oldOpens);
  assert(c.Bytes()<=c.Limit());}

 // No large files: preserve all of the old capacity, including a working set
 // between 75% and 100% of the limit (a permanently reserved quarter regresses).
 {CheckedFileCache small(4*1024*1024);const std::string contents(8192,'s');
  for(int i=0;i<420;++i){const auto path=r+"/s"+std::to_string(i);Write(path,contents);assert(small.Text(path)==contents);}
  assert(small.Bytes()>3*1024*1024&&small.Bytes()<=small.Limit());const auto before=small.opens;
  for(int i=0;i<420;++i)assert(small.Text(r+"/s"+std::to_string(i))==contents);
  assert(small.opens==before);for(int i=0;i<420;++i)unlink((r+"/s"+std::to_string(i)).c_str());}
 assert(BoundedCache<std::string>::RetainedBytes("x",std::numeric_limits<size_t>::max())==std::numeric_limits<size_t>::max());
 {BoundedCache<std::string> bound(1024);bound.Put("x",std::make_shared<std::string>("ok"),std::numeric_limits<size_t>::max());assert(bound.Bytes()==0);}
 data.clear();for(int i=0;i<800;++i)data+=prefix+"/packed"+std::to_string(i)+"\t0123456789abcdef\n";
 // Packed text follows the same bound and the member generation invalidates
 // every retained copy after any unexplained archive mutation.
 const auto jobs=r+"/jobs",job=jobs+"/n1";assert(mkdir(jobs.c_str(),0700)==0&&mkdir(job.c_str(),0700)==0);
 Write(job+"/complete.txt","complete\n");Write(job+"/angular_sources.tsv",data);PackedJobStore::Pack(job);
 const auto packedHash=c.Hash(job+"/angular_sources.tsv");assert(c.VerifiedText(job+"/angular_sources.tsv",packedHash)==data);const auto packedOpens=c.opens;
 assert(c.VerifiedText(job+"/angular_sources.tsv",packedHash)==data&&c.opens==packedOpens);
 PackedJobStore::Member m;assert(PackedJobStore::Locate(job+"/angular_sources.tsv",m));int fd=open(m.archive.c_str(),O_RDWR);assert(fd>=0);const char bad='!';assert(pwrite(fd,&bad,1,m.offset)==1);close(fd);
 Reject([&]{c.VerifiedText(job+"/angular_sources.tsv",packedHash);});assert(c.Bytes()<=c.Limit());
 unlink(m.archive.c_str());rmdir((r+"/job_packs").c_str());rmdir(jobs.c_str());
 for(int i=0;i<20;++i)unlink((r+"/m"+std::to_string(i)).c_str());unlink(p.c_str());unlink(tiny.c_str());rmdir(r.c_str());
 std::cout<<"PASS large verified text: bounded split, hit, ctime/inode edits, bad digest, eviction isolation, oversize, packed corruption\n";
}
