// Read-only, bounded-memory validation of a memory-store checkpoint.
#include <cstdio>
#include <cstdint>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <iostream>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
static const uint64_t basis=14695981039346656037ULL,prime=1099511628211ULL;
static void need(bool v,const char* msg){if(!v)throw std::runtime_error(msg);}
struct Reader{
 FILE* f;uint64_t length,position=0,hash=basis;
 Reader(FILE* file,uint64_t size):f(file),length(size){}
 void bytes(char* p,size_t n){need(position<=length&&n<=length-position,"truncated checkpoint");need(fread(p,1,n,f)==n,"read failed");position+=n;for(size_t i=0;i<n;++i)hash=(hash^static_cast<unsigned char>(p[i]))*prime;}
 uint64_t number(){char p[8];bytes(p,8);uint64_t v=0;for(int i=0;i<8;++i)v|=uint64_t(static_cast<unsigned char>(p[i]))<<(8*i);return v;}
 std::string text(size_t limit){auto n=number();need(n<=limit,"string bound exceeded");std::string s(n,'\0');if(n)bytes(&s[0],n);return s;}
};
int main(int argc,char** argv){FILE* f=nullptr;try{
 need(argc==2,"usage: check_fullauto_checkpoint FILE");int fd=open(argv[1],O_RDONLY|O_CLOEXEC|O_NOFOLLOW);need(fd>=0,"cannot open checkpoint");struct stat st{};need(fstat(fd,&st)==0&&S_ISREG(st.st_mode)&&st.st_size>=8,"invalid checkpoint file");f=fdopen(fd,"rb");need(f!=nullptr,"fdopen failed");Reader r{f,uint64_t(st.st_size)};
 need(r.text(64)=="MBS_FULLAUTO_MEMORY_STORE_V1","unsupported checkpoint schema");const auto root=r.text(4096);r.text(16*1024*1024);auto epoch=r.number(),count=r.number();need(count>=1&&count<=r.length/40,"invalid record count");std::string previous;uint64_t files=0,dirs=0;char buffer[65536];
 for(uint64_t i=0;i<count;++i){const auto path=r.text(4096);need(i==0?path.empty():(!path.empty()&&path[0]=='/'&&path>previous),"invalid path order");need(path.find("/../")==std::string::npos&&path.find("/./")==std::string::npos&&path.find('\0')==std::string::npos,"invalid path");previous=path;auto directory=r.number(),generation=r.number(),size=r.number(),expected=r.number();need(directory<=1&&generation>0&&generation<=epoch,"invalid generation");need(i!=0||directory==1,"root is not directory");uint64_t hash=basis;
  if(directory){need(size==0&&expected==basis,"invalid directory");++dirs;}
  else{++files;need(size<=r.length-r.position,"invalid payload size");for(uint64_t done=0;done<size;){size_t n=size_t(std::min<uint64_t>(sizeof(buffer),size-done));r.bytes(buffer,n);for(size_t j=0;j<n;++j)hash=(hash^static_cast<unsigned char>(buffer[j]))*prime;done+=n;}need(hash==expected,"member checksum mismatch");}
 }
 need(r.length-r.position==8,"trailing or missing data");const auto expected=r.hash;need(r.number()==expected,"checkpoint checksum mismatch");struct stat end{};need(fstat(fd,&end)==0&&end.st_size==st.st_size&&end.st_mtime==st.st_mtime,"checkpoint changed while reading");fclose(f);f=nullptr;
 std::cout<<"{\"passed\":true,\"bytes\":"<<st.st_size<<",\"files\":"<<files<<",\"directories\":"<<dirs<<",\"epoch\":"<<epoch<<",\"checkpoint_fnv1a64\":\""<<std::hex<<expected<<"\",\"read_only\":true}\n";return 0;
}catch(const std::exception& e){if(f)fclose(f);std::cerr<<e.what()<<'\n';return 1;}}
