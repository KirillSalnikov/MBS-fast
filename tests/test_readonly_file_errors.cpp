#include <cerrno>
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <cassert>
#include <iostream>
static int failures=0,openCalls=0,readCalls=0,injectedOpen=0,injectedRead=0;
static int FaultOpen(const char* path,int flags){++openCalls;if(failures){--failures;errno=injectedOpen;return -1;}return ::open(path,flags);}
static ssize_t FaultRead(int fd,void* data,size_t size){++readCalls;if(failures){--failures;errno=injectedRead;return -1;}return ::read(fd,data,size);}
#define open FaultOpen
#define read FaultRead
#include "../src/ReadOnlyFile.h"
#undef open
#undef read
static int statFailures=0,statCalls=0,injectedStat=0;
static int FaultStat(const char* path,struct stat* stamp){++statCalls;if(statFailures){--statFailures;errno=injectedStat;return -1;}return ::stat(path,stamp);}
#define stat(...) FaultStat(__VA_ARGS__)
#include "../src/CheckedFileCache.h"
#undef stat
int main(){
    injectedOpen=EIO;failures=1;errno=0;
    assert(CheckpointRead::Open("/dev/null")<0&&errno==EIO&&openCalls==1);
    injectedOpen=EACCES;failures=1;openCalls=0;
    assert(CheckpointRead::Open("/dev/null")<0&&errno==EACCES&&openCalls==1);
    injectedOpen=EINTR;failures=1;openCalls=0;
    int fd=CheckpointRead::Open("/dev/null");assert(fd>=0&&openCalls>=2);close(fd);
#ifdef O_NOATIME
    injectedOpen=EPERM;failures=1;openCalls=0;
    fd=CheckpointRead::Open("/dev/null");assert(fd>=0&&openCalls==2);close(fd);
#endif
    fd=::open("/dev/zero",O_RDONLY);assert(fd>=0);char value=1;
    injectedRead=EINTR;failures=1;readCalls=0;
    assert(CheckpointRead::Read(fd,&value,1)==1&&value==0&&readCalls==2);
    injectedRead=EIO;failures=1;readCalls=0;
    assert(CheckpointRead::Read(fd,&value,1)<0&&errno==EIO&&readCalls==1);close(fd);
    // Even when a valid packed copy exists, an open EIO must not be mistaken
    // for an absent loose member and silently replaced by that copy.
    char temporary[]="/tmp/mbs-read-errors-XXXXXX";assert(mkdtemp(temporary));const std::string root=temporary;
    assert(mkdir((root+"/jobs").c_str(),0700)==0);const auto job=root+"/jobs/j12";
    assert(mkdir(job.c_str(),0700)==0);{std::ofstream out(job+"/complete.txt");out<<"receipt\n";}
    PackedJobStore::Pack(job);assert(CheckedFileCache::Read(job+"/complete.txt")=="receipt\n");
    injectedOpen=EIO;failures=1;bool rejected=false;
    try{CheckedFileCache::Read(job+"/complete.txt");}catch(const std::exception&){rejected=true;}
    assert(rejected&&failures==0);
    for(const int error:{EIO,EACCES})for(const bool text:{false,true}){
        CheckedFileCache cache;injectedStat=error;statFailures=1;statCalls=0;rejected=false;
        try{if(text)cache.Text(job+"/complete.txt");else cache.Hash(job+"/complete.txt");}
        catch(const std::exception&){rejected=true;}
        assert(rejected&&statFailures==0&&statCalls==1);
    }
    std::cout<<"PASS read-only error injection: only EPERM fallback, EINTR retry, EIO/EACCES preserved\n";
}
