#pragma once
// Optional working-data index. SQLite is loaded only after the RAM metadata
// budget is exhausted; this is a disposable index, never a recovery checkpoint.
#include <dlfcn.h>
#include <algorithm>
#include <cstdlib>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstdint>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>
#ifndef MBS_MEMORY_CATALOG_EVENT
#define MBS_MEMORY_CATALOG_EVENT(name,path) ((void)0)
#endif

namespace FullAutoMemoryCatalog {
struct sqlite3;struct sqlite3_stmt;
inline void Require(bool okay,const std::string& message){if(!okay)throw std::runtime_error("fullauto working catalog: "+message);}
class Api {
    void* library=nullptr;
    template<class T> void Symbol(T& target,const char* name){target=reinterpret_cast<T>(dlsym(library,name));Require(target!=nullptr,std::string("missing SQLite symbol ")+name);}
public:
    int (*open_v2)(const char*,sqlite3**,int,const char*)=nullptr;
    int (*close_v2)(sqlite3*)=nullptr;
    int (*exec)(sqlite3*,const char*,int(*)(void*,int,char**,char**),void*,char**)=nullptr;
    void (*free)(void*)=nullptr;
    const char* (*errmsg)(sqlite3*)=nullptr;
    int (*prepare_v2)(sqlite3*,const char*,int,sqlite3_stmt**,const char**)=nullptr;
    int (*finalize)(sqlite3_stmt*)=nullptr;
    int (*reset)(sqlite3_stmt*)=nullptr;
    int (*clear_bindings)(sqlite3_stmt*)=nullptr;
    int (*step)(sqlite3_stmt*)=nullptr;
    int (*bind_text)(sqlite3_stmt*,int,const char*,int,void(*)(void*))=nullptr;
    int (*bind_int64)(sqlite3_stmt*,int,int64_t)=nullptr;
    const unsigned char* (*column_text)(sqlite3_stmt*,int)=nullptr;
    int (*column_bytes)(sqlite3_stmt*,int)=nullptr;
    int64_t (*column_int64)(sqlite3_stmt*,int)=nullptr;
    explicit Api(const char* explicitLibrary=nullptr){
        if(explicitLibrary)library=dlopen(explicitLibrary,RTLD_NOW|RTLD_LOCAL);
        else{library=dlopen("libsqlite3.so.0",RTLD_NOW|RTLD_LOCAL);if(!library)library=dlopen("libsqlite3.so",RTLD_NOW|RTLD_LOCAL);}
        Require(library!=nullptr,"SQLite runtime is unavailable; install libsqlite3 or increase the fullauto memory limit");
        try{
#define MBS_SQLITE_SYMBOL(name) Symbol(name,"sqlite3_" #name)
            MBS_SQLITE_SYMBOL(open_v2);MBS_SQLITE_SYMBOL(close_v2);MBS_SQLITE_SYMBOL(exec);MBS_SQLITE_SYMBOL(free);MBS_SQLITE_SYMBOL(errmsg);
            MBS_SQLITE_SYMBOL(prepare_v2);MBS_SQLITE_SYMBOL(finalize);MBS_SQLITE_SYMBOL(reset);MBS_SQLITE_SYMBOL(clear_bindings);MBS_SQLITE_SYMBOL(step);
            MBS_SQLITE_SYMBOL(bind_text);MBS_SQLITE_SYMBOL(bind_int64);MBS_SQLITE_SYMBOL(column_text);MBS_SQLITE_SYMBOL(column_bytes);MBS_SQLITE_SYMBOL(column_int64);
#undef MBS_SQLITE_SYMBOL
        }catch(...){dlclose(library);library=nullptr;throw;}
    }
    ~Api(){if(library)dlclose(library);}
    Api(const Api&)=delete;Api& operator=(const Api&)=delete;
};
inline Api& Library(){static Api library;return library;}
struct Row {bool directory=false;uint64_t size=0,digest=UINT64_C(14695981039346656037),offset=0,generation=0;};
class Catalog {
    Api& api;sqlite3* db=nullptr;std::string path;int backing=-1;bool healthy=true;
    sqlite3_stmt *lookup=nullptr,*put=nullptr,*remove=nullptr;
    void Check(int code,const char* operation)const{Require(code==0,std::string(operation)+": "+(db?api.errmsg(db):"cannot open database"));}
    void Execute(const std::string& sql){char* error=nullptr;const int code=api.exec(db,sql.c_str(),nullptr,nullptr,&error);std::string message=error?error:"";if(error)api.free(error);Require(code==0,"SQL operation failed: "+message);}
    sqlite3_stmt* Prepare(const char* sql){sqlite3_stmt* statement=nullptr;Check(api.prepare_v2(db,sql,-1,&statement,nullptr),"prepare");return statement;}
    void Text(sqlite3_stmt* statement,int field,const std::string& text)const{Check(api.bind_text(statement,field,text.data(),int(text.size()),reinterpret_cast<void(*)(void*)>(-1)),"bind path");}
    void Number(sqlite3_stmt* statement,int field,uint64_t value)const{int64_t bits;std::memcpy(&bits,&value,sizeof(bits));Check(api.bind_int64(statement,field,bits),"bind value");}
    uint64_t Number(sqlite3_stmt* statement,int column)const{const auto bits=api.column_int64(statement,column);uint64_t value;std::memcpy(&value,&bits,sizeof(value));return value;}
    Row Decode(sqlite3_stmt* statement,int first)const{Row row;const auto kind=Number(statement,first);Require(kind<=1,"invalid stored directory flag");row.directory=kind==1;
        row.size=Number(statement,first+1);row.digest=Number(statement,first+2);row.offset=Number(statement,first+3);row.generation=Number(statement,first+4);return row;}
    static std::string Parent(const std::string& name){const auto slash=name.find_last_of('/');return slash==0?"/":name.substr(0,slash);}
    void Reset(sqlite3_stmt* statement)const{api.reset(statement);api.clear_bindings(statement);}
    struct Statement {Api& api;sqlite3_stmt* value;Statement(Api& a,sqlite3_stmt* s):api(a),value(s){}~Statement(){api.finalize(value);}};
public:
    Catalog(const std::string& directory,size_t cacheBytes):api(Library()){
        std::string name=directory+"/mbs-fullauto-catalog-"+std::to_string(getpid())+"-XXXXXX";std::vector<char> buffer(name.begin(),name.end());buffer.push_back('\0');
        const int fd=mkstemp(buffer.data());Require(fd>=0,"cannot create temporary index");close(fd);path=buffer.data();
        try{
            Check(api.open_v2(path.c_str(),&db,2|4|0x8000|0x40000,nullptr),"open");
            Execute("PRAGMA page_size=4096; PRAGMA journal_mode=OFF; PRAGMA synchronous=OFF; PRAGMA locking_mode=EXCLUSIVE; PRAGMA temp_store=MEMORY; PRAGMA mmap_size=0; PRAGMA cache_spill=ON; PRAGMA cache_size=-"+std::to_string(std::max<size_t>(64,cacheBytes/1024))+";");
            Execute("CREATE TABLE records(path TEXT PRIMARY KEY,parent TEXT NOT NULL,dir INTEGER NOT NULL,size INTEGER NOT NULL,digest INTEGER NOT NULL,off INTEGER NOT NULL,generation INTEGER NOT NULL) WITHOUT ROWID; CREATE INDEX children ON records(parent,path); BEGIN EXCLUSIVE;");
            lookup=Prepare("SELECT dir,size,digest,off,generation FROM records WHERE path=?1");
            put=Prepare("INSERT OR REPLACE INTO records(path,parent,dir,size,digest,off,generation) VALUES(?1,?2,?3,?4,?5,?6,?7)");
            remove=Prepare("DELETE FROM records WHERE path=?1");
            // The connection keeps one exclusive working transaction open.
            // With rollback/WAL journals disabled it never reopens this inode.
            // Unlink only AFTER BEGIN has checked the pathname, so SIGKILL
            // automatically releases both data and metadata spill files.
            backing=::open(path.c_str(),O_RDONLY|O_CLOEXEC|O_NOFOLLOW);Require(backing>=0,"cannot retain temporary catalog descriptor");
            Require(::unlink(path.c_str())==0,"cannot unlink temporary catalog");
        }catch(...){if(lookup)api.finalize(lookup);if(put)api.finalize(put);if(remove)api.finalize(remove);if(db)api.close_v2(db);if(backing>=0)::close(backing);unlink(path.c_str());throw;}
    }
    ~Catalog(){if(lookup)api.finalize(lookup);if(put)api.finalize(put);if(remove)api.finalize(remove);if(db)api.close_v2(db);if(backing>=0)::close(backing);if(!path.empty())unlink(path.c_str());}
    Catalog(const Catalog&)=delete;Catalog& operator=(const Catalog&)=delete;
    bool Get(const std::string& key,Row& row)const{
        Reset(lookup);Text(lookup,1,key);const int result=api.step(lookup);if(result==100){row=Decode(lookup,0);Reset(lookup);return true;}
        Reset(lookup);Require(result==101,"catalog lookup failed");return false;
    }
    void Put(const std::string& key,const Row& row){try{MBS_MEMORY_CATALOG_EVENT("before_put",key);Reset(put);Text(put,1,key);Text(put,2,Parent(key));Number(put,3,row.directory?1:0);Number(put,4,row.size);Number(put,5,row.digest);Number(put,6,row.offset);Number(put,7,row.generation);
        const int result=api.step(put);Reset(put);Require(result==101,std::string("catalog write failed: ")+api.errmsg(db));}catch(...){healthy=false;throw;}}
    void Erase(const std::string& key){try{MBS_MEMORY_CATALOG_EVENT("before_erase",key);Reset(remove);Text(remove,1,key);const int result=api.step(remove);Reset(remove);Require(result==101,"catalog erase failed");}catch(...){healthy=false;throw;}}
    bool Healthy()const{return healthy;}
    template<class Visitor> void All(Visitor visitor)const{Statement statement(api,const_cast<Catalog*>(this)->Prepare("SELECT path,dir,size,digest,off,generation FROM records ORDER BY path"));Visit(statement.value,visitor);}
    template<class Visitor> void Children(const std::string& parent,bool recursive,Visitor visitor)const{
        Statement statement(api,const_cast<Catalog*>(this)->Prepare(recursive?"SELECT path,dir,size,digest,off,generation FROM records WHERE path>=?1 AND path<?2 AND path<>?3 ORDER BY path":"SELECT path,dir,size,digest,off,generation FROM records WHERE parent=?1 AND path<>?2 ORDER BY path"));
        if(recursive){const auto prefix=parent=="/"?parent:parent+"/";std::string upper=prefix;++upper.back();Text(statement.value,1,prefix);Text(statement.value,2,upper);Text(statement.value,3,parent);}
        else{Text(statement.value,1,parent);Text(statement.value,2,parent);}
        Visit(statement.value,visitor);
    }
    template<class Visitor> void Visit(sqlite3_stmt* statement,Visitor visitor)const{
        int result;while((result=api.step(statement))==100){const auto bytes=api.column_bytes(statement,0);const auto text=api.column_text(statement,0);Require(bytes>=1&&bytes<=4096&&text,"invalid stored path");
            const std::string key(reinterpret_cast<const char*>(text),size_t(bytes));Require(key.find('\0')==std::string::npos,"invalid stored path bytes");visitor(key,Decode(statement,1));}
        Require(result==101,"catalog scan failed");
    }
    uint64_t DiskBytes()const{struct stat st{};return backing>=0&&fstat(backing,&st)==0?uint64_t(st.st_size):0;}
};
}
