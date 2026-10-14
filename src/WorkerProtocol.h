#pragma once
// Private controller/worker transport. Length-delimited arguments preserve
// spaces and newlines without shell interpretation; payloads are bounded.
#include <sys/socket.h>
#include <unistd.h>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>
namespace WorkerProtocol {
inline bool Transfer(int fd,void* buffer,size_t bytes,bool sendData) {
    char* p=static_cast<char*>(buffer);
    while(bytes){ssize_t n=sendData?send(fd,p,bytes,MSG_NOSIGNAL):recv(fd,p,bytes,0);
        if(n<0&&errno==EINTR)continue;if(n<=0)return false;p+=n;bytes-=size_t(n);}
    return true;
}
inline bool Send(int fd,const std::vector<std::string>& words){
    uint32_t count=uint32_t(words.size());if(count>4096)return false;
    if(!Transfer(fd,&count,sizeof(count),true))return false;
    for(const auto& s:words){uint32_t n=uint32_t(s.size());if(s.size()>1048576)return false;
        if(!Transfer(fd,&n,sizeof(n),true)||!Transfer(fd,const_cast<char*>(s.data()),n,true))return false;}
    return true;
}
inline bool Receive(int fd,std::vector<std::string>& words){
    words.clear();uint32_t count=0;if(!Transfer(fd,&count,sizeof(count),false)||count>4096)return false;
    size_t total=0;for(uint32_t i=0;i<count;++i){uint32_t n=0;
        if(!Transfer(fd,&n,sizeof(n),false)||n>1048576||(total+=n)>16777216)return false;
        std::string s(n,'\0');if(n&&!Transfer(fd,&s[0],n,false))return false;
        if(s.find('\0')!=std::string::npos)return false;words.push_back(std::move(s));}
    return true;
}
inline bool Reply(int fd,int code){int32_t value=code;return Transfer(fd,&value,sizeof(value),true);}
inline bool Result(int fd,int& code){int32_t value=0;if(!Transfer(fd,&value,sizeof(value),false))return false;code=value;return true;}
}
