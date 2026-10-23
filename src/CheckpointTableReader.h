#pragma once
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <locale.h>
#include <stdexcept>
#include <string>
#include <utility>

// Internal checkpoint tables are decimal, whitespace-separated and have one
// header line. Parse directly from the read buffer, without constructing a
// locale-aware C++ stream for every row and every number. The C locale is local
// to this reader; changing the process locale must not change checkpoint data.
class CheckpointTableReader {
    struct NumericLocale {
        locale_t value;
        NumericLocale():value(newlocale(LC_NUMERIC_MASK,"C",nullptr)){
            if(!value)throw std::runtime_error("fullauto: cannot create checkpoint numeric locale");
        }
        ~NumericLocale(){freelocale(value);}
    };
    static locale_t Locale(){static thread_local NumericLocale numeric;return numeric.value;}
    static bool Space(char c){return c==' '||c=='\t'||c=='\r'||c=='\v'||c=='\f';}
    std::string text;
    const char *next,*end,*position=nullptr,*rowEnd=nullptr;
public:
    explicit CheckpointTableReader(std::string data):text(std::move(data)){
        end=text.data()+text.size();
        const auto header=static_cast<const char*>(std::memchr(text.data(),'\n',text.size()));
        next=header?header+1:end;
    }
    CheckpointTableReader(const CheckpointTableReader&)=delete;
    CheckpointTableReader& operator=(const CheckpointTableReader&)=delete;
    bool Next(){
        while(next<end){
            position=next;const auto newline=static_cast<const char*>(std::memchr(next,'\n',size_t(end-next)));
            rowEnd=newline?newline:end;next=newline?newline+1:end;
            if(position!=rowEnd)return true;
        }
        position=rowEnd=end;return false;
    }
    bool Number(double& value){
        if(!position)return false;
        while(position<rowEnd&&Space(*position))++position;
        if(position==rowEnd)return false;
        const char* start=position;
        if(*start=='+'||*start=='-')++start;
        if(start==rowEnd||(!(*start>='0'&&*start<='9')&&*start!='.'))return false;
        // Hex floats and locale-specific tokens are not checkpoint columns.
        if(rowEnd-start>=2&&start[0]=='0'&&(start[1]=='x'||start[1]=='X'))return false;
        char* after=nullptr;const double parsed=strtod_l(position,&after,Locale());
        if(after==position||after>rowEnd||!std::isfinite(parsed)
                ||(after<rowEnd&&!Space(*after)))return false;
        position=after;value=parsed;return true;
    }
    bool End(){
        if(!position)return false;
        while(position<rowEnd&&Space(*position))++position;
        return position==rowEnd;
    }
};
