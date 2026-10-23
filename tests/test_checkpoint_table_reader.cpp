#include "../src/CheckpointTableReader.h"
#include <array>
#include <cassert>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <vector>

int main(){
    // All serialized finite bit patterns, including signed zero, subnormals,
    // extreme exponents and values next to halfway rounding boundaries.
    std::vector<double> values={0.,-0.,1.,-1.,std::numeric_limits<double>::denorm_min(),
        std::numeric_limits<double>::min(),std::numeric_limits<double>::max(),
        std::nextafter(1.,2.),std::nextafter(1.,0.)};
    std::mt19937_64 random(20261023);
    while(values.size()<40000){const uint64_t bits=random();double v;std::memcpy(&v,&bits,sizeof v);if(std::isfinite(v))values.push_back(v);}
    for(int precision:{8,12,17}){
        std::ostringstream text;text<<"theta area all16\n"<<std::setprecision(precision);
        for(size_t i=0;i<values.size();++i)text<<values[i]<<(i%16==15?'\n':'\t');
        CheckpointTableReader reader(text.str());std::istringstream reference(text.str());std::string header;std::getline(reference,header);
        size_t i=0;while(reader.Next())for(int column=0;column<16&&i<values.size();++column,++i){
            double a,b;assert(reader.Number(a));assert(bool(reference>>b));assert(std::memcmp(&a,&b,sizeof a)==0);
            if(column==15||i+1==values.size())assert(reader.End());
        }
        assert(i==values.size());
    }
    CheckpointTableReader lines("header\r\n\n1\t-0\r\n2 .5");double v;
    assert(lines.Next()&&lines.Number(v)&&v==1&&lines.Number(v)&&std::signbit(v)&&lines.End());
    assert(lines.Next()&&lines.Number(v)&&v==2&&lines.Number(v)&&v==.5&&lines.End()&&!lines.Next());
    CheckpointTableReader truncated("header\n1\n2 3\n");assert(truncated.Next()&&truncated.Number(v)&&!truncated.Number(v));
    assert(truncated.Next()&&truncated.Number(v)&&v==2);
    for(const auto& bad:std::vector<std::string>{"nan","inf","1e9999","1e","0x1p0","+",".","1oops","1,5",std::string("1\0 2",4)}){
        CheckpointTableReader reader("header\n"+bad+"\n");assert(reader.Next());assert(!reader.Number(v));
    }
    CheckpointTableReader extra("header\n1 2 trailing\n");assert(extra.Next()&&extra.Number(v)&&extra.Number(v)&&!extra.End());
    CheckpointTableReader empty("header only");assert(!empty.Next());
    std::cout<<"PASS checkpoint reader: 120000 decimal conversions bit-exact; row boundaries, signed zero, subnormals and malformed rows checked\n";
}
