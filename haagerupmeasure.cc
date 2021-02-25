#include "itensor/all.h"

using namespace itensor;

class HaagerupSite;
using Haagerup=BasicSiteSet<HaagerupSite>;
// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br
class HaagerupSite
{
    Index s;
public:
    HaagerupSite(Index const& I) : s(I) {}
    HaagerupSite(Args const& args = Args::global()){
        auto tags=TagSet("Site,Haagerup");
        s=Index(6,tags);
    }
    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state){
        if(state=="r"){
            return s(4);
        }
        throw ITError("State "+state+" not recognized");
    }
    
    ITensor
    proj(int i) const{
        auto sP=prime(s);
        auto Op=ITensor(dag(s),sP);
        Op.set(s(i),sP(i),1);
        return Op;
    }
    
    ITensor
    FF(int i) const{
        const double zetainv=(std::sqrt(13)-3)/2;
        const double sqrtzetainv=std::sqrt(zetainv);
        auto sP=prime(s);
        auto Op=ITensor(dag(s),sP);
        Op.set(s(i),sP(i),zetainv*zetainv);
        for(int j : {4,5,6}){
            Op.set(s(i),sP(j),zetainv*sqrtzetainv);
            Op.set(s(j),sP(i),zetainv*sqrtzetainv);
        }
        for(int j: {4,5,6}){
            for(int k: {4,5,6}){
                Op.set(s(j),sP(k),zetainv);
            }
        }
        return Op;
    }
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const
    {
        if(opname=="n1"){
            return proj(1);
        }else if(opname=="na"){
            return proj(2);
        }else if(opname=="nb"){
            return proj(3);
        }else if(opname=="nr"){
            return proj(4);
        }else if(opname=="nar"){
            return proj(5);
        }else if(opname=="nbr"){
            return proj(6);
        }else if(opname=="FF1"){
            return FF(1);
        }else if(opname=="FFa"){
            return FF(2);
        }else if(opname=="FFb"){
            return FF(3);
        }
        throw ITError("Operator name "+opname+" not recognized");
    }
};


inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
const int N = 50;

void dumpEE(MPS psi){
    for(auto b=1;b<N;b++){
        psi.position(b);

        //SVD this wavefunction to get the spectrum
        //of density-matrix eigenvalues
        auto l = leftLinkIndex(psi,b);
        auto s = siteIndex(psi,b);
        auto [U,S,V] = svd(psi(b),{l,s});
        auto u = commonIndex(U,S);

        //Apply von Neumann formula
        //to the squares of the singular values
        Real SvN = 0.;
        for(auto n : range1(dim(u)))
            {
            auto Sn = elt(S,n,n);
            auto p = sqr(Sn);
            if(p > 1E-12) SvN += -p*log(p);
            }
        printfln("{%d,  %.10f},",b,SvN);
    }
}

int main(int argc,char**argv){
    auto sites = Haagerup(N);
    MPS psi(sites);
    readFromFile(argv[1],psi);
    dumpEE(psi);
    return 0;
}

