#include "itensor/all.h"

using namespace itensor;

class GoldenSite;
using Golden=BasicSiteSet<GoldenSite>;
// 1    1
// 2    t
class GoldenSite
{
    Index s;
public:
    GoldenSite(Index const& I) : s(I) {}
    GoldenSite(Args const& args = Args::global()){
        auto tags=TagSet("Site,Golden");
        s=Index(2,tags);
    }
    Index
    index() const { return s; }

    IndexVal
    state(std::string const& state){
        if(state=="t"){
            return s(2);
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
    FF() const{
        const double phiinv=(std::sqrt(5)-1)/2;
        const double sqrtphiinv=std::sqrt(phiinv);
        auto sP=prime(s);
        auto Op=ITensor(dag(s),sP);
        Op.set(s(1),sP(1),phiinv*phiinv);
        Op.set(s(1),sP(2),sqrtphiinv*phiinv);
        Op.set(s(2),sP(1),sqrtphiinv*phiinv);
        Op.set(s(2),sP(2),phiinv);
        return Op;
    }
    
    ITensor
    op(std::string const& opname,Args const& args = Args::global()) const
    {
        if(opname=="n1"){
            return proj(1);
        }else if(opname=="nt"){
            return proj(2);
        }else if(opname=="FF"){
            return FF();
        }
        throw ITError("Operator name "+opname+" not recognized");
    }
};


inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(){       
    int N = 50;
    auto sites = Golden(N);
    Real U=100*N;
    // "ferromagnetic"
    Real K=-1;
    auto ampo = AutoMPO(sites);
        
    // set up excluded pairs
    for(int j = 1; j <= N; ++j){

        ampo += U,"n1",j,"n1",mod(j+1,N);
//        ampo += U,"n1",j,"nt",mod(j+1,N);

//        ampo += U,"nt",j,"n1",mod(j+1,N);
//        ampo += U,"nt",j,"nt",mod(j+1,N);
    }

//  projectors
    for(int j = 1; j <= N; ++j){
        ampo += K,"n1",j,"nt",mod(j+1,N),"n1",mod(j+2,N);
        ampo += K,"nt",j,"FF",mod(j+1,N),"nt",mod(j+2,N);
    }

    auto H = toMPO(ampo);
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    //
    auto sweeps = Sweeps(10000);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
    auto [en,psi] = dmrg(H,InitState(sites,"t"),sweeps,{"Quiet",true});

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
    return 0;
}

