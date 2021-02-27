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

int main(){       
    auto sites = Haagerup(N);
    Real U=100*N;
    // "ferromagnetic"
    Real K=-1;
    auto ampo = AutoMPO(sites);
        
    // set up excluded pairs
    for(int j = 1; j <= N-1; ++j){

        ampo += U,"n1",j,"n1",mod(j+1,N);
        ampo += U,"n1",j,"na",mod(j+1,N);
        ampo += U,"n1",j,"nb",mod(j+1,N);
//        ampo += U,"n1",j,"nr",mod(j+1,N);
        ampo += U,"n1",j,"nar",mod(j+1,N);
        ampo += U,"n1",j,"nbr",mod(j+1,N);

        ampo += U,"na",j,"n1",mod(j+1,N);
        ampo += U,"na",j,"na",mod(j+1,N);
        ampo += U,"na",j,"nb",mod(j+1,N);
        ampo += U,"na",j,"nr",mod(j+1,N);
//        ampo += U,"na",j,"nar",mod(j+1,N);
        ampo += U,"na",j,"nbr",mod(j+1,N);

        ampo += U,"nb",j,"n1",mod(j+1,N);
        ampo += U,"nb",j,"na",mod(j+1,N);
        ampo += U,"nb",j,"nb",mod(j+1,N);
        ampo += U,"nb",j,"nr",mod(j+1,N);
        ampo += U,"nb",j,"nar",mod(j+1,N);
//        ampo += U,"nb",j,"nbr",mod(j+1,N);

//        ampo += U,"nr",j,"n1",mod(j+1,N);
        ampo += U,"nr",j,"na",mod(j+1,N);
        ampo += U,"nr",j,"nb",mod(j+1,N);
//        ampo += U,"nr",j,"nr",mod(j+1,N);
//        ampo += U,"nr",j,"nar",mod(j+1,N);
//        ampo += U,"nr",j,"nbr",mod(j+1,N);

        ampo += U,"nar",j,"n1",mod(j+1,N);
//        ampo += U,"nar",j,"na",mod(j+1,N);
        ampo += U,"nar",j,"nb",mod(j+1,N);
//        ampo += U,"nar",j,"nr",mod(j+1,N);
//        ampo += U,"nar",j,"nar",mod(j+1,N);
//        ampo += U,"nar",j,"nbr",mod(j+1,N);

        ampo += U,"nbr",j,"n1",mod(j+1,N);
        ampo += U,"nbr",j,"na",mod(j+1,N);
//        ampo += U,"nbr",j,"nb",mod(j+1,N);
//        ampo += U,"nbr",j,"nr",mod(j+1,N);
//        ampo += U,"nbr",j,"nar",mod(j+1,N);
//        ampo += U,"nbr",j,"nbr",mod(j+1,N);
    }

//  projectors
    for(int j = 1; j <= N-2; ++j){
        ampo += K,"n1",j,"nr",mod(j+1,N),"n1",mod(j+2,N);
        ampo += K,"na",j,"nar",mod(j+1,N),"na",mod(j+2,N);
        ampo += K,"nb",j,"nbr",mod(j+1,N),"nb",mod(j+2,N);
        ampo += K,"nr",j,"FF1",mod(j+1,N),"nr",mod(j+2,N);
        ampo += K,"nar",j,"FFa",mod(j+1,N),"nar",mod(j+2,N);
        ampo += K,"nbr",j,"FFb",mod(j+1,N),"nbr",mod(j+2,N);
    }

    auto H = toMPO(ampo);
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    //
    auto sweeps = Sweeps(10);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-3,1E-3,1E-4,1E-4,1E-7,1E-8;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
    auto [en,psi] = dmrg(H,MPS(InitState(sites,"r")),sweeps,{"Quiet=",true});

    writeToFile("haagerup-open-10.psi",psi);
    int s=10;
    for(;;){
        MPS psi0(sites);
        readFromFile(format("haagerup-open-%d.psi",s),psi0);
        dumpEE(psi0);
        s+=10;
        auto sw=Sweeps(10);
        sw.maxdim()=200;
        sw.cutoff()=1E-10;
        sw.niter()=2;
        sw.noise()=1E-8;
        auto [en,psi1]=dmrg(H,psi0,sw,{"Quiet",true});
        writeToFile(format("haagerup-open-%d.psi",s),psi1);
    }
    return 0;
}

