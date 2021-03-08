#include "itensor/all.h"

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}

class HaagerupSite;
using Haagerup=BasicSiteSet<HaagerupSite>;
// 1    1
// 2    a
// 3    b=a^2
// 4    r=\rho
// 5    ar
// 6    br


using  std::sqrt;
class HaagerupSite
{
    const double zetainv=(sqrt(13)-3)/2;
    const double sqrtzetainv=sqrt(zetainv);
    const double x=(2-sqrt(13))/3;
    const double z=(1+sqrt(13))/6;
    const double y1=(5-sqrt(13)-sqrt(6*(1+sqrt(13))))/12;
    const double y2=(5-sqrt(13)+sqrt(6*(1+sqrt(13))))/12;
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
    Frr(int i) const{
        auto sP=prime(s);
        auto OpL=ITensor(dag(s));
        auto OpR=ITensor(sP);
        OpL.set(s(i),sqrtzetainv);
        OpL.set(s(3+i),x);
        OpL.set(s(3+mod(i+1,3)),y1);
        OpL.set(s(3+mod(i+2,3)),y2);
        OpR.set(sP(i),sqrtzetainv);
        OpR.set(sP(3+i),x);
        OpR.set(sP(3+mod(i+1,3)),y1);
        OpR.set(sP(3+mod(i+2,3)),y2);
        return OpL*OpR;
    }

    ITensor
    Frar(int i) const{
        auto sP=prime(s);
        auto OpL=ITensor(dag(s));
        auto OpR=ITensor(sP);
        OpL.set(s(3+i),y1);
        OpL.set(s(3+mod(i+1,3)),y2);
        OpL.set(s(3+mod(i+2,3)),z);
        OpR.set(sP(3+i),y1);
        OpR.set(sP(3+mod(i+1,3)),y2);
        OpR.set(sP(3+mod(i+2,3)),z);
        return OpL*OpR;
    }

    ITensor
    Frbr(int i) const{
        auto sP=prime(s);
        auto OpL=ITensor(dag(s));
        auto OpR=ITensor(sP);
        OpL.set(s(3+i),y2);
        OpL.set(s(3+mod(i+1,3)),z);
        OpL.set(s(3+mod(i+2,3)),y1);
        OpR.set(sP(3+i),y2);
        OpR.set(sP(3+mod(i+1,3)),z);
        OpR.set(sP(3+mod(i+2,3)),y1);
        return OpL*OpR;
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
        }else if(opname=="Frr"){
            return Frr(1);
        }else if(opname=="Farar"){
            return Frr(2);
        }else if(opname=="Fbrbr"){
            return Frr(3);
        }else if(opname=="Frar"){
            return Frar(1);
        }else if(opname=="Farbr"){
            return Frar(2);
        }else if(opname=="Fbrr"){
            return Frar(3);
        }else if(opname=="Frbr"){
            return Frbr(1);
        }else if(opname=="Farr"){
            return Frbr(2);
        }else if(opname=="Fbrar"){
            return Frbr(3);
        }
        throw ITError("Operator name "+opname+" not recognized");
    }
};


const int N = 48;

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
    Real U=2;
    // "ferromagnetic"
    Real K=0;
    Real J=-1;
    auto ampo = AutoMPO(sites);
        
    // set up excluded pairs
    for(int j = 1; j <= N; ++j){

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
    for(int j = 1; j <= N; ++j){
        ampo += K,"n1",j,"nr",mod(j+1,N),"n1",mod(j+2,N);
        ampo += K,"na",j,"nar",mod(j+1,N),"na",mod(j+2,N);
        ampo += K,"nb",j,"nbr",mod(j+1,N),"nb",mod(j+2,N);
        ampo += K,"nr",j,"FF1",mod(j+1,N),"nr",mod(j+2,N);
        ampo += K,"nar",j,"FFa",mod(j+1,N),"nar",mod(j+2,N);
        ampo += K,"nbr",j,"FFb",mod(j+1,N),"nbr",mod(j+2,N);

        ampo += J,"n1",j,"nr",mod(j+1,N),"nr",mod(j+2,N);
        ampo += J,"nr",j,"nr",mod(j+1,N),"n1",mod(j+2,N);
        ampo += J,"na",j,"nar",mod(j+1,N),"nar",mod(j+2,N);
        ampo += J,"nar",j,"nar",mod(j+1,N),"na",mod(j+2,N);
        ampo += J,"nb",j,"nbr",mod(j+1,N),"nbr",mod(j+2,N);
        ampo += J,"nbr",j,"nbr",mod(j+1,N),"nb",mod(j+2,N);
        ampo += J,"nr",j,"Frr",mod(j+1,N),"nr",mod(j+2,N);
        ampo += J,"nr",j,"Frar",mod(j+1,N),"nar",mod(j+2,N);
        ampo += J,"nr",j,"Frbr",mod(j+1,N),"nbr",mod(j+2,N);
        ampo += J,"nar",j,"Farr",mod(j+1,N),"nr",mod(j+2,N);
        ampo += J,"nar",j,"Farar",mod(j+1,N),"nar",mod(j+2,N);
        ampo += J,"nar",j,"Farbr",mod(j+1,N),"nbr",mod(j+2,N);
        ampo += J,"nbr",j,"Fbrr",mod(j+1,N),"nr",mod(j+2,N);
        ampo += J,"nbr",j,"Fbrar",mod(j+1,N),"nar",mod(j+2,N);
        ampo += J,"nbr",j,"Fbrbr",mod(j+1,N),"nbr",mod(j+2,N);

        

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

    writeToFile("haagerup-R-closed-10.psi",psi);
    int s=10;
    for(;;){
        MPS psi0(sites);
        readFromFile(format("haagerup-R-closed-%d.psi",s),psi0);
        dumpEE(psi0);
        s+=10;
        auto sw=Sweeps(10);
        sw.maxdim()=200;
        sw.cutoff()=1E-10;
        sw.niter()=2;
        sw.noise()=1E-8;
        auto [en,psi1]=dmrg(H,psi0,sw,{"Quiet",true});
        writeToFile(format("haagerup-R-closed-%d.psi",s),psi1);
    }
    return 0;
}

