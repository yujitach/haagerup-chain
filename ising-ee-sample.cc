#include "itensor/all.h"

using namespace itensor;

inline int mod(int x,int N){
    if(x>N)
        return x-N;
    return x;
}
int main(){       
    int N = 100;
    auto sites = SpinHalf(N,{"ConserveQNs=",false});
    Real h=1;
    //
    // Factors of 4 and 2 are to rescale
    // spin operators into Pauli matrices
    //
    auto ampo = AutoMPO(sites);
    for(int j = 1; j <= N; ++j){
        ampo += -4,"Sz",j,"Sz",mod(j+1,N);
    } 
    for(int j = 1; j <= N; ++j){
        ampo += -2*h,"Sx",j;
    }
    auto H = toMPO(ampo);
    //
    // Set the parameters controlling the accuracy of the DMRG
    // calculation for each DMRG sweep. 
    //
    auto sweeps = Sweeps(30);
    sweeps.maxdim() = 10,20,100,100,200;
    sweeps.cutoff() = 1E-10;
    sweeps.niter() = 2;
    sweeps.noise() = 1E-7,1E-8,0.0;
    //
    // Begin the DMRG calculation
    // for the ground state
    //
    auto [en,psi] = dmrg(H,MPS(InitState(sites,"Up")),sweeps,{"Quiet=",true});

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

