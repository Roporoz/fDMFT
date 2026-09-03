#include "headers.new.h"
#include "parameters.cpp"

ofstream debug("debug.dat"), ou("main.dat"), ou_Sigma("Sigma.dat"); 



#include "solver.cpp"


int main()
{

    
    DMFT_AF_solver solver(0.); solver.init(.003);
    cout<<solver.raw_loops(1)<<"  "; cout<<(1./beta)*solver.Sf.sz()/(solver.Sf.hz_fine-solver.Sf.sz()*solver.lambda0)<<"  s="<<solver.Sf.sz()<<"  at h="<<solver.Sf.hz_fine-solver.Sf.sz()*solver.lambda0<<"\n";
    
    
//    solver.Sf.init(solver.Sr);  solver.Sf.set_g(); cout<<solver.Sf.gz[0]<<"  "<<solver.Sf.sz()<<"\n";
//    solver.Sf.hz_fine=.02; solver.Sf.set_g(); cout<<solver.Sf.gz[0]<<"  "<<solver.Sf.sz()<<"\n";
      
    cout<<solver.fine_loops(1)<<"  "; cout<<(1./beta)*solver.Sf.sz()/(solver.Sf.hz_fine-solver.Sf.sz()*solver.lambda0)<<"\n";
    
//    bath_f2r(solver.Sf, solver.Sr); bath_r2f(solver.Sr, solver.Sf); cout<<(1./beta)*solver.Sf.sz()/solver.Sf.hz_fine<<"\n";
    
    
    /*
 
    solver_raw Sr; solver_fine Sf;
    set_Iw(); for1(w, Nw) {Sf.Delta_u[w]=.1/(Iw[w]-.5)+.1/(Iw[w]+.5)+.3/(Iw[w]);  Sf.Delta_d[w]=.2/(Iw[w]-.5)+.4/(Iw[w]+.45)+.2/(Iw[w]);} Sf.mu_fine=0.1; Sf.hz_fine=-0.2;

    bath_f2r(Sf,Sr); for1(w, 2)    cout<<Sr.gu[w]<<"  "<<Sr.gd[w]<<"     "; cout<<"\n"<<flush;
    bath_r2f(Sr, Sf);for1(w, 2)    cout<<Sf.gu[w]<<"  "<<Sf.gd[w]<<"     "; cout<<"\n"<<flush;
    
    cout<<"\n";
    
    //Sf.hz_fine+=.05; 
    for1(w, Nw) Sf.Delta_u[w]+=.01/(Iw[w]-.45);
    
    Sf.set_g(); for1(w, 2)    cout<<Sf.gu[w]<<"  "<<Sf.gd[w]<<"     "; cout<<"\n"<<flush;
    bath_f2r(Sf,Sr); for1(w, 2)    cout<<Sr.gu[w]-Sf.gu[w]<<"  "<<Sr.gd[w]-Sf.gd[w]<<"     "; cout<<"\n\n"<<flush;

    
    //Sf.hz_fine+=.05;
    for1(w, Nw) Sf.Delta_u[w]+=.01/(Iw[w]-.4); 
    
    Sf.set_g(); for1(w, 2)    cout<<Sf.gu[w]<<"  "<<Sf.gd[w]<<"     "; cout<<"\n"<<flush;
    bath_f2r(Sf,Sr); for1(w, 2)    cout<<Sr.gu[w]-Sf.gu[w]<<"  "<<Sr.gd[w]-Sf.gd[w]<<"     "; cout<<"\n\n"<<flush;

    */
    return 0;
}
