complex * Iw=new complex [Nw]; double * Iw2= new double [Nw];   
bool set_Iw(){for1(w, Nw) {Iw[w]=I*(2.*w+1)*Pi/beta; Iw2[w]=-sqr((2.*w+1)*Pi/beta); }  return true;} bool set_Iw_aux=set_Iw();

ofstream ED_stream("ED.log");

class solver_raw
{
    int NFock;
    double ** H, ** Psi, ** cuPsi,** cdPsi,  *E, ** cu, **cd; long double Z;
    double adjust_bath_levels(complex * Delta_w, double * t2_bath, double * e_bath);
    double adjust_bath_couplings(complex * Delta_w, double * t2_bath,  double * e_bath);
    int n(int j, int l) {return (j>>l)%2;}
    int * nud;
    
public:
    double tolerance;
    solver_raw() ;

    double * t2u, * t2d, *eu, *ed, mu_loc, h_loc; 

    complex * gu, * gd;
    complex Gu(int w) {if (w>=0) return gu[w]; else return conj(gu[-w-1]);}
    complex Gd(int w) {if (w>=0) return gd[w]; else return conj(gd[-w-1]);}
        
    complex * sigmau, * sigmad;
    
    double s() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(gu[w]-gd[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }
    double n() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(gu[w]+gd[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }

    complex * gudh, * gddh, * gxdh;

    
    complex Gudh(int w) {if (w>=0) return gudh[w]; else return conj(gudh[-w-1]);}
    complex Gddh(int w) {if (w>=0) return gddh[w]; else return conj(gddh[-w-1]);}
    complex Gxdh(int w) {if (w>=0) return gxdh[w]; else return conj(gxdh[-w-1]);}

    
    double dszdh() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(gudh[w]-gddh[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }
    double dsxdh() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(gxdh[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }
    
    
 
    double init(complex * Delta_up, complex *Delta_down, double h_loc=0, double mu_loc=0);
    void store(complex * Du, complex * Dd, complex * Gu, complex * Gd) {for1(w, Nw) {Du[w]=0; Dd[w]=0; for1(l, Nbath_raw) {Du[w]+=t2u[l]/(Iw[w]-eu[l]); Dd[w]+=t2d[l]/(Iw[w]-ed[l]);}  Gu[w]=gu[w];Gd[w]=gd[w];}}
    
  
    
    double *** g2tau_uu, *** g2tau_du, *** g2tau_dd; //c c+ c c+
    double *gtau_u, *gtau_d;//-c c+

    complex *** gamma_du, *** gamma_ud, *** gamma_uu, *** gamma_dd;
    
    
    void set_g2();
    
    
    ~solver_raw();
};


solver_raw::solver_raw()
{
        NFock=1<<(2*Nbath_raw+2);
        H=new_double2(NFock, NFock); Psi=new_double2(NFock, NFock); cuPsi=new_double2(NFock, NFock); cdPsi=new_double2(NFock, NFock); cu=new_double2(NFock, NFock); cd=new_double2(NFock, NFock); E=new double[NFock]; nud=new int[NFock];
        
        
        t2u=new double[Nbath_raw]; t2d=new double[Nbath_raw]; eu=new double[Nbath_raw]; ed=new double[Nbath_raw]; gu=new complex [Nw]; gd=new complex [Nw]; for1(w, Nw) {gu[w]=1./Iw[w]; gd[w]=1./Iw[w];} sigmau=new complex [Nw]; sigmad=new complex [Nw]; for1(w, Nw) {sigmau[w]=0; sigmad[w]=0;} for1(l, Nbath_raw) {eu[l]=0.1*l; ed[l]=0.1*l;};
        gudh=new complex [Nw]; gddh=new complex [Nw]; gxdh=new complex [Nw];
        g2tau_uu=new_double3(Ntau4+1, Ntau4+1, Ntau4+1); g2tau_du=new_double3(Ntau4+1, Ntau4+1, Ntau4+1); g2tau_dd=new_double3(Ntau4+1, Ntau4+1, Ntau4+1);
        gtau_u=new double[Ntau4+1]; gtau_d=new double[Ntau4+1]; 
        gamma_du=new_complex3(2*Nw4, 2*Nw4, 2*Nw4); gamma_ud=new_complex3(2*Nw4, 2*Nw4, 2*Nw4); gamma_uu=new_complex3(2*Nw4, 2*Nw4, 2*Nw4); gamma_dd=new_complex3(2*Nw4, 2*Nw4, 2*Nw4);
        
}


solver_raw::~solver_raw()
{//cout<<"?\n"<<flush;return;
        for1(j, NFock) {delete [] H[j]; delete [] Psi [j]; delete [] cuPsi[j]; delete [] cdPsi[j]; delete [] cu[j]; delete [] cd[j];} delete [] H; delete [] Psi; delete [] cuPsi; delete [] cdPsi; delete [] cu; delete [] cd;  delete [] E; delete [] nud; delete [] t2u;delete [] t2d;delete [] eu;delete [] ed; delete [] gu; delete [] gd; delete [] gudh; delete [] gddh; delete [] gxdh;
        delete_double3(g2tau_uu, Ntau4+1, Ntau4+1, Ntau4+1); delete_double3(g2tau_du, Ntau4+1, Ntau4+1, Ntau4+1); delete_double3(g2tau_dd, Ntau4+1, Ntau4+1, Ntau4+1); delete [] gtau_u; delete [] gtau_d;
        delete_complex3(gamma_du, 2*Nw4, 2*Nw4, 2*Nw4); delete_complex3(gamma_ud, 2*Nw4, 2*Nw4, 2*Nw4); delete_complex3(gamma_uu, 2*Nw4, 2*Nw4, 2*Nw4); delete_complex3(gamma_dd, 2*Nw4, 2*Nw4, 2*Nw4);
}

double solver_raw::adjust_bath_levels(complex * Delta_w,   double * t2_bath, double * e_bath)
{ 
    double s, s_min=1e6; static double * e_bath_min=new double [Nbath_raw];
    static int  Ntrials=10000;  int i_min=-1;//10000
    static double de=1, *e_prev=new double [Nbath_raw]; bool t2positive;
    for1(l, Nbath_raw) e_prev[l]=e_bath[l];
 
    if (rnd()<.5) for1(j, Nbath_raw) e_bath[j]=de*real(rnd_gauss2()); else for1(j, Nbath_raw) e_bath[j]=e_prev[j]+.3*de*real(rnd_gauss2());
    s_min=adjust_bath_couplings(Delta_w, t2_bath, e_bath); for1(l, Nbath_raw) e_bath_min[l]=e_bath[l];
    
     
    for1 (i, Ntrials)
    {  
        
        if (rnd()<.5) for1(j, Nbath_raw) e_bath[j]=de*real(rnd_gauss2()); else for1(j, Nbath_raw) e_bath[j]=e_prev[j]+.3*de*real(rnd_gauss2());
         
        
        s=adjust_bath_couplings(Delta_w, t2_bath, e_bath); 
        if (!(s<0) && s<.9999*s_min) {s_min=s; for1(j, Nbath_raw) e_bath_min[j]=e_bath[j]; i_min=i;t2positive=true; for1(l, Nbath_raw) t2positive=(t2positive & t2_bath[l]>0);}        
        ED_stream<<s<<"  "<<s_min<<"\n"; 
        for1 (ii, 50)
        {
            static double de=.001, * s1=new double[Nbath_raw], ** s2=new_double2(Nbath_raw, Nbath_raw); 
            for1(j, Nbath_raw) 
            {
                e_bath[j]+=de;  
                double sp=adjust_bath_couplings(Delta_w, t2_bath, e_bath); e_bath[j]-=2.*de; 
                double sm=adjust_bath_couplings(Delta_w, t2_bath, e_bath); e_bath[j]+=de; s1[j]=(sp-sm)/(2.*de); s2[j][j]=(sp+sm-2.*s)/(de*de);
            }
            for1(j1, Nbath_raw)    for(int j2=j1+1; j2<Nbath_raw; j2++)
            {
                e_bath[j1]+=de; e_bath[j2]+=de; double sxy=adjust_bath_couplings(Delta_w, t2_bath, e_bath); e_bath[j1]-=de; e_bath[j2]-=de;
                //s2[j1][j2]=(sxy-s-(s1[j1]+s1[j2])*de-.5*(s2[j1][j1]+s2[j2][j2])*de*de )/(de*de);
                
                e_bath[j1]-=de; e_bath[j2]-=de; double sxym=adjust_bath_couplings(Delta_w, t2_bath, e_bath); e_bath[j1]+=de; e_bath[j2]+=de;
                s2[j1][j2]=(sxy+sxym-2.*s-s2[j1][j1]*de*de - s2[j2][j2]*de*de)/(2.*de*de);
                
                s2[j2][j1]=s2[j1][j2];
            }
         
            
        
            if (div_left(s2,s1,Nbath_raw))  for1(j, Nbath_raw) e_bath[j]-=s1[j];  
            double s0=s;       s=adjust_bath_couplings(Delta_w, t2_bath, e_bath);  
            if (s>=s0) break;    
            
            if (!(s<0) && s<.9999*s_min) {s_min=s; for1(j, Nbath_raw) e_bath_min[j]=e_bath[j]; i_min=i; t2positive=true; for1(l, Nbath_raw) t2positive=(t2positive & t2_bath[l]>0);} //cout<<i<<"  "<<s<<"  "; for1(j, Nbath) cout<<e_bath[j]<<"   "; cout <<"\n"; }
            ED_stream<<ii<<"   "<<s<<"  "<<s_min<<"\n"<<flush; 
        }
        
 
        
        if (s_min<1e-10 && t2positive) break;  
        if (!(s<0) && s<1e-10 && t2positive) {s_min=s; for1(j, Nbath_raw) e_bath_min[j]=e_bath[j];i_min=i; break;} 
        
        
//        if (s_min<tolerance) {cout<<i<<"\n"; break;}
    } 
        
    ED_stream<<i_min<<"/"<<Ntrials<<"   "<<de<<"    "<<s_min<<"\n";
    
    for1(j, Nbath_raw) e_bath[j]=e_bath_min[j];  
    adjust_bath_couplings(Delta_w, t2_bath, e_bath);
    if (!t2positive) {cout<<"!!! adjust_bath_levels led to unphysical bath !!!\n"<<flush; for1(l, Nbath_raw) cout<<"e="<<e_bath[l]<<"  t2="<<t2_bath[l]<<"\n";}
    
    double de_max=0; for1(j, Nbath_raw) if (de_max<abs(e_bath_min[j])) de_max=abs(e_bath_min[j]); 
    
    if (i_min>=0) Ntrials=.9*Ntrials+i_min; de=.9*de+.1*de_max; if (Ntrials<300) Ntrials=300; 
    for1(j, Nbath_raw) e_bath[j]=e_bath_min[j];        

    return adjust_bath_couplings(Delta_w, t2_bath, e_bath);
}


double solver_raw::adjust_bath_couplings(complex * Delta_w, double  * t2_bath, double * e_bath)
{
    double a, a2, a3; 
    static double ** A=new_double2(Nbath_raw, Nbath_raw);
    for1(i,Nbath_raw) for(int j=i; j<Nbath_raw; j++) {a=0; for1(w, NwED) {a+=(e_bath[i]/(sqr(e_bath[i])-Iw2[w])+e_bath[j]/(sqr(e_bath[j])-Iw2[w]))/(e_bath[i]+e_bath[j]);  //-real(1./(-Iw[w]-e_bath[i])+1./(Iw[w]-e_bath[j]))/(e_bath[i]+e_bath[j]); //real(1./((-Iw[w]-e_bath[i])*(Iw[w]-e_bath[j]))); 
        if (w+1==NwED/2) a2=a; if (w+1==3*NwED/4) a3=a;} A[i][j]=8.*a-9.*a3+2.*a2;  A[j][i]=A[i][j];}
    
    for1(i, Nbath_raw) {a=0; for1(w, NwED) {a+=real(conj(Delta_w[w])/(Iw[w]-e_bath[i])); if (w+1==NwED/2) a2=a; if (w+1==3*NwED/4) a3=a;} t2_bath[i]=8.*a-9.*a3+2.*a2; }
        
    if (!div_left(A,t2_bath,Nbath_raw)) return 100; 
    
    a=0;
    for1(w, NwED)
    {
        complex r=Delta_w[w]; for1(i, Nbath_raw) r-=t2_bath[i]/(Iw[w]-e_bath[i]);  
        a+=norm2(r);
        if (w+1==NwED/2) a2=a; if (w+1==3*NwED/4) a3=a;
    }
    return 2.*(8.*a-9.*a3+2.*a2)/beta;
}


double solver_raw::init(complex * Delta_up, complex *Delta_down, double mu_loc_,  double h_loc_)
{
 
    set_Iw(); 
    tolerance=adjust_bath_levels(Delta_up, t2u, eu)+adjust_bath_levels(Delta_down, t2d, ed); 
    //for1(l, Nbath_raw) cout<<t2u[l]<<"  "<<t2d[l]<<"   ";
    
     
    mu_loc=mu_loc_; h_loc=h_loc_;
    
    for2(i,j, NFock) {H[i][j]=0; cuPsi[i][j]=0; cdPsi[i][j]=0; cu[i][j]=0; cd[i][j]=0; }
    for1(j, NFock) {nud[j]=0; for1(l, Nbath_raw+1) nud[j]+=n(j,l)+(Nbath_raw+2)*n(j,l+Nbath_raw+1); } 
    
    
    for1(j, NFock)
    {
        H[j][j]+=-(mu_loc)*(n(j,0)+n(j, Nbath_raw+1))-(h_loc)*((n(j,0)-n(j, Nbath_raw+1)))+U*((n(j,0)-.5)*(n(j,Nbath_raw+1)-.5));
        for1(l, Nbath_raw) H[j][j]+=n(j,l+1)*eu[l]+n(j,Nbath_raw+l+2)*ed[l];
        for1(l, Nbath_raw)
        {
            int nu=0, nd=0; for1(l2,l)   {nu+=n(j,l2+1); nd+=n(j,Nbath_raw+l2+2);}
            if(n(j,0)==0  && n(j, l+1)==1)                 {int j1=j+1-(1<<(l+1));                        H[j][j1]+=sqrt(t2u[l])*(1-2*(nu%2));H[j1][j]=H[j][j1];}
            if(n(j,Nbath_raw+1)==0  && n(j, l+Nbath_raw+2)==1)   {int j1=j+(1<<(Nbath_raw+1))-(1<<(l+Nbath_raw+2));   H[j][j1]+=sqrt(t2d[l])*(1-2*(nd%2));H[j1][j]=H[j][j1];}
        }
    }
    
  
    EigenBlock(H, Psi, E, NFock, nud, sqr(Nbath_raw+2));
    
    Z=0; for1(j, NFock) Z+=expl(-beta*E[j]);  

    
    long double lnZ=logl(Z); for1(l, Nbath_raw) lnZ-=log(1+exp(-beta*eu[l]))+log(1+exp(-beta*ed[l]));

 
    for1(j, NFock) 
    {
        if (n(j,0)==1) for1(j1, NFock)  cuPsi[j1][j-1]=Psi[j1][j];
        if (n(j, Nbath_raw+1)==1) {double f=1-2*((nud[j]%(Nbath_raw+2))%2); 
            for1(j1, NFock) cdPsi[j1][j-(1<<(Nbath_raw+1))]=f*Psi[j1][j];  }
    }
        
    for1 (w, Nw) {gu[w]=0; gd[w]=0;}

    for2(j1, j2, NFock) 
    {
        if (nud[j1]==nud[j2]-1)     
        {
            double d=0; for1(j3, NFock) d+=Psi[j1][j3]*cuPsi[j2][j3]; cu[j1][j2]=d;
            double x=d*d*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); 
        
            for1(w, Nw) gu[w]+=x/(Iw[w]+E[j1]-E[j2]);
        }
        if (nud[j1]==nud[j2]-(Nbath_raw+2) )     
        {
            double d=0; for1(j3, NFock) d+=Psi[j1][j3]*cdPsi[j2][j3]; cd[j1][j2]=d;
            double x=d*d*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); 
        
            for1(w, Nw) gd[w]+=x/(Iw[w]+E[j1]-E[j2]);
        }

    }

    for1(w, Nw) {sigmau[w]=Iw[w]-1./gu[w]+h_loc+mu_loc; sigmad[w]=Iw[w]-1./gd[w]-h_loc+mu_loc; for1(l, Nbath_raw) {sigmau[w]-=t2u[l]/(Iw[w]-eu[l]); sigmad[w]-=t2d[l]/(Iw[w]-ed[l]);} }
    
    lnZ-=beta*(mu_loc);
    
    ED_stream<<lnZ<<"\n"<<flush;
 
    
    return lnZ;
}

 


void solver_raw::set_g2()
{
    double ** ExpEtau=new_double2(NFock, Ntau4+1);  
    for1 (j, NFock) for1(tau, Ntau4+1) ExpEtau[j][tau]=exp(-E[j]*tau*beta/Ntau4);
    
    
    double **** cc=new double *** [10];
    //cu cu : 0 
    //cd cu : 1
    //cu+ cu : 2
    //cd+ cu : 3
    //cu cd : 4
    //cd cd : 5 
    //cu+ cd : conj(3)
    //cd+ cd : 6
    //cu cu+ : 7
    //cd cu+ : 8
    //cu+ cu+ : conj(0)
    //cd+ cu+ : conj(4)
    //cu cd+ : conj(8)
    //cd cd+ : 9
    //cu+ cd+: conj (1)
    //cd+ cd+: conj (5)

    //tau1>tau_>tau2             cd  cd+ cu  cu+   9-7,  7-7,  9-9     perm +    
    //tau1>tau2>tau_             cd  cu  cd+ cu+   4-c4, 0-c0, 5-c5    perm -
    //tau2>tau_>tau1             cu  cd+ cd  cu+   c8-8, 7-7,  9-9     perm -
    //tau2>tau1>tau_             cu  cd  cd+ cu+   4-c4, 0-c0, 5-c5    perm +
    //tau_>tau1>tau2             cd+ cd  cu  cu+   6-7,  2-7,  6-9     perm -
    //tau_>tau2>tau1             cd+ cu  cd  cu+   3-8,  2-7,  6-9     perm +
        
    
    for1(l, 10)   {cc[l]=new double ** [NFock]; for1(i, NFock) cc[l][i]= new double * [NFock];}

    
    for1(i, NFock) 
    {   
            
        for1(j, NFock) 
        {
                if (nud[i]==nud[j]-2)            
                    {cc[0][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[0][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-1)          for1(tau, Ntau4+1) cc[0][i][j][tau]+=cu[i][k]*ExpEtau[k][tau]*cu[k][j];} else cc[0][i][j]=NULL;
                if (nud[i]==nud[j]-1-(Nbath_raw+2))    
                    {cc[1][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[1][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-1)          for1(tau, Ntau4+1) cc[1][i][j][tau]+=cd[i][k]*ExpEtau[k][tau]*cu[k][j];} else cc[1][i][j]=NULL;
                if (nud[i]==nud[j])                 
                    {cc[2][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[2][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-1)          for1(tau, Ntau4+1) cc[2][i][j][tau]+=cu[k][i]*ExpEtau[k][tau]*cu[k][j];} else cc[2][i][j]=NULL;
                if (nud[i]==nud[j]-1+(Nbath_raw+2))    
                    {cc[3][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[3][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-1)          for1(tau, Ntau4+1) cc[3][i][j][tau]+=cd[k][i]*ExpEtau[k][tau]*cu[k][j];} else cc[3][i][j]=NULL;
                if (nud[i]==nud[j]-1-(Nbath_raw+2))    
                    {cc[4][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[4][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-(Nbath_raw+2)) for1(tau, Ntau4+1) cc[4][i][j][tau]+=cu[i][k]*ExpEtau[k][tau]*cd[k][j];} else cc[4][i][j]=NULL;
                if (nud[i]==nud[j]-2*(Nbath_raw+2))    
                    {cc[5][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[5][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-(Nbath_raw+2)) for1(tau, Ntau4+1) cc[5][i][j][tau]+=cd[i][k]*ExpEtau[k][tau]*cd[k][j];}  else cc[5][i][j]=NULL;
                if (nud[i]==nud[j])                 
                    {cc[6][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[6][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]-(Nbath_raw+2)) for1(tau, Ntau4+1) cc[6][i][j][tau]+=cd[k][i]*ExpEtau[k][tau]*cd[k][j];} else cc[6][i][j]=NULL;
                if (nud[i]==nud[j])                 
                    {cc[7][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[7][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]+1)          for1(tau, Ntau4+1) cc[7][i][j][tau]+=cu[i][k]*ExpEtau[k][tau]*cu[j][k];}  else cc[7][i][j]=NULL;  
                if (nud[i]==nud[j]+1-(Nbath_raw+2))    
                    {cc[8][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[8][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]+1)          for1(tau, Ntau4+1) cc[8][i][j][tau]+=cd[i][k]*ExpEtau[k][tau]*cu[j][k];} else cc[8][i][j]=NULL;
                if (nud[i]==nud[j])                 
                    {cc[9][i][j]=new double [Ntau4+1]; for1(tau, Ntau4+1) cc[9][i][j][tau]=0; for1(k,NFock) if(nud[k]==nud[j]+(Nbath_raw+2)) for1(tau, Ntau4+1) cc[9][i][j][tau]+=cd[i][k]*ExpEtau[k][tau]*cd[j][k];} else cc[9][i][j]=NULL;
                
                
                
                
        }
    }
        
    //cout<<"calculating G2   "<<flush;
    
     
    double Z=0; for1(j, NFock) Z+=ExpEtau[j][Ntau4];
    for1(tau, Ntau4+1) {gtau_u[tau]=0; for1(j, NFock) gtau_u[tau]-=cc[7][j][j][tau]*ExpEtau[j][Ntau4-tau]; gtau_u[tau]/=Z;}
    for1(tau, Ntau4+1) {gtau_d[tau]=0; for1(j, NFock) gtau_d[tau]-=cc[9][j][j][tau]*ExpEtau[j][Ntau4-tau]; gtau_d[tau]/=Z;}    
        
    for3(tau1, tau_, tau2, Ntau4+1) {g2tau_du[tau1][tau_][tau2]=0; g2tau_uu[tau1][tau_][tau2]=0; g2tau_dd[tau1][tau_][tau2]=0;}

    
    
    
    //tau1>tau_>tau2             cd  cd+ cu  cu+   9-7,  7-7,  9-9     perm +    
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j])
        for1(tau1, Ntau4+1) for1(tau_, tau1+1) for1(tau2, tau_+1) 
            g2tau_du[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau1]*cc[9][j][i][tau1-tau_]*ExpEtau[i][tau_-tau2]*cc[7][i][j][tau2];     
        if (nud[i]==nud[j])
        for1(tau1, Ntau4+1) for1(tau_, tau1+1) for1(tau2, tau_+1) 
            g2tau_uu[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau1]*cc[7][j][i][tau1-tau_]*ExpEtau[i][tau_-tau2]*cc[7][i][j][tau2];   
        if (nud[i]==nud[j])
        for1(tau1, Ntau4+1) for1(tau_, tau1+1) for1(tau2, tau_+1) 
            g2tau_dd[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau1]*cc[9][j][i][tau1-tau_]*ExpEtau[i][tau_-tau2]*cc[9][i][j][tau2];   

        
    }

    
     
    //tau1>tau2>tau_             cd  cu  cd+ cu+   4-c4, 0-c0, 5-c5    perm -
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j]+1+(Nbath_raw+2)) 
        for1(tau1, Ntau4+1) for1(tau2, tau1+1) for1(tau_, tau2+1) 
            g2tau_du[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau1]*cc[1][j][i][tau1-tau2]*ExpEtau[i][tau2-tau_]*cc[4][j][i][tau_];     
        if (nud[i]==nud[j]+2) 
        for1(tau1, Ntau4+1) for1(tau2, tau1+1) for1(tau_, tau2+1) 
            g2tau_uu[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau1]*cc[0][j][i][tau1-tau2]*ExpEtau[i][tau2-tau_]*cc[0][j][i][tau_];   
        if (nud[i]==nud[j]+2*(Nbath_raw+2)) 
        for1(tau1, Ntau4+1) for1(tau2, tau1+1) for1(tau_, tau2+1) 
            g2tau_dd[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau1]*cc[5][j][i][tau1-tau2]*ExpEtau[i][tau2-tau_]*cc[5][j][i][tau_];   
    }
    
    
    //tau2>tau_>tau1             cu  cd+ cd  cu+   c8-8, 7-7, 9-9     perm -
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j]+1-(Nbath_raw+2)) 
        for1(tau2, Ntau4+1) for1(tau_, tau2+1) for1(tau1, tau_+1) 
            g2tau_du[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau2]*cc[8][i][j][tau2-tau_]*ExpEtau[i][tau_-tau1]*cc[8][i][j][tau1];     
        if (nud[i]==nud[j]) 
        for1(tau2, Ntau4+1) for1(tau_, tau2+1) for1(tau1, tau_+1) 
            g2tau_uu[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau2]*cc[7][j][i][tau2-tau_]*ExpEtau[i][tau_-tau1]*cc[7][i][j][tau1];   
        if (nud[i]==nud[j])
        for1(tau2, Ntau4+1) for1(tau_, tau2+1) for1(tau1, tau_+1) 
            g2tau_dd[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau2]*cc[9][j][i][tau2-tau_]*ExpEtau[i][tau_-tau1]*cc[9][i][j][tau1];   
    }
    
    
    //tau2>tau1>tau_             cu  cd  cd+ cu+   4-c4, 0-c0, 5-c5  perm +
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j]+1+(Nbath_raw+2)) 
        for1(tau2, Ntau4+1) for1(tau1, tau2+1) for1(tau_, tau1+1) 
            g2tau_du[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau2]*cc[4][j][i][tau2-tau1]*ExpEtau[i][tau1-tau_]*cc[4][j][i][tau_];     
        if (nud[i]==nud[j]+2) 
        for1(tau2, Ntau4+1) for1(tau1, tau2+1) for1(tau_, tau1+1) 
            g2tau_uu[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau2]*cc[0][j][i][tau2-tau1]*ExpEtau[i][tau1-tau_]*cc[0][j][i][tau_];   
        if (nud[i]==nud[j]+2*(Nbath_raw+2)) 
        for1(tau2, Ntau4+1) for1(tau1, tau2+1) for1(tau_, tau1+1) 
            g2tau_dd[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau2]*cc[5][j][i][tau2-tau1]*ExpEtau[i][tau1-tau_]*cc[5][j][i][tau_];   
    }
    
    
    //tau_>tau1>tau2             cd+ cd  cu  cu+   6-7,  2-7, 6-9   perm -
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j]) 
        for1(tau_, Ntau4+1) for1(tau1, tau_+1) for1(tau2, tau1+1) 
            g2tau_du[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau_]*cc[6][j][i][tau_-tau1]*ExpEtau[i][tau1-tau2]*cc[7][i][j][tau2];     
        if (nud[i]==nud[j]) 
        for1(tau_, Ntau4+1) for1(tau1, tau_+1) for1(tau2, tau1+1) 
            g2tau_uu[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau_]*cc[2][j][i][tau_-tau1]*ExpEtau[i][tau1-tau2]*cc[7][i][j][tau2];   
        if (nud[i]==nud[j]) 
        for1(tau_, Ntau4+1) for1(tau1, tau_+1) for1(tau2, tau1+1) 
            g2tau_dd[tau1][tau_][tau2]-=ExpEtau[j][Ntau4-tau_]*cc[6][j][i][tau_-tau1]*ExpEtau[i][tau1-tau2]*cc[9][i][j][tau2];   
    }

    
    //tau_>tau2>tau1             cd+ cu  cd  cu+   3-8,  2-7, 6-9  perm +
    for2(i, j, NFock)
    {
        if (nud[i]==nud[j]+1-(Nbath_raw+2)) 
        for1(tau_, Ntau4+1) for1(tau2, tau_+1) for1(tau1, tau2+1) 
            g2tau_du[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau_]*cc[3][j][i][tau_-tau2]*ExpEtau[i][tau2-tau1]*cc[8][i][j][tau1];     
        if (nud[i]==nud[j]) 
        for1(tau_, Ntau4+1) for1(tau2, tau_+1) for1(tau1, tau2+1) 
            g2tau_uu[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau_]*cc[2][j][i][tau_-tau2]*ExpEtau[i][tau2-tau1]*cc[7][i][j][tau1];   
        if (nud[i]==nud[j]) 
        for1(tau_, Ntau4+1) for1(tau2, tau_+1) for1(tau1, tau2+1) 
            g2tau_dd[tau1][tau_][tau2]+=ExpEtau[j][Ntau4-tau_]*cc[6][j][i][tau_-tau2]*ExpEtau[i][tau2-tau1]*cc[9][i][j][tau1];   
    }

 

    
    for1(tau1, Ntau4+1) for1(tau_, tau1+1) for1(tau2, tau_+1) 
    {
        g2tau_du[tau1][tau_][tau2]-=Z*gtau_u[tau2]*gtau_d[tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(gtau_u[tau2]*gtau_u[tau1-tau_]+gtau_u[tau1]*gtau_u[Ntau4+tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(gtau_d[tau2]*gtau_d[tau1-tau_]+gtau_d[tau1]*gtau_d[Ntau4+tau2-tau_]);
    }
     
    for1(tau1, Ntau4+1) for1(tau2, tau1+1) for1(tau_, tau2+1) 
    {
        g2tau_du[tau1][tau_][tau2]-=Z*gtau_u[tau2]*gtau_d[tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(gtau_u[tau2]*gtau_u[tau1-tau_]-gtau_u[tau1]*gtau_u[tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(gtau_d[tau2]*gtau_d[tau1-tau_]-gtau_d[tau1]*gtau_d[tau2-tau_]);        
    }
    
    for1(tau2, Ntau4+1) for1(tau_, tau2+1) for1(tau1, tau_+1)     
    {
        g2tau_du[tau1][tau_][tau2]-=-Z*gtau_u[tau2]*gtau_d[Ntau4+tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(-gtau_u[tau2]*gtau_u[Ntau4+tau1-tau_]-gtau_u[tau1]*gtau_u[tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(-gtau_d[tau2]*gtau_d[Ntau4+tau1-tau_]-gtau_d[tau1]*gtau_d[tau2-tau_]);        
    }

    for1(tau2, Ntau4+1) for1(tau1, tau2+1) for1(tau_, tau1+1)     
    {
        g2tau_du[tau1][tau_][tau2]-=Z*gtau_u[tau2]*gtau_d[tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(gtau_u[tau2]*gtau_u[tau1-tau_]-gtau_u[tau1]*gtau_u[tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(gtau_d[tau2]*gtau_d[tau1-tau_]-gtau_d[tau1]*gtau_d[tau2-tau_]);
    }

    for1(tau_, Ntau4+1) for1(tau1, tau_+1) for1(tau2, tau1+1) 
    {
        g2tau_du[tau1][tau_][tau2]-=-Z*gtau_u[tau2]*gtau_d[Ntau4+tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(-gtau_u[tau2]*gtau_u[Ntau4+tau1-tau_]+gtau_u[tau1]*gtau_u[Ntau4+tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(-gtau_d[tau2]*gtau_d[Ntau4+tau1-tau_]+gtau_d[tau1]*gtau_d[Ntau4+tau2-tau_]);
    }

    for1(tau_, Ntau4+1) for1(tau2, tau_+1) for1(tau1, tau2+1) 
    {
        g2tau_du[tau1][tau_][tau2]-=-Z*gtau_u[tau2]*gtau_d[Ntau4+tau1-tau_]; 
        g2tau_uu[tau1][tau_][tau2]-=Z*(-gtau_u[tau2]*gtau_u[Ntau4+tau1-tau_]+gtau_u[tau1]*gtau_u[Ntau4+tau2-tau_]);
        g2tau_dd[tau1][tau_][tau2]-=Z*(-gtau_d[tau2]*gtau_d[Ntau4+tau1-tau_]+gtau_d[tau1]*gtau_d[Ntau4+tau2-tau_]);
    }
    
    
    for2(tau, tau2, Ntau4+1)
        if (tau!=tau2) 
        {
            g2tau_du[tau][tau2][tau2]/=2.; g2tau_du[tau2][tau][tau2]/=2.; g2tau_du[tau2][tau2][tau]/=2.; 
            g2tau_uu[tau][tau2][tau2]/=2.; g2tau_uu[tau2][tau][tau2]/=2.; g2tau_uu[tau2][tau2][tau]/=2.;
            g2tau_dd[tau][tau2][tau2]/=2.; g2tau_dd[tau2][tau][tau2]/=2.; g2tau_dd[tau2][tau2][tau]/=2.;
        }
        else {g2tau_du[tau2][tau2][tau2]/=6.; g2tau_uu[tau2][tau2][tau2]/=6.; g2tau_dd[tau2][tau2][tau2]/=6.;}
     
    for3(tau1, tau_, tau2, Ntau4+1) 
    {
        g2tau_du[tau1][tau_][tau2]/=Z; g2tau_uu[tau1][tau_][tau2]/=Z; g2tau_dd[tau1][tau_][tau2]/=Z; 
//        if(abs(g2tau_du[tau1][tau_][tau2])>1e-10) cout<<tau1<<"  "<<tau_<<"  "<<tau2<<"  "<<g2tau_du[tau1][tau_][tau2]<<" ud\n";
//        if(abs(g2tau_uu[tau1][tau_][tau2])>1e-10) cout<<tau1<<"  "<<tau_<<"  "<<tau2<<"  "<<g2tau_uu[tau1][tau_][tau2]<<" uu\n";
//        if(abs(g2tau_dd[tau1][tau_][tau2])>1e-10) cout<<tau1<<"  "<<tau_<<"  "<<tau2<<"  "<<g2tau_dd[tau1][tau_][tau2]<<" dd\n";
    }
    
    for1(l,10) {for1(i, NFock) {for1(j, NFock) delete [] cc[l][i][j]; delete [] cc[l][i];} delete [] cc[l]; }

  
    
  
  
  
  
//    double *** g2tauT_uu=new_double3(Ntau4, Ntau4, Ntau4), *** g2tauT_du=new_double3(Ntau4, Ntau4, Ntau4), *** g2tauT_dd=new_double3(Ntau4, Ntau4, Ntau4);
    
//    for3(tau1, tau_, tau2, Ntau4)
//        if (tau1>=tau_) 
//            {g2tauT_uu[tau_][tau1-tau_][tau2]=g2tau_uu[tau1][tau_][tau2];         g2tauT_du[tau_][tau1-tau_][tau2]=g2tau_du[tau1][tau_][tau2];         g2tauT_dd[tau_][tau1-tau_][tau2]=g2tau_dd[tau1][tau_][tau2];}
//        else 
//            {g2tauT_uu[tau_][Ntau4+tau1-tau_][tau2]=-g2tau_uu[tau1][tau_][tau2]; g2tauT_du[tau_][Ntau4+tau1-tau_][tau2]=-g2tau_du[tau1][tau_][tau2]; g2tauT_dd[tau_][Ntau4+tau1-tau_][tau2]=-g2tau_dd[tau1][tau_][tau2];}
    

    complex *** g2tauW_uu=new_complex3(Ntau4, 2*Nw4, Ntau4), *** g2tauW_du=new_complex3(Ntau4, 2*Nw4, Ntau4), *** g2tauW_dd=new_complex3(Ntau4, 2*Nw4, Ntau4);

    for1(w_, Nw4) 
    {
        for2(tau1, tau2, Ntau4) {g2tauW_uu[tau1][Nw4+w_][tau2]=0; g2tauW_uu[tau1][Nw4-1-w_][tau2]=0; g2tauW_du[tau1][Nw4+w_][tau2]=0; g2tauW_du[tau1][Nw4-1-w_][tau2]=0; g2tauW_dd[tau1][Nw4+w_][tau2]=0; g2tauW_dd[tau1][Nw4-1-w_][tau2]=0;}
        for1(t_, Ntau4) 
        {
            complex r=(beta/Ntau4)*exp(-I*((2.*w_+1)*t_*Pi/Ntau4)), r_=(beta/Ntau4)*exp(I*((2.*w_+1)*t_*Pi/Ntau4)); 
            for2(tau1, tau2, Ntau4) 
            {
                g2tauW_uu[tau1][Nw4+w_][tau2]+=r*g2tau_uu[tau1][t_][tau2]; g2tauW_uu[tau1][Nw4-1-w_][tau2]+=r_*g2tau_uu[tau1][t_][tau2];
                g2tauW_du[tau1][Nw4+w_][tau2]+=r*g2tau_du[tau1][t_][tau2]; g2tauW_du[tau1][Nw4-1-w_][tau2]+=r_*g2tau_du[tau1][t_][tau2];
                g2tauW_dd[tau1][Nw4+w_][tau2]+=r*g2tau_dd[tau1][t_][tau2]; g2tauW_dd[tau1][Nw4-1-w_][tau2]+=r_*g2tau_dd[tau1][t_][tau2];
            }
        }
    }
    
    complex *** g2mixW_uu=new_complex3(2*Nw4, 2*Nw4, Ntau4), *** g2mixW_du=new_complex3(2*Nw4, 2*Nw4, Ntau4), *** g2mixW_dd=new_complex3(2*Nw4, 2*Nw4, Ntau4);
    
    for1(w1, Nw4)
    {
        for1(tau2, Ntau4) for1(w_, 2*Nw4) {g2mixW_uu[Nw4+w1][w_][tau2]=0;g2mixW_uu[Nw4-1-w1][w_][tau2]=0;g2mixW_du[Nw4+w1][w_][tau2]=0;g2mixW_du[Nw4-1-w1][w_][tau2]=0;g2mixW_dd[Nw4+w1][w_][tau2]=0;g2mixW_dd[Nw4-1-w1][w_][tau2]=0;}
        for1(tau1, Ntau4)
        {
            complex r=(beta/Ntau4)*exp(I*((2.*w1+1)*tau1*Pi/Ntau4)), r_=(beta/Ntau4)*exp(-I*((2.*w1+1)*tau1*Pi/Ntau4)); 
            for1(tau2, Ntau4) for1(w_, 2*Nw4) 
            {
                g2mixW_uu[Nw4+w1][w_][tau2]+=r*g2tauW_uu[tau1][w_][tau2]; g2mixW_uu[Nw4-1-w1][w_][tau2]+=r_*g2tauW_uu[tau1][w_][tau2];  
                g2mixW_du[Nw4+w1][w_][tau2]+=r*g2tauW_du[tau1][w_][tau2]; g2mixW_du[Nw4-1-w1][w_][tau2]+=r_*g2tauW_du[tau1][w_][tau2];
                g2mixW_dd[Nw4+w1][w_][tau2]+=r*g2tauW_dd[tau1][w_][tau2]; g2mixW_dd[Nw4-1-w1][w_][tau2]+=r_*g2tauW_dd[tau1][w_][tau2];
            }
        }
    }
    
    complex *** g2wW_uu=new_complex3(2*Nw4, 2*Nw4, 2*Nw4), *** g2wW_du=new_complex3(2*Nw4, 2*Nw4, 2*Nw4), *** g2wW_dd=new_complex3(2*Nw4, 2*Nw4, 2*Nw4);
    
    for1(w2, Nw4)
    {
        for1(w1, 2*Nw4) for1(w_, 2*Nw4) {g2wW_uu[w1][w_][Nw4+w2]=0;g2wW_uu[w1][w_][Nw4-1-w2]=0; g2wW_du[w1][w_][Nw4+w2]=0; g2wW_du[w1][w_][Nw4-1-w2]=0; g2wW_dd[w1][w_][Nw4+w2]=0; g2wW_dd[w1][w_][Nw4-1-w2]=0;}
        for1(tau2, Ntau4)
        {
            complex r=(beta/Ntau4)*exp(I*((2.*w2+1)*tau2*Pi/Ntau4)), r_=(beta/Ntau4)*exp(-I*((2.*w2+1)*tau2*Pi/Ntau4)); 
            for1(w1, 2*Nw4) for1(w_, 2*Nw4) 
            {
                g2wW_uu[w1][w_][Nw4+w2]+=r*g2mixW_uu[w1][w_][tau2]; g2wW_uu[w1][w_][Nw4-1-w2]+=r_*g2mixW_uu[w1][w_][tau2];  
                g2wW_du[w1][w_][Nw4+w2]+=r*g2mixW_du[w1][w_][tau2]; g2wW_du[w1][w_][Nw4-1-w2]+=r_*g2mixW_du[w1][w_][tau2];  
                g2wW_dd[w1][w_][Nw4+w2]+=r*g2mixW_dd[w1][w_][tau2]; g2wW_dd[w1][w_][Nw4-1-w2]+=r_*g2mixW_dd[w1][w_][tau2];  
            }
        }
    }
    
    
//    cout<<g2wW_du[0][0][0]<<"\n";
     
    
    for(int w1=-Nw4; w1<Nw4; w1++)  for(int w_=-Nw4; w_<Nw4; w_++) for(int w2=-Nw4; w2<Nw4; w2++)   
    {
        gamma_uu[Nw4+w1][Nw4+w_][Nw4+w2]=g2wW_uu[Nw4+w1][Nw4+w_][Nw4+w2]/(Gu(w1)*Gu(w_)*Gu(w2)*Gu(w1+w2-w_));         
        gamma_du[Nw4+w1][Nw4+w_][Nw4+w2]=g2wW_du[Nw4+w1][Nw4+w_][Nw4+w2]/(Gd(w1)*Gd(w_)*Gu(w2)*Gu(w1+w2-w_)); 
        gamma_dd[Nw4+w1][Nw4+w_][Nw4+w2]=g2wW_dd[Nw4+w1][Nw4+w_][Nw4+w2]/(Gd(w1)*Gd(w_)*Gd(w2)*Gd(w1+w2-w_)); 
        int w0=w2+w1-w_; if (w0>=-Nw4 && w0<Nw4) gamma_ud[Nw4+w2][Nw4+w0][Nw4+w1]=gamma_du[Nw4+w1][Nw4+w_][Nw4+w2];
        
        /*
        
        gamma_uu[W][Nw4+w1][Nw4+w2]=g2wW_uu[W][Nw4+w1][Nw4+w2]/(Gu(w1)*Gu(w1-W)*Gu(w2)*Gu(w2+W)); 
        gamma_du[W][Nw4+w1][Nw4+w2]=g2wW_du[W][Nw4+w1][Nw4+w2]/(Gd(w1)*Gd(w1-W)*Gu(w2)*Gu(w2+W)); 
        gamma_dd[W][Nw4+w1][Nw4+w2]=g2wW_dd[W][Nw4+w1][Nw4+w2]/(Gd(w1)*Gd(w1-W)*Gd(w2)*Gd(w2+W));

        gamma_uu[W][Nw4-1-w1][Nw4+w2]=g2wW_uu[W][Nw4-1-w1][Nw4+w2]/(Gu(-1-w1)*Gu(-1-w1-W)*Gu(w2)*Gu(w2+W)); 
        gamma_du[W][Nw4-1-w1][Nw4+w2]=g2wW_du[W][Nw4-1-w1][Nw4+w2]/(Gd(-1-w1)*Gd(-1-w1-W)*Gu(w2)*Gu(w2+W)); 
        gamma_dd[W][Nw4-1-w1][Nw4+w2]=g2wW_dd[W][Nw4-1-w1][Nw4+w2]/(Gd(-1-w1)*Gd(-1-w1-W)*Gd(w2)*Gd(w2+W)); 
       
        gamma_uu[W][Nw4+w1][Nw4-1-w2]=g2wW_uu[W][Nw4+w1][Nw4-1-w2]/(Gu(w1)*Gu(w1-W)*Gu(-1-w2)*Gu(-1-w2+W)); 
        gamma_du[W][Nw4+w1][Nw4-1-w2]=g2wW_du[W][Nw4+w1][Nw4-1-w2]/(Gd(w1)*Gd(w1-W)*Gu(-1-w2)*Gu(-1-w2+W)); 
        gamma_dd[W][Nw4+w1][Nw4-1-w2]=g2wW_dd[W][Nw4+w1][Nw4-1-w2]/(Gd(w1)*Gd(w1-W)*Gd(-1-w2)*Gd(-1-w2+W)); 
          
        gamma_uu[W][Nw4-1-w1][Nw4-1-w2]=g2wW_uu[W][Nw4-1-w1][Nw4-1-w2]/(Gu(-1-w1)*Gu(-1-w1-W)*Gu(-1-w2)*Gu(-1-w2+W)); 
        gamma_du[W][Nw4-1-w1][Nw4-1-w2]=g2wW_du[W][Nw4-1-w1][Nw4-1-w2]/(Gd(-1-w1)*Gd(-1-w1-W)*Gu(-1-w2)*Gu(-1-w2+W)); 
        gamma_dd[W][Nw4-1-w1][Nw4-1-w2]=g2wW_dd[W][Nw4-1-w1][Nw4-1-w2]/(Gd(-1-w1)*Gd(-1-w1-W)*Gd(-1-w2)*Gd(-1-w2+W)); 
       */
    }
      
    
//    for1(W,Nw4) 
//     {for(int w1=-2; w1<2; w1++) for(int w_=-2; w_<2; w_++) for(int w2=-2; w2<2; w2++) cout<<w1<<"\t"<<w_<<"\t"<<w2<<"\t"<<gamma_du[Nw4+w1][Nw4+w_][Nw4+w2]<<"\n"; cout<<"\n"<<flush;}


//cout<<gamma_du[Nw4-1][2*Nw4-1][2*Nw4-1]<<"   "<<gamma_uu[0][1][2*Nw4-1]<<"   "<<gamma_dd[0][1][2*Nw4-1]<<"\n";  
  
 //   cout<<"done\n"<<flush;   
  
    delete_double2(ExpEtau, NFock, Ntau4+1);
//    delete_double3(g2tauT_uu, Ntau4, Ntau4, Ntau4); delete_double3(g2tauT_du, Ntau4, Ntau4, Ntau4);delete_double3(g2tauT_dd, Ntau4, Ntau4, Ntau4);
    delete_complex3(g2tauW_uu, Ntau4, 2*Nw4, Ntau4); delete_complex3(g2tauW_du, Ntau4, 2*Nw4, Ntau4); delete_complex3(g2tauW_dd, Ntau4, 2*Nw4, Ntau4);
    delete_complex3(g2mixW_uu, 2*Nw4, 2*Nw4, Ntau4); delete_complex3(g2mixW_du, 2*Nw4, 2*Nw4, Ntau4); delete_complex3(g2mixW_dd, 2*Nw4, 2*Nw4, Ntau4);
    delete_complex3(g2wW_uu, 2*Nw4, 2*Nw4, 2*Nw4);   delete_complex3(g2wW_du, 2*Nw4, 2*Nw4, 2*Nw4); delete_complex3(g2wW_dd, 2*Nw4, 2*Nw4, 2*Nw4);
    
}


//=====================================================================================================================================================================


class solver_fine
{
 
    
    
    double * t2u_raw, * t2d_raw, *eu_raw, *ed_raw,   mu_raw, h_raw; 
    double * e_;
    
    
    int n(int j, int l) {return (j>>l)%2;}
    void  get_g_static(double dh, bool anomal);
    void  get_g(double t_, double e_, int flag); //flag: 0 -- up, 1 -- down, 2 -- a + x
//    void  get_g (double dmu, double dhz, double dhx, double dha);  
    
    complex * gu0, *gd0;
    complex * gu_, * gd_, * gx_, *ga_; 
    
    complex * dg0_mu, * dg0_hz, *dgz_mu, * dgz_hz, * dgx_hx, * dga_ha;   
    
    complex ** dg0_D0, ** dg0_Dz, **dgz_D0, ** dgz_Dz, ** dgx_Dx, ** dga_Da;

    void get_derivatives();
    
    
    double scalar(complex * r1, complex * r2) {complex s=0, s2=0, s3=0; for1(wn, Nw) {s+=conj(r1[wn])*r2[wn]; if (wn+1==Nw/2) s2=s; if (wn+1==3*Nw/4) s3=s; } return 2.*real(8.*s-9.*s3+2.*s2)/beta;  }
    
    complex *alphan; 
    void expand(complex * dDelta);
    
    
    
    public:
    
 
    
    solver_fine();//{t2u=new double [Nbath_raw]; t2d=new double [Nbath_raw]; eu=new double [Nbath_raw]; ed=new double [Nbath_raw]; gu=new complex[Nw]; gd=new complex[Nw]; gu0=new complex[Nw]; gd0=new complex[Nw]; guu1=new_complex2(Nbath_fine+1, Nw); gdd1=new_complex2(Nbath_fine+1, Nw); gud1=new_complex2(Nbath_fine+1, Nw);}
    ~solver_fine();// {delete [] t2u; delete [] t2d; delete [] eu; delete [] ed; delete [] gu; delete [] gd;delete [] gu0; delete [] gd0; delete_complex2(guu1, Nbath_fine+1, Nw); delete_complex2(gdd1, Nbath_fine+1, Nw); delete_complex2(gud1, Nbath_fine+1, Nw);}
    void init(solver_raw &, bool calc_g);


    complex * g0, * gz, * gu, * gd, * gx, * gy, *gax, *gay;
    complex * Delta_u, * Delta_d, * Delta_x, * Delta_y, *Delta_ax, *Delta_ay;
    double mu_fine, hz_fine, hx_fine,hy_fine, hax_fine, hay_fine;

    double sz() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(gz[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }
    double n() {double a=0, a2, a3; for1 (w, Nw) {a+=2.*real(g0[w])/beta; if (w+1==Nw/2) a2=a; if (w+1==3*Nw/4) a3=a;} return 8.*a-9.*a3+2.*a2; }

    
    void set_g_raw();
    void set_g();
};


solver_fine::solver_fine()
{
    t2u_raw=new double[Nbath_raw]; eu_raw=new double[Nbath_raw];   t2d_raw=new double[Nbath_raw]; ed_raw=new double[Nbath_raw];   
    e_=new double [Nbath_fine]; for1(ne, Nbath_fine) e_[ne]=Bath_fine_halfwidth*(-1.+ne*2./(Nbath_fine-1.));
    gu0=new complex [Nw]; gd0=new complex [Nw];
    gu_=new complex [Nw]; gd_=new complex [Nw]; gx_=new complex [Nw]; ga_=new complex [Nw];
    dg0_mu=new complex [Nw]; dg0_hz=new complex [Nw]; dgz_mu=new complex [Nw]; dgz_hz=new complex [Nw]; dgx_hx=new complex [Nw]; dga_ha=new complex [Nw];
    dg0_D0=new_complex2(Nbath_fine, Nw); dg0_Dz=new_complex2(Nbath_fine, Nw);dgz_D0=new_complex2(Nbath_fine, Nw);dgz_Dz=new_complex2(Nbath_fine, Nw);dgx_Dx=new_complex2(Nbath_fine, Nw);dga_Da=new_complex2(Nbath_fine, Nw);
    alphan=new complex [Nbath_fine];
    
    g0=new complex [Nw]; gz=new complex [Nw]; gu=new complex [Nw]; gd=new complex [Nw]; gx=new complex [Nw]; gy=new complex [Nw]; gax=new complex [Nw]; gay=new complex [Nw];
    Delta_u=new complex [Nw]; Delta_d=new complex [Nw]; Delta_x=new complex [Nw]; Delta_y=new complex [Nw]; Delta_ax=new complex [Nw]; Delta_ay=new complex [Nw];
    for1(w, Nw){Delta_u[w]=0.; Delta_d[w]=0; Delta_x[w]=0;Delta_y[w]=0;Delta_ax[w]=0;Delta_ay[w]=0;}
}

solver_fine::~solver_fine()
{
    delete [] t2u_raw; delete [] eu_raw; delete [] t2d_raw; delete [] ed_raw;
    delete [] e_;
    delete [] gu0; delete [] gd0; 
    delete [] gu_; delete [] gd_; delete [] gx_; delete [] ga_;
    delete [] dg0_mu; delete [] dg0_hz; delete [] dgz_mu; delete [] dgz_hz; delete [] dgx_hx; delete [] dga_ha;
    delete_complex2(dg0_D0, Nbath_fine, Nw); delete_complex2(dg0_Dz, Nbath_fine, Nw); delete_complex2(dgz_D0, Nbath_fine, Nw); delete_complex2(dgz_Dz, Nbath_fine, Nw); delete_complex2(dgx_Dx, Nbath_fine, Nw); delete_complex2(dga_Da, Nbath_fine, Nw);
    delete [] alphan;
    delete [] g0; delete [] gz; delete [] gu; delete [] gd; delete [] gx;delete [] gy;delete [] gax;delete [] gay; 
    delete [] Delta_u; delete [] Delta_d;delete [] Delta_x;delete [] Delta_y;delete [] Delta_ax;delete [] Delta_ay; 
}


void solver_fine::init(solver_raw & s, bool calc_g=true)
{
    for1(l, Nbath_raw) {t2u_raw[l]=s.t2u[l]; t2d_raw[l]=s.t2d[l]; eu_raw[l]=s.eu[l]; ed_raw[l]=s.ed[l];} mu_raw=s.mu_loc; h_raw=s.h_loc; 
    
    for1(w, Nw) {g0[w]=s.gu[w]+s.gd[w]; gz[w]=s.gu[w]-s.gd[w]; gu[w]=s.gu[w]; gd[w]=s.gd[w]; gx[w]=0; gy[w]=0; gax[w]=0; gay[w]=0;}     

    if (calc_g) {get_g_static(0, false); for1(w, Nw) {gu0[w]=gu_[w]; gd0[w]=gd_[w];}    get_derivatives(); }   
}


void solver_fine::get_g_static (double dh, bool anomal)
{
    set_Iw(); 
    
    int NFock1=1<<(2*Nbath_raw+2);

    static double ** H = new_double2(NFock1, NFock1), ** Psi = new_double2(NFock1, NFock1), ** cu = new_double2(NFock1, NFock1), ** cd = new_double2(NFock1, NFock1), ** cuPsi = new_double2(NFock1, NFock1), ** cdPsi = new_double2(NFock1, NFock1), *E=new double [NFock1];
    for2(i,j, NFock1) {H[i][j]=0; cuPsi[i][j]=0; cdPsi[i][j]=0; cu[i][j]=0; cd[i][j]=0; }

    
//    for1(j, NFock1) {nud[j]=0; for1(l, Nbath_raw+1) nud[j]+=n(j,l)+(Nbath_raw+2)*n(j,l+Nbath_raw+1); } 
    
    
    for1(j, NFock1)
    {
        H[j][j]+=-mu_raw*(n(j,0)+n(j, Nbath_raw+1))-h_raw*((n(j,0)-n(j, Nbath_raw+1)))+U*((n(j,0)-.5)*(n(j,Nbath_raw+1)-.5));
        
        if (anomal)
        {
            H[j][j]+=-dh*(n(j,0)-n(j, Nbath_raw+1));
            if (n(j,0)==0 && n(j,Nbath_raw+1)==0) {int j1=j+1+(1<<(Nbath_raw+1));  int nu=0; for1(l,Nbath_raw) nu+=n(j,l+1); H[j][j1]-=(1-2*(nu%2))*dh; H[j1][j]=H[j][j1];}

        }
        else
        {
            H[j][j]+=-dh*(n(j,0)+n(j, Nbath_raw+1));             
            if (n(j,0)==0 && n(j,Nbath_raw+1)==1) {int j1=j+1-(1<<(Nbath_raw+1));  int nu=0; for1(l,Nbath_raw) nu+=n(j,l+1); H[j][j1]-=(1-2*(nu%2))*dh; H[j1][j]=H[j][j1];}
        }
        
        
        for1(l, Nbath_raw) H[j][j]+=n(j,l+1)*eu_raw[l]+n(j,Nbath_raw+l+2)*ed_raw[l];
        
        for1(l, Nbath_raw)
        {
            int nu=0, nd=0; for1(l2,l)   {nu+=n(j,l2+1); nd+=n(j,Nbath_raw+l2+2);}
            if(n(j,0)==0  && n(j, l+1)==1)                 {int j1=j+1-(1<<(l+1));                        H[j][j1]+=sqrt(t2u_raw[l])*(1-2*(nu%2));H[j1][j]=H[j][j1];}
            if(n(j,Nbath_raw+1)==0  && n(j, l+Nbath_raw+2)==1)   {int j1=j+(1<<(Nbath_raw+1))-(1<<(l+Nbath_raw+2));   H[j][j1]+=sqrt(t2d_raw[l])*(1-2*(nd%2));H[j1][j]=H[j][j1];}
        }
    }
    
  
    EigenJacobi(H, Psi, E, NFock1);
    
    
    
    
    double Z=0; for1(j, NFock1) Z+=expl(-beta*E[j]);     
 
    for1(j, NFock1) 
    {
        if (n(j,0)==1) for1(j1, NFock1)  cuPsi[j1][j-1]=Psi[j1][j];
        if (n(j, Nbath_raw+1)==1) {double f=1; for1(l, Nbath_raw+1) f*=1-2*n(j, l); //1-2*((nud[j]%(Nbath_raw+2))%2); 
            for1(j1, NFock1) cdPsi[j1][j-(1<<(Nbath_raw+1))]=f*Psi[j1][j];  }
    }
    
    for1 (w, Nw) {gu_[w]=0; gd_[w]=0; gx_[w]=0; ga_[w]=0;}

    for2(j1, j2, NFock1) 
    {
            double d=0; for1(j3, NFock1) d+=Psi[j1][j3]*cuPsi[j2][j3]; cu[j1][j2]=d;
                   d=0; for1(j3, NFock1) d+=Psi[j1][j3]*cdPsi[j2][j3]; cd[j1][j2]=d;
    }
    
    for2(j1, j2, NFock1) 
    {
            double x;
            x=cu[j1][j2]*cu[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gu_[w]+=x/(Iw[w]+E[j1]-E[j2]);
            x=cd[j1][j2]*cd[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gd_[w]+=x/(Iw[w]+E[j1]-E[j2]);
            x=cu[j1][j2]*cd[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gx_[w]+=2.*x/(Iw[w]+E[j1]-E[j2]);
            x=cu[j1][j2]*cd[j2][j1]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) ga_[w]+=2.*x/(Iw[w]+E[j1]-E[j2]);            
    }

      
}



void  solver_fine::get_g(double t_, double e_, int flag) //flag: 0 -- up, 1 -- down, 2 -- a + x
{
    set_Iw(); 
    
    int NFock1=1<<(2*Nbath_raw+3);

    static double ** H = new_double2(NFock1, NFock1), ** Psi = new_double2(NFock1, NFock1), ** cu = new_double2(NFock1, NFock1), ** cd = new_double2(NFock1, NFock1), ** cuPsi = new_double2(NFock1, NFock1), ** cdPsi = new_double2(NFock1, NFock1), *E=new double [NFock1];
    for2(i,j, NFock1) {H[i][j]=0; cuPsi[i][j]=0; cdPsi[i][j]=0; cu[i][j]=0; cd[i][j]=0; }

        
    
    for1(j, NFock1)
    {
        H[j][j]+=-mu_raw*(n(j,0)+n(j, Nbath_raw+1))-h_raw*((n(j,0)-n(j, Nbath_raw+1)))+U*((n(j,0)-.5)*(n(j,Nbath_raw+1)-.5));
        
        
        
        for1(l, Nbath_raw) H[j][j]+=n(j,l+1)*eu_raw[l]+n(j,Nbath_raw+l+2)*ed_raw[l];
        
        for1(l, Nbath_raw)
        {
            int nu=0, nd=0; for1(l2,l)   {nu+=n(j,l2+1); nd+=n(j,Nbath_raw+l2+2);}
            if(n(j,0)==0  && n(j, l+1)==1)                 {int j1=j+1-(1<<(l+1));                        H[j][j1]+=sqrt(t2u_raw[l])*(1-2*(nu%2));H[j1][j]=H[j][j1];}
            if(n(j,Nbath_raw+1)==0  && n(j, l+Nbath_raw+2)==1)   {int j1=j+(1<<(Nbath_raw+1))-(1<<(l+Nbath_raw+2));   H[j][j1]+=sqrt(t2d_raw[l])*(1-2*(nd%2));H[j1][j]=H[j][j1];}
        }
        
        
        H[j][j]+=e_*n(j,2*Nbath_raw+2);
        if (flag==0)
        {
            if (n(j,0)==0 && n(j,2*Nbath_raw+2)==1)  {int j1=j+1-(1<<(2*Nbath_raw+2));                   int nu=0; for1(l,2*Nbath_raw+1) nu+=n(j,l+1);       H[j][j1]+=(1-2*(nu%2))*t_; H[j1][j]=H[j][j1];}
        }
        if (flag==1)
        {
            if (n(j,Nbath_raw+1)==0 && n(j,2*Nbath_raw+2)==1) {int j1=j+(1<<(Nbath_raw+1))-(1<<(2*Nbath_raw+2));  int nu=0; for1(l,Nbath_raw) nu+=n(j,l+Nbath_raw+2); H[j][j1]+=(1-2*(nu%2))*t_; H[j1][j]=H[j][j1];}
        }
        if (flag==2)
        {
            if (n(j,0)==0 && n(j,2*Nbath_raw+2)==1)           {int j1=j+1-(1<<(2*Nbath_raw+2));                   int nu=0; for1(l,2*Nbath_raw+1) nu+=n(j,l+1);       H[j][j1]+=(1-2*(nu%2))*t_; H[j1][j]=H[j][j1];}
            if (n(j,Nbath_raw+1)==0 && n(j,2*Nbath_raw+2)==1) {int j1=j+(1<<(Nbath_raw+1))-(1<<(2*Nbath_raw+2));  int nu=0; for1(l,Nbath_raw) nu+=n(j,l+Nbath_raw+2); H[j][j1]+=(1-2*(nu%2))*t_; H[j1][j]=H[j][j1];}   
            if (n(j,Nbath_raw+1)==1 && n(j,2*Nbath_raw+2)==1) {int j1=j-(1<<(Nbath_raw+1))-(1<<(2*Nbath_raw+2));  int nu=0; for1(l,Nbath_raw) nu+=n(j,l+Nbath_raw+2); H[j][j1]-=(1-2*(nu%2))*t_; H[j1][j]=H[j][j1];}   
        }
         
        
    }
    
  
    EigenJacobi(H, Psi, E, NFock1);
    
    
    
    
    double Z=0; for1(j, NFock1) Z+=expl(-beta*E[j]);     
 
    for1(j, NFock1) 
    {
        if (n(j,0)==1) for1(j1, NFock1)  cuPsi[j1][j-1]=Psi[j1][j];
        if (n(j, Nbath_raw+1)==1) {double f=1; for1(l, Nbath_raw+1) f*=1-2*n(j, l); //1-2*((nud[j]%(Nbath_raw+2))%2); 
            for1(j1, NFock1) cdPsi[j1][j-(1<<(Nbath_raw+1))]=f*Psi[j1][j];  }
    }
    
    for1 (w, Nw) {gu_[w]=0; gd_[w]=0; gx_[w]=0; ga_[w]=0;}

    for2(j1, j2, NFock1) 
    {
            double d=0; for1(j3, NFock1) d+=Psi[j1][j3]*cuPsi[j2][j3]; cu[j1][j2]=d;
                   d=0; for1(j3, NFock1) d+=Psi[j1][j3]*cdPsi[j2][j3]; cd[j1][j2]=d;
    }
    
    for2(j1, j2, NFock1) 
    {
            double x;
            x=cu[j1][j2]*cu[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gu_[w]+=x/(Iw[w]+E[j1]-E[j2]);
            x=cd[j1][j2]*cd[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gd_[w]+=x/(Iw[w]+E[j1]-E[j2]);
            x=cu[j1][j2]*cd[j1][j2]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) gx_[w]+=2.*x/(Iw[w]+E[j1]-E[j2]);
            x=cu[j1][j2]*cd[j2][j1]*( (expl(-beta*E[j1])+expl(-beta*E[j2]))/Z ); for1(w, Nw) ga_[w]+=2.*x/(Iw[w]+E[j1]-E[j2]);            
    }

      
}


void solver_fine::get_derivatives()
{
    double dh=.001;
    get_g_static(.5*dh, false);     for1(w, Nw) {dg0_mu[w] =(gu_[w]+gd_[w])/dh;  dgz_mu[w] =(gu_[w]-gd_[w])/dh;  dgx_hx[w] =gx_[w]/dh;}  
    get_g_static(-.5*dh, false);    for1(w, Nw) {dg0_mu[w]-=(gu_[w]+gd_[w])/dh;  dgz_mu[w]-=(gu_[w]-gd_[w])/dh;  dgx_hx[w]-=gx_[w]/dh;}  
    
    
    get_g_static(.5*dh, true);      for1(w, Nw) {dg0_hz[w] =(gu_[w]+gd_[w])/dh;  dgz_hz[w] =(gu_[w]-gd_[w])/dh;  dga_ha[w] =ga_[w]/dh;}  
    get_g_static(-.5*dh, true);     for1(w, Nw) {dg0_hz[w]-=(gu_[w]+gd_[w])/dh;  dgz_hz[w]-=(gu_[w]-gd_[w])/dh;  dga_ha[w]-=ga_[w]/dh;}  

    
    
    double t_=.01;
    for1(ne, Nbath_fine) 
    {//cout<<ne<<"\n"<<flush;
         get_g(t_, e_[ne], 0);        for1(w, Nw) {dg0_D0[ne][w] =2.*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_D0[ne][w]  =2.*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_); dg0_Dz[ne][w] =2.*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_Dz[ne][w]  =2.*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_);}        
         get_g(t_*sqrt(2), e_[ne], 0);for1(w, Nw) {dg0_D0[ne][w]-=.5*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_D0[ne][w] -=.5*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_); dg0_Dz[ne][w]-=.5*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_Dz[ne][w] -=.5*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_);}
 
         get_g(t_, e_[ne], 1);        for1(w, Nw) {dg0_D0[ne][w]+=2.*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_D0[ne][w] +=2.*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_); dg0_Dz[ne][w]-=2.*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_Dz[ne][w] -=2.*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_);}        
         get_g(t_*sqrt(2), e_[ne], 1);for1(w, Nw) {dg0_D0[ne][w]-=.5*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_D0[ne][w] -=.5*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_); dg0_Dz[ne][w]+=.5*(gu_[w] +gd_[w]-gu0[w]-gd0[w])/(-2.*t_*t_); dgz_Dz[ne][w] +=.5*(gu_[w] -gd_[w]-gu0[w]+gd0[w])/(-2.*t_*t_);}        
         
         get_g(t_, e_[ne], 2);        for1(w, Nw) {dgx_Dx[ne][w] =2.*gx_[w]/(-2.*t_*t_); dga_Da[ne][w] =2.*ga_[w]/(-2.*t_*t_);} 
         get_g(t_*sqrt(2), e_[ne], 2);for1(w, Nw) {dgx_Dx[ne][w]-=.5*gx_[w]/(-2.*t_*t_); dga_Da[ne][w]-=.5*ga_[w]/(-2.*t_*t_);}    
    }
 
    
}


void solver_fine::expand(complex *dDelta)
{
    static double  ** Pseudo1=new_double2(Nbath_fine, Nbath_fine);
    
    do_once
    {
        static double **m =new_double2(Nbath_fine, Nbath_fine), ** psik=new_double2(Nbath_fine, Nbath_fine); 
        int Nk;

        static double * e=new double [Nbath_fine]; static complex * r1=new complex [Nw], * r2=new complex [Nw];
        for2(i,j, Nbath_fine) {for1(w, Nw) {r1[w]=1./(Iw[w]-e_[i]); r2[w]=1./(Iw[w]-e_[j]);}  m[i][j]=scalar(r1,r2);}
        
  //      for1(i, Nbath_fine) {for1(j, Nbath_fine) cout<<m[i][j]<<"  "; cout<<"\n";}
        
        Nk=Eigen_leading(m, psik, e, Nbath_fine, 1e-5); //cout<<"Effective basis Nk="<<Nk<<"\n";
        
        
        for1(k, Nk)
        {
            for1(w, Nw) {r1[w]=0; for1 (j, Nbath_fine) r1[w]+=psik[k][j]/(Iw[w]-e_[j]);}
            e[k]=scalar(r1,r1);
        }
        
        
        for2(i,j, Nbath_fine) {Pseudo1[i][j]=0; for1(k, Nk) Pseudo1[i][j]+=psik[k][i]*psik[k][j]/e[k];}
    }
    end_do_once;
    
    
//    double s=0; for1(j, Nbath_fine) s+=sqr(psik[0][j]); cout<<s<<"\n";
    
    for1(j, Nbath_fine) alphan[j]=0;
    
    for1(j_, Nbath_fine)
    {
        static complex * r1=new complex [Nw]; for1(w, Nw) r1[w]=-1./(Iw[w]-e_[j_]);
        double r=scalar(r1, dDelta);
        for1(j, Nbath_fine) alphan[j]+=Pseudo1[j][j_]*r;
    }
    
     
     
//    for1(w, 5)    {complex s=0; for1(j, Nbath_fine) s+=alphan[j]/(Iw[w]-e_[j]);cout<<dDelta[w]<<"  "<<s<<"\n";}
    
//    cout<<"\n";    
    
}

void solver_fine::set_g_raw()
{
    for1(w, Nw) {g0[w]=gu0[w]+gd0[w]; gz[w]=gu0[w]-gd0[w]; gu[w]=gu0[w]; gd[w]=gd0[w]; gx[w]=0; gy[w]=0; gax[w]=0; gay[w]=0;}     
}

void solver_fine::set_g()
{
    for1(w, Nw) {g0[w]=gu0[w]+gd0[w]; gz[w]=gu0[w]-gd0[w]; gu[w]=gu0[w]; gd[w]=gd0[w]; gx[w]=0; gy[w]=0; gax[w]=0; gay[w]=0;}     
 
     
    static complex * dDelta_0=new complex [Nw], * dDelta_z=new complex [Nw];
    
    for1(w, Nw) 
    {
        dDelta_0[w]=Delta_u[w]+Delta_d[w]; dDelta_z[w]=Delta_u[w]-Delta_d[w]; 
        for1(j, Nbath_raw) {dDelta_0[w]-=t2u_raw[j]/(Iw[w]-eu_raw[j])+t2d_raw[j]/(Iw[w]-ed_raw[j]); dDelta_z[w]-=t2u_raw[j]/(Iw[w]-eu_raw[j])-t2d_raw[j]/(Iw[w]-ed_raw[j]); }
    }
    
    
    expand(dDelta_0); for1(w, Nw) {for1 (j, Nbath_fine) g0[w]+=dg0_D0[j][w]*alphan[j]; g0[w]+=dg0_mu[w]*(mu_fine-mu_raw); }
                      for1(w, Nw) {for1 (j, Nbath_fine) gz[w]+=dgz_D0[j][w]*alphan[j]; gz[w]+=dgz_mu[w]*(mu_fine-mu_raw); }
    
    
    expand(dDelta_z); for1(w, Nw) {for1 (j, Nbath_fine) g0[w]+=dg0_Dz[j][w]*alphan[j]; g0[w]+=dg0_hz[w]*(hz_fine-h_raw); }
                      for1(w, Nw) {for1 (j, Nbath_fine) gz[w]+=dgz_Dz[j][w]*alphan[j]; gz[w]+=dgz_hz[w]*(hz_fine-h_raw); }

    for1(w, Nw) {gu[w]=.5*(g0[w]+gz[w]); gd[w]=.5*(g0[w]-gz[w]);}                   
    
    expand(Delta_x);  for1(w, Nw) {for1 (j, Nbath_fine) gx[w]+=dgx_Dx[j][w]*alphan[j]; gx[w]+=dgx_hx[w]*hx_fine;  }
    expand(Delta_y);  for1(w, Nw) {for1 (j, Nbath_fine) gy[w]+=dgx_Dx[j][w]*alphan[j]; gy[w]+=dgx_hx[w]*hx_fine;  }    
    
    expand(Delta_ax); for1(w, Nw) {for1 (j, Nbath_fine) gax[w]+=dga_Da[j][w]*alphan[j];gax[w]+=dga_ha[w]*hax_fine;}
    expand(Delta_ay); for1(w, Nw) {for1 (j, Nbath_fine) gay[w]+=dga_Da[j][w]*alphan[j];gay[w]+=dga_ha[w]*hay_fine;}    
    
}    




void bath_f2r(solver_fine & Sf, solver_raw & Sr, double lambda0=0) {Sr.init(Sf.Delta_u, Sf.Delta_d,Sf.mu_fine,Sf.hz_fine-lambda0*Sf.sz() );}
void bath_r2f(solver_raw & Sr, solver_fine & Sf) {Sf.init(Sr); Sf.set_g();}



//=====================================================================================================================================================================


double eps(int kx, int ky) {return -2.*t*(cos(kx*2*Pi/Lx)+cos(ky*2*Pi/Ly))- 2.*t_1*cos(kx*2*Pi/Lx+ky*2*Pi/Ly)- 2.*t_2*cos(kx*2*Pi/Lx-ky*2*Pi/Ly)  ; }

struct DMFT_AF_solver
{
    int AF_flag;
    double lambda0;
    solver_raw Sr; solver_fine Sf;
    DMFT_AF_solver(double lambda_=0) {AF_flag=1; lambda0=lambda_;}
    ~DMFT_AF_solver(){;}
    void iter_Delta(double factor);
    void init(double hz);
    int raw_loops(double factor);
    int fine_loops(double factor);    
};



void DMFT_AF_solver::init(double hz=0)
{
        set_Iw(); for1(w, Nw) {Sf.Delta_u[w]=.1/(Iw[w]+.1*I);  Sf.Delta_u[w]=.1/(Iw[w]+.1*I);}  
        Sf.mu_fine=mu; Sf.hz_fine=hz; 
        bath_f2r(Sf,Sr); 
}


void DMFT_AF_solver::iter_Delta(double factor)
{
    
     
    for1 (w, Nw)
    {
        complex su=0, sd=0; 
        for1(kx, Lx) for1(ky, Ly/2)
        {
            complex au, ad, g00, g11, g01, g10, d;
            
            au=1./Sf.gu[w]+Sf.Delta_u[w]; ad=1./Sf.gd[w]+Sf.Delta_d[w];
            g00=-eps(kx, ky) +.5*(au+ad), g11=-eps(kx+Lx/2, ky+Ly/2)+.5*(au+ad), g01= .5*(au-ad), g10=g01;
            d=g00*g11-g01*g10;
            
            su+=(g11+g00-2.*g10)/(d*double(N)); 
            sd+=(g11+g00+2.*g10)/(d*double(N)); 
        }
        Sf.Delta_u[w]+=factor*(1./Sf.gu[w]-1./su);
        Sf.Delta_d[w]+=factor*(1./Sf.gd[w]-1./sd);                 
        if (AF_flag==0) {Sf.Delta_u[w]=.5*(Sf.Delta_u[w]+Sf.Delta_d[w]); Sf.Delta_d[w]=Sf.Delta_u[w];}
    }
    
         

}

int DMFT_AF_solver::raw_loops(double factor)
{
    complex g_=-1.;   
    for1 (i, 1000) {bath_f2r(Sf,Sr, lambda0);   Sf.init(Sr, false); iter_Delta(factor); if (abs(Sf.gu[0]-g_)>1e-4) g_=Sf.gu[0]; else return i;}
    return -1;
}

int DMFT_AF_solver::fine_loops(double factor)
{
    complex g_=-1.; Sf.init(Sr);  double h_=Sf.hz_fine; //cout<<">"<<flush;
    for1 (i, 1000) { 
        Sf.set_g();  iter_Delta(factor); Sf.hz_fine=h_-lambda0*Sf.sz();  if (abs(Sf.gu[0]-g_)>1e-6) g_=Sf.gu[0]; else {Sf.hz_fine=h_; return i;}}//Sr.h_loc-0.*lambda0*Sf.sz()
    Sf.hz_fine=h_; return -1;
}

 






