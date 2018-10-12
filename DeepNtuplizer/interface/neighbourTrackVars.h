#ifndef neighbourTrackVars_h
#define neighbourTrackVars_h



class neighbourTrackVars {
      public:
          
    double pt, eta, phi, mass, dz, dxy; //kinematics    
    double t3Dip, t3Dsip, t2Dip, t2Dsip; //impact par
    
    double dist, dsig; //distance from neighbour
    double PCA_sx, PCA_sy, PCA_sz, PCA_sxerr, PCA_syerr, PCA_szerr; //PCA on seeding track
    double PCA_tx, PCA_ty, PCA_tz, PCA_txerr, PCA_tyerr, PCA_tzerr; //PCA on neighbour 
    
    double dotprodTrack, dotprodSeed; // angle between track at PCA and its direction
    double dotprodTrackSeed2D, dotprodTrackSeed3D; //angle between tracks and neighbours
    double dotprodTrackSeedVectors2D, dotprodTrackSeedVectors3D; //angle between directions
    
    double seedPCA_pv, trackPCA_pv; //PCA-pv sidtance
    double seedMass; //seeding track mass    
    
    double PCA_JetAxis_distance, PCAPair_Jet_dotprod, PCAAxis_JetAxis_DEta, PCAAxis_JetAxis_DPhi; //PCA jAxis variables
    
    int index, seed_index; //indexing
    
    void set_PtEtaPhiMassDzDxy(double, double, double, double, double);
    void set_IPs(double, double, double, double);
    
    void set_PCAdistance ( double, double );
    void set_PCAonSeedXYZ( double, double, double,double, double, double );
    void set_PCAonTrackXYZ( double, double, double,double, double, double );
    
    void set_dotProds(double, double, double,double, double, double );
    void set_PVdistance ( double, double );
    
    void set_JetAxisVars (double, double, double, double);
    
    void set_index (int);
    void set_SeedIndex (int);      
    void setSeedMass (double);
    
    
    void set_values (double, double, double, double,  double, double, double,
    double, double, double, double,  double, double,
    double, double, double, double,  double, double,
    double, double
    );
    
     
    void set_vars (double, double, double, double,
    double, double, double, double, double
    );
    


};

inline void neighbourTrackVars::set_PtEtaPhiMassDzDxy(double pt2, double eta2, double phi2, double m2, double dz2, double dxy2){

    pt=pt2;
    eta=eta2;
    phi=phi2;
    mass=m2;
    dz=dz2; 
    dxy=dxy2; 
    
}

inline void neighbourTrackVars::set_IPs(double t2dip, double t2dsip, double t3dip, double t3dsip){

    t3Dip=t3dip;
    t3Dsip=t3dsip;
    t2Dip=t2dip;
    t2Dsip=t2dsip;

}


inline void neighbourTrackVars::set_PCAdistance(double distaaa, double dsig2){
 
    dist=distaaa;
    dsig=dsig2;
    
} 


inline void neighbourTrackVars::set_PCAonSeedXYZ( double PCA_sx2, double PCA_sy2, double PCA_sz2, double PCA_sxerr2, double PCA_syerr2, double PCA_szerr2 ){
    
    PCA_sx=PCA_sx2;
    PCA_sy=PCA_sy2;
    PCA_sz=PCA_sz2; 
    PCA_sxerr=PCA_sxerr2;
    PCA_syerr=PCA_syerr2;
    PCA_szerr=PCA_szerr2;

}


inline void neighbourTrackVars::set_PCAonTrackXYZ( double PCA_tx2, double PCA_ty2, double PCA_tz2, double PCA_txerr2, double PCA_tyerr2, double PCA_tzerr2 ){
    
    PCA_tx=PCA_tx2;
    PCA_ty=PCA_ty2;
    PCA_tz=PCA_tz2; 
    PCA_txerr=PCA_txerr2;
    PCA_tyerr=PCA_tyerr2;
    PCA_tzerr=PCA_tzerr2;

}

inline void neighbourTrackVars::set_dotProds(double dotprodTrack2, double dotprodSeed2, double t2dTS, double t3dTS, double t2dTSV, double t3dTSV) {
    
    dotprodTrack=dotprodTrack2;
    dotprodSeed=dotprodSeed2;    
    dotprodTrackSeed2D=t2dTS;
    dotprodTrackSeed3D=t3dTS;
    dotprodTrackSeedVectors2D=t2dTSV;
    dotprodTrackSeedVectors3D=t3dTSV;
  
}

inline void neighbourTrackVars::set_PVdistance(double a, double b){

    seedPCA_pv=a;
    trackPCA_pv=b;
}

inline void neighbourTrackVars::set_JetAxisVars(double jadist, double dotprod, double d_eta, double d_phi){
    PCA_JetAxis_distance=jadist;
    PCAPair_Jet_dotprod=dotprod;
    PCAAxis_JetAxis_DEta=d_eta;
    PCAAxis_JetAxis_DPhi=d_phi; 
}


inline void neighbourTrackVars::set_index (int a){

    index=a;
    
}

inline void neighbourTrackVars::set_SeedIndex (int a){

    seed_index=a;
}
     
inline void neighbourTrackVars::setSeedMass(double sm){
    
    seedMass=sm;
}


inline void neighbourTrackVars::set_values (double pt2, double eta2, double phi2, double dz2, double dxy2, double distaaa, double dsig2,
double PCA_sx2, double PCA_sy2, double PCA_sz2, double PCA_sxerr2, double PCA_syerr2, double PCA_szerr2, 
double PCA_tx2, double PCA_ty2, double PCA_tz2, double PCA_txerr2, double PCA_tyerr2, double PCA_tzerr2,
double dotprodTrack2, double dotprodSeed2) {

    pt=pt2;
    eta=eta2;
    phi=phi2; 
    dz=dz2; 
    dxy=dxy2; 
    dist=distaaa;
    dsig=dsig2;
    PCA_sx=PCA_sx2;
    PCA_sy=PCA_sy2;
    PCA_sz=PCA_sz2; 
    PCA_sxerr=PCA_sxerr2;
    PCA_syerr=PCA_syerr2;
    PCA_szerr=PCA_szerr2;
    PCA_tx=PCA_tx2;
    PCA_ty=PCA_ty2;
    PCA_tz=PCA_tz2; 
    PCA_txerr=PCA_txerr2;
    PCA_tyerr=PCA_tyerr2;
    PCA_tzerr=PCA_tzerr2;
    dotprodTrack=dotprodTrack2;
    dotprodSeed=dotprodSeed2;
   
}

inline void neighbourTrackVars::set_vars ( double m, double t2dip, double t2dsip, double t3dip, double t3dsip, double t2dTS, double t3dTS, double t2dTSV, double t3dTSV){
mass=m;
t3Dip=t3dip;
t3Dsip=t3dsip;
t2Dip=t2dip;
t2Dsip=t2dsip;
dotprodTrackSeed2D=t2dTS;
dotprodTrackSeed3D=t3dTS;
dotprodTrackSeedVectors2D=t2dTSV;
dotprodTrackSeedVectors3D=t3dTSV;
    
}


struct sortfunctionNTracks
{
    inline bool operator() (const neighbourTrackVars& struct1, const neighbourTrackVars& struct2)
    {
        return (struct1.dist < struct2.dist);
    }
};


#endif
