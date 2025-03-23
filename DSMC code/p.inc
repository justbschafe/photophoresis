C     file PARS.INC
C
C  measure of length - radius of disc a
C  measure of velocity - most probable speed v_m
C  measure of time - a/v_m
C  measure of temperature - equailibrium temperature T_0
c
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C  entrance data
      PARAMETER(
     &   igas=2, ! 1 - He3, 2 -He4, 3 - Ne, 4 - Ar, 5 - Kr, 6 - Xe ! gas species
     &   del=1.d0, ! rarefcation parameter
     &   beta=0.5d0, ! porosity
     &   T0=300.d0, ! eq. temperature
     &   Tc=1.d0, ! top temperature in T0
     &   Th=4.d0/3.d0, ! bottom temperature in T0
     &   alpt=1.0d0, ! TMAC
     &   alpn=1.0d0, ! NEAC
     &   DTM=0.01d0, ! time increment 
     &   INMt=8*10**5, ! total number of model particles
     &   NC=20, ! number of cells over disc radius
     &   NS=20, ! number of surface cells 
     &   np=60, ! number of cores
     &   INM=INMt/np, ! number of particls per core
     &   tstead=200.d0, ! time to reach steady state
     &   tstat=100 000.d0, ! time to get statistic 
     &   RR=3.d0,  ! radius of comp.domain        
     &   SRt=6.d0,  ! total length of comp.domain
     &   SR=SRt/np, ! length of comp.domain per core
     &   SC=SRt/2., ! position disc center       
     &   NCX=NC*(SR+0.001),! number of cells in longitudinal direction
     &   NCY=NC*(RR+0.001),! number of cells in radial direction
     &   jp=100*(1.+1.d-5)/DTM, ! period to print 
     &   NRW=NCY-1,
     &   NCt=NCY*NCX, ! total number of cells per core
     &   MNM=2*INM,! maximum number of model particles
     &   NNWML=MNM/3,! maximum number of new particles      
     &   INMRM=MNM/3) 
C---------------------------------------------------------------------------
      COMMON/GEN/
     &  AMEU(NCX),!number of molecules entering at each time step
     &  AMEL(2,NCY),
     &  AMRU(NCX),! reminder
     &  AMRL(2,NCY),
     &  VC(NCt),CC(NCt), ! volume of cells
     &  CCG(2,NCt),! maximum relative velocity
     &  CSC(4,NCx,NCy),! summary characteristics of cell      
     &  ANC(NCx,NCy), ! number density in n0
     &  UXC(NCx,NCy),UYC(NCx,NCy),! velocity in v_m
     &  TEMPC(NCx,NCy), ! temperature in T0
     &  CSS(2,NS), ! fluxes to surface 
     &  PP(2,MNM),PV(3,MNM),! positions and velocities
     &  ppbuf(2,NNWML), ! positions of particles im buffer
     &  pvbuf(3,NNWML), ! velocities of particles in buffer 
     &  TIME,PI,SPI,CD,CQ,CW,DENM,RW(NRW),WW(NRW),PWW(NCY),
     &  Xchi(100,900), ! deflection matrix
     &  TCS(900),  ! total cross section
     &  vm,Gv, ! most probable speed
     &  MOLNEW(NNWML),
     &  IC(2,NCt),IR(MNM), ! vector for indexation      
     &  NM,NNPCout, ! current number of particles
     &  JJJ,Nwt,nbuf,nnew

