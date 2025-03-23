C  radiometric force on dics
C  weighting factor is used
C
      INCLUDE "p.inc"
      include "mpif.h"
      call MPI_Init(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, kp, ierr)      
      call MPI_Comm_rank(MPI_COMM_WORLD, mr, ierr)
      call dsmc(mr)
      call MPI_Finalize(ierr)
      end
c
      subroutine dsmc(mr)
      INCLUDE 'p.inc'
      Nwt=0
      jj=0
      CALL INIT(mr) ! initial distibution 
      DO JJJ=1,100 000 000
        TIME=TIME+DTM
        jj=jj+1
        CALL MOVE(mr) ! free motion
        CALL INDE(mr) ! indexation
        CALL COLL  ! collisions 
        CALL SAMPLE(mr) ! samples 
        if(jj==jp) then
          jj=0 
          CALL OUT(mr)        
          if(mr==np/2-1.or.mr==np/2)print99,mr,JJJ,time,cd,cq,nm
        end if
        if(abs(time-tstead)<1.d-5.and.Nwt==1) call init(mr)
        if(time>tstat.and.Nwt==2) go to 1
      END DO
1     CALL OUT(mr)
99    FORMAT(I3,1X,I10,1x,F8.0,2(1x,d19.12),1x,i5)
      END
C
      SUBROUTINE INIT(mr)
      INCLUDE 'p.inc'
      if(Nwt==0) then
      PI=DACOS(-1.D0)
      SPI=SQRT(PI)
      au=1.66053886d-27 ! atomic unit in kg
      akb=1.380649d-23 ! Boltzmann const
      go to (201,202,203,204,205,206) igas
201     aMM= 3.0160293d0*au ! molecular mass of He3
        Gv=400.d0
        open(20,file='MATRIX/muH3.dat') ! file of viscosities
        open(12,file='MATRIX/xiHe3.dat') ! file of matrix chi
        go to 207
202     aMM=4.0026d0*au ! molecular mass of He4
        Gv=400.d0
        open(20,file='MATRIX/muH4.dat') ! file of viscosities
        open(12,file='MATRIX/xiHe4.dat') ! file of matrix chi
        go to 207
203     aMM= 20.1797d0*au ! molecular mass of Ne
       Gv=200.d0
       open(20,file='MATRIX/muNe.dat') ! file of viscosities
       open(12,file='MATRIX/xiNe.dat') ! file of matrix chi
        go to 207      
204     aMM=39.948d0*au ! molecular mass of Ar
       Gv=150.d0
       open(20,file='MATRIX/muAr.dat') ! file of viscosities
       open(12,file='MATRIX/xiAr.dat') ! file of matrix chi
       go to 207
205     aMM=83.798d0*au ! molecular mass of Kr
       Gv=100.d0
       open(20,file='MATRIX/muKr.dat') ! file of viscosities
       open(12,file='MATRIX/xiKr.dat') ! file of matrix chi
       go to 207
206     aMM=131.293d0*au ! molecular mass of Xe
       Gv=80.d0
       open(20,file='MATRIX/muXe.dat') ! file of viscosities
       open(12,file='MATRIX/xiXeC.dat') ! file of matrix chi 
207    do j=1,300
        read(20,*,end=138) Tmu,amu
        if (abs((Tmu-T0)/T0)<1.d-10) go to 137 !read the viscosity
      end do
138   print*,'wrong viscosity'
137   amu=amu*1.d-6
      close(20)
      vm=sqrt(2.d0*akb*T0/aMM) !most probable speed in m/s
      do j=1,900
        read(12,*) (Xchi(i,j),i=1,100),TCS(j)  !read cos(chi) and TCS
      end do
      close(12)
      CW=1.D0/NC ! cell size      
      do irw=1,nrw
        rw(irw)=CW*irw
        ww(irw)=(rw(irw)+CW/2.d0)/(rw(irw)-CW/2.d0) ! weight of cell
      end do
      ww(2)=ww(3)
      ww(1)=ww(2)
      PWW=1.d0
      do irw=nrw,1,-1
       do jrw=1,irw
        PWW(jrw)=PWW(jrw)/WW(irw)
       end do
      end do
      VREG=PI*RR**2*SR ! volume of comput. region
      do irw=1,nrw ! correction of comput. region
        WWW=ww(irw)-1.d0
        do irj=irw+1,nrw
          WWW=WWW*ww(irj)
        end do
        vreg=vreg+pi*rw(irw)**2*WWW*SR
      end do
      DENM=INM/VREG ! density of model molecules in equilibrium   
      DO M=1,NCt
        CCG(2,M)=rf() ! remainder       
        CCG(1,M)=4. ! C  the maximum value of relative velocity
      END DO
      MC=1
      do i=1,NCx
        x1=(i-1)*Cw-SC+mr*SR
        x2=i*Cw-SC+mr*SR
        do j=1,NCy      
         y1=(j-1)*Cw
         y2=j*Cw
         VC(MC)=PI*Cw**3*(2*j-1) ! cell volume
         do irw=1,nrw
          if(j*Cw-1.d-5<rw(irw)) VC(MC)=VC(MC)*WW(irw)
         end do
         MC=MC+1
        end do
      end do
      CC=del*AMU*DTM*1.d-20/(aMM*vm*VC*DENM) ! factor for collision number
C  number of molecules that enter at each time step
      DO Nbc=1,NCX ! upper boundary
        AMEU(Nbc)=SPI*DENM*CW*RR*DTM
        AMRU(Nbc)=rf()
      END DO
      DO Nbc=1,NCY ! left and right boundaries
        AMEL(:,Nbc)=DENM*PI*Cw**2*(2*Nbc-1)*DTM/(2.*SPI)
        do irw=1,nrw
          if((Nbc-0.5)*CW<RW(irw)) then 
            AMEL(:,Nbc)=AMEL(:,Nbc)*WW(irw)
          end if
        end do
        AMRL(1,Nbc)=rf()
        AMRL(2,Nbc)=rf()          
      END DO
c initial generation
      RNPC=rf()
      MC=1
      nm=1
      do ix=1,NCx
      do iy=1,NCy
        y1=(iy-1)*Cw
        y2=iy*Cw
        ANPC=INM*VC(MC)/Vreg+RNPC
        NNPC=ANPC
        RNPC=ANPC-NNPC
        if(ix==1.and.iy==1.and.mr==0) then
         print*,'NNPC',NNPC
        end if
        do i=1,NNPC
          PP(1,nm)=Cw*(ix-rf())-SC+mr*SR
          PP(2,nm)=sqrt(y1**2+rf()*(y2**2-y1**2))
 103      VX=8.*RF()-4. ! generation of velocity
          IF(EXP(-VX*VX)<RF()) GO TO 103
          PV(1,nm)=VX
          call RVELC(PV(2,nm),PV(3,nm),1.d0)
          nm=nm+1
        end do
        MC=MC+1
      end do
      end do
      nm=nm-1       
      end if     
      CS=0.d0
      CSC=0.d0
      CSS=0.d0
      TIME=0.d0
      Nwt=Nwt+1      
      END
C
      SUBROUTINE MOVE(mr) ! molecules are moved over the time interval DTM
      INCLUDE 'p.inc'
      include 'mpif.h'
      integer status(MPI_status_size)
      dimension  PSL(INMRM,5),PSR(INMRM,5),PRL(INMRM,5),PRR(INMRM,5)
      Nnew=0
      IFT=-1
      N=0
100   N=N+1 ! N counter of molecules
      IF(N<=NM) THEN
        AT=DTM
        IF(IFT>0) AT=rf()*DTM !time step is random fraction of DTM for entering molecules
        XI=PP(1,N) ! old coordinates
        YI=PP(2,N)
        iy=YI/Cw+1
        pw=pww(iy)
C  calculation of new coordinates
101     X= XI+PV(1,N)*AT ! new coordinates
        yd=YI+PV(2,N)*AT
        dz=PV(3,N)*AT
        Y=SQRT(yd**2+dz**2)
        if(xi*x<0.d0) then ! particle crossed the disc plane
          times=-xi/pv(1,n)
          if(times<0..or.times>AT) then
            print*,'time',time,mr
          end if
          yds=yi+pv(2,n)*times
          zs=pv(3,n)*times
          ys=sqrt(yds**2+zs**2)
          if(ys<1.d0.and.rf()>beta) then ! collision with disc happens 
           nts=ys*NS+1 
           CSS(1,nts)=CSS(1,nts)-pw*pv(1,n)
           CSS(2,nts)=CSS(2,nts)-
     &       (pv(1,n)**2+pv(2,n)**2+pv(3,n)**2)*pw
           Tw=Th
           if(xi<0.d0) Tw=Tc
           A=-alpn*Tw*LOG(RF())
           B=2.d0*Pi*RF()
           pv(1,n)=sqrt(A+2.d0*sqrt(A*(1.-alpn))*pv(1,n)*cos(B)+
     &       (1.-alpn)*pv(1,n)**2)
           if(xi<0.d0) pv(1,n)=-pv(1,n)
           A=sqrt(-alpt*(2.D0-alpt)*Tw*dlog(rf()))
           B=2.d0*Pi*rf()
           pv(2,n)=A*cos(B)+pv(2,n)*(1.D0-alpt)
           pv(3,n)=A*sin(B)+pv(3,n)*(1.D0-alpt)
           CSS(1,nts)=CSS(1,nts)+pw*pv(1,n)
           CSS(2,nts)=CSS(2,nts)+
     &       (pv(1,n)**2+pv(2,n)**2+pv(3,n)**2)*pw
           AT=AT-times
           x=pv(1,n)*AT
           yd=yd+pv(2,n)*AT
           dz=dz+pv(3,n)*AT
           y=sqrt(yd**2+dz**2)
          end if
        end if
        iyn=y/cw+1
        Prem=pw/pww(iyn)
        if(y>RR.or.x>SRt-SC.or.x<-SC.or.rf()>Prem) then ! remove particle
          PP(:,N)=PP(:,NM)
          PV(:,N)=PV(:,NM)
          NM=NM-1
          N=N-1
          go to 100
        end if
        pp(1,n)=x
        pp(2,n)=y
C rotation of velocity
        VY=PV(2,N)
        VZ=PV(3,N)
        PV(2,N)=(VY*yd+VZ*dz)/y
        PV(3,N)=(VZ*yd-VY*dz)/y
        Nmultp=int(Prem) ! number of particle to generate + 1
        if(Prem-Nmultp>rf()) Nmultp=Nmultp+1 
        if(Nmultp>1) then
          Nnewo=Nnew
          Nnew=Nnew+Nmultp-1
          do inew=1,Nmultp-1
            MOLNEW(Nnewo+inew)=N
          end do
          if(Nnew>NNWML) then
            print*,'Nnew>NNWML' 
            stop
          end if        
        end if
        GO TO 100
      ELSE
        IF(IFT<0) THEN
          IFT=1
          CALL ENTER(mr)
          N=N-1
          GO TO 100
        END IF
      END IF
      if(jjj==1) go to 103
102   do inew=1,nbuf
        nm=nm+1
        pp(:,nm)=ppbuf(:,inew) 
        pv(:,nm)=pvbuf(:,inew)        
      end do  
      if(jjj==1) go to 104
103   nbuf=nnew
      do inew=1,nbuf        
        n=molnew(inew)
        ppbuf(:,inew)=pp(:,n)
        pvbuf(:,inew)=pv(:,n)        
      end do
      if(jjj==1) go to 102      
104   N=0
      NPSL=0
      NPSR=0
      NPRL=0
      NPRR=0
150   N=N+1
      if(N<=NM) then
        if(PP(1,N)<-SC+SR*mr) then ! particles goes down
         NPSL=NPSL+1 ! number of particles going down
         if(NPSL>INMRM) then
          print*,'NPSD>INMRM'
          stop
         end if
         PSL(NPSL,1)=PP(1,N)
         PSL(NPSL,2)=PP(2,N)          
         PSL(NPSL,3)=PV(1,N)
         PSL(NPSL,4)=PV(2,N)
         PSL(NPSL,5)=PV(3,N)
         PP(:,N)=PP(:,NM)
         PV(:,N)=PV(:,NM)
         N=N-1
         NM=NM-1
         go to 150
        end if
        if(PP(1,N)>-SC+SR*(mr+1)) then ! particles goes up
         NPSR=NPSR+1 ! number of particles going up
         if(NPSR>INMRM) then
          print*,'NPSU>INMRM'
          stop
         end if
         PSR(NPSR,1)=PP(1,N)
         PSR(NPSR,2)=PP(2,N)
         PSR(NPSR,3)=PV(1,N)
         PSR(NPSR,4)=PV(2,N)
         PSR(NPSR,5)=PV(3,N)        
         PP(:,N)=PP(:,NM)
         PV(:,N)=PV(:,NM)
         N=N-1
         NM=NM-1
        end if
        go to 150 
      end if
      if (mr+1<np) then ! send information to upper domain
       call MPI_SEND(NPSR,1,MPI_integer,mr+1,1,MPI_COMM_WORLD,ierr)
       call MPI_SEND
     &   (PSR,5*INMRM,MPI_double_precision,mr+1,2,MPI_COMM_WORLD,ierr)
       call MPI_RECV(NPRR,1,MPI_integer,mr+1,3
     &   ,MPI_COMM_WORLD, status, ierr)
       call MPI_RECV(PRR,5*INMRM,MPI_double_precision,mr+1,
     &     4,MPI_COMM_WORLD, status, ierr)     
      end if 
      if (mr>0) then  ! send information to lower domain
       call MPI_RECV
     &     (NPRL,1,MPI_integer,mr-1,1,MPI_COMM_WORLD,status,ierr)
       call MPI_RECV(PRL,5*INMRM,MPI_double_precision,mr-1,2,
     &     MPI_COMM_WORLD,status,ierr)
       call MPI_SEND(NPSL,1, MPI_integer, mr-1,
     &      3, MPI_COMM_WORLD, ierr)
       call MPI_SEND(PSL,5*INMRM, MPI_double_precision, mr-1,
     &      4, MPI_COMM_WORLD, ierr)
      end if 
      do  n=1,NPRL
        PP(1,NM+N)=PRL(N,1)
        PP(2,NM+N)=PRL(N,2)        
        PV(1,NM+N)=PRL(n,3)
        PV(2,NM+N)=PRL(n,4)
        PV(3,NM+N)=PRL(n,5)
      end do
       NM=NM+NPRL
       do  n=1,NPRR
        PP(1,NM+N)=PRR(n,1)
        PP(2,NM+N)=PRR(n,2)
        PV(1,NM+N)=PRR(n,3)
        PV(2,NM+N)=PRR(n,4)
        PV(3,NM+N)=PRR(n,5)
      end do
      NM=NM+NPRR
      if(nm > mnm) then
        print*,'nm>mnm'
       stop
      end if
      END 
c
      SUBROUTINE ENTER(mr)
C  new molecules enter at boundaries
      INCLUDE 'p.inc'
      DO Nbc=1,NCx ! lateral boundary
        A=AMEU(Nbc)+AMRU(Nbc)
        M=A
        AMRU(Nbc)=A-M !M molecules enter, remainder has been reset
        DO K=1,M
          IF(NM>=MNM) THEN   
            WRITE (*,*) ' EXCESS MOLECULE LIMIT'
            STOP
          END IF
          NM=NM+1 ! number of new molecule
          PP(1,NM)=-SC+(Nbc-rf())*CW+SR*mr
          PP(2,NM)=RR-0.0001*CW
          PV(2,NM)=-SQRT(-LOG(rf()))
          CALL RVELC(PV(1,NM),PV(3,NM),1.D0) ! generates two normal velocity components
          PV(1,NM)=PV(1,NM)
        end do
      END DO
      if(mr==0.or.mr==np-1) then ! lower and upper boundary
      Nb=1
      if(mr==np-1) Nb=2
       DO Nbc=1,NCy
        A=AMEL(Nb,Nbc)+AMRL(Nb,Nbc)
        M=A
        AMRL(Nb,Nbc)=A-M !M molecules enter, remainder has been reset
        DO K=1,M
          IF(NM>=MNM) THEN   
            WRITE (*,*) ' EXCESS MOLECULE LIMIT'
            STOP
          END IF
          NM=NM+1 ! number of new molecule
          YMIN=CW*(Nbc-1)
          YMAX=CW*Nbc
          PP(2,NM)=SQRT(YMIN**2+rf()*(YMAX**2-YMIN**2))
          IF(Nb==1) THEN
            PP(1,NM)=-SC+0.0001*CW            
            CX=sqrt(-dlog(rf()))
          ELSE
            PP(1,NM)=SRt-SC-0.0001*CW                        
            CX=-sqrt(-dlog(rf()))
          END IF
          PV(1,NM)=CX
          CALL RVELC(PV(2,NM),PV(3,NM),1.D0)
        END DO
       end do
      end if
      END
c
      SUBROUTINE COLL !calculates collisions appropriate to DTM in a gas
      INCLUDE "p.inc"
      DIMENSION VRCP(3),VCCM(3)
      DO N=1,NCt
        ASEL=CCG(1,N)*IC(2,N)*(IC(2,N)-1)*CC(N)+CCG(2,N)
        NSEL=ASEL
        CCG(2,N)=ASEL-NSEL
        DO ISEL=1,NSEL
c      selects a potential collision pair and calculates the product of the M
          K=RF()*IC(2,N)+IC(1,N)+1
          L=IR(K) !first molecule L has been chosen at random in cell
100       K=RF()*IC(2,N)+IC(1,N)+1
          M=IR(K) !second molecule M is chosen at random that are in the sub-cel
          IF(L==M) GO TO 100 !choose a new second molecule if the first is again
          VR=SQRT((PV(1,L)-PV(1,M))**2+(PV(2,L)-PV(2,M))**2+
     &            (PV(3,L)-PV(3,M))**2)
          jE=dlog(1.d0+VR*vm/Gv)/dlog(1.005d0)+0.5d0 !1<jE<900
          if (jE>900) jE=900
          if (jE<1) jE=1
          if (VR*TCS(jE)>CCG(1,N)) CCG(1,N)=VR*TCS(jE) !the maximum CCG(1,N) is
          IF (((TCS(jE)*VR)/CCG(1,N))>RF()) then
            ib=100*RF()+1
            VCCM(:)=0.5d0*(PV(:,L)+PV(:,M)) !centre of mass velocity
            VRx=PV(1,L)-PV(1,M) ! Pre-colisional relative speeds.
            VRy=PV(2,L)-PV(2,M)
            VRz=PV(3,L)-PV(3,M)
            VRr=sqrt(VRy*VRy+VRz*VRz)
            EPS=2.d0*pi*RF()
            cX=Xchi(ib,jE) !cos(chi)
            sX=sqrt(1.d0-cX*cX)
            if (VRr>1.d-15) then !Calculation of post collision velocities
              VRCP(1)=VRx*cX+sX*sin(EPS)*VRr
              VRCP(2)=VRy*cX+sX*(VR*VRz*cos(EPS)-VRx*VRy*sin(EPS))/VRr
              VRCP(3)=VRz*cX-sX*(VR*VRy*cos(EPS)+VRx*VRz*sin(EPS))/VRr
            else
              VRCP(1)=VRx*cX
              VRCP(2)=VRx*sX*cos(EPS)
              VRCP(3)=VRx*sX*sin(EPS)
            end if
            PV(:,L)=VCCM(:)+0.5d0*VRCP(:)
            PV(:,M)=VCCM(:)-0.5d0*VRCP(:)
          END IF
        END DO
      END DO
      END
c
      SUBROUTINE RVELC(U,V,VMP) ! genration of velocity, cosin law
      INCLUDE 'p.inc'
      A=SQRT(-LOG(rf()))
      B=2.d0*pi*rf()
      U=A*SIN(B)*VMP
      V=A*COS(B)*VMP
      END
C
      SUBROUTINE INDE(mr) 
C  the NM molecule numbers are arranged in order of cells and,
C  within the cells, in order of the sub-cells
      INCLUDE 'p.inc'
      IC(2,:)=0
      DO N=1,NM
        MC1=int((pp(1,n)+SC-mr*SR)/Cw)
        if(abs(pp(1,n)+SC-mr*SR)<1.d-10) 
     &     MC1=int((pp(1,n)+1.d-10+SC-mr*SR)/Cw) 
        if(abs(pp(1,n)+SC-(mr+1)*SR)<1.d-10) 
     &     MC1=int((pp(1,n)-1.d-10+SC-mr*SR)/Cw) 
        MC=MC1*NCy+int(pp(2,n)/Cw)+1
        if(mc<1.or.mc>NCt) then
          print*,'MC',mr,MC,NCt,pp(1,n),pp(2,n)
          stop
        end if
        IC(2,MC)=IC(2,MC)+1
      END DO !number of molecules in the cells and sub-cells have been counted
      M=0
      DO N=1,NCt
        IC(1,N)=M
        M=M+IC(2,N)
        IC(2,N)=0
      END DO !start address has been set for the cells
      DO N=1,NM
        MC=int((pp(1,n)+SC-mr*SR)/Cw)*NCy+int(pp(2,n)/Cw)+1
        IC(2,MC)=IC(2,MC)+1
        K=IC(1,MC)+IC(2,MC)
        IR(K)=N
C  the molecule number N has been set in the cross-reference array
      END DO
      END
c
      SUBROUTINE SAMPLE(mr)!  sample the molecules in the flow.
      INCLUDE 'p.inc'
      do i=1,NM  
        ix=(PP(1,i)+SC-mr*SR)/Cw+1 ! Cartesian coordinates
        if(ix==0) ix=1
        if(ix==Ncx+1) ix=NCx
        iy=PP(2,i)/Cw+1
        if(ix<0.or.ix>NCx+1) then
           print*,'ix',ix,NCx,pp(1,i)
           stop
        end if
        if(iy<1.or.iy>NCy) then
           print*,'iy',iy,NCy
           stop
        end if
        CSC(1,ix,iy)=CSC(1,ix,iy)+1.D0
        CSC(2,ix,iy)=CSC(2,ix,iy)+PV(1,i)
        CSC(3,ix,iy)=CSC(3,ix,iy)+PV(2,i)          
        CSC(4,ix,iy)=CSC(4,ix,iy)+PV(1,i)**2+PV(2,i)**2+PV(3,i)**2
      end do
      END
C
      SUBROUTINE OUT(mr)! output a progressive set of results
      INCLUDE 'p.inc'
      character*13 fn
1002  format('F00',i1,'.dat')
1003  format('F0' ,i2,'.dat')
1008  format('F'  ,i3,'.dat')
      if(mr<10) write(fn,1002) mr
      if(mr>=10.and.mr<100)write(fn,1003) mr
      if(mr>=100) write(fn,1008) mr      
      OPEN (11,FILE=fn)
      Uxc=0.
      Uyc=0.
      Tempc=0.
      MC=1
      do ix=1,NCx
      do iy=1,NCy
        AnC(ix,iy)=CSC(1,ix,iy)*DTM/(TIME*DENM*VC(MC)) ! density in n0
        MC=MC+1
        if(CSC(1,ix,iy)>0.) then
         UxC(ix,iy)=CSC(2,ix,iy)/CSC(1,ix,iy) ! velocity components in vm
         UyC(ix,iy)=CSC(3,ix,iy)/CSC(1,ix,iy)
         TempC(ix,iy)=2./3.*(CSC(4,ix,iy)/CSC(1,ix,iy)- ! temperature in T0
     &    (UXC(ix,iy)**2+UYC(ix,iy)**2))
        end if
        x=(ix-0.5)*SR/NCx-SC+mr*SR
        y=(iy-0.5)*RR/NCy
        if(iy==1) y=0
        WRITE (11,9904) x,y,AnC(ix,iy),
     &   UxC(ix,iy),UyC(ix,iy),TempC(ix,iy),ic(2,mc-1)
      end do
      write(11,*)      
      end do
      close(11)
9904  FORMAT(2(2X,F10.4),4(2X,E12.5),1x,i10,1x,e12.5)
      if(mr==np/2-1.or.mr==np/2) then
       if(mr==np/2-1) open(4,file='R1.dat')
       if(mr==np/2)   open(4,file='R2.dat')       
       CD=0.d0
       CQ=0.d0
       do its=1,NS
         CD=CD+CSS(1,its)
         CQ=CQ+CSS(2,its)
       end do
       CD=2.*CD/(time*DENM*pi*(Th-Tc))
       CQ=CQ/(time*DENM*pi*(Th-Tc))       
       WRITE(4,9902) DT,del,DTM,NC,SRt,RR,SC,time,CD,CQ
       DO its=1,NS
        y0=float(its-1)/NS
        y1=float(its)/NS
        CCS=PI*(y1**2-y0**2)
        CP=2.*CSS(1,its)/(TIME*DENM*CCS)
        CH=   CSS(2,its)/(TIME*DENM*CCS)
        write(4,9903) (its-0.5)/NS,CP,CH
       END DO
       close(4)
      end if
9902  FORMAT('#',1X,F5.2,2X,F6.3,3X,F6.4,1X,I4,
     &       3(1x,F5.1),3x,F7.0,2(3X,F10.7))
9903  FORMAT(1x,F10.5,2x,6(2X,E12.5)) 
      END
C
      DOUBLE PRECISION FUNCTION RF() ! random number generator
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      data IX1,IX2,IX3,IX4,IX5,IX6/1500419,1400159,1364,1528,1,3/
      RR1=1.0/FLOAT(IX1)
      RR2=1.0/FLOAT(IX2)
      IX5=MOD(IX5*IX3,IX1)
      IX6=MOD(IX6*IX4,IX2)
      RF=RR1*IX5+RR2*IX6
      IF(RF.GE.1.0)RF=RF-1.0
      END
      
