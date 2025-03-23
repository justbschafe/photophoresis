      INCLUDE 'p.inc'
      dimension CP(Ns),Cf(Ns),Ch(Ns)
      character*13 fn
      pi=acos(-1.d0)
      open(200,file='F.dat')
      do mr=0,np-1 ! merge Files
       if(mr<10) write(fn,1002) mr
       if(mr>=10.and.mr<100) write(fn,1003) mr
       if(mr>=100) write(fn,1008) mr       
       OPEN (100,FILE=fn)
       do ix=1,NCx
        do iy=1,NCy
         read (100,9904) x,y,a,b,c,t,i
         write(200,9904) x,y,a,b,c,t,i
        end do
        read(100,*) 
        write(200,*)       
       end do
       close(100)
      end do
      close(200)
c
      ch=0.
      Cd=0.
      OPEN (101,FILE='R1.dat')
      OPEN (102,FILE='R2.dat')      
      read(101,99012) time,CD1,CQ1
      read(102,99012) time,CD2,CQ2
99012 FORMAT(49x,F7.0,2(3x,F10.7))
      CD=CD1+CD2
      open(601,file='d.lst')
      read(601,9907) NPC
      open(600,file='Res.dat')
      WRITE(600,9901)
      WRITE(600,9902) del,beta,DTM,NPC,NC,SRt,RR,time
      WRITE(600,9908) CD,CQ1,CQ2      
      write(600,9905)
       do its=1,Ns
        read(101,9903) rs,cp1,ch1
        read(102,9903) rs,cp2,ch2
        write(600,9906) rs,(-1+beta-cp1)/(Th-Tc),(cp2-1+beta)/(Th-Tc),
     &    ch1/(Th-Tc),ch2/(Th-Tc)
      END DO
      close(600)
9901  FORMAT('#',2X,'del',3x,'beta',6x,'DTM',3X,'NPC',
     &       2X,'NC',3X,'SR',4X,'RR',2x,'time')
9902  FORMAT('#',1X,F5.1,2x,f4.2,3X,F6.4,2X,I3,1X,I3,
     &       2(1x,F5.1),1X,F7.0)
9908  format('#',70(1h-)/'#  F=',F10.7,1x,',  Qt=',F10.7,1x,
     &     ',  Qb=',F10.7)
9903  FORMAT(1x,F10.5,2x,4(2X,E12.5)) 
9904  FORMAT(2(2X,F10.4),4(2X,E12.5),1x,i10,1x,e12.5)
9905  format('#',70(1h-)/'# local pressure and heat flux'/
     & '#',4x,'r',8x,'p top',9x,'p bot',9x,'h top',
     & 9x,'h bot')    
9906  FORMAT(1x,F7.3,2x,4(2X,E12.5)) 
9907  format(14x,i3)
1002  format('F00',i1,'.dat')
1003  format('F0' ,i2,'.dat')
1008  format('F' ,i3,'.dat')
      end
      
