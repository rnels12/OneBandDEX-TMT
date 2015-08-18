!
        program one_band
! 
!     DMF code for one-band model based on a code 
!     previously developed by Juana Moreno 
!     to calculate finite temperature properties. 
!     This code calculates average & typical DOS at T=0 given 
!     doping, filling and coupling Jc.
!     Store DOS, chemical potential, and
!     Sigma(w) for all temperatures.
!     Energy in units of the hoping "t".
!     developed by R. Nelson

        use myconstants
        use mymatrices
!
        implicit none
!
         character*80 checkfile,Sigmafile0,checkfor
	 integer i,k,n,q
	 real*8 im,re
         real*8 prevmue(0:15),aux1,aux2
         real*8 falpos
         external falpos
!        
        ii = (0.0d0,1.0d0)
        zero = (0.d0,0.d0)
        dzero=0.d0
        one = (1.0d0,0.d0)
        done=1.d0
        pi = 2.d0*dasin(1.d0)

        open(unit=10, file='in_dat')
        read(10,*) phase        !0=PM, otherwise FM
        read(10,*) Xi           !doping
        read(10,*) nf           !filling value
        read(10,*) Jc	   !coupling
        read(10,*) Vc	   !spin-dependent impurity potential
        read(10,*) nw           !# of real frequencies
        read(10,*) nep           !# of frequencies for bare-dos
        read(10,*) wmin,wmax    !energy interval
        read(10,*) S            !impurity Spin
        read(10,*) H	           !initial magnetic field
        read(10,*) DH	   !increments on magnetic field
!  For smaller real frequencies more k-points are used on Green_loc integral
        read(10,*) nncxint      !how many intervals
        read(10,*) (Ncxint(i),i=1,nncxint)	!# k-points in each int.
        read(10,*) ntheta	   !# points on the angle theta integr.
        read(10,*) Tolrns       !Convergence when "error"<  Tolrns
        read(10,*) delta        !delta= broadening of quasiparticle peak
        read(10,*) errorfilling !Maximum error in filling for 1st filling
        read(10,*) acurfill     !Accuracy on filling equation
        read(10,*) damp         !linear mixing factor
        read(10,*) no,Sigmafile0	   !if no.ne.0 read Sigma from Sigmafile0
        close(unit=10)

        Ncx = Ncxint(1) 
        call allocation

!   Set up a table for angles theta
        call gauleg(dzero,pi,theta,wtheta,ntheta)
!
        if(Vc.lt.0.d0) then
           write(checkfor,400)
        else
           write(checkfor,600)
        end if
        write(checkfile,checkfor) xi,nf,Jc,Vc,phase
        open(unit=10,file=checkfile,status='unknown',position='append') 

        write(10,200) Xi,Nf,Jc,Vc,nw
        write(10,210) wmin,wmax,S
        write(10,230) H,DH
        write(10,240) (Ncxint(i),i=1,nncxint),ntheta
        write(10,250) Tolrns,delta,errorfilling
        
!  Define the different frequency intervals for the k-integration
!  Use the bandwidth for COSINE DISPERSION
        do i=1,nncxint-1
           wnint(i)=(10**(i-1))*12.d0
        enddo
        write(10,270) (wnint(i),i=1,nncxint-1)

!    If "no".ne.0 -> read the selfenergy function: sigma
        if(no.ne.0) then
           open(unit=67,file=Sigmafile0,status='unknown')
           read(67,*) mue,aux1,aux1,aux1
           do n=1,nw
              do i = 1,2
                 read(67,*) q, re, im
                 Sigmaf(n,1,i) = dcmplx(re,im)   
              end do !i
           end do !n
           close(67) 
           H=dzero
        else    
!   mue=Fermi energy of the free electron gas COSINE DISPERSION
           mue=(3.d0*pi*pi*nf)**(2.d0/3.d0)
           Sigmaf=-ii*delta
        end if !no
!
!       Set up a table for real freq GF.
        do n = 1,nw
           w(n) = ((nw-n)*wmin+wmax*(n-1))/real(nw-1)
        end do
        diffw=w(2)-w(1)

!       Set up a table for freq od dos0.
!        do n = 1,nep
!           ep(n) = ((nep-n)*wmin+wmax*(n-1))/real(nep-1)
!        end do
!        d_ep=ep(2)-ep(1)
!        call dos_init

        call get_dos0
        call scf
        call deallocate
        
        close(10)

	stop
 100     format(e15.9,1x,e15.9)
 150     format(e15.9,1x,e15.9,1x,e15.9,1x,e15.9)
 400     format(60h('check.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))
 600     format(60h('check.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))
 200     format('# Xi=',e15.9,' Nf=',e15.9,' Jc=',e15.9,' Vc=',e15.9,' nw=',i6)
 210     format('# wmin=',e15.9,' wmax=',e15.9,' S=',e15.9)
 230     format('# H=',e15.9,' DH=',e15.9)
 240     format('# Ncx=',i4,1x,i4,1x,i4,1x,i4,1x,i4,' ntheta=',i6)
 250     format('# Tolrns=',e15.9,' delta=',e15.9,' errorfilling=',e15.9)
 270     format('#  wnint=',e15.9,1x,e15.9,1x,e15.9,1x,e15.9)

      end program one_band
