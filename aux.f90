!
subroutine allocation
!
  use myconstants
  use mymatrices
  implicit none
  integer :: info,infot
  
  infot=0
  allocate (Kpoint(1:Ncx),stat=info)
  infot=infot+info
  allocate (theta(1:ntheta),stat=info)
  infot=infot+info
  allocate (wkpoint(1:Ncx),stat=info)
  infot=infot+info
  allocate (wtheta(1:ntheta),stat=info)
  infot=infot+info
  allocate (w(1:nw),stat=info)   
  infot=infot+info
  allocate (ep(1:nep),stat=info)   
  infot=infot+info
  allocate (Sigmaf(1:nw,2,2),stat=info) 
  infot=infot+info
  allocate (Greenf(1:nw,2),stat=info)
  infot=infot+info
  allocate (green_inv0f(1:nw,2),stat=info)
  infot=infot+info
  allocate (Gef(1:nw,2),stat=info)
  infot=infot+info
  allocate (Gif(1:nw,2),stat=info)
  infot=infot+info 
  allocate (dos0(1:nep),stat=info)
  infot=infot+info 

  if (infot.ne.0) then 
     write(*,*) 'Allocation request denied'
     stop
  endif
  
  return  
end subroutine allocation

       
subroutine deallocate
!    This routine deallocates all the arrays
  use myconstants
  use mymatrices
  implicit none
  
  deallocate(Kpoint) 
  deallocate(theta)
  deallocate(wkpoint)
  deallocate(wtheta)
  deallocate(w)
  deallocate(Sigmaf)    
  deallocate(Greenf)  
  deallocate(green_inv0f)   
  deallocate(Gef) 
  deallocate(Gif) 
  
  return
end subroutine deallocate

!************************************************************
subroutine get_dos0
  use myconstants
  use mymatrices

  implicit none
  integer is,n
  real*8 r1
  character*80 dosinitfor,dosinit

  Ncx = Ncxint(1)
  write(dosinitfor,617)int(log10(Ncx*1.d0))+1,int(log10(nep*1.d0))+1
  write(dosinit,dosinitfor)Ncx,nep
  open(unit=14, file=dosinit,status='old',action='read',iostat=is)
!  open(unit=14, file='p_rho.bin',status='old',action='read',iostat=is)
  if(is.ne.0) then
    write(*,*)'file doesnt exist'
     stop
  end if

  r1=0.d0
  do n=1,nep
     read(14,*)ep(n),dos0(n)
     r1=r1+dos0(n)
  end do
  d_ep=ep(2)-ep(1)
  write(*,*)'totdos0=',r1*d_ep
  close(14)

617 format(14h('dosin.Ncx',ii1,8h'.nep',ii1,1h))
  return
end subroutine get_dos0

!************************************************************
subroutine dos_init
  use myconstants
  use mymatrices

  implicit none
  integer i,j,k,n
  real*8 eng,weight,r1,width
  character*80 dosinitfor,dosinit
  complex*16 z1

  if(Vc.lt.0.d0) then
     write(dosinitfor,417)
  else
     write(dosinitfor,617)
  end if
  write(dosinit,dosinitfor)xi,nf,Jc,Vc,phase

  Ncx = Ncxint(1)
  call intcosdisp(kpoint,wkpoint,Nf,Ncx)

!============binning method
!  dos0=0.d0
!  do i=1,Ncx
!     do j=1,Ncx
!        do k=1,Ncx
!           eng=dcos(kpoint(i))+dcos(kpoint(j))+ dcos(kpoint(k))
!!           eng=6.d0-2.d0*eng
!           eng=-2.d0*eng
!           weight=wkpoint(i)*wkpoint(j)*wkpoint(k)
!  
!           n=int((eng-wmin)/d_ep)+1
!           if(n .le. nep) then
!              dos0(n)=dos0(n)+weight
!           else
!              write(*,*)'energy out of range'
!              stop
!           endif
!        end do   ! j
!     end do   ! k
!  end do   ! i
!  dos0(:)=dos0(:)/real(d_ep*pi**3)
!  r1=sum(dos0(:))*d_ep


!======================broading method
  dos0=0.d0;r1=0.d0; width=1.5d0
  do n=1,nep
     if(abs( ep(n) ) .le. width) then
        z1=zero
        do i=1,Ncx
           do j=1,Ncx
              do k=1,Ncx
                 eng=dcos(kpoint(i))+dcos(kpoint(j))+ dcos(kpoint(k))
                 !              eng=6.d0-2.d0*eng
                 eng=-0.5d0*eng
                 weight=wkpoint(i)*wkpoint(j)*wkpoint(k)
                 z1=z1+weight/(ep(n)-eng+ii*5.d-4) 
              end do   ! j
           end do   ! k
        end do   ! i
        dos0(n)=-dimag(z1); r1=r1+dos0(n)
     end if
  end do
  dos0=dos0/pi**4; r1=r1*d_ep/pi**4

!=====================semicircular dos
!  dos0=0.d0; r1=0.d0; width=1.5d0
!  do n=1,nep
!     if(abs( ep(n) ) .le. width) then
!        dos0(n)=2.d0/(pi*width)*sqrt( 1.d0 - (ep(n)/width)**2 )
!        r1=r1+dos0(n)
!     end if
!  end do
!ri=r1*d_ep

  write(6,*) 'Norm of dos0=',r1
  open(unit=17,file=dosinit,status='unknown') 
  do n=1,nep
     write(17,*)ep(n),dos0(n)
  end do
  close(17)

217 format(3(e15.9,1x))
417 format(60h('dosin.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))
617 format(60h('dosin.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))
  return
end subroutine dos_init
!***************************************************************************

!***************************************************************************
subroutine adjust_mue(efer)
  use myconstants
  use mymatrices
  implicit none
  integer n
  real*8 aux,aux2,chp0,chp1,filling,fill0
  real*8 efer,falpos
  external falpos,filling

!    Adjust chemical potential with the updated DOS
  aux=dzero
  do n=1,nw-1
     aux2=aux
     aux=aux-0.5d0*diffw*(dimag(Gif(n,1))+dimag(Gif(n+1,1)))-       &
          &          0.5d0*diffw*(dimag(Gif(n,2))+dimag(Gif(n+1,2)))
     if (aux.gt.(pi*Nf)) then
        efer=w(n)+(pi*Nf-aux2)*diffw/(aux-aux2)
        exit
     endif
  end do
  chp0=efer; fill0=filling(chp0)  
  if(fill0.lt.dzero) then
     do
        chp1=chp0+5.0d0*diffw
        if(filling(chp1).gt.dzero) exit
        chp0=chp1
     enddo
  else
     do
        chp1=chp0
        chp0=chp0-5.0*diffw
        if(filling(chp0).lt.dzero) exit
     end do
  endif

  !  Solve for the chemical potential
  mue=falpos(filling,chp0,chp1,acurfill)

  return
end subroutine adjust_mue
!***************************************************************************

!***************************************************************************
subroutine intgpoints(kpoint,wkpoint,Ni,Ncx)
  implicit none
  integer Ncx,Ncx8,l
  real*8 Ni,kpoint(ncx),wkpoint(ncx),aux,pi
  
  pi = 2.d0*dasin(1.d0)
  Ncx8=Ncx/8
  if (dabs((dble(Ncx)/8.d0)-Ncx8) > 0.d0) then
     write(*,*) 'Ncx must be multiple of 8' 
     stop
  else
     aux=(3.d0*(pi**2)*Ni)**(1.d0/3.d0)
     if (aux == 0.d0) then
        call gauleg(-pi,pi,kpoint,wkpoint,8*ncx8)
     else
        call gauleg(-pi,-aux,kpoint,wkpoint,ncx8)
        call gauleg(-aux,aux,kpoint(ncx8+1),wkpoint(ncx8+1),6*ncx8)
        do l=1,ncx8
           kpoint(7*ncx8+l)=-kpoint(ncx8-l+1)
           wkpoint(7*ncx8+l)=wkpoint(ncx8-l+1)
        enddo
     end if
  endif
  return
end subroutine intgpoints
!***************************************************************************


!***************************************************************************
subroutine intcosdisp(kpoint,wkpoint,Ni,Ncx)
  implicit none
  integer Ncx,Ncxhalf
  real*8 kpoint(ncx),wkpoint(ncx),Ni,aux,pi
  ! 
  pi = 2.d0*dasin(1.d0)
  Ncxhalf=Ncx/2
  if (dabs((dble(Ncx)/2.d0)-Ncxhalf) > 0.d0) then
     write(*,*) 'Ncx must be multiple of 2' 
     stop
  endif
!   aux=Fermi momentum of the free electron gas 
  aux=(3.d0*pi*pi*ni)**(1.d0/3.d0)
  if (aux == 0.d0) then
     call gauleg(0.d0,pi,kpoint,wkpoint,Ncx)
!     call gauleg(-pi,pi,kpoint,wkpoint,Ncx)
  else
     call gauleg(0.d0,aux,kpoint,wkpoint,ncxhalf)
     call gauleg(aux,pi,kpoint(ncxhalf+1),wkpoint(ncxhalf+1),ncxhalf)
!     call gauleg(-pi,aux,kpoint,wkpoint,ncxhalf)
!     call gauleg(aux,pi,kpoint(ncxhalf+1),wkpoint(ncxhalf+1),ncxhalf)
  end if
  return
end subroutine intcosdisp
!***************************************************************************


SUBROUTINE gauleg(x1,x2,x,weight,n)
  implicit none
  INTEGER n
  REAL*8 x1,x2,x(n),weight(n)
  DOUBLE PRECISION EPS
  PARAMETER (EPS=3.d-14)
  INTEGER i,j,m
  DOUBLE PRECISION p1,p2,p3,pp,xl,xm,z,z1
  m=(n+1)/2
  xm=0.5d0*(x2+x1)
  xl=0.5d0*(x2-x1)
  do 12 i=1,m
     z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
1    continue
     p1=1.d0
     p2=0.d0
     do 11 j=1,n
        p3=p2
        p2=p1
        p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
11      continue
        pp=n*(z*p1-p2)/(z*z-1.d0)
        z1=z
        z=z1-p1/pp
        if(abs(z-z1).gt.EPS)goto 1
        x(i)=xm-xl*z
        x(n+1-i)=xm+xl*z
        weight(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
        weight(n+1-i)=weight(i)
12      continue
        return
     END
!
!
     real*8 function filling(chpot)
       use myconstants
       use mymatrices

       implicit none
       integer n,isp
       real*8 aux,fermi1,fermi2,chpot,eng

       aux=dzero
       do n=1,nw-1
          eng=(w(n)-chpot)
          if (eng.gt.dzero) then
             fermi1=0.d0
          else
             fermi1=1.d0
          endif
          eng=(w(n+1)-chpot)
          if (eng.gt.dzero) then
             fermi2=0.d0
          else
             fermi2=1.d0
          endif
          do isp=1,2
             aux=aux-(fermi1*dimag(Gif(n,isp))+fermi2*dimag(Gif(n+1,isp)))
          enddo
       end do
       filling=(aux*0.5d0*diffw/pi)-nf
       return 
     end function filling

      real*8 FUNCTION falpos(func,x1,x2,xacc)
      INTEGER*4 MAXIT
      REAL*8 x1,x2,xacc,func
      EXTERNAL func
      PARAMETER (MAXIT=5000)
      INTEGER*4 j
      REAL*8 del,dx,f,fh,fl,swap,xh,xl

      fl=func(x1)
      fh=func(x2)
!      if(fl*fh.gt.0.) pause 'falpos: root must be bracketed'
      if(fl*fh.gt.0.) then 
         write(6,*) 'falpos: root must be bracketed'
         write(6,*) 'fl',fl,x1,'fh',fh,x2
         stop
      endif
      if(fl.lt.0.)then
        xl=x1
        xh=x2
      else
        xl=x2
        xh=x1
        swap=fl
        fl=fh
        fh=swap
      endif

      dx=xh-xl
      do 11 j=1,MAXIT
        falpos =xl+dx*fl/(fl-fh)
        f=func(falpos)
        if(f.lt.0.) then
          del=xl-falpos
          xl=falpos
          fl=f
        else
          del=xh-falpos
          xh=falpos
          fh=f
        endif
        dx=xh-xl
        if(abs(del).lt.xacc.or.f.eq.0.) return
11    continue

      write(6,*) 'falpos: exceeding maximum iterations'
      write(6,*) falpos,'falpos'
      END
