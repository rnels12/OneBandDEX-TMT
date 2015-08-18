subroutine scf

  use myconstants
  use mymatrices
  implicit none
  integer n,m,iter,jj,i,j,k,l,irealconv,iter_max
  real*8 aux,aux2,aux3,weight,costheta,zeroplus,toten
  real*8 sqbt2,Ni,rmstotal,rms,epsilon,efermi,Gif_r(nw,2)
  real*8 dos_av(nw,2),dos_typ(nw,2)
  real*8 eng,engren,rmstotmax,Vcc
  real*8 mauxr(2),mauxi(2),maux2r(2),maux2i(2)
  complex*16 maux(2),maux2(2),determ,determ2,caux
  character*80 sigfor,dosfor,errfor,dostypfor
  character*80 Sigmafile,DOSfile,errfile,dostypfile
  parameter(epsilon=1.d-2,zeroplus=1.d-6)
!        

  iter_max=200

  Vcc = 0.d0
  if(Vc.lt.0.d0) then
     write(sigfor,500)
     write(dosfor,510)
     write(dostypfor,511)
     write(errfor,512)
  else
     write(sigfor,600)
     write(dosfor,610)
     write(dostypfor,611)
     write(errfor,612)
  end if
  write(Sigmafile,sigfor) xi,nf,Jc,Vc,phase
  write(DOSfile,dosfor) xi,nf,Jc,Vc,phase
  write(dostypfile,dostypfor) xi,nf,Jc,Vc,phase
  write(errfile,errfor) xi,nf,Jc,Vc,phase
  open(unit=11,file=Sigmafile,status='unknown',position='append') 
  open(unit=12,file=DOSfile,status='unknown',position='append') 
  open(unit=13,file=dostypfile,status='unknown',position='append') 
  open(unit=33,file=errfile,status='unknown',position='append') 
  open(unit=31,file='Gi.txt',status='unknown') 
  write(10,140)mue

  write(*,*)'scf=',sum(dos0(:))*d_ep

!    Convergence loop:If (irealconv=1) calculation converged 
  irealconv=0
  do 20 iter = 1,iter_max	! max number of iterations = 1000
     Greenf=zero; Gif=zero
     do n=1,nw
        !coarsed grained GF: G(wn)=int_{ep} d_ep * dos(ep)/[w-ep-Sigma(w)]   
        maux(:)=w(n)-Sigmaf(n,1,:) 
        mauxr=dreal(maux)
        mauxi=dimag(maux)
        do m=1,nep
           aux=ep(1)+(real(m)-0.5d0)*d_ep
           aux2=ep(1)+(real(m)+0.5d0)*d_ep
           maux2r(:)=-0.5d0 * dlog( ( (mauxr(:)-aux2)**2 + mauxi(:)**2 )/ &
                ( (mauxr(:)-aux)**2 + mauxi(:)**2 ) )
           maux2i(:)= datan( (aux-mauxr(:))/dabs(mauxi(:)) ) - &
                datan( (aux2-mauxr(:))/dabs(mauxi(:)) )
           Greenf(n,:)=Greenf(n,:) + dos0(m)*(maux2r(:) + ii*maux2i(:))
        end do !m

        !mean field inversed:[G~0]^-1=[G]^-1 + Sigma
        green_inv0f(n,:) = (one/Greenf(n,:)) + Sigmaf(n,1,:)
     end do   !n

     ! G~i(w):
     if(phase .eq. 0) then
        do 45 j=1,ntheta    !Loop over angle theta
           weight=0.5d0*wtheta(j)*dsin(theta(j))
           costheta=0.5d0*Jc*S*dcos(theta(j)) 
           sqbt2=0.d0 !0.25d0*Jc*Jc*S*S*(done-dcos(theta(j))*dcos(theta(j)))           
           do 55 n = 1,nw
              maux(1) = green_inv0f(n,1) - Vc + costheta
              maux(2) = green_inv0f(n,2) - Vc - costheta
              determ=maux(1)*maux(2)-sqbt2
              Gef(n,1)=maux(2)/determ
              Gef(n,2)=maux(1)/determ
              !the spin orientation average:
              Gif(n,:) = Gif(n,:) + weight*Gef(n,:)
55         end do  !n
45      end do  !j=1,ntheta
     elseif(phase .eq. 1) then
        costheta=0.5d0*Jc*S
        do 56 n = 1,nw
           maux(1) = green_inv0f(n,1) - Vc + costheta
           maux(2) = green_inv0f(n,2) - Vc - costheta
           Gif(n,:)=one/maux(:)
56      end do  !n
     else
        costheta=0.5d0*Jc*S;maux=zero
        do n=1,nw
           do jj=1,2
              maux(1) = green_inv0f(n,1) - Vc + (-done)**jj * costheta
              maux(2) = green_inv0f(n,2) - Vc - (-done)**jj * costheta
              Gif(n,:)=Gif(n,:) + one/maux(:)
           end do !jj
              Gif(n,:)=Gif(n,:)*0.5d0
        end do !n
     end if

     do n=1,nw
        maux2(:)= green_inv0f(n,:) - Vcc
        maux(:) = one/maux2(:)
        mauxr(:)=-dimag(Gif(n,:))/pi;maux2r(:)=-dimag(maux(:))/pi
        do jj=1,2
           if(mauxr(jj).lt.0.d0)  mauxr(jj) = mauxr(jj) + epsilon !zeroplus
           if(maux2r(jj).lt.0.d0) maux2r(jj)= maux2r(jj)+ epsilon !zeroplus
        end do !jj
        mauxi(:)=xi*log(mauxr(:)) + (1.d0-xi)*log(maux2r(:))
        dos_typ(n,:)=exp(mauxi(:))
     end do

!     if(irealconv .eq. 1 .or. iter .eq. iter_max) then
!        aux3=0.d0
        do n=1,nw
           maux2(:)= green_inv0f(n,:) - Vcc
           maux(:) = one/maux2(:)
           mauxr(:)=-dimag(Gif(n,:))/pi;maux2r(:)=-dimag(maux(:))/pi
           dos_av(n,:)=xi*mauxr(:)+(1.d0-xi)*maux2r(:)
!           do jj=1,2
!              if(dos_av(n,jj).lt.0.d0) dos_av(n,jj)=dos_av(n,jj)+ zeroplus
!           end do !jj
!              aux3=aux3+dos_av(n,1)*diffw
        end do !n
     if(irealconv .eq. 1 .or. iter .eq. iter_max) then
        write(*,*)iter,sum(dos_av(:,1))*diffw,sum(dos_av(:,2))*diffw
        exit
     end if

     do jj=1,2
!        call hilbert_w_logcorr(Gif(:,jj),Gif_r(:,jj),dos_typ(:,jj))
        call hilbert(Gif(:,jj),Gif_r(:,jj),dos_typ(:,jj))
     end do

     aux=sum( abs(Gif(:,:)) )
     if(aux .ne. aux) then
        write(*,*)'green function NaN, please check parameters';stop
     endif

     do n=1,nw   ! average on configurations
        Sigmaf(n,2,:)= green_inv0f(n,:) - (one/Gif(n,:))
        do jj=1,2
           if(abs( dimag( Sigmaf(n,2,jj) ) ) .lt. epsilon) &
                Sigmaf(n,2,jj)= dreal(Sigmaf(n,2,jj))-ii*epsilon
        end do
     end do  ! n


!******************************Test-box
!     if(iter .eq. 80) then
!        do n=1,nw
!           write(31,'(10(f10.4,1X))')w(n),( dreal(Gif(n,j)),dimag(Gif(n,j)), &
!                dos_typ(n,j), & 
!                dreal(sigmaf(n,2,j)),dimag(sigmaf(n,2,j)),j=1,2)
!           if(n .eq. nw) stop
!        end do
!     end if
!***********************************************************************

!     call adjust_mue(efermi)

!       Check difference between new iteration and previous one
     rmstotmax=dzero!;jj=-1
!     do n=1,nw
        aux=dzero
        jj=int(nw/4)
        do i=1,2
!           rms=abs(Sigmaf(n,1,i) - Sigmaf(n,2,i))
           rms=abs(Sigmaf(jj,1,i) - Sigmaf(jj,2,i))
           aux=aux+rms
        enddo
        if (aux.gt.rmstotmax) then
           rmstotmax=aux!;jj=n
        endif
!     enddo !n
     write(33,*)iter,rmstotmax,Tolrns,jj;flush(33)
     if (rmstotmax .lt. Tolrns) irealconv=1
     
     !Calculate the DOS at the edge 
     eng=-0.5d0*diffw*dimag(Gif(nw,1)+Gif(nw,2))
     aux=-0.5d0*diffw*dimag(Gif(1,1)+Gif(1,2))      
     write(10,215) iter,mue,efermi,eng,aux,rmstotmax; call flush(10)

     !Linear mixing: Sigma_new= damp*Sigma_old + (1-damp)*Sigma_out
     do n=1,nw
        Sigmaf(n,1,:)= damp*Sigmaf(n,1,:)+(done-damp)*Sigmaf(n,2,:)
     enddo


     if (irealconv .eq. 1)  then
        write(10,*) 'Real freq. code converged after ',iter, ' iterations!' 
!        exit
!        goto 125
     endif
20 end do  ! iter: End of convergence loop

write(10,230)  mue,nf
!write(11,*) mue,Ni
do n=1,nw
   !Storage real selfenergies
   write(11,*) w(n),(dreal(Sigmaf(n,1,i)), dimag(Sigmaf(n,1,i)),i=1,2)

   !Storage DOS_average
   write(12,100) w(n),dos_av(n,1),-dos_av(n,2)

   !Storage DOS_typical
   write(13,100) w(n),dos_typ(n,1),-dos_typ(n,2)
end do !n

!storage total-energy
toten=0.d0
aux = 0.5d0*(-dimag(Gif(1,1))/pi - dimag(Gif(1,2))/pi)*w(1)
do n=2,nw
   if(w(n).le.mue) then
      aux2 = 0.5d0*(-dimag(Gif(n,1))/pi - dimag(Gif(n,2))/pi)*w(n)
      toten=toten+aux+aux2
      aux=aux2
   end if
end do
write(10,145)toten*diffw

close(11);close(12);close(13)

        
100 format(3(e15.9,1x))
140 format('mue=',e15.9)
145 format('total energy=',e15.6)
215 format(i5,' mue=',e15.9,' efermi=',e15.9,' DOS_edge=',e15.9,  &
         &         1x,e15.9,' error=',e15.9)
230 format(e15.9,1x,e15.9,'  mu')
500 format(60h('Sigma.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))
510 format(60h('DOSAV.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))
511 format(60h('DOSty.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))
512 format(60h('error.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f5.2,'.phase_',i1))

600 format(60h('Sigma.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))
610 format(60h('DOSAV.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))
611 format(60h('DOSty.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))
612 format(60h('error.x',f4.2,'.nf',f5.3,'.J',f4.2,'.V',f4.2,'.phase_',i1))

return
end subroutine scf


subroutine hilbert_w_logcorr(Gt,rlsigma,rhosigma)    	    	
  use myconstants
  use mymatrices
  
  implicit none
  integer :: i,j,n
  integer :: info
  real*8 :: rhosigma(nw),rlsigma(nw)
  real*8 :: woj,dspj,resgj,spi,dspi,r1,dwi,woji1,woji2
  real*8, allocatable  :: rlsg(:)
  complex*16 :: Gt(nw)
!	************************************************************

  allocate(rlsg(0:nw+1),stat=info)
  
  do j=2,nw
     woj=.5D0*(w(j-1)+w(j))
     dspj=rhosigma(j)-rhosigma(j-1)
     resgj=0.0D0
     
     do i=1,j-2
        spi=rhosigma(i)
        dspi=rhosigma(i+1)-rhosigma(i)
        dwi=w(i+1)-w(i)
        woji1=w(i)-woj
        woji2=w(i+1)-woj
        r1=dlog(woji2/woji1)

        resgj=resgj-(spi*r1 + dspi )
        resgj=resgj-(dspi/dwi)*(woj -w(i))*r1
     end do

        !	 skip the interval (j-1) to j

     do i=j,nw-1               
        spi=rhosigma(i)
        dspi=rhosigma(i+1)-rhosigma(i)
        dwi=w(i+1)-w(i)
        woji1=w(i)-woj
        woji2=w(i+1)-woj
        r1=dlog(woji2/woji1)
        
        resgj=resgj-(spi*r1 + dspi )
        resgj=resgj-(dspi/dwi)*(woj -w(i))*r1
     end do

     resgj=resgj - dspj
     rlsg(j)=resgj
  end do


  rlsg(nw+1)=rlsg(nw)
  rlsg(1)=rlsg(2)
    
  do i=1,nw
     rlsigma(i)=0.5d0*(rlsg(i)+rlsg(i+1))
  end do

!c	Now load the real and imaginary part into one
  Gt = zero
  do n=1,nw
     Gt(n)= rlsigma(n) - ii*pi*(rhosigma(n))
  end do
  
  deallocate(rlsg)
  return
end subroutine Hilbert_w_logcorr




subroutine hilbert(Gt,rlsigma,rhosigma)    	    	
  use myconstants
  use mymatrices
  
  implicit none
  integer :: i,j,n
  integer :: info
  real*8 :: rhosigma(nw),rlsigma(nw)
  real*8 :: woj,dspj,resgj,spi,dspi,r1,dwi,woji1,woji2
  complex*16 :: Gt(nw)
!	************************************************************
  
  rlsigma=0.d0
  do j=1,nw   
     do i=1,j-1 
        rlsigma(j)=rlsigma(j)+rhosigma(i)*diffw/(w(j)-w(i))
     end do

     do i=j+1,nw
        rlsigma(j)=rlsigma(j)+rhosigma(i)*diffw/(w(j)-w(i))
     end do
  end do

!c	Now load the real and imaginary part into one
  Gt = zero
  do n=1,nw
     Gt(n)= rlsigma(n) - ii*pi*(rhosigma(n))
  end do

  return
end subroutine hilbert
