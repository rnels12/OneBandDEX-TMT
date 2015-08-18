      module myconstants
      integer  Ncx,nw,nep,ntheta
      integer  nncxint,Ncxint(10),no,phase
      real*8 pi,dzero,done,nf,wnint(10),mue,wmax,wmin,diffw,d_ep
      real*8 errorfilling,S,xi,Jc,Vc,H,DH,damp,Tolrns,delta,acurfill
      complex*16 ii,zero,one
      end module myconstants
! 
      module mymatrices
      complex*16, allocatable :: Greenf(:,:),Gef(:,:),Gif(:,:)
      complex*16, allocatable :: Sigmaf(:,:,:),green_inv0f(:,:)
      real*8, allocatable :: kpoint(:),theta(:),wkpoint(:),wtheta(:)  
      real*8, allocatable :: w(:),ep(:),dos0(:)
      end module mymatrices     
!     






