!gfortran -I../LibDualzn64 -o ExJ ExJacobian.f90 -L../LibDualzn64 -ldualzn

!module with example of functions
module function_mod
  use dualzn_mod
  implicit none
  private

  public :: fstest, fvectest

contains
  
  function fstest(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r

    type(dualzn) :: fr
    type(dualzn) :: x,y,z

    x = r(1); y = r(2); z = r(3)

    fr = sin(x*y*z) + cos(x*y*z)
  end function fstest


  !f = fvectest(r) is a function f:D^m ---> Dn 
  !
  function fvectest(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r
    type(dualzn), allocatable, dimension(:) :: fr

    type(dualzn) :: f1,f2,f3
    type(dualzn) :: x,y,z,w

    x = r(1); y = r(2); z = r(3); w = r(4)

    f1 = sin(x*y*z*w)
    f2 = cos(x*y*z*w)*sqrt(w/y - x/z)
    f3 = sin(log(x*y*z*w))

    allocate(fr(3))
    fr = [f1,f2,f3]
  end function fvectest

end module function_mod


!-----------------------------------------------------------------------
!main program
program main
  use precision_mod
  use dualzn_mod
  use diff_mod
  use function_mod
  implicit none

  integer, parameter :: nf =3, mq = 4
  complex(prec), dimension(mq) :: q, vec
  complex(prec), dimension(nf,mq) :: Jmat
  complex(prec), parameter  :: ii = (0,1)
  integer :: i

  complex(prec), dimension(nf) :: JV

  complex(prec) :: aux
  complex(prec), dimension(3) :: qR3,GV
  real(prec), dimension(3) :: eiR3
  complex(prec), dimension(3,3) :: Hmat


  vec = [1.0_prec,2.0_prec,3.0_prec,4.0_prec]
  q = vec/10.0_prec + ii

  Jmat = Jacobian(fvectest, q , nf)
  do i=1,nf
     write(*,*) Jmat(i,:)
     print*,'.....'
  end do

  print*,'====xv===='
  JV = matmul(Jmat,vec)
  write(*,*) JV
  print*,'.....'
  JV = d1fvector(fvectest,vec,q,nf)- matmul(Jmat,vec)
  write(*,*) JV
  print*,'====xv===='
  
  
  eiR3 = [1.0_prec,0.0_prec,0.0_prec]
  qR3 = q(1:3)

  aux = d2fscalar(fstest,cmplx(eiR3,kind=prec),cmplx(eiR3,kind=prec),qR3)
  write(*,*) aux

  aux = d2fscalar(fstest,cmplx(eiR3,kind=prec),qR3)
  write(*,*) aux


  print*,'====xxx===='
  Hmat = Hessian(fstest,qR3)
  do i=1,3
     write(*,*) Hmat(i,:)
     print*,'.....'
  end do
  
  print*,'====xxxx===='
  print*,d1fscalar(fstest,cmplx([1,0,0],kind=prec),qR3)
  print*,d1fscalar(fstest,cmplx([0,1,0],kind=prec),qR3)
  print*,d1fscalar(fstest,cmplx([0,0,1],kind=prec),qR3)
  print*,'.......'
  GV = gradient(fstest,qR3)
  do i=1,3
     write(*,*) GV(i)
  end do  
  print*,'====xxxx===='
  

end program main
