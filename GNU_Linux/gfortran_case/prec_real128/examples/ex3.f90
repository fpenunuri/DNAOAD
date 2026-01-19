!gfortran -I../LibDualzn128 -o ExDO ex3.f90 -L../LibDualzn128 -ldualzn

!module with example of functions
module function_mod
  use dualzn_mod
  implicit none
  private

  public :: fstest, fvectest

contains
  !Example of scalar function f = sin(x*y*z) + cos(x*y*z)
  function fstest(r) result(fr)
    type(dualzn), intent(in), dimension(:) :: r
    type(dualzn) :: fr
    type(dualzn) :: x,y,z

    x = r(1); y = r(2); z = r(3)
    fr = sin(x*y*z) + cos(x*y*z)
  end function fstest

  !Example of vector function f = [f1,f2,f3]
  !f = fvectest(r) is a function f:D^m ---> Dn 
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
  complex(prec), parameter :: ii = (0,1)
  complex(prec), dimension(mq) :: q, vec
  complex(prec), dimension(nf,mq) :: Jmat
  complex(prec), dimension(nf) :: JV
  complex(prec), dimension(3) :: GV
  complex(prec), dimension(3,3) :: Hmat
  integer :: i
  
  vec = [1.0_prec,2.0_prec,3.0_prec,4.0_prec]
  q = vec/10.0_prec + ii

  print*,"Jv using matmul"
  Jmat = Jacobian(fvectest, q , nf)
  JV = matmul(Jmat,vec)  
  do i=1,nf
     write(*,*) JV(i)
  end do
  write(*,*)
  
  print*,"Jv using vector directional derivative"
  JV = d1fvector(fvectest,vec,q,nf)
  do i=1,nf
     write(*,*) JV(i)
  end do
  write(*,*)
  
  print*,"---Hessian matrix---"
  Hmat = Hessian(fstest,q(1:3))
  do i=1,3
     write(*,"(A,i0)") "row:",i
     write(*,*) Hmat(i,:)
  end do
  write(*,*)
  
  print*,"---Gradient---"
  GV = gradient(fstest,q(1:3))
  do i=1,3
     write(*,*) GV(i)
  end do   
end program main
