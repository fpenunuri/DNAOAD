!ifx -I../LibDualzn64 -o e2 ejemplo2.f90 -L../LibDualzn64 -ldualzn
program main
  use precision_mod
  use dualzn_mod
  implicit none

  complex(prec) :: z0
  type(dualzn) :: r, fval
  integer :: k
  real :: t1,t2

  call set_order(15) !we set the order to work with

  !since a dualzn numbers is an allocatable entity, do not forget to
  !initialize it
  r = 0 !<--- initializing r to 0 
  r%f(0) = (1.1_prec,0.0_prec)
  r%f(1) = 1 !since we want to differentiate, r = r0 +1*eps_1
  !all the other components are 0 as r was initialized to 0  

  call cpu_time(t1)
  fval = ftest(r)
  call cpu_time(t2)
  
  !Computing the derivatives, from the 0th derivative up to the
  !order-th derivative.
  print*,"derivatives"
  do k=0,order
     write(*,"(i0,a,f0.1,a,e17.10)") k,"-th derivative at x = ", &
          real(r%f(0)),":",real(fval%f(k))
  end do

  print*,"elapsed time (s):",t2-t1  

contains
  function ftest(x) result(fr)
    type(dualzn), intent(in) :: x
    type(dualzn) :: fr
    integer :: k

    !nested function f(x) = sin(x) * exp(-x^2), f(f(...(f(x))...))
    !applied 1000 times
    fr = sin(x)*exp(-x*x)
    do k=1, 1000-1
       fr = sin(fr)*exp(-fr*fr)
    end do
  end function ftest
end program main


