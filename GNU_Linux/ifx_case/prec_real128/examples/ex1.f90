!ifx -I../LibDualzn128 -o e1 ex1.f90 -L../LibDualzn128 -ldualzn
program main
  use precision_mod
  use dualzn_mod
  implicit none

  complex(prec) :: z0
  type(dualzn) :: r, fval
  integer :: k

  call set_order(5) !before anything we set the 'order' to work with
  
  z0 = (1.1_prec,2.2_prec)
  r = 0 !we initialize the dual number to 0

  !we set the 0-th and 1-th components. If dual numbers are used to
  !calculate D^n f(z0) then r must be of the form r = r0 + 1*eps_1 
  r%f(0) = z0
  r%f(1) = 1 

  fval = sin(r)**log(r*r)

  !writing the derivatives, from the 0th derivative up to the
  !order-th derivative.
  print*,"derivatives"
  do k=0,order
     print*,fval%f(k)
  end do

end program main

  
  
  
