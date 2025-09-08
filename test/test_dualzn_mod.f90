module test_dualzn_mod
  use testdrive, only: error_type, unittest_type, new_unittest, check
  use dualzn_mod
  use precision_mod, only: dp => prec
  implicit none
  private

  public :: collect_dualzn_mod

contains
  subroutine collect_dualzn_mod(testsuite)
    type(unittest_type), allocatable, intent(out) :: testsuite(:)
    testsuite = [new_unittest("add", test_add)]
  end subroutine collect_dualzn_mod

  subroutine test_add(error)
    type(error_type), allocatable, intent(out) :: error

    integer, parameter :: order = 2
    type(dualzn) :: x, y

    call set_order(order)
    call initialize_dualzn(x)
    call initialize_dualzn(y)

    x%f(0) = 2.0_dp
    x%f(1) = 1.0_dp

    ! Add two dual numbers
    y = x + x
    call check(error, y%f(0), (4.0_dp, 0.0_dp))
    call check(error, y%f(1), (2.0_dp, 0.0_dp))
    call check(error, y%f(2), (0.0_dp, 0.0_dp))

  end subroutine test_add

end module test_dualzn_mod
