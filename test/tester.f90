program tester
  use, intrinsic :: iso_fortran_env, only: error_unit
  use testdrive, only: run_testsuite
  use test_dualzn_mod, only: collect_dualzn_mod
  implicit none

  integer :: stat
  stat = 0

  call run_testsuite(collect_dualzn_mod, error_unit, stat)

  if (stat > 0) then
    write(error_unit, '(i0, 1x, "test(s) failed!")') stat
  end if

end program tester
