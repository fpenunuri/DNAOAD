module precision_mod
  use iso_fortran_env, only: int32, int64, real32, real64, real128
  implicit none

  !integer, parameter :: prec = real64
  integer, parameter :: prec = real128
end module precision_mod
