!Dual numbers to order n of complex components
!This can be used to compute first, second,... nth order derivatives
!F. Pe~nu~nuri
!UADY, Merida Yucatan Mexico
!2024
!

module dualzn_mod
  use precision_mod
  implicit none

  private
  !---------------------------------------------------------------------
  !Some Module variables
  !Default order, can be modified with set_order
  integer, public  :: order = 1
  real(prec), public, parameter :: Pi = 4.0_prec*atan(1.0_prec)
  !---------------------------------------------------------------------

  !dual number definition to any order
  type, public :: dualzn
     complex(prec), allocatable, dimension(:) :: f
  end type dualzn

  public :: set_order, initialize_dualzn, f_part
  public :: binomial, BellY, Dnd
  public :: itodn, realtodn, cmplxtodn, Mset_fpart

  public :: inv, sin, cos, tan, exp, log, sqrt, asin, acos, atan, asinh
  public :: acosh, atanh, sinh, cosh, tanh, absx, atan2
  public :: conjg, sum, product, matmul

  public :: assignment (=)
  public :: operator(*), operator(/), operator(+), operator(-)
  public :: operator(**), operator(==),  operator(/=)

  !equal assignment
  interface assignment (=)
     module procedure igualz16_dzn   !dualzn <--- complex16
     module procedure igualz8_dzn    !dualzn <--- complex8
     module procedure igualz4_dzn    !dualzn <--- complex4
     module procedure igualr16_dzn   !dualzn <--- real16
     module procedure igualr8_dzn    !dualzn <--- real8
     module procedure igualr4_dzn    !dualzn <--- real4
     module procedure iguali8_dzn    !dualzn <--- integer8
     module procedure iguali4_dzn    !dualzn <--- integer4
  end interface assignment (=)
  !---------------------------------------------------------------------

  !functions to change from integer, real and complex numbers to dual
  interface itodn
     module procedure i4todn
     module procedure i8todn
  end interface itodn

  interface realtodn
     module procedure r4todn
     module procedure r8todn
     module procedure r16todn
  end interface realtodn

  interface cmplxtodn
     module procedure c4todn
     module procedure c8todn
     module procedure c16todn
  end interface cmplxtodn
  !---------------------------------------------------------------------

  !Logical equal operator 
  interface operator (==) 
     module procedure eq_dzn
  end interface operator (==)
  !---------------------------------------------------------------------

  !Logical not equal operator 
  interface operator (/=) 
     module procedure noteq_dzn
  end interface operator (/=)
  !---------------------------------------------------------------------

  !overloaded operators
  interface operator(*)
     module procedure timesdX
     module procedure timesXd
  end interface operator(*)

  interface operator(/)
     module procedure divdX
     module procedure divXd
  end interface operator(/)

  interface operator (+)
     module procedure masd      !unary
     module procedure sumadX
     module procedure sumaXd
  end interface operator (+)

  interface operator (-)
     module procedure menosd    !unary
     module procedure restadX
     module procedure restaXd
  end interface operator (-)

  interface operator(**)
     module procedure powerdX
     module procedure powerXd
  end interface operator(**)
  !---------------------------------------------------------------------
  !some matrix and vector operations
  interface sum
     module procedure sumR2dzn   !sum(dual_rank2,dir)
     module procedure sumR20dzn  !sum(dual_rank2)
     module procedure sumR1dzn   !sum(dual_rank1)
  end interface sum
  !---------------------------------------------------------------------

  !product for dual 
  interface product
     module procedure prodR2dzn   !product(dual_rank2,dir)
     module procedure prodR20dzn  !product(dual_rank2)
     module procedure prodR1dzn   !product(dual_rank1)
  end interface product
  !---------------------------------------------------------------------

  interface matmul
     module procedure MtimesdX
     module procedure Mtimesc128d 
     module procedure Mtimesc64d
     module procedure Mtimesc32d
     module procedure Mtimesr128d 
     module procedure Mtimesr64d
     module procedure Mtimesr32d
     module procedure Mtimesi64d
     module procedure Mtimesi32d     
  end interface matmul
  !---------------------------------------------------------------------

  !overloaded functions
  interface conjg
     module procedure conjg_dzn
  end interface conjg

  interface atan2
     module procedure atan2d_
  end interface atan2

  interface tanh
     module procedure tanhd
  end interface tanh

  interface cosh
     module procedure coshd
  end interface cosh

  interface sinh
     module procedure sinhd
  end interface sinh

  interface atanh
     module procedure atanhd
  end interface atanh

  interface acosh
     module procedure acoshd
  end interface acosh

  interface asinh
     module procedure asinhd
  end interface asinh

  interface atan
     module procedure atand_
     ! Fortran 2008 standard for 2 argument atan:
     module procedure atan2d_
  end interface atan

  interface acos
     module procedure acosd_
  end interface acos

  interface asin
     module procedure asind_
  end interface asin

  interface log
     module procedure logd
  end interface log

  interface exp
     module procedure expd
  end interface exp

  interface sqrt
     module procedure sqrtd
  end interface sqrt

  interface sin
     module procedure sind_
  end interface sin

  interface cos
     module procedure cosd_
  end interface cos

  interface tan
     module procedure tand_
  end interface tan
  !---------------------------------------------------------------------

  !interface to define a function of complex variable returning a dual
  !variable
  !f:z --> dual
  abstract interface
     pure function funzdual(z_val) result(f_result)
       use precision_mod
       import :: dualzn
       complex(prec), intent(in) :: z_val
       type(dualzn) :: f_result
     end function funzdual
  end interface
  !---------------------------------------------------------------------

  !Functions
contains
  !to set the order for dual numbers
  subroutine set_order(new_order)
    integer, intent(in) :: new_order
    if(new_order<1) then
       stop "order must be greater than 0"
    end if
    order = new_order
  end subroutine set_order
  !---------------------------------------------------------------------

  !Assignment, equal operator
  !dualzn <--- complex16
  elemental subroutine igualz16_dzn(A, z)
    type(dualzn), intent(out) :: A
    complex(real128), intent(in) :: z

    A = cmplxtodn(z)
  end subroutine igualz16_dzn

  !dualzn <--- complex8
  elemental subroutine igualz8_dzn(A, z)
    type(dualzn), intent(out) :: A
    complex(real64), intent(in) :: z

    A = cmplxtodn(z)
  end subroutine igualz8_dzn

  !dualzn <--- complex4
  elemental subroutine igualz4_dzn(A, z)
    type(dualzn), intent(out) :: A
    complex(real32), intent(in) :: z

    A = cmplxtodn(z)
  end subroutine igualz4_dzn
  !---------------------------------------------------------------------

  !dualzn <--- real16
  elemental subroutine igualr16_dzn(A, x)
    type(dualzn), intent(out) :: A
    real(real128), intent(in) :: x

    A = realtodn(x)
  end subroutine igualr16_dzn

  !dualzn <--- real8
  elemental subroutine igualr8_dzn(A, x)
    type(dualzn), intent(out) :: A
    real(real64), intent(in) :: x

    A = realtodn(x)
  end subroutine igualr8_dzn

  !dualzn <--- real4
  elemental subroutine igualr4_dzn(A, x)
    type(dualzn), intent(out) :: A
    real(real32), intent(in) :: x

    A = realtodn(x)
  end subroutine igualr4_dzn
  !---------------------------------------------------------------------

  !dualzn <--- integer8
  elemental subroutine iguali8_dzn(A, ix)
    type(dualzn), intent(out) :: A
    integer(int64), intent(in) :: ix

    A = itodn(ix)
  end subroutine iguali8_dzn

  !dualzn <--- integer4
  elemental subroutine iguali4_dzn(A, ix)
    type(dualzn), intent(out) :: A
    integer(int32), intent(in) :: ix

    A = itodn(ix)
  end subroutine iguali4_dzn
  !---------------------------------------------------------------------

  !to initialize a dual number to zero
  elemental subroutine initialize_dualzn(zdn)
    type(dualzn), intent(out) :: zdn

    allocate(zdn%f(0:order)) !notice, it will be order+1 components
    zdn%f=(0.0_prec,0.0_prec)
  end subroutine initialize_dualzn
  !---------------------------------------------------------------------

  !to extract the f%(k) part of a dual number or array of duals
  elemental function f_part(x,k) result(fr)
    type(dualzn), intent(in) :: x
    integer, intent(in) :: k
    complex(prec) :: fr

    fr = x%f(k)
  end function f_part
  !---------------------------------------------------------------------

  !Logical equal operator
  elemental function eq_dzn(lhs, rhs) result(fr)
    type (dualzn), intent(in) :: lhs, rhs
    logical :: fr
    logical :: eqfk
    integer :: k

    fr = .true.
    do k=0,order
       eqfk = lhs%f(k) == rhs%f(k)
       if(.not. eqfk) then
          fr = .false.
          exit
       end if
    end do
  end function eq_dzn
  !---------------------------------------------------------------------
  !Logical not equal operator
  elemental function noteq_dzn(lhs, rhs) result(f_res)
    type (dualzn), intent(in) :: lhs, rhs
    logical :: f_res

    f_res = .not.(lhs == rhs)
  end function noteq_dzn
  !---------------------------------------------------------------------
  elemental function c16todn(x) result(fr)
    complex(real128), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=cmplx(x,kind=prec)
  end function c16todn

  elemental function c8todn(x) result(fr)
    complex(real64), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=cmplx(x,kind=prec)
  end function c8todn

  elemental function c4todn(x) result(fr)
    complex(real32), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=cmplx(x,kind=prec)
  end function c4todn
  !---------------------------------------------------------------------

  elemental function r16todn(x) result(fr)
    real(real128), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=real(x,kind=prec)
  end function r16todn

  elemental function r8todn(x) result(fr)
    real(real64), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=real(x,kind=prec)
  end function r8todn

  elemental function r4todn(x) result(fr)
    real(real32), intent(in) :: x
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=real(x,kind=prec)
  end function r4todn
  !---------------------------------------------------------------------

  elemental function i4todn(ix) result(fr)
    integer(int32), intent(in) :: ix
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=real(ix,kind=prec)
  end function i4todn

  elemental function i8todn(ix) result(fr)
    integer(int64), intent(in) :: ix
    type(dualzn) :: fr

    call initialize_dualzn(fr)
    fr%f(0)=real(ix,kind=prec)
  end function i8todn
  !---------------------------------------------------------------------

  !dual + class
  elemental function sumadX(B,X) result(fr)
    class(*), intent(in) :: X
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (X)
    type is (dualzn)
       Xd = X
    type is(complex(kind=real128))
       Xd = cmplxtodn(X)
    type is (complex(kind=real64))
       Xd =cmplxtodn(X)
    type is (complex(kind=real32))
       Xd =cmplxtodn(X)
    type is (real(kind=real128))
       Xd = realtodn(X)
    type is (real(kind=real64))
       Xd = realtodn(X)
    type is (real(kind=real32))
       Xd = realtodn(X)   
    type is (integer(kind=int64))
       Xd = itodn(X)
    type is (integer(kind=int32))
       Xd = itodn(X)
    end select
    fr = sumad(B,Xd)
  end function sumadX
  !---------------------------------------------------------------------

  !class + dual
  elemental function sumaXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (XX)
    type is (dualzn)
       Xd = XX
    type is(complex(kind=real128))
       Xd = cmplxtodn(XX)
    type is (complex(kind=real64))
       Xd =cmplxtodn(XX)
    type is (complex(kind=real32))
       Xd =cmplxtodn(XX)
    type is (real(kind=real128))
       Xd = realtodn(XX)
    type is (real(kind=real64))
       Xd = realtodn(XX)
    type is (real(kind=real32))
       Xd = realtodn(XX)
    type is (integer(kind=int64))
       Xd = itodn(XX)
    type is (integer(kind=int32))
       Xd = itodn(XX)
    end select
    fr = sumad(Xd, BB)
  end function sumaXd
  !---------------------------------------------------------------------

  !A+B
  elemental function sumad(A,B) result(fr)
    type(dualzn), intent(in) :: A,B
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order))
    do k=0,order
       fr%f(k) = A%f(k) + B%f(k)
    end do
  end function sumad
  !---------------------------------------------------------------------

  !+dualzn (unary)
  elemental function masd(A) result(fr)
    type(dualzn), intent(in) :: A
    type(dualzn) :: fr

    fr = A
  end function masd
  !---------------------------------------------------------------------

  !dual - class
  elemental function restadX(B,X) result(fr)
    class(*), intent(in) :: X
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (X)
    type is (dualzn)
       Xd = X
    type is(complex(kind=real128))
       Xd = cmplxtodn(X)
    type is (complex(kind=real64))
       Xd =cmplxtodn(X)
    type is (complex(kind=real32))
       Xd =cmplxtodn(X)
    type is (real(kind=real128))
       Xd = realtodn(X)
    type is (real(kind=real64))
       Xd = realtodn(X)
    type is (real(kind=real32))
       Xd = realtodn(X)       
    type is (integer(kind=int64))
       Xd = itodn(X)
    type is (integer(kind=int32))
       Xd = itodn(X)
    end select
    fr = restad(B,Xd)
  end function restadX
  !---------------------------------------------------------------------

  !class - dual
  elemental function restaXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (XX)
    type is (dualzn)
       Xd = XX
    type is(complex(kind=real128))
       Xd = cmplxtodn(XX)
    type is (complex(kind=real64))
       Xd =cmplxtodn(XX)
    type is (complex(kind=real32))
       Xd =cmplxtodn(XX)
    type is (real(kind=real128))
       Xd = realtodn(XX)
    type is (real(kind=real64))
       Xd = realtodn(XX)
    type is (real(kind=real32))
       Xd = realtodn(XX)   
    type is (integer(kind=int64))
       Xd = itodn(XX)
    type is (integer(kind=int32))
       Xd = itodn(XX)
    end select
    fr = restad(Xd, BB)
  end function restaXd
  !---------------------------------------------------------------------

  !-dualzn (unary)
  elemental function menosd(A) result(fr)
    type(dualzn), intent(in) :: A
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order))
    do k=0,order
       fr%f(k) = -A%f(k)
    end do
  end function menosd

  !A-B
  elemental function restad(A,B) result(fr)
    type(dualzn), intent(in) :: A,B
    type(dualzn) :: fr

    fr = -B+A
  end function restad
  !---------------------------------------------------------------------

  !dual**class
  elemental function powerdX(B,X) result(fr)
    class(*), intent(in) :: X
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr
    type(dualzn) :: Xd
    real(real64) :: Xr8

    select type (X)
    type is (dualzn)
       Xd = X
       fr = powerd(B,Xd)

    type is(complex(kind=real128))
       Xd = cmplxtodn(X)
       fr = powerd(B,Xd)

    type is (complex(kind=real64))
       Xd =cmplxtodn(X)
       fr = powerd(B,Xd)

    type is (complex(kind=real32))
       Xd =cmplxtodn(X)
       fr = powerd(B,Xd)

    type is (real(kind=real128))
       fr = power_dr16(B,X)

    type is (real(kind=real64))
       fr = power_dr8(B,X)

    type is (real(kind=real32))
       Xr8 = X
       fr = power_dr8(B,Xr8)

    type is (integer(kind=int64))
       fr = power_dint8(B,X)

    type is (integer(kind=int32))
       fr = power_dint(B,X)
    end select

  end function powerdX
  !---------------------------------------------------------------------

  !class**dual
  elemental function powerXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (XX)
    type is (dualzn)
       Xd = XX
    type is(complex(kind=real128))
       Xd = cmplxtodn(XX)
    type is (complex(kind=real64))
       Xd =cmplxtodn(XX)
    type is (complex(kind=real32))
       Xd =cmplxtodn(XX)
    type is (real(kind=real128))
       Xd = realtodn(XX)
    type is (real(kind=real64))
       Xd = realtodn(XX)
    type is (real(kind=real32))
       Xd = realtodn(XX)   
    type is (integer(kind=int64))
       Xd = itodn(XX)
    type is (integer(kind=int32))
       Xd = itodn(XX)
    end select
    fr = powerd(Xd, BB)
  end function powerXd
  !---------------------------------------------------------------------

  function MtimesdX(A,X) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A
    class(*), intent(in), dimension(:,:) :: X
    type(dualzn), dimension(size(A,1),size(X,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    select type (X)
    type is (dualzn)
       Xd = X
    type is(complex(kind=real128))
       Xd = X
    type is(complex(kind=real64))
       Xd = X
    type is(complex(kind=real32))
       Xd = X
    type is(real(kind=real128))
       Xd = X
    type is(real(kind=real64))
       Xd = X
    type is(real(kind=real32))
       Xd = X
    type is(integer(kind=int64))
       Xd = X
    type is(integer(kind=int32))
       Xd = X         
    end select
    fr = Mtimesd(A,Xd)
    !fr = Mtimesd_LbnzRule(A,Xd)
  end function MtimesdX
  !---------------------------------------------------------------------

  !cmplx*dual
  function Mtimesc128d(X,B) result(fr)
    complex(real128), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesc128d


  function Mtimesc64d(X,B) result(fr)
    complex(real64), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesc64d

  function Mtimesc32d(X,B) result(fr)
    complex(real32), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesc32d
  !---------------------------------------------------------------------

  !real*dual
  function Mtimesr128d(X,B) result(fr)
    real(real128), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesr128d

  function Mtimesr64d(X,B) result(fr)
    real(real64), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesr64d

  function Mtimesr32d(X,B) result(fr)
    real(real32), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesr32d
  !---------------------------------------------------------------------

  function Mtimesi64d(X,B) result(fr)
    integer(real64), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesi64d

  function Mtimesi32d(X,B) result(fr)
    integer(real32), intent(in), dimension(:,:) :: X
    type(dualzn), intent(in), dimension(:,:) :: B
    type(dualzn), dimension(size(X,1),size(B,2)) :: fr
    type(dualzn), dimension(size(X,1),size(X,2)) :: Xd

    Xd = X

    fr = Mtimesd(Xd,B)
    !fr = Mtimesd_LbnzRule(Xd,B)
  end function Mtimesi32d
  !---------------------------------------------------------------------

  !direct implementation of dual matrix multiplication
  function Mtimesd(A,B) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A, B
    type(dualzn), dimension(size(A,1),size(B,2)) :: fr
    integer :: i, j, k, m, n, p, q

    m = size(A,1)
    n = size(A,2)
    p = size(B,1)
    q = size(B,2)

    if (n /= p) then
       print*,"Error: Matrix dimensions do not align for multiplication"
       return
    end if

    call initialize_dualzn(fr)

    do i = 1, m
       do j = 1, q
          do k = 1, n
             fr(i,j) = fr(i,j) + A(i,k)*B(k,j)
          end do
       end do
    end do
  end function Mtimesd
  !...................

  !Dual matrix multiplication following the Leibniz rule.
  !This can be useful for implementing dual matrix multiplication 
  !using the BLAS or MKL libraries. (Not implemented in this module).
  function Mtimesd_LbnzRule(A,B) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A, B
    type(dualzn), dimension(size(A,1),size(B,2)) :: fr
    integer :: k
    complex(prec), dimension(size(A,1),size(B,2)) :: ABk

    call initialize_dualzn(fr)
    do k=0,order
       ABk = Mtimesdzn(A,B,k)
       call Mset_fpart(k,ABk,fr)
    end do
  end function Mtimesd_LbnzRule

  function Mtimesdzn(A,B,k) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: A, B
    integer, intent(in) :: k     
    complex(prec), dimension(size(A,1),size(B,2)) :: fr
    complex(prec), dimension(size(A,1),size(A,2)) :: Aizaux
    complex(prec), dimension(size(B,1),size(B,2)) :: Bkmizaux
    integer :: i

    fr=0.0_prec
    do i=0,k
       Aizaux = f_part(A,i)
       Bkmizaux = f_part(B,k-i)
       fr=fr+binomial(k,i)*matmul(Aizaux,Bkmizaux)
    end do
  end function Mtimesdzn

  !to set the dual-k component to cm in matrix A
  subroutine Mset_fpart(k,cm,A)
    integer, intent(in) :: k
    complex(prec),intent(in), dimension(:,:) :: cm
    type(dualzn), intent(inout), dimension(size(cm,1),size(cm,1)) :: A
    integer :: i,j

    do i=1,size(cm,1)
       do j=1,size(cm,2)
          A(i,j)%f(k) = cm(i,j)
       end do
    end do
  end subroutine Mset_fpart
  !---------------------------------------------------------------------

  !A**B
  elemental function powerd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr
    type(dualzn) :: intdual
    integer :: iaux
    logical :: BIQ

    iaux = int(real(B%f(0)))
    intdual = itodn(iaux)

    BIQ = intdual == B

    if(BIQ) then
       fr = A**iaux
    else
       fr = expd(B*logd(A))
    end if
  end function powerd

  !A**n (n integer)
  elemental function power_dint(A,n) result(fr)
    type(dualzn), intent(in) :: A
    integer, intent(in) :: n
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order))
    fr%f=(0.0_prec,0.0_prec)

    if(A==itodn(1)) then
       fr=A
       return
    elseif(A==itodn(0) .and. n>0)then
       fr=A
       return
    end if

    !0**0 ---> 1
    if(n==0) then
       fr%f(0) = (1.0_prec,0.0_prec)
    elseif(n>=1) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,n
          fr = fr*A
       end do
    elseif(n<0) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,-n
          fr = fr*A
       end do
       fr = inv(fr)
    end if
  end function power_dint
  !---------------------------------------------------------------------

  !A**n (n integer)
  elemental function power_dint8(A,n) result(fr)
    type(dualzn), intent(in) :: A
    integer(int64), intent(in) :: n
    type(dualzn) :: fr
    integer(int64) :: k

    allocate(fr%f(0:order))
    fr%f=(0.0_prec,0.0_prec)

    if(A==itodn(1)) then
       fr=A
       return
    elseif(A==itodn(0) .and. n>0)then
       fr=A
       return
    end if

    !0**0 ---> 1
    if(n==0) then
       fr%f(0) = (1.0_prec,0.0_prec)
    elseif(n>=1) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,n
          fr = fr*A
       end do
    elseif(n<0) then
       fr%f(0) = (1.0_prec,0.0_prec)
       do k=1,-n
          fr = fr*A
       end do
       fr = inv(fr)
    end if
  end function power_dint8
  !---------------------------------------------------------------------

  !A**x16
  elemental function power_dr16 (A,x16) result(fr)
    type(dualzn), intent(in) :: A
    real(real128), intent(in) :: x16
    type(dualzn) :: fr
    type(dualzn) :: x16dn
    integer :: iaux
    logical :: x16IQ

    iaux = int(x16)

    x16IQ = real(iaux,kind=real128) == x16

    if(x16IQ) then
       fr = A**iaux
    else
       x16dn = realtodn(x16)
       fr = A**x16dn
    end if
  end function power_dr16

  !A**x8
  elemental function power_dr8 (A,x8) result(fr)
    type(dualzn), intent(in) :: A
    real(real64), intent(in) :: x8
    type(dualzn) :: fr
    type(dualzn) :: x8dn
    integer :: iaux
    logical :: x8IQ

    iaux = int(x8)

    x8IQ = real(iaux,kind=real64) == x8 

    if(x8IQ) then
       fr = A**iaux
    else
       x8dn = realtodn(x8)
       fr = A**x8dn
    end if
  end function power_dr8
  !---------------------------------------------------------------------

  !chain rule, D^n (f(g))
  pure function Dnd(fc,gdual,n) result(dnfc)
    procedure(funzdual) :: fc
    type(dualzn), intent(in) :: gdual
    integer, intent(in) :: n
    complex(prec) :: dnfc
    type(dualzn) :: fvd
    complex(prec) :: g0, suma
    complex(prec), allocatable, dimension(:) :: xvg
    integer :: k, j

    g0 = gdual%f(0)
    fvd = fc(g0)
    if(n==0) then
       dnfc = fvd%f(0)
    else
       suma = 0
       do k=1,n
          allocate(xvg(1:n-k+1))
          do j=1,n-k+1
             xvg(j) = gdual%f(j)
          end do
          suma = suma + fvd%f(k)*BellY(n,k,xvg)
          deallocate(xvg)
       end do
       dnfc = suma
    end if
  end function Dnd
  !---------------------------------------------------------------------

  !dual/class
  elemental function divdX(B,X) result(fr)
    class(*), intent(in) :: X
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (X)
    type is (dualzn)
       Xd = X
    type is(complex(kind=real128))
       Xd = cmplxtodn(X)
    type is (complex(kind=real64))
       Xd =cmplxtodn(X)
    type is (complex(kind=real32))
       Xd =cmplxtodn(X)
    type is (real(kind=real128))
       Xd = realtodn(X)
    type is (real(kind=real64))
       Xd = realtodn(X)
    type is (real(kind=real32))
       Xd = realtodn(X)       
    type is (integer(kind=int64))
       Xd = itodn(X)
    type is (integer(kind=int32))
       Xd = itodn(X)
    end select
    fr = divd(B,Xd)
  end function divdX
  !---------------------------------------------------------------------

  !class/dual
  elemental function divXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (XX)
    type is (dualzn)
       Xd = XX
    type is(complex(kind=real128))
       Xd = cmplxtodn(XX)
    type is (complex(kind=real64))
       Xd =cmplxtodn(XX)
    type is (complex(kind=real32))
       Xd =cmplxtodn(XX)       
    type is (real(kind=real128))
       Xd = realtodn(XX)
    type is (real(kind=real64))
       Xd = realtodn(XX)
    type is (real(kind=real32))
       Xd = realtodn(XX)   
    type is (integer(kind=int64))
       Xd = itodn(XX)
    type is (integer(kind=int32))
       Xd = itodn(XX)
    end select
    fr = divd(Xd, BB)
  end function divXd
  !---------------------------------------------------------------------

  !A/B
  elemental function divd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr

    fr = A*inv(B)
  end function divd
  !---------------------------------------------------------------------

  !dual*class
  elemental function timesdX(B,X) result(fr)
    class(*), intent(in) :: X
    type(dualzn), intent(in) :: B
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (X)
    type is (dualzn)
       Xd = X
    type is(complex(kind=real128))
       Xd = cmplxtodn(X)
    type is (complex(kind=real64))
       Xd =cmplxtodn(X)
    type is (complex(kind=real32))
       Xd =cmplxtodn(X)
    type is (real(kind=real128))
       Xd = realtodn(X)
    type is (real(kind=real64))
       Xd = realtodn(X)
    type is (real(kind=real32))
       Xd = realtodn(X)         
    type is (integer(kind=int64))
       Xd = itodn(X)
    type is (integer(kind=int32))
       Xd = itodn(X)
    end select
    fr = timesd(B,Xd)
  end function timesdX
  !---------------------------------------------------------------------

  !class*dual
  elemental function timesXd(XX,BB) result(fr)
    class(*), intent(in) :: XX
    type(dualzn), intent(in) :: BB
    type(dualzn) :: fr
    type(dualzn) :: Xd

    select type (XX)
    type is (dualzn)
       Xd = XX
    type is(complex(kind=real128))
       Xd = cmplxtodn(XX)
    type is (complex(kind=real64))
       Xd =cmplxtodn(XX)
    type is (complex(kind=real32))
       Xd =cmplxtodn(XX)
    type is (real(kind=real128))
       Xd = realtodn(XX)
    type is (real(kind=real64))
       Xd = realtodn(XX)
    type is (real(kind=real32))
       Xd = realtodn(XX)         
    type is (integer(kind=int64))
       Xd = itodn(XX)
    type is (integer(kind=int32))
       Xd = itodn(XX)
    end select
    fr = timesd(Xd, BB)
  end function timesXd
  !---------------------------------------------------------------------

  !A*B
  elemental function timesd(A,B) result(fr)
    type(dualzn), intent(in) :: A, B
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order))
    do k=0,order
       fr%f(k)=timesdzn(A,B,k)
    end do
  end function timesd

  pure function timesdzn(A,B,k) result(fr)
    type(dualzn), intent(in) :: A, B
    integer, intent(in) :: k     
    complex(prec) :: fr
    integer :: i

    fr=0.0_prec
    do i=0,k
       fr=fr+binomial(k,i)*A%f(i)*B%f(k-i)
    end do
  end function timesdzn
  !---------------------------------------------------------------------

  !inverse function
  elemental function inv(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(invzdn,g,k)
    end do
  end function inv

  pure function invzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k,signo

    allocate(fr%f(0:order)) 

    signo=1
    do k=0,order
       fr%f(k)=signo*gamma(real(k+1,kind=prec))/(z**(k+1))
       signo=-signo
    end do
  end function invzdn
  !---------------------------------------------------------------------

  !acos function
  elemental function acosd_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(acoszdn,g,k)
    end do
  end function acosd_

  pure function acoszdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, sqz
    type(dualzn) :: sqzdn, zdn

    allocate(fr%f(0:order))
    fr%f(0) = acos(z)
    if (order == 0) return

    fr%f(1) = -1.0_prec/sqrt(1.0_prec - z**2)
    if (order == 1) return

    call initialize_dualzn(zdn)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    sqz = sqrt(1.0_prec-z*z)
    sqzdn = sqrt(realtodn(1.0_prec) - zdn*zdn)
    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*sqzdn%f(i)*fr%f(k-i)/sqz
       end do
       fr%f(k) = -sumterm
    end do
  end function acoszdn
  !---------------------------------------------------------------------

  !asin function
  elemental function asind_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(asinzdn,g,k)
    end do
  end function asind_

  pure function asinzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, sqz
    type(dualzn) :: sqzdn, zdn

    allocate(fr%f(0:order))
    fr%f(0) = asin(z)
    if (order == 0) return

    fr%f(1) = 1.0_prec/sqrt(1.0_prec - z**2)
    if (order == 1) return

    call initialize_dualzn(zdn)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    sqz = sqrt(1.0_prec-z*z)
    sqzdn = sqrt(realtodn(1.0_prec) - zdn*zdn)
    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*sqzdn%f(i)*fr%f(k-i)/sqz
       end do
       fr%f(k) = -sumterm
    end do
  end function asinzdn
  !---------------------------------------------------------------------

  !atan function
  elemental function atand_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(atanzdn,g,k)
    end do
  end function atand_

  pure function atanzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn

    allocate(fr%f(0:order))
    fr%f(0) = atan(z)
    fr%f(1) = 1.0_prec/(1.0_prec + z**2)

    den = 1.0_prec + z*z

    call initialize_dualzn(auxdn)
    auxdn%f(0) = 1.0_prec + z*z
    auxdn%f(1) = 2.0_prec * z
    auxdn%f(2) = 2.0_prec

    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function atanzdn
  !---------------------------------------------------------------------

  !asinh function
  elemental function asinhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(asinhzdn,g,k)
    end do
  end function asinhd

  pure function asinhzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn, zdn

    allocate(fr%f(0:order))
    fr%f(0) = asinh(z)
    if (order == 0) return

    fr%f(1) = 1.0_prec/sqrt(1.0_prec + z**2)
    if (order == 1) return

    den = sqrt(1.0_prec + z*z)

    call initialize_dualzn(zdn)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    auxdn = sqrt(realtodn(1.0_prec) + zdn*zdn)

    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function asinhzdn
  !---------------------------------------------------------------------

  !acosh function
  elemental function acoshd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(acoshzdn,g,k)
    end do
  end function acoshd

  pure function acoshzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn, zdn, onedn

    allocate(fr%f(0:order))
    fr%f(0) = acosh(z)
    if (order == 0) return

    !do not 'simplify' they are complex
    den = sqrt(z - 1.0_prec) * sqrt(1.0_prec + z)
    fr%f(1) = 1.0_prec/den
    if (order == 1) return

    call initialize_dualzn(zdn)
    zdn%f(0) = z
    zdn%f(1) = 1.0_prec

    onedn =  realtodn(1.0_prec)
    auxdn = sqrt(zdn - onedn) * sqrt(zdn + onedn)
    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function acoshzdn
  !---------------------------------------------------------------------

  !atanh function
  elemental function atanhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(atanhzdn,g,k)
    end do
  end function atanhd

  pure function atanhzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, i
    complex(prec) :: sumterm, den
    type(dualzn) :: auxdn

    allocate(fr%f(0:order))
    fr%f(0) = atanh(z)

    !do not 'simplify' they are complex
    den = 1.0_prec - z*z
    fr%f(1) = 1.0_prec/den

    call initialize_dualzn(auxdn)
    auxdn%f(0) = 1.0_prec - z*z
    auxdn%f(1) = -2.0_prec * z
    auxdn%f(2) = -2.0_prec

    do k=2,order
       sumterm = 0.0_prec
       do i=1,k-1
          sumterm = sumterm + binomial(k-1,i)*auxdn%f(i)*fr%f(k-i)/den
       end do
       fr%f(k) = -sumterm
    end do
  end function atanhzdn
  !---------------------------------------------------------------------

  !atan2d function
  elemental function atan2d_(y,x) result(fr)
    type(dualzn), intent(in) :: y, x
    type(dualzn) :: fr
    complex(prec) :: x0, y0

    fr = atan(y/x)
    x0 = x%f(0)
    y0 = y%f(0)
    fr%f(0) = atan2_z(y0,x0)
  end function atan2d_

  !Atan2 for complex arguments
  elemental function atan2_z(zy,zx) result (f_res)
    complex(prec), intent(in) :: zy, zx
    complex(prec) :: f_res
    complex(prec), parameter :: ii = cmplx(0,1,prec)
    complex(prec) :: num, den, divnd, t1, t2
    real(prec) :: x1, x2, y1, y2
    complex(prec) :: x1c, x2c, y1c, y2c

    x1 = real(zx,kind=prec)
    x2 = aimag(zx)

    y1 = real(zy,kind=prec)
    y2 = aimag(zy)

    x1c = x1
    x2c = x2
    y1c = y1
    y2c = y2

    num = x1c + ii*x2c + ii*(y1c + ii*y2c)
    den = sqrt((x1c + ii*x2c)**2 + (y1c + ii*y2c)**2)
    divnd = num/den;
    t1 = atan2(aimag(divnd),real(divnd,kind=prec))
    t2 = ii*log(sqrt((x2c + y1c)**2 + (x1c - y2c)**2)/((2.0_prec*x1c*&
         x2c +  2.0_prec*y1c*y2c)**2 + (x1c**2 - x2c**2 + y1c**2 -   &
         y2c**2)**2)**0.25_prec)

    f_res = t1 - t2
  end function atan2_z
  !---------------------------------------------------------------------

  !sinh function
  elemental function sinhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    fr = (exp(g)-exp(-g))
    do k=0,order
       fr%f(k) = 0.5_prec*fr%f(k)
    end do
  end function sinhd
  !---------------------------------------------------------------------

  !cosh function
  elemental function coshd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    fr = (exp(g)+exp(-g))
    do k=0,order
       fr%f(k) = 0.5_prec*fr%f(k)
    end do
  end function coshd
  !---------------------------------------------------------------------

  !tanh function
  elemental function tanhd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = (exp(g) - exp(-g))/(exp(g) + exp(-g))
  end function tanhd
  !---------------------------------------------------------------------

  !absx = sqrt(z*z) is not sqrt(z*conjg(z))
  elemental function absx(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    complex(prec) :: g0
    integer :: k

    g0 = g%f(0)

    fr = g
    do k=0,order
       fr%f(k) = g%f(k) * g0/sqrt(g0*g0)
    end do
  end function absx
  !---------------------------------------------------------------------

  !conjg
  !notice tat the conjugation operation is not differentiable. In the
  !below definitions we mean (df)* not d(f*)
  elemental function conjg_dzn(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    call initialize_dualzn(fr)
    do k=0, order
       fr%f(k) = conjg(g%f(k))
    end do
  end function conjg_dzn
  !--------------------------------------------------------------------

  !sum(AR2,dir) for rank2 arrays
  !the result is given as array of rank 1
  function sumR2dzn(AR2,dir) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    integer, intent(in) :: dir
    type(dualzn), dimension(size(AR2,2/dir)) :: fr
    complex(prec), dimension(size(AR2,2/dir)) :: sAk
    complex(prec), dimension(size(AR2,1),size(AR2,2)) :: Ak
    integer :: k,j, dimfr

    call initialize_dualzn(fr)
    dimfr = size(AR2,2/dir)

    do k=0,order
       Ak = f_part(AR2,k)
       sAk = sum(Ak,dir)
       do j=1,dimfr
          fr(j)%f(k) = sAk(j)
       end do
    end do
  end function sumR2dzn
  !---------------------------------------------------------------------

  !sum for Rank 2 array
  function sumR20dzn(AR2) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    type(dualzn) :: fr
    complex(prec), dimension(size(AR2,1),size(AR2,2)) :: Azaux
    integer :: k

    call initialize_dualzn(fr)
    do k=0,order
       Azaux = f_part(AR2,k)
       fr%f(k) = sum(Azaux)
    end do
  end function sumR20dzn
  !---------------------------------------------------------------------

  !sum for Rank 1 array
  function sumR1dzn(AR1) result(fr)
    type(dualzn), intent(in), dimension(:) :: AR1
    type(dualzn) :: fr
    complex(prec), dimension(size(AR1)) :: Azaux
    integer :: k

    call initialize_dualzn(fr)
    do k=0,order
       Azaux = f_part(AR1,k)
       fr%f(k) = sum(Azaux)
    end do
  end function sumR1dzn
  !---------------------------------------------------------------------

  !product(AR2,dir) for Rank 2 array
  function prodR2dzn(AR2,dir) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    integer, intent(in) :: dir
    type(dualzn), dimension(size(AR2,2/dir)) :: fr
    type(dualzn), dimension(size(AR2,dir)) :: vkdir
    integer :: k

    if(dir==1) then
       do k = 1, size(AR2,2)
          vkdir = AR2(:,k)
          fr(k) = prodR1dzn(vkdir)
       end do
    else if(dir==2) then
       do k = 1, size(AR2,1)
          vkdir = AR2(k,:)
          fr(k) = prodR1dzn(vkdir)
       end do
    else 
       stop 'use 1 (2) to collapse rows (columns) in product function'
    end if
  end function prodR2dzn

  !product for Rank 2 array
  function prodR20dzn(AR2) result(fr)
    type(dualzn), intent(in), dimension(:,:) :: AR2
    type(dualzn) :: fr
    integer :: k, m

    m=size(AR2,1)
    fr = 1.0_prec
    do k=1,m
       fr = fr*prodR1dzn(AR2(k,:))
    end do
  end function prodR20dzn

  !product for Rank 1 array
  function  prodR1dzn(x) result(fr)
    type(dualzn), intent(in), dimension(:) :: x
    type(dualzn) :: fr
    integer :: k

    fr = 1.0_prec
    do k=1,size(x)
       fr = fr*x(k)
    end do
  end function prodR1dzn
  !---------------------------------------------------------------------

  !log function
  elemental function logd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(logzdn,g,k)
    end do
  end function logd

  pure function logzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k, signo

    allocate(fr%f(0:order))
    fr%f(0) = log(z)

    signo = 1
    do k=1,order
       fr%f(k) = signo*gamma(real(k,kind=prec))/(z**k)
       signo = -signo
    end do
  end function logzdn
  !---------------------------------------------------------------------

  !exp function
  elemental function expd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(expzdn,g,k)
    end do
  end function expd

  pure function expzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = exp(z)
    end do
  end function expzdn
  !---------------------------------------------------------------------

  !sin function
  elemental function sind_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = Dnd(sinzdn,g,k)
    end do
  end function sind_

  pure function sinzdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = sin(z + k*Pi/2)
    end do
  end function sinzdn
  !---------------------------------------------------------------------

  !cos function
  elemental function cosd_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order))

    do k=0,order
       fr%f(k) = Dnd(coszdn,g,k)
    end do
  end function cosd_

  pure function coszdn(z) result(fr)
    complex(prec), intent(in) :: z
    type(dualzn) :: fr
    integer :: k

    allocate(fr%f(0:order)) 

    do k=0,order
       fr%f(k) = cos(z + k*Pi/2)
    end do
  end function coszdn
  !---------------------------------------------------------------------

  !tan function
  elemental function tand_(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = sind_(g)*inv(cosd_(g))
  end function tand_
  !---------------------------------------------------------------------

  !sqrt function
  elemental function sqrtd(g) result(fr)
    type(dualzn), intent(in) :: g
    type(dualzn) :: fr

    fr = g**(0.5_prec)
  end function sqrtd
  !---------------------------------------------------------------------

  !Combinations
  !m!/((m-n)! * n!) 
  pure function binomial(m, n) result(binom)
    integer, intent(in) :: m, n
    real(prec) :: binom
    integer :: j

    if (n == 0 .or. n == m) then
       binom = 1.0
    else
       binom = 1.0
       do j = 1, n
          binom = binom*(m-j+1)/j
       end do
    endif
  end function binomial

  pure function BellY(n, k, x) result(result_value)
    integer, intent(in) :: n, k
    complex(prec), dimension(:), intent(in) :: x
    complex(prec) :: result_value
    complex(prec), dimension(n+1,k+1) :: dp
    complex(prec) :: sum_val
    complex(prec), dimension(:), allocatable :: newx
    integer :: nn, kk, ii, LX, nx

    dp = 0.0
    dp(1, 1) = 1.0  

    do nn = 1, n
       dp(nn+1, 1) = 0.0
    end do

    do kk = 1, k
       dp(1, kk+1) = 0.0
    end do

    !special cases
    if (n == 0 .and. k == 0) then
       result_value = 1.0
       return
    elseif (size(x) == 0) then
       result_value = 0.0
       return
    end if

    !Main loop to compute BellY[n, k, x]
    do nn = 1, n
       do kk = 1, k
          LX = size(x)
          nx = max(nn - kk + 1, LX)
          if (nx > 0) then
             allocate(newx(nx))
             newx = 0.0
             newx(1:LX) = x
             sum_val = 0.0
             do ii = 0, nn - kk
                sum_val = sum_val + binomial(nn - 1, ii) * &
                     newx(ii + 1)*dp(nn - ii, kk)
             end do
             dp(nn + 1, kk + 1) = sum_val
             deallocate(newx)
          end if
       end do
    end do

    result_value = dp(n + 1, k + 1)
  end function BellY
end module dualzn_mod
