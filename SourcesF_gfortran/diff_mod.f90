module diff_mod
  use precision_mod !Depending on how to link, it will set prec to be
  !real64 or real128
  use dualzn_mod
  implicit none

  private
  public :: d1fscalar, d2fscalar, d1fvector
  public :: Jacobian, Hessian, gradient

  !Abstract interfaces for functions passed as arguments to other
  !functions.
  !
  !Interface for a scalar dual function f: D^m --> D (similar to
  !f: R^m --> R)
  abstract interface
     function fsdual(xd) result(frsd)
       use dualzn_mod
       type(dualzn), intent(in), dimension(:) :: xd
       type(dualzn) :: frsd
     end function fsdual
  end interface

  !Interface for a vector dual function f: D^m --> D^n
  interface
     function fvecdual(xd) result(frd)
       use dualzn_mod
       type(dualzn), intent(in), dimension(:)  :: xd
       type(dualzn), allocatable, dimension(:) :: frd
     end function fvecdual
  end interface

  interface d2fscalar
     module procedure d2fscalarvv
     module procedure d2fscalaruv       
  end interface d2fscalar
  !---------------------------------------------------------------------

contains
  !Hessian operator
  function Hessian(fsd,qcmplx) result(Hmat)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: qcmplx    
    complex(prec), dimension(size(qcmplx),size(qcmplx)) :: Hmat

    complex(prec), dimension(size(qcmplx)) :: ei, ej
    integer :: i,j,m

    m = size(qcmplx)

    Hmat = 0
    do i=1,m
       ei = 0
       ei(i) = 1
       Hmat(i,i) = d2fscalarvv(fsd,ei,qcmplx)
    end do

    do i=1, m
       ei = 0
       ei(i) = 1
       do j=i+1,m
          ej = 0
          ej(j) = 1
          Hmat(i,j) = d2fscalaruv(fsd,ei,ej,qcmplx)
          Hmat(j,i) = Hmat(i,j)
       end do
    end do
  end function Hessian

  !d2fscalaruv(f,u,v,q) gives the second-order directional derivative
  !of the scalar function f: D^m --> D, where f = f(q), along vectors u
  !and v, evaluated at point q
  function d2fscalaruv(fsd,x,y,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: x, y, q
    complex(prec) :: fr

    if(all(x == y)) then
       fr = d2fscalarvv(fsd,x,q)
    else
       fr = 0.5_prec*(d2fscalarvv(fsd,x + y,q) - d2fscalarvv(fsd,x,q) -&
            d2fscalarvv(fsd,y,q))
    end if
  end function d2fscalaruv

  !Second-order directional derivative of a scalar function along vector
  !v, evaluated at point q
  function d2fscalarvv(fsd,v,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: v, q
    complex(prec) :: fr
    type(dualzn) :: eps1

    call set_order(2)   
    eps1 = 0
    eps1%f(1) = 1

    fr = f_part(fsd(add_vecdzn(cmplxtodn(q),eps1*v)),2)  
  end function d2fscalarvv

  !Jacobian operator: To optimize efficiency, we include the parameter
  !'n', representing the dimension of 'fvecd'.
  function Jacobian(fvecd,qcmplx,n) result(Jmat)
    procedure(fvecdual) :: fvecd
    complex(prec), intent(in), dimension(:) :: qcmplx
    integer, intent(in) :: n
    complex(prec), dimension(n,size(qcmplx)) :: Jmat

    complex(prec), dimension(size(qcmplx),n) :: TrpJmat
    complex(prec), dimension(size(qcmplx))      :: ei
    integer :: i

    do i = 1,size(qcmplx)
       ei = 0
       ei(i) = 1      
       TrpJmat(i,:) = d1fvector(fvecd,ei,qcmplx,n)
    end do

    Jmat = transpose(TrpJmat)
  end function Jacobian

  !First-order directional derivative of a vector function along vector
  !v, evaluated at point q. To optimize efficiency, we include the
  !parameter 'n', representing the dimension of 'fvecd'.
  !
  function d1fvector(fvecd,v,q,n) result(fr)
    procedure(fvecdual) :: fvecd
    complex(prec), intent(in), dimension(:) :: v, q
    integer, intent(in) :: n
    complex(prec), dimension(n) :: fr  
    type(dualzn) :: eps1
    type(dualzn), allocatable, dimension(:) ::  auxv
    integer :: k

    call set_order(1)
    eps1 = 0
    eps1%f(1) = 1

    allocate(auxv(n))
    auxv = fvecd(add_vecdzn(cmplxtodn(q),eps1*v))
    fr = f_part(auxv,1)

    do k=1,n 
      deallocate(auxv(k)%f)
    end do 
    deallocate(auxv)

  end function d1fvector

  function gradient(fsd,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: q
    complex(prec), dimension(size(q)) :: fr
    complex(prec), dimension(size(q)) :: ei
    integer :: i

    do i = 1,size(q)
       ei = 0
       ei(i) = 1      
       fr(i) = d1fscalar(fsd,ei,q)
    end do
  end function gradient

  !First-order directional derivative of a scalar function along vector
  !v, evaluated at point q
  function d1fscalar(fsd,v,q) result(fr)
    procedure(fsdual) :: fsd
    complex(prec), intent(in), dimension(:) :: v, q
    complex(prec) :: fr  
    type(dualzn) :: eps1

    call set_order(1)
    eps1 = 0
    eps1%f(1) = 1

    fr = f_part(fsd(add_vecdzn(cmplxtodn(q),eps1*v)),1)  
  end function d1fscalar

  function add_vecdzn(x,y) result(fr)
    type(dualzn), intent(in), dimension(:) :: x, y
    type(dualzn), dimension(size(x)) :: fr
    
    integer :: k
    
    do k=1, size(x)
      fr(k)=x(k) + y(k)
    end do
    end function add_vecdzn

end module diff_mod
