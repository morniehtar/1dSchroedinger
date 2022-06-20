!Поиск корней
!Нормировка в.ф.
!Отрисовка в.ф.

module config
    implicit none
    public

    !Integration precision (also affects smoothness)
    integer, parameter :: prec = 10000 ! 100000 is precise and fast
    !Draw smoothness
    integer, parameter :: dots = 1000
    !Threads count
    !integer, parameter :: thr = 4
!---------------------------------------------------------------
    !Effective infinity
    real(8), parameter :: xinf = 10.d0 ! U(x)~0 for x=10.d0
    !Cross-linking point
    real(8), parameter :: xcrs = -0.456d0
!---------------------------------------------------------------
    !Energy effective infinity
    real(8), parameter :: einf = 10.d0 ! >0
    !Energy effective nought
    real(8), parameter :: enou = 4.d-3 ! >0; 4.d-3 for U1
    !Machine epsilon (root search limiter)
    real(8), parameter :: eps = 1.d-8 ! 1.d-16 is machine epsilon
    !Energy search precision
    integer, parameter :: eprec = 60
!---------------------------------------------------------------
end module config

module potent
    implicit none
    public

    abstract interface
        pure real(8) function ufct(x)
        real(8), intent(in) :: x
        end function ufct
    end interface


    procedure(U0), pointer :: uptr => U0

    real(8), parameter :: u = 1.d1
    !real(8), parameter :: u = 2.d0
    !real(8), parameter :: u = 1.1d0


contains

    pure real(8) function U0(x)
        implicit none
        real(8), intent(in) :: x
        U0 = - u/cosh(x)
    end function U0

    pure real(8) function U1(x)
        implicit none
        real(8), intent(in) :: x
        U1 = - u/(cosh(x)**2)
    end function U1

    pure real(8) function Y(r)
        implicit none
        real(8), intent(in) :: r
        real(8), parameter :: g = 2.d0
        integer, parameter :: l = 0
        Y = - g*exp(-r)/r + 2.d-1*l*(l+1)/(r**2)
    end function Y

end module potent

program main
    use config
    use potent
    implicit none

    real(8) :: erg, stp
    real(8) :: bt, tp, estep

    integer :: i, k

    abstract interface
        real(8) function cnd(nrg)
        real(8), intent(in) :: nrg
        end function cnd
    end interface

    procedure(cnd), pointer :: cndPtr

    !vvvvvvvvvvvvvvv
    cndPtr => cndWron ! cndLog doesn't work at all with finding roots (why?!)
    !^^^^^^^^^^^^^^^

    stp = (einf-enou)/dots
    open(unit = 1, file = "cndPlot.dat")
    do i = 1, dots+1
        erg = -einf + (i-1)*stp
        write(1, *) erg, cndPtr(erg)
    end do
    close(1)
    print *, "Cnd(erg) is drawn"

    k=0
    estep = (einf-enou)/eprec
    do i = 1, eprec
        bt = -einf + (i-1)*estep
        tp = -einf + i*estep
        call getRoot(cndPtr, bt, tp, k)
    end do

    !print *, cndWron(-1.5d0)

contains

    real(8) function cndLog(nrg)
        implicit none

        real(8), intent(in) :: nrg
        real(8) :: yleft, zleft, yright, zright

        call integL(nrg, yleft, zleft)
        call integR(nrg, yright, zright)

        cndLog = zright/yright - zleft/yleft

    end function cndLog

    real(8) function cndWron(nrg)
        implicit none

        real(8), intent(in) :: nrg
        real(8) :: yleft, zleft, yright, zright

        call integL(nrg, yleft, zleft)
        call integR(nrg, yright, zright)

        cndWron = yleft*zright - zleft*yright

    end function cndWron

    !y'=func(x, y, z); y = psi
    !z'=gunc(x, y, z); zdx = dpsi

    pure real(8) function func(x, y, z, nrg)
        implicit none
        real(8), intent(in) :: x, y, z, nrg
        func = z
    end function func

    pure real(8) function gunc(x, y, z, nrg)
        implicit none
        real(8), intent(in) :: x, y, z, nrg
        !ħ=m=1
        gunc = 2.d0*(uptr(x)-nrg)*y
    end function gunc

    subroutine integL(nrg, yleft, zleft)
        use config
        implicit none

        real(8), intent(inout) :: yleft, zleft
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i

        fptr => func
        gptr => gunc

        yleft = exp(-sqrt(-2.d0*nrg)*xinf)
        zleft = sqrt(-2.d0*nrg)*exp(-sqrt(-2.d0*nrg)*xinf)

        step = (xcrs + xinf) / prec

        !open(unit = 1, file = "rkLeft.dat")
        do i = 1, prec+1
            x = -xinf + (i-1)*step
            !write(1,*) x, yleft, zleft
            k(1) = fptr(x, yleft, zleft, nrg)
            l(1) = gptr(x, yleft, zleft, nrg)

            k(2) = fptr(x+step/2, yleft+k(1)*step/2, zleft+k(1)*step/2, nrg)
            l(2) = gptr(x+step/2, yleft+l(1)*step/2, zleft+l(1)*step/2, nrg)

            k(3) = fptr(x+step/2, yleft+k(2)*step/2, zleft+k(2)*step/2, nrg)
            l(3) = gptr(x+step/2, yleft+l(2)*step/2, zleft+l(2)*step/2, nrg)

            k(4) = fptr(x+step, yleft+k(3)*step, zleft+k(3)*step, nrg)
            l(4) = gptr(x+step, yleft+l(3)*step, zleft+l(3)*step, nrg)

            yleft = yleft + step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zleft = zleft + step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do
        !close(1)

    end subroutine integL

    subroutine integR(nrg, yright, zright)
        use config
        implicit none

        real(8), intent(inout) :: yright, zright
        real(8), intent(in) :: nrg

        abstract interface
            real(8) function fct(x, y, z, nrg)
            real(8), intent(in) :: x, y, z, nrg
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i

        fptr => func
        gptr => gunc

        yright = exp(-sqrt(-2.d0*nrg)*xinf)
        zright = -sqrt(-2.d0*nrg)*exp(-sqrt(-2.d0*nrg)*xinf)

        step = (xinf - xcrs) / prec

        !open(unit = 1, file = "rkRight.dat")
        do i = 1, prec+1
            x = xinf - (i-1)*step
            !write(1,*) x, yright, zright
            k(1) = fptr(x, yright, zright, nrg)
            l(1) = gptr(x, yright, zright, nrg)

            k(2) = fptr(x-step/2, yright-k(1)*step/2, zright-k(1)*step/2, nrg)
            l(2) = gptr(x-step/2, yright-l(1)*step/2, zright-l(1)*step/2, nrg)

            k(3) = fptr(x-step/2, yright-k(2)*step/2, zright-k(2)*step/2, nrg)
            l(3) = gptr(x-step/2, yright-l(2)*step/2, zright-l(2)*step/2, nrg)

            k(4) = fptr(x-step, yright-k(3)*step, zright-k(3)*step, nrg)
            l(4) = gptr(x-step, yright-l(3)*step, zright-l(3)*step, nrg)

            yright = yright - step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zright = zright - step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do
        !close(1)

    end subroutine integR

    subroutine getRoot(fptr, bt, tp, k)
        implicit none

        abstract interface
            real(8) function fct(x)
            real(8), intent(in) :: x
            end function fct
        end interface

        procedure(fct), pointer, intent(in) :: fptr
        real(8), intent(inout) :: bt, tp
        integer, intent(inout) :: k
        real(8) :: md
        integer :: j

        if (fptr(bt)*fptr(tp).gt.0) return
        k = k+1

        md = (tp + bt) / 2
        j = 0
        do while ((abs(fptr(md)).gt.eps).and.((tp-bt).gt.eps))
            md = (tp + bt) / 2
            if (fptr(bt) * fptr(md) .gt. 0) then
                bt = md
            else if (fptr(bt) * fptr(md) .lt. 0) then
                tp = md
            end if
            j = j + 1
        end do

        print *, "(k =", k, ") E =", md, j , "steps"
    end subroutine getRoot

end program main
