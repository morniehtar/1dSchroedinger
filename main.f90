module config
    implicit none
    public

    !Integration precision (also affects smoothness)
    integer, parameter :: dots = 100000 ! 100000 is precise and fast
    !Threads count
    !integer, parameter :: thr = 4
!---------------------------------------------------------------
    !Effective infinity
    real(8), parameter :: xinf = 10.d0 ! Epsilon is 1.d-16, U(x)~0 for x=1.d1
    !Crosslinking point
    real(8), parameter :: xcrs = -0.5d0
    !Energy seeking break limit
    real(8), parameter :: einf = 10.d0
    !Machine epsilon (bisection limiter)
    real(8), parameter :: eps = 1.d-9
    !
    real(8), parameter :: estep = 1.d0
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

    !Current U(x) function
    procedure(ufct), pointer :: uptr => U0

contains

    pure real(8) function U0(x)
        implicit none
        real(8), intent(in) :: x
        U0 = - 1.1d0/cosh(x)
    end function U0

    pure real(8) function U1(x)
        implicit none
        real(8), intent(in) :: x
        U1 = - 1.1d0/(cosh(x)**2)
    end function U1

    pure real(8) function U2(r)
        implicit none
        real(8), intent(in) :: r
        real(8), parameter :: g = 2.d0
        U2 = - g*exp(-r)/r
    end function U2

end module potent

program main
    use config
    use potent
    implicit none

    procedure(getLog), pointer :: ptr
    real(8) :: yleft, zleft, yright, zright
    integer :: i, lim
    real(8) :: eBt, eTp, eMid, stBot, stMid, stE
    real(8) :: erg

    !erg = -8.0004364013671960d0
    !erg = -4.5006561279297079d0
    !erg = -2.0004974365234567d0
    !erg = -0.50028381347658124d0
    !erg = -0.592d0

    lim = einf/estep
    ptr => getLog

    do i = 1, lim+1
        eBt = -einf + (i-1)*estep
        eTp = -einf + i*estep
        call bisection(ptr, eBt, eTp)
    end do

    !print *, "Log derivatives:", getLog(yleft, zleft, yright, zright)
    !print *, "Wronsky:", getWron(yleft, zleft, yright, zright)


contains

    subroutine bisection(fptr, a, b)
        !Please, ensure that the function is monotoneous, otherwise
        !the criteria f(a)*f(b)<0 doesn't really make that much sense
        !The a and b aren't initial guesses, those are the borders
        implicit none

        abstract interface
            real(8) function ftst(x)
            real(8), intent(in) :: x
            end function
        end interface

        procedure(ftst), pointer, intent(in) :: fptr
        real(8), intent(in) :: a, b
        real(8) :: bt, md, tp
        integer :: j

        if (fptr(a) * fptr(b) .gt. 0) then
            print *, "Bisection error: root criteria isn't met"
            return
        end if

        bt = a
        tp = b
        md = (tp + bt) / 2
        j = 0

        do while ((abs(fptr(md)) .gt. eps).and.((tp-bt) .gt. eps))
        !do while (abs(fptr(md)) .gt. eps)
            md = (tp + bt) / 2
            if (fptr(bt) * fptr(md) .gt. 0) then
                bt = md
            else if (fptr(bt) * fptr(md) .lt. 0) then
                tp = md
            end if
            j = j + 1
        end do

        print *, "Bisection gives x =", md, "as a root, in", j , "steps"

    end subroutine bisection

    pure real(8) function getWron(yleft, zleft, yright, zright)
        implicit none
        real(8), intent(in) :: yleft, zleft, yright, zright
        getWron = yleft*zright-zleft*zright
    end function getWron

    real(8) function getLog(param)
        implicit none

        abstract interface
            pure real(8) function fct(x, y, z, param)
            real(8), intent(in) :: x, y, z, param
            end function fct
        end interface

        procedure(fct), pointer :: fptr, gptr

        real(8), intent(in) :: param
        real(8) :: yleft, zleft, yright, zright

        fptr => func
        gptr => gunc

        call integL(fptr, gptr, -xinf, xcrs, param)
        call integR(fptr, gptr, xcrs, xinf, param)

        getLog = zright/yright - zleft/yleft
    end function getLog

    !y'=func(x, y, z); y ~ psi
    !z'=gunc(x, y, z); z ~ dpsi/dx

    pure real(8) function func(x, y, z, param)
        implicit none
        real(8), intent(in) :: x, y, z, param
        func = z
    end function func

    pure real(8) function gunc(x, y, z, param)
        implicit none
        real(8), intent(in) :: x, y, z, param

        gunc = 2.d0*(uptr(x)-param)*y
        !ħ=m=1
    end function gunc

    subroutine integL(fptr, gptr, st, ed, param)
        use config
        implicit none

        abstract interface
            pure real(8) function fct(x, y, z, param)
            real(8), intent(in) :: x, y, z, param
            end function fct
        end interface

        procedure(fct), pointer, intent(in) :: fptr, gptr
        real(8), intent(in) :: st, ed, param
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i, j
        real(8) :: yleft, zleft

        yleft = exp(-sqrt(-2.d0*param)*xinf)
        zleft = sqrt(-2.d0*param)*exp(-sqrt(-2.d0*param)*xinf)

        step = abs(ed - st) / dots

        !open(unit = 1, file = "rkLeft.dat")
        do i = 1, dots+1
            x = st + (i-1)*step
            !write(1,*) x, yleft, zleft
            k(1) = fptr(x, yleft, zleft, param)
            l(1) = gptr(x, yleft, zleft, param)

            k(2) = fptr(x+step/2, yleft+k(1)*step/2, zleft+k(1)*step/2, param)
            l(2) = gptr(x+step/2, yleft+l(1)*step/2, zleft+l(1)*step/2, param)

            k(3) = fptr(x+step/2, yleft+k(2)*step/2, zleft+k(2)*step/2, param)
            l(3) = gptr(x+step/2, yleft+l(2)*step/2, zleft+l(2)*step/2, param)

            k(4) = fptr(x+step, yleft+k(3)*step, zleft+k(3)*step, param)
            l(4) = gptr(x+step, yleft+l(3)*step, zleft+l(3)*step, param)

            yleft = yleft + step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zleft = zleft + step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do
        !close(1)

    end subroutine integL

    subroutine integR(fptr, gptr, st, ed, param)
        use config
        implicit none

        abstract interface
            pure real(8) function fct(x, y, z, param)
            real(8), intent(in) :: x, y, z, param
            end function fct
        end interface

        procedure(fct), pointer, intent(in) :: fptr, gptr
        real(8), intent(in) :: st, ed, param
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i, j
        real(8) :: yright, zright

        yright = exp(-sqrt(-2.d0*param)*xinf)
        zright = sqrt(-2.d0*param)*exp(-sqrt(-2.d0*param)*xinf)

        step = abs(ed - st) / dots

        !open(unit = 1, file = "rkRight.dat")
        do i = 1, dots+1
            x = ed - (i-1)*step
            !write(1,*) x, yright, zright
            k(1) = fptr(x, yright, zright, param)
            l(1) = gptr(x, yright, zright, param)

            k(2) = fptr(x-step/2, yright-k(1)*step/2, zright-k(1)*step/2, param)
            l(2) = gptr(x-step/2, yright-l(1)*step/2, zright-l(1)*step/2, param)

            k(3) = fptr(x-step/2, yright-k(2)*step/2, zright-k(2)*step/2, param)
            l(3) = gptr(x-step/2, yright-l(2)*step/2, zright-l(2)*step/2, param)

            k(4) = fptr(x-step, yright-k(3)*step, zright-k(3)*step, param)
            l(4) = gptr(x-step, yright-l(3)*step, zright-l(3)*step, param)

            yright = yright - step*(k(1)+2*k(2)+2*k(3)+k(4))/6
            zright = zright - step*(l(1)+2*l(2)+2*l(3)+l(4))/6

        end do
        !close(1)

    end subroutine integR

end program main

