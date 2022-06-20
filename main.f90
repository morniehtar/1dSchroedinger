module config
    implicit none
    public

    !Integration precision (also affects smoothness)
    integer, parameter :: dots = 1000 ! 100000 is precise and fast
!---------------------------------------------------------------
    !Effective infinity
    real(8), parameter :: xinf = 10.d0 ! Epsilon is 1.d-16, U(x)~0 for x=10.d0
    !Cross-linking point
    real(8), parameter :: xcrs = -0.5d0
    !Energy seeking break limit
    real(8), parameter :: einf = 10.d0
    !Machine epsilon (bisection limiter)
    real(8), parameter :: eps = 1.d-9
    !Energy search width
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
    procedure(U0), pointer :: uptr => U1

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

    real(8) :: yleft, zleft, yright, zright
    real(8) :: erg

    erg = -0.592d0

    call integL(erg, yleft, zleft)
    call integR(erg, yright, zright)

    print *, getLog(yleft, zleft, yright, zright)

contains

    pure real(8) function getWron(yleft, zleft, yright, zright)
        implicit none
        real(8), intent(in) :: yleft, zleft, yright, zright
        getWron = yleft*zright - zleft*zright
    end function getWron

    pure real(8) function getLog(yleft, zleft, yright, zright)
        implicit none
        real(8), intent(in) :: yleft, zleft, yright, zright
        getLog = zright/yright - zleft/yleft
    end function getLog

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

        procedure(gunc), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i, j

        fptr => func
        gptr => gunc

        yleft = exp(-sqrt(-2.d0*nrg)*xinf)
        zleft = sqrt(-2.d0*nrg)*exp(-sqrt(-2.d0*nrg)*xinf)

        step = (xcrs + xinf) / dots

        !open(unit = 1, file = "rkLeft.dat")
        do i = 1, dots+1
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

        procedure(gunc), pointer :: fptr, gptr
        real(8), dimension(4) :: k, l
        real(8) :: x, step
        integer :: i, j

        fptr => func
        gptr => gunc

        yright = exp(-sqrt(-2.d0*nrg)*xinf)
        zright = sqrt(-2.d0*nrg)*exp(-sqrt(-2.d0*nrg)*xinf)

        step = (xinf - xcrs) / dots

        !open(unit = 1, file = "rkRight.dat")
        do i = 1, dots+1
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

end program
