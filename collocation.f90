!<--! This is a module for collocation method. -->
module collocation
    implicit none
contains

    ! compute the Legendre polynomial
    elemental function legendrep(n, x) result(y)
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: x
        real :: y

        integer :: k
        real :: t1, t2, alpha

        if (n < 0) then
            y = 0.0
        else ! n >= 0
            ! k = 0
            t1 = 0.0   ! p_{-2}(x)
            t2 = 0.0   ! p_{-1}(x)
            y = 1.0     ! p_{0}(x)

            do k = 1, n
                t1 = t2
                t2 = y
                alpha = 1.0 / k
                y = (2.0 - alpha) * x * t2 - (1.0 - alpha) * t1
            end do
        end if
    end function legendrep

    ! compute the derivative of the Legendre polynomial
    elemental function dlegendrep(n, x) result(dy)
        implicit none
        integer, intent(in) :: n
        real, intent(in) :: x
        real :: y, dy
        
        integer :: k, n_mod_2
        real :: t1, t2, alpha

        if (n < 0) then
            y = 0.0
            dy = 0.0
        else if (n == 0) then
            y = 1.0
            dy = 0.0
        else ! n > 0
            n_mod_2 = mod(n,2)

            ! k = 0
            t1 = 0.0   ! p_{-2}(x)
            t2 = 0.0   ! p_{-1}(x)
            y = 1.0     ! p_{0}(x)

            if (0 /= n_mod_2) then
                dy = y
            else
                dy = 0.0
            end if

            do k = 1, n
                t1 = t2
                t2 = y
                alpha = 1.0 / k
                y = (2.0 - alpha) * x * t2 - (1.0 - alpha) * t1

                if (mod(k,2) /= n_mod_2) then
                    dy = dy + (2*k+1) * y
                end if
            end do
        end if
    end function dlegendrep

    ! compute the diagonal mass matrix of the Legendre polynomials
    pure subroutine legendre_mass_diag(n, ans)
        implicit none
        integer, intent(in) :: n
        real, dimension(0:n), intent(out) :: ans
        integer :: i
        do i = 0, n
            ans(i) = 1.0 / (i + 0.5)
        end do
    end subroutine legendre_mass_diag

    ! compute the differentiation matrix of the Legendre polynomials
    pure subroutine legendre_diff_mat(n, ans)
        implicit none
        integer, intent(in) :: n
        real, dimension(0:n, 0:n), intent(out) :: ans
        integer :: i, j
        do j = 0, n
            do i = 0, j-1
                if (mod(i+j,2)/=0) then
                    ans(i,j) = 2*i+1
                else
                    ans(i,j) = 0.0
                end if
            end do
            ans(j:n, j) = 0.0
        end do
    end subroutine legendre_diff_mat

    ! compute the n-point Gauss-Legendre quadrature nodes and weights on [-1,1]
    subroutine gauss_legendre(n, x, w)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: x, w
        select case(n)
        case (:0)
            error stop 'gauss_legendre: n must be positive'
        case (1)
            x = 0.0
            w = 2.0
        case (2)
            x(2) = 0.5773502691896257645091487805019574556476 ! sqrt(1.0/3)
            x(1) = -x(2)

            w(1) = 1.0
            w(2) = w(1)
        case (3)
            x(3) = 0.7745966692414833770358530799564799221666 ! sqrt(0.6)
            x(2) = 0.0
            x(1) = -x(3)

            w(1) = 0.5555555555555555555555555555555555555556 ! 5.0/9
            w(2) = 0.8888888888888888888888888888888888888889 ! 8.0/9
            w(3) = w(1)
        case (4)
            x(4) = 0.8611363115940525752239464888928095050957 ! sqrt((3.0 + 2*sqrt(1.2))/7)
            x(3) = 0.3399810435848562648026657591032446872006 ! sqrt((3.0 - 2*sqrt(1.2))/7)
            x(2) = -x(3)
            x(1) = -x(4)

            w(1) = 0.3478548451374538573730639492219994072353 ! 0.5 - sqrt(30.0)/36
            w(2) = 0.6521451548625461426269360507780005927647 ! 0.5 + sqrt(30.0)/36
            w(3) = w(2)
            w(4) = w(1)
        case (5)
            x(5) = 0.9061798459386639927976268782993929651257 ! sqrt((35.0 + 2*sqrt(70.0))/63)
            x(4) = 0.5384693101056830910363144207002088049673 ! sqrt((35.0 - 2*sqrt(70.0))/63)
            x(3) = 0.0
            x(2) = -x(4)
            x(1) = -x(5)

            w(1) = 0.2369268850561890875142640407199173626433 ! (322.0 - 13*sqrt(70.0))/900
            w(2) = 0.4786286704993664680412915148356381929123 ! (322.0 + 13*sqrt(70.0))/900
            w(3) = 0.5688888888888888888888888888888888888889 ! 128.0/225
            w(4) = w(2)
            w(5) = w(1)
        case default
            error stop 'gauss_legendre: n is too large'
        end select
    end subroutine gauss_legendre

    ! compute the n-point Gauss-Lobatto quadrature nodes and weights on [-1,1]
    subroutine gauss_lobatto(n, x, w)
        implicit none
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: x, w
        select case(n)
        case (:0)
            error stop 'gauss_lobatto: n must be positive'
        case (1)
            ! use the Gauss-Legendre quadrature instead
            print *, 'gauss_lobatto: Warning: Gauss-Lobatto quadrature with n=1 is not recommended'
            x = 0.0
            w = 2.0
        case (2)
            x(2) = 1.0
            x(1) = -x(2)

            w(1) = 1.0
            w(2) = w(1)
        case (3)
            x(3) = 1.0
            x(2) = 0.0
            x(1) = -x(3)

            w(1) = 0.3333333333333333333333333333333333333333 ! 1.0/3
            w(2) = 1.3333333333333333333333333333333333333333 ! 4.0/3
            w(3) = w(1)
        case (4)
            x(4) = 1.0
            x(3) = 0.4472135954999579392818347337462552470881 ! sqrt(0.2)
            x(2) = -x(3)
            x(1) = -x(4)

            w(1) = 0.1666666666666666666666666666666666666667 ! 1.0/6
            w(2) = 0.8333333333333333333333333333333333333333 ! 5.0/6
            w(3) = w(2)
            w(4) = w(1)
        case (5)
            x(5) = 1.0
            x(4) = 0.6546536707079771437982924562468583555692 ! sqrt(3.0/7)
            x(3) = 0.0
            x(2) = -x(4)
            x(1) = -x(5)

            w(1) = 0.1
            w(2) = 0.5444444444444444444444444444444444444444 ! 49.0/90
            w(3) = 0.7111111111111111111111111111111111111111 ! 32.0/45
            w(4) = w(2)
            w(5) = w(1)
        case (6)
            x(6) = 1.0
            x(5) = 0.7650553239294646928510029739593381503657 ! sqrt((7.0 + 2*sqrt(7.0))/21)
            x(4) = 0.2852315164806450963141509940408790719190 ! sqrt((7.0 - 2*sqrt(7.0))/21)
            x(3) = -x(4)
            x(2) = -x(5)
            x(1) = -x(6)

            w(1) = 0.06666666666666666666666666666666666666667 ! 1.0/15
            w(2) = 0.3784749562978469803166128082120246524763 ! (14.0 - sqrt(7.0))/30
            w(3) = 0.5548583770354863530167205251213086808570 ! (14.0 + sqrt(7.0))/30
            w(4) = w(3)
            w(5) = w(2)
            w(6) = w(1)
        case (7)
            x(7) = 1.0
            x(6) = 0.8302238962785669298720322139674651395872 ! sqrt((15.0 + 2*sqrt(15.0))/33)
            x(5) = 0.4688487934707142138037718819087663294056 ! sqrt((15.0 - 2*sqrt(15.0))/33)
            x(4) = 0.0
            x(3) = -x(5)
            x(2) = -x(6)
            x(1) = -x(7)

            w(1) = 0.04761904761904761904761904761904761904762 ! 1.0/21
            w(2) = 0.2768260473615659480107004062900662934976 ! (124.0 - 7*sqrt(15.0))/350
            w(3) = 0.4317453812098626234178710222813622779309 ! (124.0 + 7*sqrt(15.0))/350
            w(4) = 0.4876190476190476190476190476190476190476 ! 256.0/525
            w(5) = w(3)
            w(6) = w(2)
            w(7) = w(1)
        case default
            error stop 'gauss_lobatto: n is too large'
        end select
    end subroutine gauss_lobatto

    ! compute the collocation matrices for the Gauss-Lobatto quadrature
    ! n: order of collocation (n+1 collocation nodes)
    ! qx: collocation points
    ! qw: collocation weights
    ! M: mass matrix
    ! Minv: inverse of the mass matrix
    ! D: differentiation matrix
    ! Tl: left boundary transformation matrix
    ! Tr: right boundary transformation matrix
    ! Pmn: from Legendre modal to Lobatto nodal
    ! Pnm: from Lobatto nodal to Legendre modal
    subroutine Gauss_Lobatto_collocation_matrices(n, qx, qw, M, Minv, D, Tl, Tr, Pmn, Pnm)
        implicit none
        integer, intent(in) :: n
        real, dimension(0:n), intent(out) :: qx, qw
        real, dimension(0:n), intent(out) :: Tl, Tr
        real, dimension(0:n, 0:n), intent(out) :: M, Minv, D, Pmn, Pnm
        call gauss_lobatto(n, qx, qw)
        Tl(:) = 0.0
        Tr(:) = 0.0
        Tl(0) = 1.0
        Tr(n) = 1.0
        select case(n)
        case (:-1)
            error stop 'Gauss_Lobatto_collocation_matrices: n must be non-negatigve!'
        case (0)
            M(0,0) = 2.0
            Minv(0,0) = 0.5
            D(0,0) = 0.0
            Pmn(0,0) = 1.0
            Pnm(0,0) = 1.0
        case (1)
            M = transpose(reshape([ &
                2.0/3, 1.0/3, &
                1.0/3, 2.0/3], shape(M)))
            Minv = transpose(reshape([ &
                2.0, -1.0, &
                -1.0, 2.0], shape(Minv)))
            D = transpose(reshape([ &
                -0.5, 0.5, &
                -0.5, 0.5], shape(D)))
            Pmn = transpose(reshape([ &
                1.0, -1.0, &
                1.0, 1.0], shape(Pmn)))
            Pnm = transpose(reshape([ &
                0.5, 0.5, &
                -0.5, 0.5], shape(Pnm)))
        case (2)
            M = transpose(reshape([ &
                4.0/15, 2.0/15, -1.0/15, &
                2.0/15, 16.0/15, 2.0/15, &
                -1.0/15, 2.0/15, 4.0/15], shape(M)))
            Minv = transpose(reshape([ &
                4.5, -0.75, 1.5, &
                -0.75, 1.125, -0.75, &
                1.5, -0.75, 4.5], shape(Minv)))
            D = transpose(reshape([ &
                -1.5, 2.0, -0.5, &
                -0.5, 0.0, 0.5, &
                0.5, -2.0, 1.5], shape(D)))
            Pmn = transpose(reshape([ &
                1.0, -1.0, 1.0, &
                1.0, 0.0, -0.5, &
                1.0, 1.0, 1.0], shape(Pmn)))
            Pnm = transpose(reshape([ &
                1.0/6, 2.0/3, 1.0/6, &
                -0.5, 0.0, 0.5, &
                1.0/3, -2.0/3, 1.0/3], shape(Pnm)))
        case (3)
            M = reshape([ &
                1.0/7, sqrt(5.0)/42, -sqrt(5.0)/42, 1.0/42, &
                sqrt(5.0)/42, 5.0/7, 5.0/42, -sqrt(5.0)/42, &
                -sqrt(5.0)/42, 5.0/42, 5.0/7, sqrt(5.0)/42, &
                1.0/42, -sqrt(5.0)/42, sqrt(5.0)/42, 1.0/7], shape(M))
            Minv = transpose(reshape([ &
                8.0, -2.0/sqrt(5.0), 2.0/sqrt(5.0), -2.0, &
                -2.0/sqrt(5.0), 1.6, -0.4, 2.0/sqrt(5.0), &
                2.0/sqrt(5.0), -0.4, 1.6, -2.0/sqrt(5.0), &
                -2.0, 2.0/sqrt(5.0), -2.0/sqrt(5.0), 8.0], shape(Minv)))
            D = transpose(reshape([ &
                -3.0, 1.25*(1.0 + sqrt(5.0)), -1.25*(-1.0 + sqrt(5.0)), 0.5, &
                0.25*(-1.0 - sqrt(5.0)), 0.0, 0.5*sqrt(5.0), 0.25*(1.0 - sqrt(5.0)), &
                0.25*(-1.0 + sqrt(5.0)), -0.5*sqrt(5.0), 0.0, 0.25*(1.0 + sqrt(5.0)), &
                -0.5, 1.25*(-1.0 + sqrt(5.0)), -1.25*(1.0 + sqrt(5.0)), 3.0], shape(D)))
            Pmn = transpose(reshape([ &
                1.0, -1.0, 1.0, -1.0, &
                1.0, -1.0/sqrt(5.0), -0.2, 1.0/sqrt(5.0), &
                1.0, 1.0/sqrt(5.0), -0.2, -1.0/sqrt(5.0), &
                1.0, 1.0, 1.0, 1.0], shape(Pmn)))
            Pnm = transpose(reshape([ &
                1.0/12, 5.0/12, 5.0/12, 1.0/12, &
                -0.25, -0.25*sqrt(5.0), 0.25*sqrt(5.0), 0.25, &
                5.0/12, -5.0/12, -5.0/12, 5.0/12, &
                -0.25, 0.25*sqrt(5.0), -0.25*sqrt(5.0), 0.25], shape(Pnm)))
            return
        case (4)
            M = transpose(reshape([ &
                4.0/45, 7.0/270, -4.0/135, 7.0/270, -1.0/90, &
                7.0/270, 196.0/405, 28.0/405, -49.0/810, 7.0/270, &
                -4.0/135, 28.0/405, 256.0/405, 28.0/405, -4.0/135, &
                7.0/270, -49.0/810, 28.0/405, 196.0/405, 7.0/270, &
                -1.0/90, 7.0/270, -4.0/135, 7.0/270, 4.0/45], shape(M)))
            Minv = transpose(reshape([ &
                12.5, -15.0/14, 0.9375, -15.0/14, 2.5, &
                -15.0/14, 225.0/98, -45.0/112, 45.0/98, -15.0/14, &
                0.9375, -45.0/112, 225.0/128, -45.0/112, 0.9375, &
                -15.0/14, 45.0/98, -45.0/112, 225.0/98, -15.0/14, &
                2.5, -15.0/14, 0.9375, -15.0/14, 12.5], shape(Minv)))
            D = transpose(reshape([ &
                -5.0, (7.0/12)*(7.0 + sqrt(21.0)), -8.0/3, -(7.0/12)*(-7.0 + sqrt(21.0)), -0.5, &
                -(3.0/28)*(7.0 + sqrt(21.0)), 0.0, 8.0/sqrt(21.0), -0.5*sqrt(7.0/3), -(3.0/28)*(-7.0 + sqrt(21.0)), &
                0.375, -0.875*sqrt(7.0/3), 0.0, 0.875*sqrt(7.0/3), -0.375, &
                (3.0/28)*(-7.0 + sqrt(21.0)), 0.5*sqrt(7.0/3), -8.0/sqrt(21.0), 0.0, (3.0/28)*(7.0 + sqrt(21.0)), &
                0.5, (7.0/12)*(-7.0 + sqrt(21.0)), 8.0/3, -(7.0/12)*(7.0 + sqrt(21.0)), 5.0], shape(D)))
            Pmn = transpose(reshape([ &
                1.0, -1.0, 1.0, -1.0, 1.0, &
                1.0, -sqrt(3.0/7), 1.0/7, (3.0/7)*sqrt(3.0/7), -3.0/7, &
                1.0, 0.0, -0.5, 0.0, 0.375, &
                1.0, sqrt(3.0/7), 1.0/7, -(3.0/7)*sqrt(3.0/7), -3.0/7, &
                1.0, 1.0, 1.0, 1.0, 1.0], shape(Pmn)))
            Pnm = transpose(reshape([ &
                0.05, 49.0/180, 16.0/45, 49.0/180, 0.05, &
                -0.15, -0.35*sqrt(7.0/3), 0.0, 0.35*sqrt(7.0/3), 0.15, &
                0.25, 7.0/36, -8.0/9, 7.0/36, 0.25, &
                -0.35, 0.35*sqrt(7.0/3), 0.0, -0.35*sqrt(7.0/3), 0.35, &
                0.2, -7.0/15, 8.0/15, -7.0/15, 0.2], shape(Pnm)))
            return
        case default
            error stop 'Gauss_Lobatto_collocation_matrices: n is too large!'
            ! Actually, we can still implement the general case here, just without using the explicit expression for each coefficient. 
        end select
    end subroutine Gauss_Lobatto_collocation_matrices

end module collocation