program main
    use collocation
    implicit none
    integer, parameter :: n = 4
    real, dimension(0:n) :: qx, qw
    real, dimension(0:n) :: Tl, Tr
    real, dimension(0:n, 0:n) :: M, Minv, D, Pmn, Pnm
    call Gauss_Lobatto_collocation_matrices(n, qx, qw, M, Minv, D, Tl, Tr, Pmn, Pnm)
    print *, matmul(M, Minv)
end program main