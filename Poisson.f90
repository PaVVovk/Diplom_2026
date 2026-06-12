    !Решение уравнения Пуассона
    subroutine solve_poisson_trapezoidal(n_e_m1,n_i_m1,E_r_m1)
    use variables
    implicit none
    real(dp), intent(in) :: n_e_m1(0:M), n_i_m1(0:M)
    real(dp), intent(out) :: E_r_m1(0:M-1)

    integer :: k
    real(dp) :: F(0:M-1)
    E_r_m1(0) = const*(n_i_m1(0) - n_e_m1(0))*r_half_k(0)
    F(0) = r_half_k(0)*E_r_m1(0)
    do k = 1, M-1
        F(k) = F(k-1) + const*(n_i_m1(k) - n_e_m1(k))*h_half_k(k)*(r_half_k(k)+r_half_k(k-1))
        E_r_m1(k) = F(k)/r_half_k(k)
    end do
    end subroutine solve_poisson_trapezoidal
