    !Решение уравнения Пуассона
    subroutine solve_poisson_trapezoidal(n_e_m1,n_i_m1,E_r_m1)
    use variables
    implicit none
    real(dp), intent(in) :: n_e_m1(0:M), n_i_m1(0:M)
    real(dp), intent(out) :: E_r_m1(0:M)

    integer :: k
    real(dp) :: F(1:M)
    E_r_m1(0) = 0
!SQ    E_r_m1(1) = const*h_k(0)*( (n_i_m1(1)-n_e_m1(1)) + (n_i_m1(0) - n_e_m1(0)) )
    E_r_m1(1) = const*h_k(0)*( (n_i_m1(1)-n_e_m1(1)) + (n_i_m1(0) - n_e_m1(0)) )/2.0_dp
     F(1) = r_k(1)*E_r_m1(1)
    do k = 2, M
        F(k) = F(k-1) + const*h_k(k-1)*(r_k(k-1)*(n_i_m1(k-1)-n_e_m1(k-1))+ &
               r_k(k)*(n_i_m1(k) - n_e_m1(k)))
        E_r_m1(k) = F(k)/r_k(k)
    end do
    end subroutine solve_poisson_trapezoidal
