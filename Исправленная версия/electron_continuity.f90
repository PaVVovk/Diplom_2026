    !Прогонка для электронов

    subroutine solve_electron_continuity(tau,n_e_m1,n_i_im,E_r_im)
    use variables
    implicit none
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: n_i_im(0:M), E_r_im(0:M)
    real(dp), intent(out) :: n_e_m1(0:M)

    integer :: k
    real(dp) :: A(0:M), B(0:M), C(0:M), F(0:M)
    real(dp) :: alpha(0:M+1), beta(0:M+1)
    real(dp) :: y(0:M)

    real(dp) :: a_k, b_k, u_k, w_k, sitau
!----------------------------------------------------------------------------------------
    sitau = 1.0_dp/(sigma*tau)
!
    B(0) = 4.0_dp*D_e/(h_k(0)**2) + 2.0_dp*k_e*E_r_im(1)/h_k(0)
    C(0) = sitau - nu_ion + beta_ei * n_i_im(0) + 4.0_dp * D_e / (h_k(0)**2)
    F(0) = n_e_i(0)/tau + sigm1 * (B(0)*n_e_i(1) - (C(0) - sitau)*n_e_i(0))


    do k = 1, M-1
        a_k = a_km(k); b_k = b_km(k); u_k = u_km(k); w_k = w_km(k)
        A(k) = a_k * (D_e/h_k(k-1) - k_e*E_r_im(k-1)/2.0_dp)
        B(k) = b_k * (D_e/h_k(k) + k_e*E_r_im(k+1)/2.0_dp)
        C(k) = sitau - nu_ion + beta_ei * n_i_im(k) + w_k * D_e - u_k * k_e * E_r_im(k)
        F(k) = n_e_i(k)/tau + sigm1 * (A(k)*n_e_i(k-1) - (C(k) - sitau)*n_e_i(k) + B(k)*n_e_i(k+1))
    end do

    !A(M) = gamma_e*l_e/h_k(M-1)
    !C(M) = 1.0_dp + gamma_e * l_e / h_k(M-1)
    A(M) = 0.0_dp
    C(M) = 1.0_dp
    B(M) = 0.0_dp
    F(M) = 0.0_dp

    ! Метод прогонки

    ! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
    !Сохранить знаменатель в отдельную перем
        alpha(k+1) = B(k) / (C(k) - alpha(k) * A(k))
    end do

    do k = 1, M
        beta(k+1) = (beta(k) * A(k) + F(k)) / (C(k) - alpha(k) * A(k))
    end do
    ! Обратный ход
    y(M) = beta(M+1)
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
    end do

    do k = 0, M
        n_e_m1(k) = y(k) / sigma
    end do

    end subroutine solve_electron_continuity
