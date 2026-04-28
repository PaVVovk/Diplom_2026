!Прогонка для ионов
subroutine solve_ion_continuity(tau,n_i_m1,n_e_im,E_r_im)
    use variables
    implicit none
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: n_e_im(0:M), E_r_im(0:M)
    real(dp), intent(out) :: n_i_m1(0:M)

    integer :: k
    real(dp) :: A(1:M), B(0:M), C(0:M), F(0:M)
    real(dp) :: alpha(1:M), beta(1:M)
    real(dp) :: y(0:M)

    real(dp) :: a_k, b_k, u_k, w_k, sitau, ttt
!---------------------------------------------
    sitau = 1.0_dp/(sigma*tau)
    B(0) = 4.0_dp*D_i/(h_k(0)**2) - 2.0_dp*k_i*E_r_im(1)/h_k(0)
    C(0) = 1.0_dp/(sigma*tau) + beta_ei*n_e_im(0) + 4.0_dp*D_i/(h_k(0)**2)
    F(0) = n_i_i(0)/tau + nu_ion*n_e_im(0) + &
           sigm1*(B(0)*n_i_i(1) - (C(0) - sitau)*n_i_i(0))
    do k = 1, M-1
        a_k = a_km(k); b_k = b_km(k); u_k = u_km(k); w_k = w_km(k)
        A(k) = a_k*(D_i/h_k(k-1) + k_i*E_r_im(k-1)/2)
        B(k) = b_k*(D_i/h_k(k)) !- k_i*E_r_im(k+1)/2)
        C(k) = sitau + beta_ei*n_e_im(k) + w_k*D_i + b_k*k_i*E_r_im(k) !u_k*k_i*E_r_im(k)
        F(k) = n_i_i(k)/tau + nu_ion*n_e_im(k) + &
               sigm1*(A(k)*n_i_i(k-1) - &
               (C(k) - sitau)*n_i_i(k) + B(k)*n_i_i(k+1))
    end do
    !A(M) = gamma_i * l_i / h_k(M-1)
    !C(M) = 1.0_dp + gamma_i * l_i / h_k(M-1)
    A(M) = 0.0_dp
    C(M) = 1.0_dp
    B(M) = 0.0_dp
    F(M) = 0.0_dp

! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
        ttt = C(k) - alpha(k) * A(k)
        alpha(k+1) = B(k) / ttt
        beta(k+1) = (beta(k) * A(k) + F(k)) / ttt
    end do
    ttt = (beta(M) * A(M) + F(M)) / (C(M) - alpha(M) * A(M)) !beta(M+1)
    ! Обратный ход
    y(M) = ttt; n_i_m1(M) = y(M) / sigma
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
        n_i_m1(k) = y(k) / sigma
    end do

    end subroutine solve_ion_continuity
