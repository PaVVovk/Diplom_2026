!Прогонка для ионов
subroutine solve_ion_continuity(tau,n_i_m1,n_e_im,E_r_im)
    use variables
    use functions
    implicit none
    real(dp), intent(in) :: tau
    real(dp), intent(in) :: n_e_im(0:M), E_r_im(0:M-1)
    real(dp), intent(out) :: n_i_m1(0:M)

    integer :: k
    real(dp) :: A(1:M), B(0:M), C(0:M), F(0:M)
    real(dp) :: alpha(1:M), beta(1:M)
    real(dp) :: y(0:M), eps_i(0:M-1)

    real(dp) :: a_k, b_k, w_k, sitau, denom, v_i, l_i
!---------------------------------------------
    sitau = 1.0_dp/(sigma*tau)
    eps_i = epsilon_i(E_r_im)
    
    v_i = sqrt(3*8.31_dp*T_gas/mu)
    l_i = 3.0_dp*D_i/v_i
    
    B(0) = 4.0_dp*D_i/(h_k(0)**2) - 2.0_dp*(1-eps_i(0))*k_i*E_r_im(0)/r_half_k(0)
    C(0) = 1.0_dp/(sigma*tau) + beta_ei*n_e_im(0) + 4.0_dp*D_i/(h_k(0)**2) + 2.0_dp*eps_i(0)*k_i*E_r_im(0)/r_half_k(0)
    F(0) = n_i_i(0)/tau + nu_ion*n_e_im(0) + &
           sigm1*(B(0)*n_i_i(1) - (C(0) - sitau)*n_i_i(0))
    do k = 1, M-1
        a_k = a_km(k); b_k = b_km(k); w_k = w_km(k)
        A(k) = a_k*(D_i/h_k(k-1) + (1-eps_i(k-1))*k_i*E_r_im(k-1))
        B(k) = b_k*(D_i/h_k(k) - eps_i(k)*k_i*E_r_im(k))
        C(k) = sitau + beta_ei*n_e_im(k) + w_k*D_i + k_i*(b_k*(1-eps_i(k))*E_r_im(k) - a_k*eps_i(k-1)*E_r_im(k-1))
        F(k) = n_i_i(k)/tau + nu_ion*n_e_im(k) + &
               sigm1*(A(k)*n_i_i(k-1) - &
               (C(k) - sitau)*n_i_i(k) + B(k)*n_i_i(k+1))
    end do
    
    if (bc_mode == 1) then
		A(M) = gamma_i*l_i/h_k(M-1)
		C(M) = 1.0_dp + gamma_i * l_i / h_k(M-1)
	else
		A(M) = 0.0_dp
		C(M) = 1.0_dp
    end if
    B(M) = 0.0_dp
    F(M) = 0.0_dp

! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
        denom = C(k) - alpha(k) * A(k)
        alpha(k+1) = B(k) / denom
        beta(k+1) = (beta(k) * A(k) + F(k)) / denom
    end do
    denom = (beta(M) * A(M) + F(M)) / (C(M) - alpha(M) * A(M))
    ! Обратный ход
    y(M) = denom; n_i_m1(M) = y(M) / sigma
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
        n_i_m1(k) = y(k) / sigma
    end do

    end subroutine solve_ion_continuity
