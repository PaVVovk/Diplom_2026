    !Итерации
    subroutine solve_iterations(tau,repeat_flag, big_data)
    use variables
    use functions
    use bolsig_data
    implicit none
    ! Входные параметры
    logical, intent(inout) :: repeat_flag
    real(dp), intent(inout) :: tau
    ! Выходные параметры - передаются через module
    !Локальные переменные
    integer :: k, j
    real(dp) :: error_ne, error_ni, EN_m, EN_im
    real(dp) :: n_e_m(0:M), n_i_m(0:M), E_r_m(0:M-1)
    real(dp) :: n_e_im(0:M), n_i_im(0:M), E_r_im(0:M-1)
    real(dp) :: n_e_m1(0:M), n_i_m1(0:M), E_r_m1(0:M-1)
    real(dp) :: max_error_e, max_error_i, ne_max, ni_max
    real(dp) :: D_e_i, D_i_i, k_e_i, k_i_i, nu_ion_i, E_mean_i 
    type(bolsig_data_t), intent(in) :: big_data
	
	D_e_i = D_e; D_i_i = D_i
	k_e_i = k_e; k_i_i = k_i
	nu_ion_i = nu_ion
	E_mean_i = E_mean
     	
    n_e_m = n_e_i
    n_i_m = n_i_i
    E_r_m = E_r_i


    do k = 1, iter
        max_error_e = 0.0_dp; max_error_i = 0.0_dp
        do j = 0, M - 1 
            n_i_im(j) = sigma * n_i_m(j) + sigm1 * n_i_i(j)
            E_r_im(j) = sigma * E_r_m(j) + sigm1 * E_r_i(j)
        end do
        n_i_im(M) = sigma * n_i_m(M) + sigm1 * n_i_i(M)
        call solve_electron_continuity(tau,n_e_m1,n_i_im,E_r_im) 
        do j = 0, M
            error_ne = abs(n_e_m1(j) - n_e_m(j))
            max_error_e= max(error_ne, max_error_e)
            n_e_im(j) = sigma * n_e_m1(j) + sigm1 * n_e_i(j)        
        end do
        ne_max = maxval(abs(n_e_m))
        max_error_e = max_error_e/ne_max
        call solve_ion_continuity(tau,n_i_m1, n_e_im, E_r_im)
        ni_max = maxval(abs(n_i_m))
        do j = 0, M
            error_ni = abs(n_i_m1(j) - n_i_m(j))
            max_error_i = max(error_ni, max_error_i)
        end do
        max_error_i = max_error_i/ni_max
        call solve_poisson_trapezoidal(n_e_m1, n_i_m1, E_r_m1)
        
        if (mode == 0) then
			I = n_e_m1(0)*r_half_k(0)**2
			do j = 1, M - 1
				I = I + n_e_m1(j)*h_half_k(j)*(r_half_k(j)+r_half_k(j-1))
			end do
				I = e*k_e*pi*R*I
				EN_m = U_D/(L+I)
				EN_m = EN_m/N * 1.0E21_dp
				EN_im = sigma * EN_m + sigm1 * EN
				D_e = extrapolate_diff_e(EN_im, big_data) 
				k_e = extrapolate_mob_e(EN_im, big_data)
				k_i = (0.286_dp + 0.669_dp*EXP(-EN_im/179.5_dp) + 0.679_dp*EXP(-EN_im/1305_dp))*1.0d-4
				k_i = k_i*N_L/N 
				D_i = k_i*(k_b*T_gas/e)
				nu_ion = extrapolate_freq(EN_im, big_data)
				E_mean = extrapolate_mean_e(EN_im, big_data)
		end if
        
		do j = 0, M-1
            n_e_m(j) = n_e_m1(j)
            n_i_m(j) = n_i_m1(j)
            E_r_m(j) = E_r_m1(j)
        end do
        n_e_m(M) = n_e_m1(M)
        n_i_m(M) = n_i_m1(M)
		  
        if (max_error_e + max_error_i < tol) then
            n_e_final = n_e_m1
            n_i_final = n_i_m1
            E_r_final = E_r_m1
            
            if (mode == 0) EN = EN_m
            
            if (k<=k_maxit) tau = tau*1.1_dp
            repeat_flag = .false.
            return
        end if
    end do
    tau = tau/2.0_dp
    repeat_flag = .true.
    if (mode == 0) then
		D_e = D_e_i; D_i = D_i_i
		k_e = k_e_i; k_i = k_i_i
		nu_ion = nu_ion_i
		E_mean = E_mean_i
    end if
    end subroutine solve_iterations
