program concentration
    use variables
    use bolsig_data
    use functions
    implicit none
    integer :: k, step_count, t_counter, bolsig_index
    integer, parameter :: T_count = 50
    real(dp) :: r0, delta, h_first, tau, E_all, k_ion
    real(dp) :: t_max,t,t0,dt,D_amb,nu_amb
    real(dp) :: t_out(0:T_count)
    logical :: repeat_flag
    type(bolsig_data_t) :: big_data
    
    ch = char(9)
    big_data = read_bolsig_file('collected_dataf.dat')
    
    ch_nt0 = '_Data\'
    k_maxit = 3; iter = 8; tol = 1.0E-6_dp
    beta_ei = 1.0d-9*1.0d-6  !fav 1 cm^3/s -> m^3/s
    const = 1.602176634E-19_dp/(2.0_dp*epsilon0)  !from 4\pi e = e/\epsilon_0
    k_ion = 2.517E-18_dp
    sigma = 0.5_dp
    
    ! Ввод параметров
    M = 101     !fav
    print *, 'Enter the number of intervals M. Now M = ',M     !fav
    read *, M
    r0 = 1.0_dp     !fav
    print *, 'Enter the radius of the tube r0. Now r0 = ',r0 ,' cm'    !fav
    read *, r0; r0 = r0/100.0_dp
    delta = -0.02_dp     !fav
    print *, 'Enter the DELTA thickening parameter. Now delta = ',delta     !fav
    read *, delta
    p = 100.0_dp     !fav
    print *, 'Enter the pressure in pascals. Now p = ',p     !fav
    read *, p
    mode = 0
    print *, 'Choose computation mode: 0 - connected power and resistor, 1 - constant E/N. Now mode is ', mode
    read *, mode
    bc_mode = 1
    print *, 'Choose boundary conditions: 0 - zero, 1 - effective. Now conditions are ', bc_mode
    read *, bc_mode
    
    
! Выделение памяти
    allocate(r_k(0:M));         allocate(r_half_k(0:M-1))
    allocate(h_k(0:M-1));       allocate(h_half_k(1:M-1))
    allocate(n_e_i(0:M));       allocate(n_i_i(0:M))
    allocate(E_r_i(0:M-1));       allocate(n_e_final(0:M))
    allocate(n_i_final(0:M));   allocate(E_r_final(0:M-1))
    allocate(a_km(1:M-1));      allocate(b_km(1:M-1))
    allocate(w_km(1:M-1))
    !Значения коэффициентов
    N_L = 1.01d5/(k_b*300.0_dp)     !fav
    N = p/(k_b*T_gas)
    print *, 'N = ', N
    
    I = 0.05_dp
    U_D = 900.0_dp
    R = 1.0e4_dp
    L = 0.1_dp
    bolsig_index = 51
    
    if (mode == 1) then
		EN = big_data%E_over_N(bolsig_index)
		E_all = EN * N / 1e21
		nu_ion = big_data%ionization_freq_N(bolsig_index)*N
		D_e = big_data%diffusion_N(bolsig_index)/N
		k_e = big_data%mobility_N(bolsig_index)/N
		E_mean = big_data%mean_energy(bolsig_index)*e
    else
		print *, 'EN_MAX = ', U_D/(L*N)*1.0E21_dp
		E_all = (U_D-I*R)/L
		EN = E_all/N*1e21
		nu_ion =  extrapolate_freq(EN,big_data)
		D_e = extrapolate_diff_e(EN,big_data)
		k_e = extrapolate_mob_e(EN,big_data)
		E_mean = extrapolate_mean_e(EN, big_data)
    end if
    
    k_i = (0.286_dp + 0.669_dp*EXP(-EN/179.5_dp) + 0.679_dp*EXP(-EN/1305_dp))*1.0d-4
    k_i = k_i*N_L/N
	D_i = k_i*(k_b*T_gas/e)
    ne0 = 2.4048_dp*I/(2*e*pi*k_e*E_all*r0**2*bessel_j1(2.4048_dp))
    
    	
    
    D_amb = (D_e*k_i + D_i*k_e)/(k_e+k_i); nu_amb = D_amb*(2.4048/r0)**2
    print *, 'EN = ', EN
    print *, 'ne0 = ', ne0
    print *, 'k_i,D_i = ', k_i,D_i 
    print *, 'k_e,D_e = ', k_e,D_e
    print *, 'nu_ion,b_ei = ',nu_ion,beta_ei
    print *, 'D_amb,nu_amb =',D_amb,nu_amb
    print *, 'Press any button to begin'
    READ (*,*)
!
    
!
! Генерация основной сетки
    r_k(0) = 0.0_dp
    if(delta == 0) then
        h_first = r0/M
    else
        h_first = (delta * r0) / ((1.0_dp + delta)**M - 1.0_dp)
    end if
    do k = 0, M-1
        h_k(k) = h_first * (1.0_dp + delta)**(k)
        r_k(k+1) = r_k(k) + h_k(k)
        r_half_k(k) = (r_k(k) + r_k(k+1)) / 2.0_dp
    end do
    do k = 1, M-1
        h_half_k(k) = (h_k(k-1) + h_k(k)) / 2.0_dp
        a_km(k) = 2*(r_half_k(k-1)/(r_half_k(k-1) + r_half_k(k)))/h_half_k(k)
        b_km(k) = 2*(r_half_k(k)/(r_half_k(k-1) + r_half_k(k)))/h_half_k(k)
        w_km(k) = b_km(k) / h_k(k) + a_km(k) / h_k(k-1)
    end do
!
! Вывод результатов
    print *, '========================================'
    print *, 'Grid parameters:'
    print *, 'M =', M
    print *, 'r0 =', r0
    print *, 'delta =', delta
    print *, 'Initial step h_first =', h_first
    print *, 'Examination: r(M) =', r_k(M), 'must be =', r0
    print *, 'Min step h_k:', minval(h_k)
    if (minval(h_k) < 1.0e-30) then
        print *, 'ATTENTION: very small step!'
    end if
    print *, '========================================'

    !Заполняем массив t_out
    t_max = 20_dp/nu_amb
    dt =  10.0_dp ** 0.1_dp            !fav
    t0 = t_max / dt**T_count            !fav
    t_out(0) = t0
    do k = 1, T_count
        t_out(k) = t_out(k-1) * dt                      
    end do
    tau = minval(h_k**2/(4.0_dp*D_e))
    call ambipolar(D_amb, nu_amb)
    
    t = 0
    step_count = 0
    t_counter = 0
    
    call show_parameters(t_counter, t)
    
    do while(t < t_max)
        repeat_flag = .false.
        step_count = step_count + 1

        if (mod(step_count, 10000) == 0) then
            print *, 'Step:', step_count, 't =', t, 'tau =', tau, 'EN =', EN 
        end if
        call solve_iterations(tau,repeat_flag,big_data)
        do while(repeat_flag .eqv. .true.)
            call solve_iterations(tau,repeat_flag,big_data)
        end do
        
        n_e_i = n_e_final
        n_i_i = n_i_final
        E_r_i = E_r_final
        
        t = t + tau
          
        if (t > t_out(t_counter-1)) then
            call show_parameters(t_counter, t)
        end if
    end do


    print *, '========================================'
    do k = 0, M, M  !fav
       print *, 'Electron concentration n_e =', n_e_i(k), ' 1/m^3'
    end do
    do k = 0, M, M  !fav
       print *, 'Ion concentration n_i =', n_i_i(k), ' 1/m^3'
    end do
    do k = 0, M-1, M-1  !fav
       print *, 'Field E_r = ', E_r_i(k), ' V/m'
    end do
    print *, 'Field E_z = ', EN, ' Td'
    deallocate(r_k, r_half_k, h_k, h_half_k, n_e_i, n_i_i, &
               E_r_i, n_e_final, n_i_final, E_r_final)
    deallocate(a_km,b_km,w_km)
end program concentration
