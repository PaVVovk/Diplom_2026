module variables
    integer, parameter :: dp = kind(1.0d0)
    integer :: M,k_maxit,iter
    real(dp) :: sigma, sigm1, const,tol,T_e
    real(dp), allocatable :: r_k(:), r_half_k(:), h_k(:), h_half_k(:)
    real(dp) :: l_e,l_i,D_e,D_i,k_e,k_i,nu_ion,p,N,beta_ei,N_L
    real(dp), parameter :: k_b = 1.380649E-23_dp
    real(dp), parameter :: e = 1.602176634E-19_dp, epsilon0 = 8.85E-12_dp
    real(dp), parameter :: gamma_e = 0.7104_dp
    real(dp), parameter :: gamma_i = 0.7104_dp
    real(dp), parameter :: T_gas = 300.0_dp
    real(dp), allocatable :: n_e_i(:), n_i_i(:), E_r_i(:)
    real(dp), allocatable :: n_e_final(:), n_i_final(:), E_r_final(:)
    real(dp), allocatable :: a_km(:),b_km(:),u_km(:),w_km(:)
    character ch    !fav
    character*6 ch_nt0

end module variables
!----------------------------------------------------------------
program concentration
    use variables
    use bolsig_data
    implicit none
    integer :: k, step_count, t_counter
    integer, parameter :: T_count = 50
    real(dp) :: r0, delta, h_first, tau, E_all, k_ion
    real(dp) :: t_max,t,t0,dt,EN,ne0,D_amb,nu_amb
    real(dp) :: t_out(0:T_count)
    logical :: repeat_flag
    type(bolsig_data_t) :: big_data
    
    big_data = read_bolsig_file('collected_data.dat')
    EN = big_data%E_over_N(1)
    ch_nt0 = '_Data\'
    k_maxit = 3; iter = 8; tol = 1.0E-6_dp
    beta_ei = 1.0d-8*1.0d-6  !fav 1 cm^3/s -> m^3/s
    ch = char(9)    !fav tabulation
!    const = 2.0_dp*3.141592653589793_dp*1.602176634E-19_dp  !from 4\pi e
    const = 1.602176634E-19_dp/(2.0_dp*epsilon0)  !from 4\pi e = e/\epsilon_0
    EN = 100.0_dp !fav E/N in Td = 1.0e-17 V cm^2 = 1.0e-21 V m^2
    EN = 40.37_dp; T_e = 6.087_dp*2.0_dp/3.0_dp
    k_e = 9.006E+23_dp; D_e = 7.419E+24_dp
    k_ion = 2.517E-18_dp
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
    !print *, 'Enter a numeric parameter sigma'
    !read *, sigma
    sigma = 0.5_dp
    p = 100.0_dp     !fav
    print *, 'Enter the pressure in pascals. Now p = ',p     !fav
    read *, p
! Выделение памяти
    allocate(r_k(0:M));         allocate(r_half_k(0:M-1))
    allocate(h_k(0:M-1));       allocate(h_half_k(1:M-1))
    allocate(n_e_i(0:M));       allocate(n_i_i(0:M))
    allocate(E_r_i(0:M));       allocate(n_e_final(0:M))
    allocate(n_i_final(0:M));   allocate(E_r_final(0:M))
    allocate(a_km(1:M-1));      allocate(b_km(1:M-1))
    allocate(u_km(1:M-1));      allocate(w_km(1:M-1))
    !Значения коэффициентов
    N_L = 1.01d5/(k_b*300.0_dp)     !fav
    N = p/(k_b*T_gas)
    print *, 'N = ', N
    nu_ion = k_ion *N!   !*1.0d6 !!!0.8093E-16_dp
    E_all = 1.0d-21*N*EN
    k_e = k_e/N !0.7890E+24_dp/N
!    k_i = (0.286_dp + 0.669_dp*EXP(-E_all/(N*179.5_dp)) + 0.679_dp*EXP(-E_all/(N*1305_dp)))*1e-4
    k_i = (0.286_dp + 0.669_dp*EXP(-EN/179.5_dp) + 0.679_dp*EXP(-EN/1305_dp))*1.0d-4    !fav
!fav   k_i is for normal conditions. k_i in cm^2/(V s) = 1.0e-4 m^2/(V s)
    k_i = k_i*N_L/N  !fav
    D_e = D_e/N    !0.6932d25/N
    D_i = k_i*(k_b*T_gas/e) !fav
    D_amb = (D_e*k_i + D_i*k_e)/(k_e+k_i); nu_amb = 5.78_dp*D_amb/r0**2
    l_e = 1.0_dp/(N*1E-20_dp)
    l_i = 1.0_dp/(N*1E-20_dp)
    ne0 = nu_ion/beta_ei
    n_i_i = ne0
    n_e_i = n_i_i
    print *, 'k_i,D_i,n_i_i = ', k_i,D_i,n_i_i(1)  !tmp
    print *, 'k_e,D_e,n_e_i = ', k_e,D_e,n_e_i(1)  !tmp
    print *, 'nu_ion,b_ei = ',nu_ion,beta_ei
    print *, 'D_amb,nu_amb =',D_amb,nu_amb
    print *, 'Press any button to begin'
    READ (*,*)
!
    E_r_i = 0.0_dp
    t = 0
    t_max = 1.0d-4
    step_count = 0
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
!       r_half_k(M) = r_k(M)    !fav
    do k = 1, M-1
        h_half_k(k) = (h_k(k-1) + h_k(k)) / 2.0_dp
        a_km(k) = (r_half_k(k-1)/r_k(k))/h_half_k(k)
        b_km(k) = (r_half_k(k)/r_k(k))/h_half_k(k)!(r_half_k(k+1)/r_k(k))/h_half_k(k)
        u_km(k) = (b_km(k) - a_km(k)) / 2.0_dp
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
    print *, '========================================'
    write (36,*) 'k',ch,'Main mesh (r_k):',ch,'grid steps h_k',ch,'Aux. grid (r_{k+1/2}',ch,'h_half_k', &
            ch,'a_k',ch,'b_k',ch,'u_k',ch,'w_k' !fav
    k = 0
    write (36,'(I4, 3(a, 1p,E14.6))') k,ch,r_k(k),ch,h_k(k),ch,r_half_k(k)
    do k = 1, M-1
        write (36,'(I4, 9(a, 1p,E14.6))') k,ch,r_k(k),ch,h_k(k),ch,r_half_k(k),ch,h_half_k(k),  &
           ch,a_km(k),ch,b_km(k),ch,u_km(k),ch,w_km(k)!fav
    end do
    k = M
    write (36,'(I4, (a, 1p,E14.6))') k,ch,r_k(k)

    close (36)

go to 123   !fav
    print *, '========================================'
    print *, 'Main grid steps h_k:'
    do k = 0, M-1
        print '(I4, F12.6)', k, h_k(k)
    end do

    print *, '========================================'
    print *, 'Auxiliary grid (r_{k+1/2}):'
    do k = 0, M-1
        print '(I4, F12.6)', k, r_half_k(k)
    end do

    print *, '========================================'
    print *, 'Auxiliary grid steps h_half_k:'
    do k = 1, M-1
        print '(I4, F12.6)', k, h_half_k(k)
    end do

    ! Проверка последнего узла
123 continue
    print *, '========================================'
    print *, 'Examination: r(M) =', r_k(M), 'must be =', r0

    print *, 'Min step h_k:', minval(h_k)
    if (minval(h_k) < 1.0e-30) then
        print *, 'ATTENTION: very small step!'
    end if

    !Заполняем массив t_out
    dt =  10.0_dp ** 0.1_dp            !fav
    t0 = t_max / dt**T_count            !fav
    t_out(0) = t0
    write (32,'(I4, a, 1p,E14.6)'), 0,ch,t_out(0)  !fav
    do k = 1, T_count
        t_out(k) = t_out(k-1) * dt                      !fav
        ! fav print *, t_out(k)
         write (32,'(I4, a, 1p,E14.6)'), k,ch,t_out(k)  !fav
    end do

    tau = minval(h_k**2/(4.0_dp*D_e))
    call ambipolar(D_amb, nu_amb, ne0)
    !do k = 0, M
        !print *, n_e_i(k), ch, n_i_i(k)
    !end do
    t_counter = 0
    call show_parameters(t_counter, t)
    t_counter = t_counter+1
    do while(t < t_max)
        repeat_flag = .false.
        step_count = step_count + 1

        if (mod(step_count, 10000) == 0) then
            print *, 'Step:', step_count, 't =', t, 'tau =', tau
        end if
        call solve_iterations(tau,repeat_flag)
        do while(repeat_flag .eqv. .true.)
            call solve_iterations(tau,repeat_flag)
        end do

        if (t > t_out(t_counter-1)) then
            call show_parameters(t_counter, t)
            t_counter = t_counter+1
        end if
        !print *, 'n_e_final = ', n_e_final(M/2)
        !print *, 'n_i_final = ', n_i_final(M/2)
        !print *, 'E_r_final = ', E_r_final(M/2)
        n_e_i = n_e_final
        n_i_i = n_i_final
        E_r_i = E_r_final
        t = t + tau
    end do
! fav go to 124


    print *, '========================================'
    do k = 0, M, M  !fav
       print *, 'Electron concentration n_e =', n_e_i(k)
    end do
    do k = 0, M, M  !fav
       print *, 'Ion concentration n_i =', n_i_i(k)
    end do
    do k = 0, M, M  !fav
       print *, 'Field E_r =', E_r_i(k)
    end do
124 continue
    deallocate(r_k, r_half_k, h_k, h_half_k, n_e_i, n_i_i, &
               E_r_i, n_e_final, n_i_final, E_r_final)
    deallocate(a_km,b_km,u_km,w_km)
end program concentration
