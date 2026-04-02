program concentration
implicit none 
integer, parameter :: dp = kind(1.0d0)
    integer :: M, k, step_count
    real(dp) :: r0, delta, h_first, sigma, tau 
    real(dp), allocatable :: r_k(:), r_half_k(:), h_k(:), h_half_k(:)
    real(dp), allocatable :: n_e_i(:), n_i_i(:), E_r_i(:)
    real(dp), allocatable :: n_e_final(:), n_i_final(:), E_r_final(:)
    real(dp) :: gamma_e, l_e, gamma_i, l_i, D_e, D_i, k_e, k_i, nu_ion, beta_ei
    logical :: repeat_flag
    real(dp) :: e, p, N, k_b, T_gas, t_max, t
    ! Ввод параметров
    print *, 'Enter the number of intervals M:'
    read *, M
    print *, 'Enter the radius of the tube r0:'
    read *, r0
    print *, 'Enter the DELTA thickening parameter:'
    read *, delta
    print *, 'Enter a numeric parameter sigma'
    read *, sigma 
    print *, 'Enter the pressure in pascals'
    read *, p
    
    ! Выделение памяти
    allocate(r_k(0:M))
    allocate(r_half_k(0:M-1))
    allocate(h_k(0:M-1))
    allocate(h_half_k(1:M-1))
    allocate(n_e_i(0:M))
    allocate(n_i_i(0:M))
    allocate(E_r_i(0:M))
    allocate(n_e_final(0:M))
    allocate(n_i_final(0:M))
    allocate(E_r_final(0:M))
    
    !Элементарный заря
    
    e = 1.602176634E-19_dp
    
    !Значения коэффициентов
    gamma_e = 0.7104_dp
    gamma_i = 0.7104_dp
    T_gas = 300_dp
    k_b = 1.380649E-23_dp
    N = p/(k_b*T_gas)
    print *, 'N = ', N
    nu_ion = 0.8093E-16_dp*N
    !beta_ei = 10E-19_dp*(T_gas/300)**(-4.5_dp)*N
    beta_ei = 10E-12
    k_e = 0.7890E+24_dp/N
    k_i = 0.286_dp + 0.669_dp*EXP(-E/(N*179.5_dp)) + 0.679_dp*EXP(-E/(N*1305_dp))
    D_e = 0.6932E+25_dp/N
    D_i = k_i*T/e
    l_e = 1/(N*10E-20_dp)
    l_i = 1/(N*10E-20_dp)
    n_i_i = nu_ion/beta_ei
    n_e_i = n_i_i
    E_r_i = 0.0_dp
    t = 0
    t_max = 10E-4
    step_count = 0
  
    
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
    end do
      ! Вывод результатов
    print *, '========================================'
    print *, 'Grid parameters:'
    print *, 'M =', M
    print *, 'r0 =', r0
    print *, 'delta =', delta
    print *, 'Initial step h_first =', h_first
    print *, '========================================'
    print *, 'Main mesh (r_k):'
    do k = 0, M
        print '(I4, F12.6)', k, r_k(k)
    end do
    
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
    print *, '========================================'
    print *, 'Examination: r(M) =', r_k(M), 'must be =', r0
    
    print *, 'Min step h_k:', minval(h_k)
    if (minval(h_k) < 1e-30) then
        print *, 'ATTENTION: very small step!'
    end if
    
    tau = minval(h_k**2/(4*D_e))  
    do while(t < t_max)
        repeat_flag = .false.
        step_count = step_count + 1
    
        if (mod(step_count, 10) == 0) then
            print *, 'Step:', step_count, 't =', t, 'tau =', tau
        end if
        call solve_iterations(M, r_k, r_half_k, sigma, tau, gamma_e, l_e, gamma_i, l_i, &
                           D_e, D_i, k_e, k_i, nu_ion, beta_ei, &
                           n_e_i, n_i_i, E_r_i, &
                           n_e_final, n_i_final, E_r_final, &
                           repeat_flag, h_half_k, h_k)
        do while(repeat_flag .eqv. .true.)
            call solve_iterations(M, r_k, r_half_k, sigma, tau, gamma_e, l_e, gamma_i, l_i, &
                               D_e, D_i, k_e, k_i, nu_ion, beta_ei, &
                               n_e_i, n_i_i, E_r_i, &
                               n_e_final, n_i_final, E_r_final, &
                               repeat_flag,  h_half_k, h_k)
        end do
        !print *, 'n_e_final = ', n_e_final(M/2)
        !print *, 'n_i_final = ', n_i_final(M/2)
        !print *, 'E_r_final = ', E_r_final(M/2)
        n_e_i = n_e_final 
        n_i_i = n_i_final 
        E_r_i = E_r_final
        t = t + tau
    end do                                          
    

    
    print *, '========================================'
    do k = 0, M
       print *, 'Electron concentration n_e =', n_e_i(k)
    end do
    do k = 0, M
       print *, 'Ion concentration n_i =', n_i_i(k)
    end do
    do k = 0, M
       print *, 'Field E_r =', E_r_i(k)
    end do
    
    deallocate(r_k, r_half_k, h_k, h_half_k, n_e_i, n_i_i, &
               E_r_i, n_e_final, n_i_final, E_r_final)
contains 
    
    !Прогонка для электронов
    
    subroutine solve_electron_continuity(M, sigma, tau, h_k, h_half_k, r_k, &
                                      r_half_k, gamma_e, l_e, n_e_m1, D_e, & 
                                      nu_ion, beta_ei, k_e, n_i_im, &
                                      n_e_i, E_r_im)
    implicit none
    integer, intent(in) :: M
    real(dp), intent(in) :: sigma, tau
    real(dp), intent(in) :: h_k(0:M-1), h_half_k(1:M-1), r_k(0:M), r_half_k(0:M-1)
    real(dp), intent(in) :: n_e_i(0:M)
    real(dp), intent(in) :: n_i_im(0:M), E_r_im(0:M)
    real(dp), intent(in) :: gamma_e, l_e, D_e, nu_ion, beta_ei, k_e
    real(dp), intent(out) :: n_e_m1(0:M)
    
    integer :: k
    real(dp) :: A(0:M), B(0:M), C(0:M), F(0:M) 
    real(dp) :: alpha(0:M+1), beta(0:M+1)       
    real(dp) :: y(0:M)   
  
    real(dp) :: a_k, b_k, u_k, w_k
    
   
    B(0) = 4.0_dp*D_e/(h_k(0)**2) + 2.0_dp*k_e*E_r_im(1)/h_k(0)
    C(0) = 1.0_dp/(sigma*tau) - nu_ion + beta_ei * n_i_im(0) + 4.0_dp * D_e / (h_k(0)**2)
    F(0) = n_e_i(0)/tau + (1.0_dp - sigma) * (B(0)*n_e_i(1) - (C(0) - 1.0_dp/(sigma*tau))*n_e_i(0))

    
    do k = 1, M-1
        a_k = (r_half_k(k-1)/r_k(k)) * (1.0_dp/h_half_k(k))  
        b_k = (r_half_k(k+1)/r_k(k)) * (1.0_dp/h_half_k(k))    
        u_k = (b_k - a_k) / 2.0_dp
        w_k = b_k / h_k(k) + a_k / h_k(k-1)
        A(k) = a_k * (D_e/h_k(k-1) - k_e*E_r_im(k-1)/2.0_dp)
        B(k) = b_k * (D_e/h_k(k) + k_e*E_r_im(k+1)/2.0_dp)
        C(k) = 1.0_dp/(sigma*tau) - nu_ion + beta_ei * n_i_im(k) + w_k * D_e - u_k * k_e * E_r_im(k)
        F(k) = n_e_i(k)/tau + (1.0_dp - sigma) * (A(k)*n_e_i(k-1) - (C(k) - 1.0_dp/(sigma*tau))*n_e_i(k) + B(k)*n_e_i(k+1))
    end do
    	
    A(M) = gamma_e*l_e/h_k(M-1)
    C(M) = 1.0_dp + gamma_e * l_e / h_k(M-1)
    B(M) = 0.0_dp  ! B_M не используется
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
    
    !Прогонка для ионов 
		
    subroutine solve_ion_continuity(M, sigma, tau, h_k, h_half_k, r_k, &
                                      r_half_k, gamma_i, n_i_i, l_i, n_i_m1, D_i, & 
                                      nu_ion, n_e_im, E_r_im,beta_ei, k_i)
    integer, intent(in) :: M
    real(dp), intent(in) :: sigma, tau
    real(dp), intent(in) :: h_k(0:M-1), h_half_k(1:M-1), r_k(0:M), r_half_k(0:M-1)
    real(dp), intent(in) :: n_e_im(0:M), E_r_im(0:M), n_i_i(0:M)
    real(dp), intent(in) :: gamma_i, l_i, D_i, nu_ion, beta_ei, k_i
    real(dp), intent(out) :: n_i_m1(0:M)
    
    integer :: k
    real(dp) :: A(0:M), B(0:M), C(0:M), F(0:M) 
    real(dp) :: alpha(0:M+1), beta(0:M+1)       
    real(dp) :: y(0:M)   
  
    real(dp) :: a_k, b_k, u_k, w_k
    
    
    B(0) = 4.0_dp*D_i/(h_k(0)**2) - 2.0_dp*k_i*E_r_im(1)/h_k(0)
    C(0) = 1.0_dp/(sigma*tau) + beta_ei*n_e_im(0) + 4.0_dp*D_i/(h_k(0)**2)
    F(0) = n_i_i(0)/tau + nu_ion*n_e_im(0) + &
           (1-sigma)*(B(0)*n_i_i(1) - (C(0) - 1/(sigma*tau))*n_i_i(0))
    do k = 1, M-1
        a_k = (r_half_k(k-1)/r_k(k)) * (1.0_dp/h_half_k(k))  
        b_k = (r_half_k(k+1)/r_k(k)) * (1.0_dp/h_half_k(k))    
        u_k = (b_k - a_k) / 2.0_dp
        w_k = b_k / h_k(k) + a_k / h_k(k-1)
        A(k) = a_k*(D_i/h_k(k-1) + k_i*E_r_im(k-1)/2)
        B(k) = b_k*(D_i/h_k(k) - k_i*E_r_im(k+1)/2)
        C(k) = 1/(sigma*tau) + beta_ei*n_e_im(k) + w_k*D_i + u_k*k_i*E_r_im(k)
        F(k) = n_i_i(k)/tau + nu_ion*n_e_im(k) + &
               (1 - sigma)*(A(k)*n_i_i(k-1) - &
               (C(k) - 1/(sigma* tau))*n_i_i(k) + B(k)*n_i_i(k+1))
    end do
    A(M) = gamma_i * l_i / h_k(M-1)
    C(M) = 1.0_dp + gamma_i * l_i / h_k(M-1)
    B(M) = 0.0_dp
    F(M) = 0.0_dp   
    
    ! Прямой ход
    alpha(1) = B(0) / C(0)
    beta(1) = F(0) / C(0)
    do k = 1, M-1
        alpha(k+1) = B(k) / (C(k) - alpha(k) * A(k))
    end do
    
    do k =1, M
        beta(k+1) = (beta(k) * A(k) + F(k)) / (C(k) - alpha(k) * A(k))
    end do
    ! Обратный ход 
    y(M) = beta(M+1)
    do k = M-1, 0, -1
        y(k) = alpha(k+1) * y(k+1) + beta(k+1)
    end do
    
    do k = 0, M
        n_i_m1(k) = y(k) / sigma
    end do
       
    end subroutine solve_ion_continuity 
    
    !Решение уравнения Пуассона
    subroutine solve_poisson_trapezoidal(M, r_k, h_k, n_e_m1, n_i_m1, E_r_m1)
    implicit none
    integer, intent(in) :: M
    real(dp), intent(in) :: r_k(0:M), h_k(0:M-1)
    real(dp), intent(in) :: n_e_m1(0:M), n_i_m1(0:M)
    real(dp), intent(out) :: E_r_m1(0:M)
    
    integer :: k
    real(dp) :: F(1:M), const   
    const = 2.0_dp*3.141592653589793_dp*1.602176634E-19_dp
    E_r_m1(0) = 0
    E_r_m1(1) = const*h_k(0)*0.5_dp*((n_i_m1(1)-n_e_m1(1)) + & 
                 (n_i_m1(0) - n_e_m1(0)))
    F(1) = r_k(1)*E_r_m1(1)
    do k = 2, M
        F(k) = F(k-1) + const*h_k(k-1)*(r_k(k-1)*(n_i_m1(k-1)-n_e_m1(k-1))+ &
               r_k(k)*(n_i_m1(k) - n_e_m1(k)))
        E_r_m1(k) = F(k)/r_k(k)        
    end do  
    end subroutine solve_poisson_trapezoidal
    
    !Итерации 
    
    subroutine solve_iterations(M, r_k, r_half_k, sigma, tau, gamma_e, l_e, gamma_i, l_i, &
                           D_e, D_i, k_e, k_i, nu_ion, beta_ei, &
                           n_e_i, n_i_i, E_r_i, &
                           n_e_final, n_i_final, E_r_final, &
                           repeat_flag, h_half_k, h_k) 
    
    ! Входные параметры
    
    integer, intent(in) :: M
    real(dp), intent(in) :: sigma, r_k(0:M), r_half_k(0:M-1), h_half_k(1:M-1), h_k(0:M-1)
    real(dp), intent(in) :: gamma_e, l_e, gamma_i, l_i
    real(dp), intent(in) :: D_e, D_i, k_e, k_i, nu_ion, beta_ei
    real(dp), intent(in) :: n_e_i(0:M), n_i_i(0:M), E_r_i(0:M)
    logical, intent(inout) :: repeat_flag
    real(dp), intent(inout) :: tau
    
    ! Выходные параметры
    
    real(dp), intent(out) :: n_e_final(0:M), n_i_final(0:M), E_r_final(0:M)
    
    !Локальные переменные 
    
    integer :: iter, k, j
    real(dp) :: error_ne, error_ni, error_E
    real(dp) :: n_e_m(0:M), n_i_m(0:M), E_r_m(0:M)
    real(dp) :: n_e_im(0:M), n_i_im(0:M), E_r_im(0:M)
    real(dp) :: n_e_m1(0:M), n_i_m1(0:M), E_r_m1(0:M)
    real(dp) :: max_error, tol
     
    n_e_m = n_e_i
    n_i_m = n_i_i
    E_r_m = E_r_i
    
    tol = 10E-8
    iter = 10
    
    
    do k = 1, iter
        max_error = 0.0_dp
        do j = 0, M
            n_e_im(j) = sigma * n_e_m(j) + (1.0_dp - sigma) * n_e_i(j)
            n_i_im(j) = sigma * n_i_m(j) + (1.0_dp - sigma) * n_i_i(j)
            E_r_im(j) = sigma * E_r_m(j) + (1.0_dp - sigma) * E_r_i(j)  
        end do
        call solve_electron_continuity(M, sigma, tau, h_k, h_half_k, r_k, &
                                      r_half_k, gamma_e, l_e, n_e_m1, D_e, & 
                                      nu_ion, beta_ei, k_e, n_i_im, &
                                      n_e_i, E_r_im)  
  
        call solve_ion_continuity(M, sigma, tau, h_k, h_half_k, r_k, &
                                      r_half_k, gamma_i, n_i_i, l_i, n_i_m1, D_i, & 
                                      nu_ion, n_e_im, E_r_im,beta_ei, k_i) 
        call solve_poisson_trapezoidal(M, r_k, h_k, n_e_m1, n_i_m1, E_r_m1)
         
         do j = 0, M
            error_ne = abs(n_e_m1(j) - n_e_m(j)) / max(abs(n_e_m1(j)), 1.0_dp)
            error_ni = abs(n_i_m1(j) - n_i_m(j)) / max(abs(n_i_m1(j)), 1.0_dp)
            error_E =  abs(E_r_m1(j) - E_r_m(j)) / max(abs(E_r_m1(j)), 1.0_dp)
            max_error = max(error_ne, error_ni, error_E, max_error)
        end do
       ! print *, 'n_e_m1= ', n_e_m1(M/2)
        !print *, 'n_i_m1= ', n_i_m1(M/2)
        !print *, 'E_r_m1= ', E_r_m1(M/2)
       ! print *, 'error_ne= ', error_ne
        !print *, 'error_ni= ', error_ni
        !print *, 'error_E= ',  error_E
        !print *, 'max_error= ', max_error
        do j = 0, M
            n_e_m(j) = n_e_m1(j)
            n_i_m(j) = n_i_m1(j)
            E_r_m(j) = E_r_m1(j)
        end do
        if (max_error < tol) then
            n_e_final = n_e_m1 
            n_i_final = n_i_m1
            E_r_final = E_r_m1
            exit  
        end if     
        if (k == iter) then
            tau = tau/2 
            print *, 'tau = ', tau
            repeat_flag = .true.
        else 
            repeat_flag = .false.    
        end if  
    end do
    end subroutine solve_iterations               
end program concentration
