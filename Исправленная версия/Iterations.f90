    !Итерации
    subroutine solve_iterations(tau,repeat_flag)
    use variables
    implicit none
    ! Входные параметры
    logical, intent(inout) :: repeat_flag
    real(dp), intent(inout) :: tau
    ! Выходные параметры - передаются через module
    !Локальные переменные
    integer :: k, j
    real(dp) :: error_ne, error_ni, error_E
    real(dp) :: n_e_m(0:M), n_i_m(0:M), E_r_m(0:M)
    real(dp) :: n_e_im(0:M), n_i_im(0:M), E_r_im(0:M)
    real(dp) :: n_e_m1(0:M), n_i_m1(0:M), E_r_m1(0:M)
    real(dp) :: max_error_e,max_error_i,max_err_E, ne_max,ni_max,E_max

    n_e_m = n_e_i
    n_i_m = n_i_i
    E_r_m = E_r_i

!f    tol = 1.0E-5
!f    iter = 10

    do k = 1, iter
        max_error_e = 0.0_dp; max_error_i = 0.0_dp; max_err_E = 0.0_dp
        do j = 0, M
            n_i_im(j) = sigma * n_i_m(j) + sigm1 * n_i_i(j)
            E_r_im(j) = sigma * E_r_m(j) + sigm1 * E_r_i(j)
        end do
        call solve_electron_continuity(tau,n_e_m1,n_i_im,E_r_im)
        do j = 0, M
            error_ne = abs(n_e_m1(j) - n_e_m(j))! / ne_max
            max_error_e= max(error_ne, max_error_e)
            n_e_im(j) = sigma * n_e_m1(j) + sigm1 * n_e_i(j)        !fav
        end do
        ne_max = maxval(abs(n_e_m))
        max_error_e = max_error_e/ne_max
!n_e_im(j) = sigma * n_e_m(j) + sigm1 * n_e_i(j)
        call solve_ion_continuity(tau,n_i_m1, n_e_im, E_r_im)
        ni_max = maxval(abs(n_i_m))
        do j = 0, M
            error_ni = abs(n_i_m1(j) - n_i_m(j))! / ni_max
!            error_E =  abs(E_r_m1(j) - E_r_m(j))! / E_max
            max_error_i = max(error_ni, max_error_i)!, error_E)!)
!            n_i_im(j) = sigma * n_i_m1(j) + sigm1 * n_i_i(j)        !fav
        end do
        max_error_i = max_error_i/ni_max
        call solve_poisson_trapezoidal(n_e_m1, n_i_m1, E_r_m1)
!f  E_max = maxval(abs(E_r_m1))
!f           do j = 0, M
!f            error_E =  abs(E_r_m1(j) - E_r_m(j))! / E_max
!f            max_err_E = max(error_E, max_err_E)!, error_E)!)
!f        end do
       do j = 0, M
            n_e_m(j) = n_e_m1(j)
            n_i_m(j) = n_i_m1(j)
            E_r_m(j) = E_r_m1(j)
        end do
        if (max_error_e + max_error_i < tol) then
            n_e_final = n_e_m1
            n_i_final = n_i_m1
            E_r_final = E_r_m1
            if (k<=k_maxit) tau = tau*1.1_dp
            repeat_flag = .false.
            exit
        end if
        if (k == iter) then
            tau = tau/2.0_dp
            print *, 'tau = ', tau  !; READ(*,*)
            repeat_flag = .true.
            if (tau < 1.0e-30) then
                print *, 'ERRROR: tau is too low'
            end if 
        end if
    end do
    end subroutine solve_iterations
