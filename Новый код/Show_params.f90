subroutine show_parameters(t_counter, time)
use variables
use functions
implicit none
    integer, intent(in) :: t_counter
    real(dp), intent(in) :: time
    real(dp) :: potential(0:M), j_e(0:M), j_i(0:M), difference_e(1:M-1), difference_i(1:M-1)
    integer :: k
    real(dp) :: eps_e(0:M-1), eps_i(0:M-1)
    character(len=20) :: time_name
    character*4 name
!
    write (time_name, '(F20.18)') time
    write (name, '(I4.4)') t_counter
    open(11, file = ch_nt0//'output_'//name//'.txt')!!!!, status = 'new')
!        open(11, file = 'potential_t='//trim(time_name)//'.txt', status = 'new')
!        open(12, file = 'j_e_t='//trim(time_name)//'.txt', status = 'new')
!        open(13, file = 'j_i_t='//trim(time_name)//'.txt', status = 'new')
        !open(14, file = 'n_e_t='//trim(time_name)//'.txt', status = 'new')
        !open(15, file = 'n_i_t='//trim(time_name)//'.txt', status = 'new')

        potential(0) = 0

        do k = 1, M
            potential(k) = potential(k-1) + (E_r_i(k) + E_r_i(k-1)) * h_k(k-1) / 2.0_dp
        end do

        j_e(0) = 0.0_dp
        j_i(0) = 0.0_dp
        eps_e = epsilon_e(E_r_i)
        eps_i = epsilon_i(E_r_i)
        do k = 1, M 
            j_e(k-1) = r_half_k(k-1)*(-(D_e*n_e_i(k) - D_e*n_e_i(k-1))/h_k(k-1) - &
            (eps_e(k-1)*n_e_i(k) + (1 - eps_e(k-1))*n_e_i(k-1))*k_e*E_r_i(k-1))
            j_i(k-1) = r_half_k(k-1)*(-(D_e*n_i_i(k) - D_e*n_i_i(k-1))/h_k(k-1) + &
            (eps_i(k-1)*n_i_i(k) + (1 - eps_i(k-1))*n_i_i(k-1))*k_i*E_r_i(k-1))
        end do

        !j_e(M) = -D_e * (n_e(M) - n_e(M-1)) / h_k(M-1) - k_e * n_e(M) * E_r(M)
        !j_i(M) = -D_i * (n_i(M) - n_i(M-1)) / h_k(M-1) + k_i * n_i(M) * E_r(M)
        write (11, *) 'r_k',ch,'n_e(',ch,'n_i(',ch,'potential',ch,'E_r',ch,'j_e',ch,'j_i'
        do k = 0, M - 1
        write (11, 11) r_k(k),ch,n_e_i(k),ch,n_i_i(k),ch,potential(k),ch,E_r_i(k),ch,j_e(k),ch,j_i(k)
!            write (11, '(2ES20.10)') r_k(k), potential(k)
!            write (12, '(2ES20.10)') r_k(k), j_e(k)
!            write (13, '(2ES20.10)') r_k(k), j_i(k)
        end do
        ! Линейной экстраполяцией по двум последним точкам вспомогательной сетки около стенки найти потоки на стенку.
        write (11, 11) r_k(M),ch,n_e_i(M),ch,n_i_i(M),ch,potential(M),ch,E_r_i(M)
11 format (1p,E15.6,6(a,E15.6))
        !do k = 0, M
            !write (14, '(2ES20.10)') r_k(k), n_e(k)
            !write (15, '(2ES20.10)') r_k(k), n_i(k)
        !end do

        close(11)
!        close(12)
!        close(13)
        !close(14)
        !close(15)

    end subroutine show_parameters
