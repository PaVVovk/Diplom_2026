subroutine show_parameters(t_counter, time)
use variables
use functions
implicit none
    integer, intent(inout) :: t_counter
    real(dp), intent(in) :: time
    real(dp) :: potential(0:M), j_e(0:M), j_i(0:M)
    integer :: k
    real(dp) :: eps_e(0:M-1), eps_i(0:M-1)
    character*4 name

    write (name, '(I4.4)') t_counter
    open(11, file = ch_nt0//'output_'//name//'.txt')
	
	write(11, *) 'Time: t = ', time, ' s'
	
        potential(0) = 0

        do k = 0, M-1 
            potential(k+1) = potential(k) + E_r_i(k) * h_k(k)
        end do

        eps_e = epsilon_e(E_r_i)
        eps_i = epsilon_i(E_r_i)
        do k = 1, M 
            j_e(k-1) = -e*(-D_e*(n_e_i(k) - n_e_i(k-1))/h_k(k-1) - &
            (eps_e(k-1)*n_e_i(k) + (1 - eps_e(k-1))*n_e_i(k-1))*k_e*E_r_i(k-1))
            j_i(k-1) = e*(-D_i*(n_i_i(k) - n_i_i(k-1))/h_k(k-1) + &
            (eps_i(k-1)*n_i_i(k) + (1 - eps_i(k-1))*n_i_i(k-1))*k_i*E_r_i(k-1))
        end do 
        
        write (11, *) 'r_k',ch,'n_e',ch,'n_i',ch,'potential',ch,'r_half_k',ch,'E_r',ch,'j_e',ch,'j_i'
        write (11, 11) r_k(0),ch,n_e_i(0),ch,n_i_i(0),ch,potential(0),ch,0.0_dp,ch,0.0_dp,ch,0.0_dp,ch,0.0_dp
        do k = 1, M
			write (11, 11) r_k(k),ch,n_e_i(k),ch,n_i_i(k),ch,potential(k),ch,r_half_k(k-1),ch,E_r_i(k-1),ch,j_e(k-1),ch,j_i(k-1)
        end do
        11 format (1p,E15.6,7(a,E15.6))
        
        close(11)
		t_counter = t_counter + 1
		
    end subroutine show_parameters

