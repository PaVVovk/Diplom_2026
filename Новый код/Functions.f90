    module functions
        implicit none
    contains
    function df_dr(f1, f2, delta)
		integer, parameter :: dp = kind(1.0d0)      !fav
		real(dp) :: f1, f2, delta, df_dr
        df_dr = (f2 - f1)/delta
        !df_dr = (delta2*(f2 - f1) / delta1  + delta1*(f3 - f2) / delta2) / ( delta1 + delta2)
    end function df_dr
    
    function epsilon_e(E_r) result(eps_e)
		use variables
		integer :: k
		real(dp), intent(in) :: E_r(0:M-1)
		real(dp) :: eps_e(0:M-1)
 		do k = 0, M-1 
		    if (-k_e*E_r(k) < 0) then
		        eps_e(k) = 1.0_dp
		    else 
		        eps_e(k) = 0.0_dp
		    end if
		end do     
	end function epsilon_e	     
    
    function epsilon_i(E_r) result(eps_i)
        use variables
        integer :: k
        real(dp), intent(in) :: E_r(0:M-1)
		real(dp) :: eps_i(0:M-1)
        do k = 0, M-1
            if (k_i*E_r(k) < 0) then
		        eps_i(k) = 1.0_dp
		    else 
		        eps_i(k) = 0.0_dp
		    end if
		end do 
	end function epsilon_i
    end module functions	
