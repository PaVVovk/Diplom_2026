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
	
	function extrapolate_mean_e(E_in, b_data) result(res)
	    use variables
	    use bolsig_data
	    integer :: n_max, j
	    real(dp) :: x1,x2,y1,y2,a,b,E_in,res
	    type(bolsig_data_t), intent(in) :: b_data
	    n_max = b_data%n_points
	    if (E_in < b_data%E_over_N(1)) then
			res = 2.397_dp * e
			return
		end if
	    do j = 1, n_max - 1
	        if (b_data%E_over_N(j) <= E_in .AND. E_in <= b_data%E_over_N(j+1)) then
	            x2 = b_data%E_over_N(j+1)
				x1 = b_data%E_over_N(j)
				y2 = b_data%mean_energy(j+1)*e
				y1 = b_data%mean_energy(j)*e
				b = (y2 - y1) / (x2 - x1)
				a = y1-b*x1
				res = a + b*E_in
				exit
			end if 	
	    end do 
	end function extrapolate_mean_e
	
	function extrapolate_diff_e(E_in, b_data) result(res)
	    use variables
	    use bolsig_data
	    integer :: n_max, j
	    real(dp) :: x1,x2,y1,y2,a,b,E_in,res
	    type(bolsig_data_t), intent(in) :: b_data
	    n_max = b_data%n_points
	    if (E_in < b_data%E_over_N(1)) then
			res = 0.3050e25_dp/N
			return
		end if
	    do j = 1, n_max - 1
	        if (b_data%E_over_N(j) <= E_in .AND. E_in <= b_data%E_over_N(j+1)) then
	            x2 = b_data%E_over_N(j+1)
				x1 = b_data%E_over_N(j)
				y2 = b_data%diffusion_N(j+1)/N
				y1 = b_data%diffusion_N(j)/N
				b = (y2 - y1) / (x2 - x1)
				a = y1-b*x1
				res = a + b*E_in
				exit
			end if 	
	    end do 
	end function extrapolate_diff_e
	
	function extrapolate_mob_e(E_in, b_data) result(res)
	    use variables
	    use bolsig_data
	    integer :: n_max, j
	    real(dp) :: x1,x2,y1,y2,a,b,E_in, res
	    type(bolsig_data_t), intent(in) :: b_data
	    n_max = b_data%n_points
	    if (E_in < b_data%E_over_N(1)) then
			res = 0.1174e27_dp/N
		    return
	    end if
	    do j = 1, n_max - 1
	        if (b_data%E_over_N(j) <= E_in .AND. E_in <= b_data%E_over_N(j+1)) then
	            x2 = b_data%E_over_N(j+1)
				x1 = b_data%E_over_N(j)
				y2 = b_data%mobility_N(j+1)/N
				y1 = b_data%mobility_N(j)/N
				b = (y2 - y1) / (x2 - x1)
				a = y1-b*x1
				res = a + b*E_in 
				exit
			end if 	
	    end do
	end function extrapolate_mob_e 
	
	function extrapolate_freq(E_in, b_data) result(res)
	    use variables
	    use bolsig_data
	    integer :: n_max, j
	    real(dp) :: x1,x2,y1,y2,a,b,E_in,res
	    type(bolsig_data_t), intent(in) :: b_data
	    n_max = b_data%n_points
	    if (E_in < b_data%E_over_N(1)) then
			res = 0.0_dp
			return
	    end if
	    do j = 1, n_max - 1
	        if (b_data%E_over_N(j) <= E_in .AND. E_in <= b_data%E_over_N(j+1)) then
	            x2 = b_data%E_over_N(j+1)
				x1 = b_data%E_over_N(j)
				y2 = log(b_data%ionization_freq_N(j+1)*N)
				y1 = log(b_data%ionization_freq_N(j)*N)
				b = (y2 - y1) / (1/x2 - 1/x1)
				a = y2 - b/x2
				res = exp(a+b/E_in) 
				exit
			end if 	
		end do	
	end function extrapolate_freq 
	  
    end module functions	
