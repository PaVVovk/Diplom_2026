    function df_dr(f1, f2, delta)
     integer, parameter :: dp = kind(1.0d0)      !fav
       real(dp) :: f1, f2, delta, df_dr
        df_dr = (f2 - f1)/delta
        !df_dr = (delta2*(f2 - f1) / delta1  + delta1*(f3 - f2) / delta2) / ( delta1 + delta2)
    end function df_dr
