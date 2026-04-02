module text_transfer
contains

function ru_windos(string, win_dos)
	character(*), intent(in) :: string
	character (len(string)) :: ru_windos, ru_temp
	logical, intent(in) :: win_dos
	integer dos_win_code, dif
	ru_temp=string
	l=len_trim(ru_temp)
	do i=1,l
		dos_win_code=iachar(ru_temp(i:i))
		dif=0
		if (win_dos) then
			select case (dos_win_code)
			case (#C0:#EF)
				dif=-#40
			case (#F0:#FF)
				dif=-#10
			case (#A8)
				dif=#48
			case (#B8)
				dif=#39
			case (#B9)
				dif=#43
			end select
		else
			select case(dos_win_code)
			case(#80:#AF)
				dif=#40
			case(#E0:#EF)
				dif=#10
			case (#F0)
				dif=-#48
			case (#F1)
				dif=-#39
			case (#FC)
				dif=-#43
			end select
		endif
		if (dif/=0) ru_temp(i:i)=char(dos_win_code+dif)
	end do
	ru_windos=ru_temp
end function ru_windos

function ru_win(string)
	character(*), intent(in) :: string
	character (len(string)) :: ru_win, ru_temp
	integer dos_win_code, dif
	ru_temp=string
	l=len_trim(ru_temp)
	do i=1,l
		dos_win_code=iachar(ru_temp(i:i))
		dif=0
		select case(dos_win_code)
		case(#80:#AF)
			dif=#40
		case(#E0:#EF)
			dif=#10
		case (#F0)
			dif=-#48
		case (#F1)
			dif=-#39
		case (#FC)
			dif=-#43
		end select
		if (dif/=0) ru_temp(i:i)=char(dos_win_code+dif)
	end do
	ru_win=ru_temp
end function ru_win

function ru_dos(string)
	character(*), intent(in) :: string
	character (len(string)) :: ru_dos, ru_temp
	integer dos_win_code, dif
	ru_temp=string
	l=len_trim(ru_temp)
	do i=1,l
		dos_win_code=iachar(ru_temp(i:i))
		dif=0
		select case (dos_win_code)
		case (#C0:#EF)
			dif=#40
		case (#F0:#FF)
			dif=#10
		case (#A8)
			dif=-#48
		case (#B8)
			dif=-#39
		case (#B9)
			dif=-#43
		end select
		if (dif/=0) ru_temp(i:i)=char(dos_win_code-dif)
	end do
	ru_dos=ru_temp
end function ru_dos

function up_string(string)
	character(*), intent(in) :: string
	character (len(string)) :: up_string
	integer dos_win_code, dif
	up_string=string
	l=len_trim(string)
	do i=1,l
		dos_win_code=iachar(up_string(i:i))
		dif=0
		select case (dos_win_code)
		case (#61:#7A,#A0:#AF)
			dif=#20
		case(#E0:#EF)
			dif=#50
		end select
		if (dif/=0) up_string(i:i)=char(dos_win_code-dif)
	enddo
end function up_string

function low_string(string)
	character(*), intent(in) :: string
	character (len(string)) :: low_string
	integer dos_win_code, dif
	low_string=string
	l=len_trim(string)
	do i=1,l
		dos_win_code=iachar(low_string(i:i))
		dif=0
		select case (dos_win_code)
		case (#41:#5A,#80:#8F)
			dif=#20
		case(#90:#9F)
			dif=#50
		end select
		if (dif/=0) low_string(i:i)=char(dos_win_code+dif)
	enddo
end function low_string

function str_cmp(string1, string2)
	integer str_cmp
	character(*), intent(in) :: string1, string2
	character (len(string1)) :: s1
	character (len(string2)) :: s2
	s1=up_string(string1)
	s2=up_string(string2)
	if (s1==s2) then
		str_cmp=0
	elseif (s1>s2) then
		str_cmp=1
	else
		str_cmp=-1
	endif

end function str_cmp


end module text_transfer
