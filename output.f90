! Подпрограммы для вывода файлов

    subroutine create_folder(folder_name)
        character(len=*), intent(in) :: folder_name
        character(len=10), parameter :: MKDIR = 'mkdir '
        integer :: folder_status

        folder_status = system(MKDIR // folder_name)
        if (folder_status /= 0) then
            print *, 'mkdir: failed to create folder!'
        end if
    end subroutine

    character(len=40) function data_string()
        use functions
        integer :: dt(8)
        call date_and_time(values=dt)
        !character(len=2) :: day, month, minute, second
        write (data_string, "(A,A,A,A,A,I0,A,A,A,A,A,A)") "Data_", &
                            str(dt(3)), "-", str(dt(2)), "-", dt(1), "_", &
                            str(dt(5)), "-", str(dt(6)), "-", str(dt(7))
    end function data_string

