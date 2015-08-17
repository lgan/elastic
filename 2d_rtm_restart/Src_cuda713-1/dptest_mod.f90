! dot product test module
module dptest_mod

  implicit none

  contains

  subroutine dptest(oper, n_model, n_data, dot1, dot2)
    integer,            intent(in)  :: n_model, n_data
    real, dimension(2), intent(out) :: dot1, dot2
    interface
      function oper(adj, add, model, data) result(stat)
        logical, intent(in)     :: adj, add
        real*8, dimension(:)   :: model, data
        integer                 :: stat
      end function
    end interface 
    real*8, dimension(n_model)          :: model1, model2
    real*8, dimension(n_data)           :: data1, data2
    integer                           :: stat

    call random_number(model1)
    call random_number(data2)

! dot product test with add = false
    stat = oper(.false., .false., model1, data1)
    dot1(1) = dot_product(data1, data2)
    stat = oper(.true., .false., model2, data2)
    dot1(2) = dot_product(model1, model2)

! dot product test with add = true
    stat = oper(.false., .true., model1, data1)
    dot2(1) = dot_product(data1, data2)
    stat = oper(.true., .true., model2, data2)
    dot2(2) = dot_product(model1, model2)

  end subroutine

end module
