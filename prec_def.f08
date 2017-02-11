module prec_def
  implicit none
  save
  integer, parameter:: sp = kind( 1.0 ), &
       dp = SELECTED_REAL_kind( 2 * precision( 1._sp ) ), &
       qp = SELECTED_REAL_kind( 2 * precision( 1._dp ) ), &
       pr = sp, & ! precision in use
       tpr = dp ! precision for time variables: 
  ! it must be larger than pr in order to have a smaller truncation error
  ! with respect to the other variables 
end module prec_def

