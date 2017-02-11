module constants
  use prec_def

  implicit none
  save
  integer, parameter :: NOD = 225, N_TERM = 125 ! NOD = number of disks
  integer, parameter :: MAX_ITER = 75000, nspeed = 30 ! max number of iterations

  real( pr ), parameter :: PI = 3.14159265358979323846264338327950288419716939937510_pr, ETA_M = PI/4
end module constants
