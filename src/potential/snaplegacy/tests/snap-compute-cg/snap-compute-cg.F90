! DEBUG_SNAP_CG_PRG: compile or not program
! => this allow to call compute_cg_vector from another program
#ifdef DEBUG_SNAP_CG_PRG
program dummy
   implicit none

   character(len=3) :: arg
   double precision :: j_max
   integer :: ntype

   call get_command_argument(1, arg)
   write(*,*) "snap-compute-cg-f90: 1st arg ", arg
   read(arg,*) j_max
   write(*,*) "snap-compute-cg-f90: j_max ", j_max

   call get_command_argument(2, arg)
   write(*,*) "snap-compute-cg-f90: 2d  arg ", arg
   read(arg,*) ntype
   write(*,*) "snap-compute-cg-f90: ntype ", ntype

   call compute_cg_vector(j_max, ntype)
end program dummy
#endif

function cg(j1,m1,j2,m2,j,m,ntype)
   implicit none

   integer,intent(in) :: j1,j2,j,m1,m2,m,ntype
   double precision :: cg
   double precision ::  wigner3j

   cg = (-1.d0)**((j1-j2+m)/ntype) * dsqrt(2.d0*dble(j)/dble(ntype)+1.d0) * wigner3j(j1,m1,j2,m2,j,-m,ntype)

end function cg


function wigner3j(j1,m1,j2,m2,j,m,ntype)
   implicit none

   integer,intent(in) :: j1,j2,j,m1,m2,m, ntype
   double precision ::  wigner3j
   integer :: t,tmin,tmax
   double precision :: factorial, triangle,coef,sum_coef

   coef = (-1.d0)**((j1-j2-m)/ntype) * dsqrt(factorial((j1+m1)/ntype)*factorial((j1-m1)/ntype)* &
                                                     factorial((j2+m2)/ntype)*factorial((j2-m2)/ntype)* &
                                                     factorial((j+m)/ntype)*factorial((j-m)/ntype))

   triangle = dsqrt(factorial((j1+j2-j)/ntype)*factorial((j1-j2+j)/ntype)* &
               factorial((-j1+j2+j)/ntype)/factorial((j1+j2+j)/ntype+1))

   sum_coef=0.d0
   tmin = max((j2-j-m1)/ntype,(j1+m2-j)/ntype,0)
   tmax = min((j1+j2-j)/ntype,(j1-m1)/ntype,(j2+m2)/ntype)
   do t=tmin,tmax
      sum_coef = sum_coef + (-1.d0)**t / (                   &
      factorial(t)*factorial((j1+j2-j)/ntype-t)*factorial((j1-m1)/ntype-t)* &
                   factorial((j2+m2)/ntype-t)*factorial((j-j2+m1)/ntype+t)*factorial((j-j1-m2)/ntype+t) )
   enddo

   wigner3j = coef * triangle * sum_coef

end function wigner3j


function factorial (n) result(res)

   !use ml_in_ndm_module, only : rangml
   implicit none

   integer,intent(in) :: n
   double precision :: res
   integer :: i

   !if ((n<0).and.(rangml==0)) then
   !   stop "Argument of factorial function should be positive"
   !elseif (n==0) then
   !   res=1
   !endif

   res=1.d0
   do i=2,n
      res = res*i
   enddo

end function factorial




subroutine compute_cg_vector(j_max,ntype)
!use ml_in_ndm_module, only : rangml, cg_vector

implicit none
double precision, intent(in) :: j_max
integer, intent(in) ::  ntype
double precision, allocatable, dimension(:,:,:,:,:,:) :: cg_vector
double precision :: cg

integer :: j, j1, j2, m, m1,m2, j0_max, j1_max, j2_max

!if ((ntype==1) .and. dabs(int(j_max)-j_max).gt.0.1d0) then
!      if (rangml==0) then
!        write(*,*) 'Probably you intend to use some SO3 or bi-SO3 -  related descriptor'
!        write(*,*) 'in that case j_max can be only integer, now j_max', j_max
!      end if
!      stop
!end if

open(unit=12, file="cg.ref", action="write")

j0_max=int(j_max*ntype)
j1_max=int(j_max*ntype)
j2_max=int(j_max*ntype)

if (allocated(cg_vector)) deallocate(cg_vector) ; allocate(cg_vector( 0:j1_max,-j1_max:j1_max, &
                                                                         0:j2_max,-j2_max:j2_max, &
                                                                         0:j0_max,-j0_max:j0_max ))
 cg_vector(:,:,:,:,:,:)=0.d0
do j1=0,j1_max
   do j2=0,j2_max
      do j =0,j0_max
         do m1=-j1, j1,ntype
            do m2=-j2, j2,ntype
               do m=-j, j,ntype
                   if (j1+j2-j <0) cycle
                   if (j-j1+j2 <0) cycle
                   if (j-j2+j1 <0) cycle
                   if (mod(j+j1+j2,ntype)==1) cycle
                   if (m-m1-m2/=0) then
!                      write(12,'(6i5, e20.10)') j1,j2,j,m1,m2, m, cg_vector(j1,m1,j2,m2,j,m) 
!                      write(6,'(a, 6i5, e20.10)') "f90:",j1,j2,j,m1,m2,m,cg_vector(j1,m1,j2,m2,j,m)
                      cycle
                   endif
                   cg_vector(j1,m1,j2,m2,j,m) =  cg(j1,m1,j2,m2,j,m, ntype)
!                   write(12,'(6i5, e20.10)') j1,j2,j,  m1,m2, m, cg_vector(j1,m1,j2,m2,j,m)
                   write(12,'(6i5, e20.10)') j,j1,j2,m,m1,m2, cg_vector(j1,m1,j2,m2,j,m)
#ifdef DEBUG_SNAP_CG
! Debug traces added in cmake logs: convenient for debugging failed tests
!                   write(6,'(a, 6i5, e20.10)') "f90:",j1,j2,j,  m1,m2, m, cg_vector(j1,m1,j2,m2,j,m)
                   write(6,'(a, 6i5, e20.10)') "f90:",j,j1,j2,m,m1,m2, cg_vector(j1,m1,j2,m2,j,m)
#endif
               end do
            end do
         end do
      end do
   end do
end do
!stop 'test cg in compute cg vector'
return
end subroutine  compute_cg_vector
