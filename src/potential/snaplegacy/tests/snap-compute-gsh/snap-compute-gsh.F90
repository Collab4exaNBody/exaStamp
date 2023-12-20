! DEBUG_SNAP_GSH_PRG: compile or not program
! => this allow to call spherical_4d from another program
#ifdef DEBUG_SNAP_GSH_PRG
program dummy
   implicit none

   character(len=7) :: arg
   real*8 :: r, rcut
   real*8, dimension(3) :: x
   integer :: jmax, j, m1, m2
   complex*16, dimension(:,:,:), allocatable   :: Umm
   complex*16, dimension(:,:,:,:), allocatable :: dUmm

   call get_command_argument(1, arg)
   write(*,*) "snap-compute-gsh-f90: 1st arg ", arg
   read(arg,*) x(1)
   write(*,*) "snap-compute-gsh-f90: x ", x(1)

   call get_command_argument(2, arg)
   write(*,*) "snap-compute-gsh-f90: 2d  arg ", arg
   read(arg,*) x(2)
   write(*,*) "snap-compute-gsh-f90: y ", x(2)

   call get_command_argument(3, arg)
   write(*,*) "snap-compute-gsh-f90: 3d  arg ", arg
   read(arg,*) x(3)
   write(*,*) "snap-compute-gsh-f90: z ", x(3)

   r = dsqrt(x(1)*x(1) + x(2)*x(2) + x(3)*x(3))

   call get_command_argument(4, arg)
   write(*,*) "snap-compute-gsh-f90: 4th arg ", arg
   read(arg,*) rcut
   write(*,*) "snap-compute-gsh-f90: rcut ", rcut

   call get_command_argument(5, arg)
   write(*,*) "snap-compute-gsh-f90: 5th arg ", arg
   read(arg,*) jmax
   write(*,*) "snap-compute-gsh-f90: jmax ", jmax
   jmax = 2*jmax ! Needed for spherical_4d j=2J (j f90 == J Bartok).

   allocate( Umm(-jmax:jmax,-jmax:jmax,0:jmax))
   allocate(dUmm(-jmax:jmax,-jmax:jmax,0:jmax,3))

   call spherical_4d(x, r, Umm, dUmm, rcut, jmax)

   open(unit=12, file="gsh.ref", action="write")
   open(unit=13, file="dgsh.ref", action="write")
   do j=0,jmax
      do m1=-j,j,2
         do m2=-j,j,2
            write(12,'(3i5, 2e20.10)') j,m2,m1, real(Umm(m2,m1,j)),imag(Umm(m2,m1,j))
            write(13,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,1)),imag(dUmm(m2,m1,j,1))
            write(13,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,2)),imag(dUmm(m2,m1,j,2))
            write(13,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,3)),imag(dUmm(m2,m1,j,3))
#ifdef DEBUG_SNAP_GSH
! Debug traces added in cmake logs: convenient for debugging failed tests
            write(6,'(3i5, 2e20.10)') j,m2,m1, real(Umm(m2,m1,j)),imag(Umm(m2,m1,j))
            write(6,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,1)),imag(dUmm(m2,m1,j,1))
            write(6,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,2)),imag(dUmm(m2,m1,j,2))
            write(6,'(3i5, 2e20.10)') j,m2,m1, real(dUmm(m2,m1,j,3)),imag(dUmm(m2,m1,j,3))
#endif
         end do
      end do
   end do
end program dummy
#endif

   subroutine spherical_4d(x, rr, Umm, dUmm, r_cut, jj_max)
   ! Input:
   !    x(3) -> x, y, z ;
   !    rr=sqrt(x**2 + y**2 + z**2)
   ! Output:
   !    Spherical 4D and derivatives
   !    Umm, dUmm
      !use ml_in_ndm_module, only : r0,jj_max, one_pi,r_cut
      implicit none
!      real(kind=kind(1.d0)),dimension(3),intent(in) :: x
!      real(kind=kind(1.d0)),intent(in) :: rr
!      double precision,dimension(3),intent(in) :: x
!      double precision,intent(in) :: rr
      real*8,dimension(3),intent(in) :: x
      real*8,intent(in) :: rr
      double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max), intent(out) :: Umm
      double complex,dimension(-jj_max:jj_max,-jj_max:jj_max,0:jj_max,3), intent(out):: dUmm
!      double precision :: r_cut
      real*8, intent(in)     :: r_cut
      integer, intent(in) :: jj_max

!      double precision :: r0, one_pi
      real*8 :: r0, one_pi
      integer :: j,m1,m2
      !real(kind=kind(1.d0)) :: z0,  th0, l0, l0i, dz, dil0
      !real(kind=kind(1.d0)),   dimension(3) :: drr
!      double precision :: z0,  th0, l0, l0i, dz, dil0
!      double precision,   dimension(3) :: drr
      real*8 :: z0,  th0, l0, l0i, dz, dil0
      real*8,   dimension(3) :: drr
      double complex :: cplx_l0i, cplx_i, z_plus,z_minus,x_plus,x_minus, tjm1, tjm2
      double complex,dimension(3) :: dz_plus,dz_minus,dx_plus,dx_minus
      real*8 :: const
      one_pi = 3.14159265359
!      const=2.d-2
      cplx_i=cmplx(0.d0,1.d0, kind=kind(1.d0))
      r0 = r_cut/(one_pi-0.02)
      th0=rr/r0
      l0= rr/dsin(th0)
      l0i= dsin(th0)/rr
      z0=rr/dtan(th0)
!      write(*,*) "theta0:", th0
!      write(*,*) "l0: ", l0
!      write(*,*) "r0: ", r0

      ! z_+/- is cos(theta0) +/- i sin(theta0) * cos(theta) ; z+=u*;z-=u
      z_plus= cmplx(z0, x(3), kind=kind(1.d0))/l0 
      z_minus=cmplx(z0,-x(3), kind=kind(1.d0))/l0 
      ! x +/- is sin(theta) * sin(theta0) * exp(+/- i \phi) x+= v*exp(i phi) ; x- = v*exp(-i phi)
      x_plus= cmplx(x(1), x(2), kind=kind(1.d0))/l0
      x_minus=cmplx(x(1),-x(2), kind=kind(1.d0))/l0


      !Here are the derivatives ...
      drr(1:3)=x(1:3)/rr

      !error dz  = (1.d0/tan(th0) - th0/tan(th0)**2)
      dz  = (1.d0/dtan(th0) - th0/dsin(th0)**2)
      dil0= (dcos(th0)/r0 - l0i)/rr
      cplx_l0i=cmplx(0.d0,l0i, kind=kind(1.d0))

      dz_plus(1:3)= ( cmplx(z0,x(3), kind=kind(1.d0))*dil0 + dz*l0i )*drr(1:3)
      dz_plus(3) = dz_plus(3) + cplx_l0i

      dz_minus(1:3)= ( cmplx(z0,-x(3), kind=kind(1.d0))*dil0 + dz*l0i )*drr(1:3)
      dz_minus(3) = dz_minus(3) - cplx_l0i


      dx_plus(1:3) = cmplx(x(1),x(2), kind=kind(1.d0))*dil0*drr(1:3)
      dx_plus(1) = dx_plus(1) + cmplx(l0i,0.d0, kind=kind(1.d0))
      dx_plus(2) = dx_plus(2) + cplx_l0i

      dx_minus(1:3) = cmplx(x(1),-x(2), kind=kind(1.d0))*dil0*drr(1:3)
      dx_minus(1) = dx_minus(1) + cmplx(l0i,0.d0, kind=kind(1.d0))
      dx_minus(2) = dx_minus(2) - cplx_l0i

      write(6,*) "rr", rr
      write(6,*) "rcut", r_cut
      write(6,*) "r0", r0
      write(6,*) "th0", th0
      write(6,*) "l0", l0
      write(6,*) "l0i", l0i
      write(6,*) "z0", z0
      write(6,*) "z_plus", z_plus
      write(6,*) "z_minus", z_minus
      write(6,*) "x_plus", x_plus
      write(6,*) "x_minus", x_minus
      write(6,*) "drr", drr(1), drr(2), drr(3)
      write(6,*) "dz", dz
      write(6,*) "dz_plus(1)", dz_plus(1)
      write(6,*) "dz_plus(2)", dz_plus(2)
      write(6,*) "dz_plus(3)", dz_plus(3)
      write(6,*) "dz_minus(1)", dz_minus(1)
      write(6,*) "dz_minus(2)", dz_minus(2)
      write(6,*) "dz_minus(3)", dz_minus(3)

      write(6,*) "dx_plus(1)", dx_plus(1)
      write(6,*) "dx_plus(2)", dx_plus(2)
      write(6,*) "dx_plus(3)", dx_plus(3)
      write(6,*) "dx_minus(1)", dx_minus(1)
      write(6,*) "dx_minus(2)", dx_minus(2)
      write(6,*) "dx_minus(3)", dx_minus(3)

      Umm(:,:,:)=(0.d0,0.d0)
      Umm(0,0,0)=(1.d0,0.d0)
      do j=1,jj_max
        do m1=-j,j,2
          ! what is that baby ?
          Umm(m1,m1,j)=(1.d0,0.d0)
        end do
      end do

      dUmm(:,:,:,:)=(0.d0,0.d0)

      do j=1,jj_max
         do m1=-j,j,2
            do m2=-j,j,2
               if (m1.eq.j) then
                 if (m2.eq.-j) then
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
                 else if (m2.eq.j) then
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
                 else
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
                                    -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
                 end if
               else if (m1.eq.-j) then
                 if (m2.eq.j) then
                   tjm1=cplx_i*dsqrt(dble(j+m2)/dble(j-m1))
                   Umm(m2,m1,j) =  - tjm1*x_minus*Umm(m2-1,m1+1,j-1)
                   !write(*,*) 'ddd', tjm1, x_minus, Umm(m2-1,m1+1,j-1)
                   dUmm(m2,m1,j,:) =  - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))
                 else if (m2.eq.-j) then
                   tjm2= dsqrt(dble(j-m2)/dble(j-m1))
                   Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1)
                   !if (j==1) write(*,*)  'ddd', Umm(m2,m1,j)
                   dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:))
                 else

                   tjm1=cplx_i*dsqrt(dble(j+m2)/dble(j-m1))
                   tjm2= dsqrt(dble(j-m2)/dble(j-m1))
                   Umm(m2,m1,j) = tjm2*z_plus*Umm(m2+1,m1+1,j-1) - tjm1*x_minus*Umm(m2-1,m1+1,j-1)

                   dUmm(m2,m1,j,:) = tjm2*(dz_plus(:)*Umm(m2+1,m1+1,j-1)+z_plus*dUmm(m2+1,m1+1,j-1,:)) &
                                 - tjm1*(dx_minus(:)*Umm(m2-1,m1+1,j-1)+x_minus*dUmm(m2-1,m1+1,j-1,:))

                 end if
               else
                 if (m2.eq.-j) then
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = - tjm2*( dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:) )
                 else if (m2.eq.j) then
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1)+z_minus*dUmm(m2-1,m1-1,j-1,:))
                 else
                   tjm1=dsqrt(dble(j+m2)/dble(j+m1))
                   tjm2=cplx_i*dsqrt(dble(j-m2)/dble(j+m1))
                   Umm(m2,m1,j) = tjm1*z_minus*Umm(m2-1,m1-1,j-1) - tjm2*x_plus*Umm(m2+1,m1-1,j-1)
                   dUmm(m2,m1,j,:) = tjm1*(dz_minus(:)*Umm(m2-1,m1-1,j-1) + z_minus*dUmm(m2-1,m1-1,j-1,:))   &
                                    -tjm2*(dx_plus(:)*Umm(m2+1,m1-1,j-1)+x_plus*dUmm(m2+1,m1-1,j-1,:))
                 end if
               end if
               !debug
               !write(*,'(3i4,2G15.6)') j, m2,m1, Umm(m2,m1,j)
            enddo
         enddo
      enddo
   
!   Write(6,*) "Premiers coefs : "
!   Write(6,*) "U(1,1,1)=",z_minus
!   Write(6,*) "U(1,-1,1)=", -cplx_i*x_plus
!   Write(6,*) "U(1,1,-1)=", -cplx_i*x_minus
!   Write(6,*) "U(1,-1,-1)=", z_plus
   return
   end subroutine spherical_4d
