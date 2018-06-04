!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
module fockbias_subs
!------------------------------------------------------------------------------!
!
!    This module controls the aplication of charge biases that are directly
! applied through the fock matrix. It is implementation specific and consists
! of the following subroutines:
!
! fockbias_setup0: This subroutine sets up the shape of the bias.
!
! fockbias_setorb: This subroutine sets up the atomic charges.
!
! fockbias_setmat: This subroutine sets up the bias matrix (full amplitude).
!                ( changes if atomic positions change )
!
! fockbias_loads: subroutine that reads the atomic biases from list in a
!                 an external file (and calls setorb)
!
! fockbias_apply: Takes the fock matrix and adds the fockbias term (with the
!                 corresponding time shape). This subroutine can take 
!                 either a complex or a real fock input (this is kind of 
!                 useless, fock is always real, but it is implemented this way.
!
! REFERENCE: Physical Review B 74, 155112, 2006
!
!------------------------------------------------------------------------------!
   implicit none

   interface fockbias_apply
      module procedure fockbias_apply_d
      module procedure fockbias_apply_c
      module procedure fockbias_apply_z
   end interface fockbias_apply

   contains
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setup0( activate, do_shape, grow_start, fall_start, setamp )

   use fockbias_data, only: fockbias_is_active, fockbias_is_shaped &
                         &, fockbias_timegrow , fockbias_timefall  &
                         &, fockbias_timeamp0
   logical, intent(in) :: activate
   logical, intent(in) :: do_shape
   real*8 , intent(in) :: grow_start
   real*8 , intent(in) :: fall_start
   real*8 , intent(in) :: setamp

   fockbias_is_active = activate
   fockbias_is_shaped = do_shape
   fockbias_timegrow  = grow_start
   fockbias_timefall  = fall_start
   fockbias_timeamp0  = setamp
end subroutine fockbias_setup0
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setorb( qweight_of_atom, atom_of_orb )

   use fockbias_data, only: fockbias_is_active, fockbias_orbqw

   real*8 , intent(in) :: qweight_of_atom(:)
   integer, intent(in) :: atom_of_orb(:)
   integer             :: Nbasis
   integer             :: nn

   Nbasis = size(atom_of_orb)
   allocate( fockbias_orbqw(Nbasis) )
   do nn = 1, Nbasis
      fockbias_orbqw(nn) = qweight_of_atom( atom_of_orb(nn) )
   end do

end subroutine fockbias_setorb
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_setmat(sqsmat)
   use fockbias_data, only: fockbias_is_active, fockbias_orbqw, fockbias_matrix

   real*8 , intent(in) :: sqsmat(:,:)
   real*8              :: newterm
   integer             :: Nbasis
   integer             :: ii, jj, kk

   if (.not.fockbias_is_active) return

   Nbasis = size( fockbias_orbqw ) !this is M

   allocate( fockbias_matrix(Nbasis,Nbasis) )
   do jj = 1, Nbasis
   do ii = 1, Nbasis
      fockbias_matrix(ii,jj) = 0.0d0
      do kk = 1, Nbasis
         newterm = fockbias_orbqw(kk) * sqsmat(ii,kk) * sqsmat(kk,jj)
         fockbias_matrix(ii,jj) = fockbias_matrix(ii,jj) + newterm
      end do
   end do
   end do

end subroutine fockbias_setmat
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_loads( Natom, atom_of_orb, file_unit_in, file_name_in )
   use fockbias_data, only: fockbias_is_active, fockbias_readfile

   integer         , intent(in)           :: Natom
   integer         , intent(in)           :: atom_of_orb(:)
   integer         , intent(in), optional :: file_unit_in
   character(len=*), intent(in), optional :: file_name_in
   integer                                :: file_unit
   real*8          , allocatable          :: qweight_of_atom(:)
   integer                                :: nn, ios


   fockbias_readfile = file_name_in
   file_unit = file_unit_in

   open( file=fockbias_readfile, unit=file_unit, iostat=ios )

!  Now read all atomic weights and set it up.
   allocate( qweight_of_atom(Natom) )
   do nn = 1, Natom
      read( unit=file_unit, fmt=*, iostat=ios ) qweight_of_atom(nn)
   end do
   call fockbias_setorb( qweight_of_atom, atom_of_orb )


end subroutine fockbias_loads
!
!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
subroutine fockbias_apply_d( timepos, fockmat )
   !here is where timegrow, timefall and timeamp are used
   use fockbias_data, only: fockbias_is_active, fockbias_is_shaped &
                         &, fockbias_timegrow , fockbias_timefall  &
                         &, fockbias_timeamp0 , fockbias_matrix
   implicit none
   real*8 , intent(in)    :: timepos
   real*8 , intent(inout) :: fockmat(:,:)

   real*8  :: time_shape, exparg
   integer :: Nbasis
   integer :: ii, jj

   Nbasis = size( fockbias_matrix, 1 )

   time_shape = 1.0d0

   if ( fockbias_is_shaped ) then
   ! this applies a time dependence to the bias of the form:
   !        ___________________________
   !       /                           \
   !      /                             \
   !   /                                   \
   !_/      t_grow              t_fall       \_
   ! where the constant part in the middle is delimited,
   ! between t_grow and t_fall, the growth and fall are gaussian

      if ( timepos <= fockbias_timegrow ) then
         exparg = (timepos - fockbias_timegrow) / fockbias_timeamp0
         exparg = (-1.0d0) * exparg * exparg
         time_shape = time_shape * dexp( exparg )
      end if

      if ( timepos >= fockbias_timefall ) then
         exparg = (timepos - fockbias_timefall) / fockbias_timeamp0
         exparg = (-1.0d0) * exparg * exparg
         time_shape = time_shape * dexp( exparg )
      end if

   end if

   do jj = 1, Nbasis
   do ii = 1, Nbasis
      fockmat(ii,jj) = fockmat(ii,jj) + time_shape * fockbias_matrix(ii,jj)
   enddo
   enddo

end subroutine fockbias_apply_d
