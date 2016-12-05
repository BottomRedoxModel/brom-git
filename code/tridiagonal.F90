!#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: mtridiagonal --- solving the system\label{sec:tridiagonal}
!
! !INTERFACE:
   MODULE mtridiagonal
!
! !DESCRIPTION:
!
!  Solves a linear system of equations with a tridiagonal matrix
!  using Gaussian elimination.
!
   use fabm_types, only: rk  
! !PUBLIC MEMBER FUNCTIONS:
   public init_tridiagonal, tridiagonal, clean_tridiagonal
!
! !PUBLIC DATA MEMBERS:
   real(kind=rk), dimension(:), allocatable     :: au,bu,cu,du
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!   
! Phil Wallhead 24/02/2016: Commented out #include"cppdefs.h" statement on line 1   
!                           Replaced real(kind=rk) declarations with real(kind=rk) and use mkind statement to avoid use of cpps   
!
!EOP
!
!  private data members
   real(kind=rk), private, dimension(:), allocatable  ::  ru,qu
!
!-----------------------------------------------------------------------

   contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Allocate memory
!
! !INTERFACE:
   subroutine init_tridiagonal(N)
!
! !DESCRIPTION:
!  This routines allocates memory necessary to perform the Gaussian
!  elimination.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: rc
!
!-----------------------------------------------------------------------
!BOC
   write(*,*) '   ', 'init_tridiagonal'
   allocate(au(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating au)'
   au = 0.0_rk

   allocate(bu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating bu)'
   bu = 0.0_rk

   allocate(cu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating cu)'
   cu = 0.0_rk

   allocate(du(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating du)'
   du = 0.0_rk

   allocate(ru(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating ru)'
   ru = 0.0_rk

   allocate(qu(0:N),stat=rc)
   if (rc /= 0) stop 'init_tridiagonal: Error allocating qu)'
   qu = 0.0_rk

   return
   end subroutine init_tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: Simplified Gaussian elimination
!
! !INTERFACE:
   subroutine tridiagonal(N,fi,lt,value)
!
! !DESCRIPTION:
! A linear equation with tridiagonal matrix structure is solved here. The main
! diagonal is stored on {\tt bu}, the upper diagonal on {\tt au}, and the
! lower diagonal on {\tt cu}, the right hand side is stored on {\tt du}.
! The method used here is the simplified Gauss elimination, also called
! \emph{Thomas algorithm}.
!
! !USES:
   IMPLICIT NONE
!
! !INPUT PARAMETERS:
   integer, intent(in)                 :: N,fi,lt
!
! !OUTPUT PARAMETERS:
   real(kind=rk)                            :: value(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Hans Burchard & Karsten Bolding
!
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
!
!-----------------------------------------------------------------------
!BOC
   ru(lt)=au(lt)/bu(lt)
   qu(lt)=du(lt)/bu(lt)

   do i=lt-1,fi+1,-1
      ru(i)=au(i)/(bu(i)-cu(i)*ru(i+1))
      qu(i)=(du(i)-cu(i)*qu(i+1))/(bu(i)-cu(i)*ru(i+1))
   end do

   qu(fi)=(du(fi)-cu(fi)*qu(fi+1))/(bu(fi)-cu(fi)*ru(fi+1))

   value(fi)=qu(fi)
   do i=fi+1,lt
      value(i)=qu(i)-ru(i)*value(i-1)
   end do


   return
   end subroutine tridiagonal
!EOC

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: De-allocate memory
!
! !INTERFACE:
   subroutine clean_tridiagonal()
!
! !DESCRIPTION:
!  De-allocates memory allocated in init\_tridiagonal.
!
! !USES:
   IMPLICIT NONE
!
! !REVISION HISTORY:
!  Original author(s): Karsten Bolding
!
!EOP
!-----------------------------------------------------------------------
!BOC
   write(*,*) '   ', 'clean_tridiagonal'

   if (allocated(au)) deallocate(au)

   if (allocated(bu)) deallocate(bu)

   if (allocated(cu)) deallocate(cu)

   if (allocated(du)) deallocate(du)

   if (allocated(ru)) deallocate(ru)

   if (allocated(qu)) deallocate(qu)

   return
   end subroutine clean_tridiagonal
!EOC
!-----------------------------------------------------------------------

   end module mtridiagonal

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
