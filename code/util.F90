!#include"cppdefs.h"
!-----------------------------------------------------------------------
!BOP
!
! !MODULE: util --- parameters and interfaces for utilities\label{sec:utils}
!
! !INTERFACE:
   MODULE util
!
! !DESCRIPTION:
! This module is an encapsulation of a number of parameters used by different
! routines found in the {\tt util} directory. It should make it easier to
! read the code, since finding a line like
!
! \vspace{7mm}
!
! {\tt if (method.eq.UPSTREAM) then ...}
!
! \vspace{7mm}
!
! in a subroutine for advection methods tells you more than reading only
!
! \vspace{7mm}
!
! {\tt if (method.eq.1) then ...}

!
! \vspace{7mm}

!
! !USES:
   IMPLICIT NONE
!
! !DEFINED PARAMETERS:

!  type of advection scheme
   integer,parameter                    :: UPSTREAM       = 1
   integer,parameter                    :: P1             = 2
   integer,parameter                    :: P2             = 3
   integer,parameter                    :: Superbee       = 4
   integer,parameter                    :: MUSCL          = 5
   integer,parameter                    :: P2_PDM         = 6

!  boundary condition type
!  for diffusion scheme
   integer,parameter                    :: Dirichlet      = 0
   integer,parameter                    :: Neumann        = 1

!  boundary condition type
!  for advection schemes
   integer,parameter                    :: flux           = 1
   integer,parameter                    :: value          = 2
   integer,parameter                    :: oneSided       = 3
   integer,parameter                    :: zeroDivergence = 4

!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
! Phil Wallhead 24/02/2016: Commented out #include"cppdefs.h" statement on line 1   
!
!EOP


   end module util

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
