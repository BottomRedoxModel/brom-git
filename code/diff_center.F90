!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diffusion schemes --- grid centers\label{sec:diffusionMean}
!
! !INTERFACE:
   subroutine diff_center(N,dt,cnpar,posconc,h,Vc,Af,Bcup,Bcdw, &
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,Y)
!
! !DESCRIPTION:
! This subroutine solves the one-dimensional diffusion equation
! including source terms,
!  \begin{equation}
!   \label{YdiffCenter}
!    \partder{Y}{t}
!    = \partder{}{z} \left( \nu_Y \partder{Y}{z} \right)
!    - \frac{1}{\tau_R}(Y-Y_{obs})
!    + Y L_{\text{sour}} + Q_{\text{sour}}
!    \comma
!  \end{equation}
! for al variables defined at the centers of the grid cells, and
! a diffusion coefficient $\nu_Y$ defined at the faces.
! Relaxation with time scale $\tau_R$ towards observed values
! $Y_{\text{obs}}$ is possible. $L_{\text{sour}}$ specifies a
! linear source term, and $Q_{\text{sour}}$ a constant source term.
! Central differences are used to discretize the problem
! as discussed in \sect{SectionNumericsMean}. The diffusion term,
! the linear source term, and the linear part arising from the
! relaxation term are treated
! with an implicit method, whereas the constant source term is treated
! fully explicit.
!
! The input parameters {\tt Bcup} and {\tt Bcdw} specify the type
! of the upper and lower boundary conditions, which can be either
! Dirichlet or Neumann-type. {\tt Bcup} and {\tt Bcdw} must have integer
! values corresponding to the parameters {\tt Dirichlet} and {\tt Neumann}
! defined in the module {\tt util}, see \sect{sec:utils}.
! {\tt Yup} and {\tt Ydw} are the values of the boundary conditions at
! the surface and the bottom. Depending on the values of {\tt Bcup} and
! {\tt Bcdw}, they represent either fluxes or prescribed values.
! The integer {\tt posconc} indicates if a quantity is
! non-negative by definition ({\tt posconc}=1, such as for concentrations)
! or not ({\tt posconc}=0). For {\tt posconc}=1 and negative
! boundary fluxes, the source term linearisation according to
! \cite{Patankar80} is applied.
!
! Note that fluxes \emph{entering} a boundary cell are counted positive
! by convention. The lower and upper position for prescribing these fluxes
! are located at the lowest und uppermost grid faces with index "0" and
! index "N", respectively. If values are prescribed, they are located at
! the centers with index "1" and index "N", respectively.
!
! !USES:
   use util,          only  : Dirichlet, Neumann
   use mtridiagonal
   use fabm_types, only: rk

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: N

!  time step (s)
   real(rk), intent(in)                :: dt

!  "implicitness" parameter
   real(rk), intent(in)                :: cnpar

!  1: non-negative concentration, 0: else
   integer, intent(in)                 :: posconc

!  layer thickness (m)
   real(rk), intent(in)                :: h(0:N)

!  hypsograph at grid centre
   real(rk), intent(in)                :: Vc(0:N)

!  hypsograph at grid face
   real(rk), intent(in)                :: Af(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real(rk), intent(in)                :: Yup

!  value of lower BC
   real(rk), intent(in)                :: Ydw

!  diffusivity of Y
   real(rk), intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
   real(rk), intent(in)                :: Lsour(0:N)

!  constant source term
!  (treated explicitly)
   real(rk), intent(in)                :: Qsour(0:N)

!  relaxation time (s)
   real(rk), intent(in)                :: Taur(0:N)

!  observed value of Y
   real(rk), intent(in)                :: Yobs(0:N)
!
! !INPUT/OUTPUT PARAMETERS:
   real(rk)                            :: Y(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  Phil Wallhead 19/05/2016: Removed #include"cppdefs.h" statement on line 1   
!                            Replaced real(rk) declarations with real(rk) and introduced 'use fabm_types' statement to avoid use of cpps
!                            Replaced _ZERO_ and _ONE_ with 0.0_rk and 1.0_rk respectively
!                            Replaced FATAL with "write(*,*) 'FATAL ERROR: '," 
!EOP
!
! !LOCAL VARIABLES:
   integer                   :: i
   real(rk)                  :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  set up matrix
   do i=2,N-1
      c     = 2.0d0*dt*Af(i)*nuY(i)    /(h(i)+h(i+1))/Vc(i)
      a     = 2.0d0*dt*Af(i-1)*nuY(i-1)/(h(i)+h(i-1))/Vc(i)
      l     =     dt*Lsour(i)

      cu(i) =-cnpar*c
      au(i) =-cnpar*a
      bu(i) = 1.0_rk + cnpar*(a + c) - l
      du(i) = (1.0_rk - (1.0_rk-cnpar)*(a + c))*Y(i)                  &
            + (1.0_rk - cnpar)*( a*Y(i-1) + c*Y(i+1) ) + dt*Qsour(i)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = 2.0d0*dt*Af(N-1)*nuY(N-1)/(h(N)+h(N-1))/Vc(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*a
      if (posconc .eq. 1 .and. Yup.lt.0.0_rk) then ! Patankar (1980) trick
         bu(N) =  1.0_rk - au(N) - l  - dt*Yup/Y(N)/h(N)
         du(N) = Y(N) + dt*Qsour(N)   &
               + (1.0_rk - cnpar)*a*(Y(N-1)-Y(N))
      else
         bu(N) =  1.0_rk - au(N) - l
         du(N) = Y(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1.0_rk - cnpar)*a*(Y(N-1)-Y(N))
      end if
   case(Dirichlet)
      au(N) = 0.0_rk
      bu(N) = 1.0_rk
      du(N) = Yup
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for upper boundary'
      stop  'diff_center.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = 2.0d0*dt*Af(1)*nuY(1)/(h(1)+h(2))/Vc(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*c
      if (posconc.eq.1 .and. Ydw.lt.0.0_rk) then ! Patankar (1980) trick
         bu(1) = 1.0_rk - cu(1) - l - dt*Ydw/Y(1)/h(1)
         du(1) = Y(1) + dt*(Qsour(1))   &
               + (1.0_rk - cnpar)*c*(Y(2)-Y(1))
      else
         bu(1) = 1.0_rk - cu(1) - l
         du(1) = Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
               + (1.0_rk - cnpar)*c*(Y(2)-Y(1))
      end if
   case(Dirichlet)
      cu(1) = 0.0_rk
      bu(1) = 1.0_rk
      du(1) = Ydw
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for lower boundary'
      stop  'diff_center.F90'
   end select

!  relaxation to observed value
   if (minval(Taur).lt.1.d10) then
      do i=1,N
         bu(i)=bu(i)+dt/Taur(i)
         du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
   end if

!  solve linear system
   call tridiagonal(N,1,N,Y)

   return
    end subroutine diff_center
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

    
    

    
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diffusion schemes --- grid centers\label{sec:diffusionMean}
!
! !INTERFACE:
   subroutine diff_center1(N,dt,cnpar,posconc,h,Bcup,Bcdw, &
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,pF,Y)
!
! !DESCRIPTION:
! This subroutine solves the one-dimensional diffusion equation
! including source terms,
!  \begin{equation}
!   \label{YdiffCenter}
!    \partder{Y}{t}
!    = \partder{}{z} \left( \nu_Y \partder{pY}{z} \right)
!    - \frac{1}{\tau_R}(Y-Y_{obs})
!    + Y L_{\text{sour}} + Q_{\text{sour}}
!    \comma
!  \end{equation}
! for al variables defined at the centers of the grid cells, and
! a diffusion coefficient $\nu_Y$ defined at the faces.
! Relaxation with time scale $\tau_R$ towards observed values
! $Y_{\text{obs}}$ is possible. $L_{\text{sour}}$ specifies a
! linear source term, and $Q_{\text{sour}}$ a constant source term.
! Central differences are used to discretize the problem
! as discussed in \sect{SectionNumericsMean}. The diffusion term,
! the linear source term, and the linear part arising from the
! relaxation term are treated
! with an implicit method, whereas the constant source term is treated
! fully explicit.
! p is a porosity factor to account any necessary changes of units from 
! those of Y (assumed to be [mass per unit total volume]) to either 
! [mass per unit volume pore water] or [mass per unit volume solids]
! when calculating Fickian gradients in the sediments 
! (see Berner, 1980; Boudreau, 1997).
! p = 1 in the water column
! p = 1/phi for solutes in the sediments
! p = 1/(1-phi) for solids in the sediments
! Note that $\nu_Y$ is assumed to already include any area restriction
! due to porosity (phi for solutes, (1-phi) for solids).
! A further volume factor is not needed to account for porosity
! because this subroutine assumes that [mass per unit total volume]
! is the modelling currency for all concentrations.
!
! The input parameters {\tt Bcup} and {\tt Bcdw} specify the type
! of the upper and lower boundary conditions, which can be either
! Dirichlet or Neumann-type. {\tt Bcup} and {\tt Bcdw} must have integer
! values corresponding to the parameters {\tt Dirichlet} and {\tt Neumann}
! defined in the module {\tt util}, see \sect{sec:utils}.
! {\tt Yup} and {\tt Ydw} are the values of the boundary conditions at
! the surface and the bottom. Depending on the values of {\tt Bcup} and
! {\tt Bcdw}, they represent either fluxes or prescribed values.
! The integer {\tt posconc} indicates if a quantity is
! non-negative by definition ({\tt posconc}=1, such as for concentrations)
! or not ({\tt posconc}=0). For {\tt posconc}=1 and negative
! boundary fluxes, the source term linearisation according to
! \cite{Patankar80} is applied.
!
! Note that fluxes \emph{entering} a boundary cell are counted positive
! by convention. The lower and upper position for prescribing these fluxes
! are located at the lowest und uppermost grid faces with index "0" and
! index "N", respectively. If values are prescribed, they are located at
! the centers with index "1" and index "N", respectively.
!
! !USES:
   use util,          only  : Dirichlet, Neumann
   use mtridiagonal
   use fabm_types, only: rk

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: N

!  time step (s)
   real(rk), intent(in)                :: dt

!  "implicitness" parameter
   real(rk), intent(in)                :: cnpar

!  1: non-negative concentration, 0: else
   integer, intent(in)                 :: posconc

!  layer thickness (m)
   real(rk), intent(in)                :: h(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real(rk), intent(in)                :: Yup

!  value of lower BC
   real(rk), intent(in)                :: Ydw

!  diffusivity of Y
   real(rk), intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
   real(rk), intent(in)                :: Lsour(0:N)

!  constant source term
!  (treated explicitly)
   real(rk), intent(in)                :: Qsour(0:N)

!  relaxation time (s)
   real(rk), intent(in)                :: Taur(0:N)

!  observed value of Y
   real(rk), intent(in)                :: Yobs(0:N)
   
!  porosity factor for Fickian gradient
   real(rk), intent(in)                :: pF(0:N)   

! !INPUT/OUTPUT PARAMETERS:
   real(rk)                            :: Y(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  Phil Wallhead 24/02/2016: Commented out #include"cppdefs.h" statement on line 1   
!                            Replaced real(rk) declarations with real(rk) and introduced 'use fabm_types' statement to avoid use of cpps
!                            Replaced _ZERO_ and _ONE_ with 0.0_rk and 1.0_rk respectively
!                            Replaced FATAL with "write(*,*) 'FATAL ERROR: '," 
!  Phil Wallhead 18/05/2016: Introduced porosity factors pF to convert units for calculating Fickian gradients in the sediments
!   
!EOP
!
! !LOCAL VARIABLES:
   integer                :: i
   real(rk)               :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  set up matrix
   do i=2,N-1
      c     = 2.0d0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0d0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     =     dt*Lsour(i)

      cu(i) =-cnpar*pF(i+1)*c
      au(i) =-cnpar*pF(i-1)*a
      bu(i) = 1.0_rk + cnpar*pF(i)*(a + c) - l
      du(i) = (1.0_rk - (1.0_rk-cnpar)*pF(i)*(a + c))*Y(i)                  &
            + (1.0_rk - cnpar)*( pF(i-1)*a*Y(i-1) + pF(i+1)*c*Y(i+1) ) + dt*Qsour(i)
   !  Note: all terms of coefficients involving cnpar or (1-cnpar) derive from the Fickian gradient
   !        and therefore need to be multiplied by the appropriate porosity factors (pF)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = 2.0d0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*pF(N-1)*a
      if (posconc .eq. 1 .and. Yup.lt.0.0_rk) then ! Patankar (1980) trick
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l  - dt*Yup/Y(N)/h(N)
         du(N) = Y(N) + dt*Qsour(N)   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      else
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l
         du(N) = Y(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      end if
   case(Dirichlet)
      au(N) = 0.0_rk
      bu(N) = 1.0_rk
      du(N) = Yup
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for upper boundary'
      stop  'diff_center.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = 2.0d0*dt*nuY(1)/(h(1)+h(2))/h(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*pF(2)*c
      if (posconc.eq.1 .and. Ydw.lt.0.0_rk) then ! Patankar (1980) trick
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l - dt*Ydw/Y(1)/h(1)
         du(1) = Y(1) + dt*(Qsour(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      else
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l
         du(1) = Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      end if
   case(Dirichlet)
      cu(1) = 0.0_rk
      bu(1) = 1.0_rk
      du(1) = Ydw
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for lower boundary'
      stop  'diff_center.F90'
   end select

!  relaxation to observed value
   if (minval(Taur).lt.1.d10) then
      do i=1,N
         bu(i)=bu(i)+dt/Taur(i)
         du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
   end if

!  solve linear system
   call tridiagonal(N,1,N,Y)

   return
    end subroutine diff_center1
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------

    
    

    
!-----------------------------------------------------------------------
!BOP
!
! !ROUTINE: Diffusion schemes --- grid centers\label{sec:diffusionMean}
!
! !INTERFACE:
   subroutine diff_center2(N,dt,cnpar,posconc,h,Bcup,Bcdw, &
                          Yup,Ydw,nuY,Lsour,Qsour,Taur,Yobs,pF, &
                          pFSWIup, pFSWIdw, i_sed_top, Y)
!
! !DESCRIPTION:
! This subroutine solves the one-dimensional diffusion equation
! including source terms,
!  \begin{equation}
!   \label{YdiffCenter}
!    \partder{Y}{t}
!    = \partder{}{z} \left( \nu_Y \partder{pY}{z} \right)
!    - \frac{1}{\tau_R}(Y-Y_{obs})
!    + Y L_{\text{sour}} + Q_{\text{sour}}
!    \comma
!  \end{equation}
! for al variables defined at the centers of the grid cells, and
! a diffusion coefficient $\nu_Y$ defined at the faces.
! Relaxation with time scale $\tau_R$ towards observed values
! $Y_{\text{obs}}$ is possible. $L_{\text{sour}}$ specifies a
! linear source term, and $Q_{\text{sour}}$ a constant source term.
! Central differences are used to discretize the problem
! as discussed in \sect{SectionNumericsMean}. The diffusion term,
! the linear source term, and the linear part arising from the
! relaxation term are treated
! with an implicit method, whereas the constant source term is treated
! fully explicit.
! p is a porosity factor to account any necessary changes of units from 
! those of Y (assumed to be [mass per unit total volume]) to either 
! [mass per unit volume pore water] or [mass per unit volume solids]
! when calculating Fickian gradients in the sediments 
! (see Berner, 1980; Boudreau, 1997).
! p = 1 in the water column
! p = 1/phi for solutes in the sediments
! p = 1/(1-phi) for solids in the sediments
! Note that $\nu_Y$ is assumed to already include any area restriction
! due to porosity (phi for solutes, (1-phi) for solids).
! A further volume factor is not needed to account for porosity
! because this subroutine assumes that [mass per unit total volume]
! is the modelling currency for all concentrations.
! This routine also makes special cases of the cells neighbouring the
! sediment-water interface, allowing for interphase mixing.
!
! The input parameters {\tt Bcup} and {\tt Bcdw} specify the type
! of the upper and lower boundary conditions, which can be either
! Dirichlet or Neumann-type. {\tt Bcup} and {\tt Bcdw} must have integer
! values corresponding to the parameters {\tt Dirichlet} and {\tt Neumann}
! defined in the module {\tt util}, see \sect{sec:utils}.
! {\tt Yup} and {\tt Ydw} are the values of the boundary conditions at
! the surface and the bottom. Depending on the values of {\tt Bcup} and
! {\tt Bcdw}, they represent either fluxes or prescribed values.
! The integer {\tt posconc} indicates if a quantity is
! non-negative by definition ({\tt posconc}=1, such as for concentrations)
! or not ({\tt posconc}=0). For {\tt posconc}=1 and negative
! boundary fluxes, the source term linearisation according to
! \cite{Patankar80} is applied.
!
! Note that fluxes \emph{entering} a boundary cell are counted positive
! by convention. The lower and upper position for prescribing these fluxes
! are located at the lowest und uppermost grid faces with index "0" and
! index "N", respectively. If values are prescribed, they are located at
! the centers with index "1" and index "N", respectively.
!
! !USES:
   use util,          only  : Dirichlet, Neumann
   use mtridiagonal
   use fabm_types, only: rk

   IMPLICIT NONE
!
! !INPUT PARAMETERS:

!  number of vertical layers
   integer,  intent(in)                :: N

!  time step (s)
   real(rk), intent(in)                :: dt

!  "implicitness" parameter
   real(rk), intent(in)                :: cnpar

!  1: non-negative concentration, 0: else
   integer, intent(in)                 :: posconc

!  layer thickness (m)
   real(rk), intent(in)                :: h(0:N)

!  type of upper BC
   integer,  intent(in)                :: Bcup

!  type of lower BC
   integer,  intent(in)                :: Bcdw

!  value of upper BC
   real(rk), intent(in)                :: Yup

!  value of lower BC
   real(rk), intent(in)                :: Ydw

!  diffusivity of Y
   real(rk), intent(in)                :: nuY(0:N)

!  linear source term
!  (treated implicitly)
   real(rk), intent(in)                :: Lsour(0:N)

!  constant source term
!  (treated explicitly)
   real(rk), intent(in)                :: Qsour(0:N)

!  relaxation time (s)
   real(rk), intent(in)                :: Taur(0:N)

!  observed value of Y
   real(rk), intent(in)                :: Yobs(0:N)
   
!  porosity factor for Fickian gradient
   real(rk), intent(in)                :: pF(0:N)   
   
!  porosity factor for cell just above the sediment-water interface
   real(rk), intent(in)                :: pFSWIup
   
!  porosity factor for cell just below the sediment-water interface
   real(rk), intent(in)                :: pFSWIdw   

!  index to identify the cell just below the sediment-water interface
   integer, intent(in)                 :: i_sed_top
   
! !INPUT/OUTPUT PARAMETERS:
   real(rk)                            :: Y(0:N)
!
! !REVISION HISTORY:
!  Original author(s): Lars Umlauf
!
!  Phil Wallhead 24/02/2016: Commented out #include"cppdefs.h" statement on line 1   
!                            Replaced real(rk) declarations with real(rk) and introduced 'use fabm_types' statement to avoid use of cpps
!                            Replaced _ZERO_ and _ONE_ with 0.0_rk and 1.0_rk respectively
!                            Replaced FATAL with "write(*,*) 'FATAL ERROR: '," 
!  Phil Wallhead 18/05/2016: Introduced porosity factors pF to convert units for calculating Fickian gradients in the sediments
!  Phil Wallhead 13/06/2016: Introduced special case porosity factors pFSWIup, pFSWIdw and index i_sed_top to allow for interphase mixing across SWI
!   
!EOP
!
! !LOCAL VARIABLES:
   integer                :: i
   real(rk)               :: a,c,l
!
!-----------------------------------------------------------------------
!BOC
!
!  set up matrix
   do i=2,N-1
      c     = 2.0d0*dt*nuY(i)  /(h(i)+h(i+1))/h(i)
      a     = 2.0d0*dt*nuY(i-1)/(h(i)+h(i-1))/h(i)
      l     =     dt*Lsour(i)

      if (i.eq.i_sed_top) then        !Special treatment for top cell of sediments (PWA)
         !pF(i+1) -> pFSWIup, pF(i)*c -> pFSWIdw*c
         cu(i) =-cnpar*pFSWIup*c
         au(i) =-cnpar*pF(i-1)*a
         bu(i) = 1.0_rk + cnpar*(pF(i)*a + pFSWIdw*c) - l
         du(i) = (1.0_rk - (1.0_rk-cnpar)*(pF(i)*a + pFSWIdw*c))*Y(i)                  &
               + (1.0_rk - cnpar)*( pF(i-1)*a*Y(i-1) + pFSWIup*c*Y(i+1) ) + dt*Qsour(i)          
      else if (i.eq.i_sed_top+1) then !Special treatment for bottom cell of water column (PWA)
         !pF(i-1) -> pFSWIdw, pF(i)*a -> pFSWIup*a
         cu(i) =-cnpar*pF(i+1)*c
         au(i) =-cnpar*pFSWIdw*a
         bu(i) = 1.0_rk + cnpar*(pFSWIup*a + pF(i)*c) - l
         du(i) = (1.0_rk - (1.0_rk-cnpar)*(pFSWIup*a + pF(i)*c))*Y(i)                  &
               + (1.0_rk - cnpar)*( pFSWIdw*a*Y(i-1) + pF(i+1)*c*Y(i+1) ) + dt*Qsour(i)              
      else
         cu(i) =-cnpar*pF(i+1)*c
         au(i) =-cnpar*pF(i-1)*a
         bu(i) = 1.0_rk + cnpar*pF(i)*(a + c) - l
         du(i) = (1.0_rk - (1.0_rk-cnpar)*pF(i)*(a + c))*Y(i)                  &
               + (1.0_rk - cnpar)*( pF(i-1)*a*Y(i-1) + pF(i+1)*c*Y(i+1) ) + dt*Qsour(i)
      end if
      !  Note: all terms of coefficients involving cnpar or (1-cnpar) derive from the Fickian gradient
      !        and therefore need to be multiplied by the appropriate porosity factors (pF)
   end do

!   set up upper boundary condition
   select case(Bcup)
   case(Neumann)
      a     = 2.0d0*dt*nuY(N-1)/(h(N)+h(N-1))/h(N)
      l     = dt*Lsour(N)

      au(N) =-cnpar*pF(N-1)*a
      if (posconc .eq. 1 .and. Yup.lt.0.0_rk) then ! Patankar (1980) trick
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l  - dt*Yup/Y(N)/h(N)
         du(N) = Y(N) + dt*Qsour(N)   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      else
         bu(N) =  1.0_rk + cnpar*pF(N)*a - l
         du(N) = Y(N) + dt*(Qsour(N)+Yup/h(N))   &
               + (1.0_rk - cnpar)*a*( pF(N-1)*Y(N-1) - pF(N)*Y(N) )
      end if
   case(Dirichlet)
      au(N) = 0.0_rk
      bu(N) = 1.0_rk
      du(N) = Yup
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for upper boundary'
      stop  'diff_center.F90'
   end select

!   set up lower boundary condition
   select case(Bcdw)
   case(Neumann)
      c     = 2.0d0*dt*nuY(1)/(h(1)+h(2))/h(1)
      l     = dt*Lsour(1)

      cu(1) =-cnpar*pF(2)*c
      if (posconc.eq.1 .and. Ydw.lt.0.0_rk) then ! Patankar (1980) trick
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l - dt*Ydw/Y(1)/h(1)
         du(1) = Y(1) + dt*(Qsour(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      else
         bu(1) = 1.0_rk + cnpar*pF(1)*c - l
         du(1) = Y(1) + dt*(Qsour(1)+Ydw/h(1))   &
               + (1.0_rk - cnpar)*c*( pF(2)*Y(2) - pF(1)*Y(1) )
      end if
   case(Dirichlet)
      cu(1) = 0.0_rk
      bu(1) = 1.0_rk
      du(1) = Ydw
   case default
      write(*,*) 'FATAL ERROR: ', 'invalid boundary condition type for lower boundary'
      stop  'diff_center.F90'
   end select

!  relaxation to observed value
   if (minval(Taur).lt.1.d10) then
      do i=1,N
         bu(i)=bu(i)+dt/Taur(i)
         du(i)=du(i)+dt/Taur(i)*Yobs(i)
      end do
   end if

!  solve linear system
   call tridiagonal(N,1,N,Y)

   return
   end subroutine diff_center2
!EOC

!-----------------------------------------------------------------------
! Copyright by the GOTM-team under the GNU Public License - www.gnu.org
!-----------------------------------------------------------------------
