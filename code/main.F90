! This file is part of Bottom RedOx Model (BROM, v.1.1).
! BROM is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free
! Software Foundation (https://www.gnu.org/licenses/gpl.html).
! It is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE. A copy of the license is provided in
! the COPYING file at the root of the BROM distribution.
!-----------------------------------------------------------------------
! Original author(s): Evgeniy Yakushev, Shamil Yakubov, 
!                     Elizaveta Protsenko, Phil Wallhead
!-----------------------------------------------------------------------
    
    
    program main
    
    use brom_transport, only: init_brom_transport, do_brom_transport, clear_brom_transport

    !initializing, writing from gotm data and fabm.yaml included
    call init_brom_transport()
    !main cycle
    call do_brom_transport()
    !clear all
    !call clear_brom_transport()

    end program main
