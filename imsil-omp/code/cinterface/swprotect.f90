module Cinterface_Swprotect_mod
!    use iso_c_binding
    implicit none

    interface
!       bool initLMXLicense();
        function initLMXLicense() bind(C, name="initLMXLicense")
            use iso_c_binding
            logical(c_bool) :: initLMXLicense
        end function initLMXLicense

!       void shutdownLMXLicense();
        subroutine shutdownLMXLicense() bind(C, name="shutdownLMXLicense")
            use iso_c_binding
        end subroutine shutdownLMXLicense

    end interface

    public :: initGTSLicense
    public :: shutdownGTSLicense

    contains
    
    subroutine initGTSLicense( )
        if (.not.initLMXLicense( )) then
            error stop 'Failed to obtain license "prosimMonteCarloImplant".'
        end if
    end subroutine initGTSLicense

    subroutine shutdownGTSLicense( )
        call shutdownLMXLicense( )
    end subroutine shutdownGTSLicense

end module Cinterface_Swprotect_mod
