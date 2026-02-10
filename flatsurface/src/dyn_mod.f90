
module Dyn_mod
use CelType_mod
use basicfun_mod
use point_mod
implicit none

    contains
    subroutine Gradually_Vary_Vol(C,k,t,p)
        type(Cel)::C(:)
        integer,intent(in):: k,t
        real,intent(in)::p
            if(abs(C(k)%vol_total_zero-C(k)%vol_rel).gt.0.000001)then
                C(k)%vol_total_zero=C(k)%vol_total_zero+(C(k)%vol_rel- &
                C(k)%vol_total_zero)*p
            end if
            if(abs(C(k)%vol_total_zero-C(k)%vol_rel).lt.0.1)C(k)%vol_total_zero &
            =C(k)%vol_rel
    end subroutine Gradually_Vary_Vol
    
    subroutine Gradually_Vary_Area(C,k,t,p)
        type(Cel)::C(:)
        integer,intent(in):: k,t
        real,intent(in)::p
            if(abs(C(k)%area_tot_zero-C(k)%area_rel).gt.0.1.and.mod(t,C(k)%N_V).eq.0)then
                C(k)%area_tot_zero=C(k)%area_tot_zero+ p*500 !(C(k)%area_rel- &
                !C(k)%area_tot_zero)*p
            end if
            if(abs(C(k)%area_tot_zero-C(k)%area_rel).lt.0.1)C(k)%area_tot_zero &
            =C(k)%area_rel
    end subroutine Gradually_Vary_Area
    
end module Dyn_mod