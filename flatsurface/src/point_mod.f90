
module point_mod
implicit none
type PtCell

 real*8 rv_move(1:3)
 real*8 PV_Force_Point(1:3), PM_Force_Point(1:3), MS_Force_Point(1:3)
 real*8 Rep_Force_Point(1:3)
 real*8 Darea_Point(1:20,1:3),Dvol_Point(1:3)
 integer V_move

end type PtCell
end module point_mod
