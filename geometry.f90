! Created by egor on 2/7/25.

module geometry

    use, intrinsic :: iso_c_binding, only: c_double, c_ptr
    implicit none

    type, bind(C) :: geometry_props
        real(c_double) :: area, kappa
    end type geometry_props

    interface

        function initialize_constant_mean_curvature_surface(a, b, c) result(d) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            real(c_double), intent(in), value :: a, b, c
            type(c_ptr) :: d
        end function initialize_constant_mean_curvature_surface

        subroutine delete_constant_mean_curvature_surface(a) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
        end subroutine delete_constant_mean_curvature_surface

        function get_liquid_props(a, b) result(c) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none

            type, bind(C) :: geometry_props
                real(c_double) :: area, kappa
            end type geometry_props

            type(c_ptr), intent(in), value :: a
            real(c_double), intent(in), value :: b
            type(geometry_props) :: c
        end function get_liquid_props

        function get_max_liquid_volume(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_max_liquid_volume

        function get_filling_angle(a, b) result(c) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double), intent(in), value :: b
            real(c_double) :: c
        end function get_filling_angle

        function get_neck_volume(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_neck_volume

        function get_r_part(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_r_part

    end interface

end module geometry