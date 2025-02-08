! Created by egor on 2/7/25.

module geometry

    ! c_double is a C-compatible double precision floating point number
    ! c_ptr is an opaque pointer that is owned by C code
    use, intrinsic :: iso_c_binding, only: c_double, c_ptr
    implicit none

    ! Struct returned by get_liquid_props
    ! contains surface area and curvature
    type, bind(C) :: geometry_props
        real(c_double) :: area, kappa
    end type geometry_props

    interface

        ! Allocates and initializes a surface object on the heap
        ! Arguments:
        !   c_double - contact angle, deg
        !   c_double - neck filling angle, deg
        !   c_double - primary particle radius, m
        ! Return value:
        !   c_ptr - pointer to the surface object
        function initialize_constant_mean_curvature_surface(a, b, c) result(d) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            real(c_double), intent(in), value :: a, b, c
            type(c_ptr) :: d
        end function initialize_constant_mean_curvature_surface

        ! Deletes the surface object (memory freed) and sets the pointer to NULL
        ! Arguments:
        !   c_ptr - pointer to the surface object
        subroutine delete_constant_mean_curvature_surface(a) bind(C)
            use, intrinsic :: iso_c_binding, only: c_ptr
            implicit none
            type(c_ptr), intent(in) :: a
        end subroutine delete_constant_mean_curvature_surface

        ! Returns a geometry_props structure (with surface area and curvature)
        ! Arguments:
        !   c_ptr - pointer to the surface objet
        !   c_double - total volume of coating in the gap, m^3
        ! Return value:
        !   geometry_props - structure with area and curvature
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

        ! Returns the maximum liquid volume that can be accommodated by the geometry model
        ! Arguments:
        !   c_ptr - pointer to the surface object
        ! Return value:
        !   c_double - maximum liquid volume, m^3
        function get_max_liquid_volume(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_max_liquid_volume

        ! Returns the filling angle as a function of coating volume
        ! Arguments:
        !   c_ptr - pointer to the surface object
        !   c_double - total condensate volume in the gap, m^3
        ! Return value:
        !   c_double - filling angle, deg
        function get_filling_angle(a, b) result(c) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double), intent(in), value :: b
            real(c_double) :: c
        end function get_filling_angle

        ! Returns the volume of the neck
        ! Arguments:
        !   c_ptr - pointer to the surface object
        ! Return value:
        !   c_double - neck volume, m^3
        function get_neck_volume(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_neck_volume

        ! Returns the primary particle radius
        ! Arguments:
        !   c_ptr - pointer to the surface object
        ! Return value:
        !   c_double - primary particle radius, m
        function get_r_part(a) result(b) bind(C)
            use, intrinsic :: iso_c_binding, only: c_double, c_ptr
            implicit none
            type(c_ptr), intent(in), value :: a
            real(c_double) :: b
        end function get_r_part

    end interface

end module geometry