! This program validates that the geometry model reproduces the expected results
! for a simple cylindrical coating (90deg filling angle, 0deg contact angle)

program main

    use geometry
    use, intrinsic :: iso_c_binding, only: c_double, c_ptr
    implicit none

    real(c_double) :: PI
    real(c_double) :: r_part, contact_angle, neck_filling_angle, neck_volume, max_volume
    real(c_double) :: liquid_volume, target_area, target_kappa, condensate_volume, filling_angle, target_filling_angle
    type(c_ptr) :: surface_ptr
    type(geometry_props) :: liquid_props

    PI = 4.0_c_double * atan(1.0_c_double)
    contact_angle = 0.0_c_double
    neck_filling_angle = 10.0_c_double
    r_part = 0.1_c_double
    liquid_volume = 2.0_c_double / 3.0_c_double * PI * r_part ** 3.0_c_double
    target_area =  4.0_c_double * PI * r_part ** 2.0_c_double
    target_kappa = 0.5_c_double / r_part
    target_filling_angle = 90.0_c_double

    surface_ptr = initialize_constant_mean_curvature_surface(contact_angle, neck_filling_angle, r_part)

    neck_volume = get_neck_volume(surface_ptr)
    condensate_volume = liquid_volume - neck_volume

    liquid_props = get_liquid_props(surface_ptr, condensate_volume)
    filling_angle = get_filling_angle(surface_ptr, condensate_volume)
    max_volume = get_max_liquid_volume(surface_ptr)

    print *, "Area", target_area, liquid_props%area
    print *, "Curvature", target_kappa, liquid_props%kappa
    print *, "Filling angle", target_filling_angle, filling_angle
    print *, "Max volume", liquid_volume, max_volume

    call delete_constant_mean_curvature_surface(surface_ptr)

end program main