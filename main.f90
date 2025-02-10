! Created by egor on 2/7/25.

program main

    use geometry
    use, intrinsic :: iso_c_binding, only: c_double, c_ptr
    implicit none

    ! Parameters needed to initialize the geometry object
    real(c_double) :: contact_angle, neck_filling_angle, r_part
    ! When the geometry object is initialized, a C pointer is returned to Fortran
    type(c_ptr) :: constant_mean_curvature_surface
    ! A structure that stores volume and curvature
    type(geometry_props) :: props
    ! A variable that contains the total volume of condensate in the gap
    real(c_double) :: condensate_volume
    ! A variable that contains the maximum coating volume allowed by the capillary condensation model
    real(c_double) :: max_coating_volume

    contact_angle = 0.0         ! Contact angle of the coating material (degrees)
    neck_filling_angle = 10.0   ! Filling angle of the neck between monomers (degrees)
    r_part = 1.0                ! Radius of primary particles
    condensate_volume = 0.1     ! Total volume of condensate

    ! Initialize the constant mean curvature surface object
    constant_mean_curvature_surface = initialize_constant_mean_curvature_surface(contact_angle, neck_filling_angle, r_part)

    ! The simulation loop that uses the geometry model starts here...

    ! Compute surface area and curvature as a function of condensate volume
    props = get_liquid_props(constant_mean_curvature_surface, condensate_volume)
    ! Get the maximum coating volume that can be accommodated by the model
    max_coating_volume = get_max_liquid_volume(constant_mean_curvature_surface)

    print *, 'Surface area:', props%area, 'Curvature:', props%kappa, 'Max volume:', max_coating_volume

    ! The simulation loop ends here and the geometry model is no longer needed. Clean up...

    ! Delete the constant mean curvature surface object when it's no longer needed
    ! IMPORTANT to avoid memory leaks
    call delete_constant_mean_curvature_surface(constant_mean_curvature_surface)

end program main