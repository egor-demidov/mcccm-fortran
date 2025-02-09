# Fortran binding for capillary meniscus geometry

This is a Fortran binding for the part of the 
[MCCCM C++ library](https://github.com/egor-demidov/mcccm)
that computes the area and curvature of the capillary
meniscus as a function of condensate volume.

## API reference

The Fortran module `geometry` defines the following functions:

### initialize_constant_mean_curvature_surface()

Function; Allocates and initializes a surface object on the heap

#### Arguments

- `c_double` - contact angle, deg
- `c_double` - neck filling angle, deg
- `c_double` - primary particle radius, m

#### Return value

- `c_ptr` - pointer to the surface object

### delete_constant_mean_curvature_surface()

Subroutine; Deletes the surface object (memory is freed) and sets the pointer to `nullptr`

> [!WARNING]
> This function must be called as soon the geometry model object is no longer needed to avoid memory leaks

#### Arguments

- `c_ptr` - pointer to the surface object

### get_liquid_props()

Function; Returns a geometry_props structure (with surface area and curvature)

#### Arguments

- `c_ptr` - pointer to the surface objet
- `c_double` - total volume of coating in the gap, m^3

#### Return value

- `geometry_props` - structure with area and curvature

### get_max_liquid_volume()

Function; Returns the maximum liquid volume that can be accommodated by the geometry model

#### Arguments

- `c_ptr` - pointer to the surface object

#### Return value

- `c_double` - maximum liquid volume, m^3

### get_filling_angle()

Function; Returns the filling angle as a function of coating volume

#### Arguments

- `c_ptr` - pointer to the surface object
- `c_double` - total condensate volume in the gap, m^3

#### Return value

- `c_double` - filling angle, deg

### get_neck_volume()

Function; Returns the volume of the neck

#### Arguments

- `c_ptr` - pointer to the surface object

#### Return value

- `c_double` - neck volume, m^3

### get_r_part()

Function; Returns the primary particle radius

#### Arguments

- `c_ptr` - pointer to the surface object

#### Return value

- `c_double` - primary particle radius, m

## Example program

```fortran
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

    ! The simulation loop that uses the geometry model will go here...

    ! Compute surface area and curvature as a function of condensate volume
    props = get_liquid_props(constant_mean_curvature_surface, condensate_volume)
    ! Get the maximum coating volume that can be accommodated by the model
    max_coating_volume = get_max_liquid_volume(constant_mean_curvature_surface)

    print *, 'Surface area:', props%area, 'Curvature:', props%kappa, 'Max volume:', max_coating_volume

    ! The simulation loop has completed and the geometry model is no longer needed. Clean up...

    ! Delete the constant mean curvature surface object when it's no longer needed
    ! IMPORTANT to avoid memory leaks
    call delete_constant_mean_curvature_surface(constant_mean_curvature_surface)

end program main
```

## Compiling the example program manually

To manually compile the example program on Linux with the GNU compiler, run the following commands from
the project root directory:

```bash
# Create a build directory for an out-of-tree build
mkdir -p build
cd build

# Compile the Fortran sources
gfortran -c ../geometry.f90
gfortran -c ../main.f90
# Compile the C++ sources
g++ -c ../geometry.cpp -o geometry_cxx.o    # custom object file name to prevent name overlaps
g++ -c ../generated_geometry.cpp
g++ -c ../geometry_c_binding.cpp
# Link the object files into an executable (-lgfortran is required to link the Fortran standard library)
#g++ *.o -lgfortran -o main
# Alternatively, use gfortran (-lstdc++ is required to link the C++ standard library)
gfortran *.o -lstdc++ -o main

# Return to source directory
cd ..
```
