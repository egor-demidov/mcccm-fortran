# Fortran binding for capillary meniscus geometry

This is a Fortran binding to the part of the 
[MCCCM C++ library](https://github.com/egor-demidov/mcccm)
that computes the volume and curvature of the capillary
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
