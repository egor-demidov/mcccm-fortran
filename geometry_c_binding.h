//
// Created by egor on 2/7/25.
//

#ifndef MCCCM_FORTRAN_GEOMETRY_C_BINDING_H
#define MCCCM_FORTRAN_GEOMETRY_C_BINDING_H

#ifdef __cplusplus
#define EXTERNC extern "C"
#else
#define EXTERNC
#endif

struct GeometryProps {
    double area, kappa;
};

EXTERNC void * initialize_constant_mean_curvature_surface(double contact_angle, double neck_filling_angle, double r_part);
EXTERNC void delete_constant_mean_curvature_surface(void ** constant_mean_curvature_surface);
EXTERNC GeometryProps get_liquid_props(void * constant_mean_curvature_surface, double condensate_volume);
EXTERNC double get_max_liquid_volume(void * constant_mean_curvature_surface);
EXTERNC double get_filling_angle(void * constant_mean_curvature_surface, double condensate_volume);
EXTERNC double get_neck_volume(void * constant_mean_curvature_surface);
EXTERNC double get_r_part(void * constant_mean_curvature_surface);

#endif //MCCCM_FORTRAN_GEOMETRY_C_BINDING_H
