//
// Created by egor on 2/7/25.
//

#include "geometry.h"
#include "geometry_c_binding.h"

inline geometry_interfaces::ConstantMeanCurvatureSurface * to_geom_ptr(void * void_ptr) {
    return static_cast<geometry_interfaces::ConstantMeanCurvatureSurface *>(void_ptr);
}

void * initialize_constant_mean_curvature_surface(double contact_angle, double neck_filling_angle, double r_part) {
    return new geometry_interfaces::ConstantMeanCurvatureSurface(contact_angle, neck_filling_angle, r_part);
}

void delete_constant_mean_curvature_surface(void ** constant_mean_curvature_surface) {
    delete to_geom_ptr(*constant_mean_curvature_surface);
    *constant_mean_curvature_surface = nullptr;
}

GeometryProps get_liquid_props(void * constant_mean_curvature_surface, double condensate_volume) {
    auto geom_ptr = to_geom_ptr(constant_mean_curvature_surface);
    auto [area, kappa] = geom_ptr->get_liquid_props(condensate_volume);
    return {area, kappa};
}

double get_max_liquid_volume(void * constant_mean_curvature_surface) {
    auto geom_ptr = to_geom_ptr(constant_mean_curvature_surface);
    return geom_ptr->get_max_liquid_volume();
}

double get_filling_angle(void * constant_mean_curvature_surface, double condensate_volume) {
    auto geom_ptr = to_geom_ptr(constant_mean_curvature_surface);
    return geom_ptr->get_filling_angle(condensate_volume);
}

double get_neck_volume(void * constant_mean_curvature_surface) {
    auto geom_ptr = to_geom_ptr(constant_mean_curvature_surface);
    return geom_ptr->get_neck_volume();
}

double get_r_part(void * constant_mean_curvature_surface) {
    auto geom_ptr = to_geom_ptr(constant_mean_curvature_surface);
    return geom_ptr->get_r_part();
}
