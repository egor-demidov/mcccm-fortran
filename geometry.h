//
// Created by Egor Demidov on 12/12/24.
//

#ifndef WASM_CAPILLARY_CONDENSATION_GEOMETRY_H
#define WASM_CAPILLARY_CONDENSATION_GEOMETRY_H

#include <array>
#include <exception>
#include <numeric>
#include <limits>

template <unsigned long num_fa, unsigned long num_ca>
struct Geometry {

    struct CASeries {
        const double ca;
        const std::array<double, num_fa> filling_angle, volume, area, kappa;
    };

    const double begin_fa, end_fa;
    const double begin_ca, end_ca;
    const std::array<CASeries, num_ca> ca_sets;

    static constexpr unsigned long n_fa = num_fa;
    static constexpr unsigned long n_ca = num_ca;
};

#include "generated_geometry.h"

class IndexOutOfBoundsException : public std::exception {
public:
    const char * what() const noexcept override {
        return "Index passed to a getter is out of bounds";
    }
};

class InvalidContactAngleException : public std::exception {
public:
    const char * what() const noexcept override {
        return "Invalid contact angle provided to GeometryInterpolator";
    }
};

class VolumeOutOfBoundsException : public std::exception {
public:
    const char * what() const noexcept override {
        return "Volume of condensate is out of bounds";
    }
};

class FillingAngleOutOfBoundsException : public std::exception {
public:
    const char * what() const noexcept override {
        return "Filling angle of condensate is out of bounds";
    }
};

struct LinearInterpolation {
    unsigned long index_lo, index_hi;
    double weight_lo, weight_hi;
};

class GeometryInterpolator {
public:
    explicit GeometryInterpolator(double contact_angle);

    double get_filling_angle(unsigned long index) const;
    double get_volume(unsigned long index) const;
    double get_area(unsigned long index) const;
    double get_kappa(unsigned long index) const;

    std::pair<double, LinearInterpolation> volume_to_filling_angle(double volume) const;
    double interpolate_area(LinearInterpolation const & filling_angle_interpolation) const;
    double interpolate_kappa(LinearInterpolation const & filling_angle_interpolation) const;
    LinearInterpolation get_filling_angle_interpolation(double filling_angle) const;
    double interpolate_volume(double filling_angle) const;
    double interpolate_area(double filling_angle) const;
    double interpolate_kappa(double filling_angle) const;

private:
    const double ca, dca;
    LinearInterpolation ca_interpolation;
};

namespace geometry_interfaces {

    struct GeometryProps {
        double area, kappa;
    };

    class ConstantMeanCurvatureSurface {
    public:
        ConstantMeanCurvatureSurface(double contact_angle, double neck_filling_angle, double r_part);

        // These methods are required for any geometry interface
        GeometryProps get_liquid_props(double condensate_volume) const;
        double get_max_liquid_volume() const;

        // This method is specific to the capillary condensation interface
        double get_filling_angle(double condensate_volume) const;
        double get_neck_volume() const;
        double get_r_part() const;

    private:
        const double contact_angle, neck_filling_angle, r_part;
        const GeometryInterpolator liquid_interpolator;
        double neck_volume;
        double max_liquid_volume;
    };

    class SphericalSurface {
    public:
        explicit SphericalSurface(double r_part);

        // These methods are required for any geometry interface
        GeometryProps get_liquid_props(double condensate_volume) const;
        double get_max_liquid_volume() const;

        // This method is specific to the spherical condensation interface
        double get_core_volume() const;
        double get_equivalent_radius(double condensate_volume) const;
        double get_r_part() const;

    private:
        const double r_part;
        const double core_volume;
        static constexpr double max_liquid_volume = std::numeric_limits<double>::infinity();
    };

}

#endif //WASM_CAPILLARY_CONDENSATION_GEOMETRY_H
