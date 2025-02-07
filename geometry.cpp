//
// Created by egor on 12/14/24.
//

#include "geometry.h"

#include <fmt/format.h>

GeometryInterpolator::GeometryInterpolator(double contact_angle)
    : ca{contact_angle}
    , dca{(generated_geometry.end_ca - generated_geometry.begin_ca) / static_cast<double>(Geometry_t::n_ca - 1)}
{
    ca_interpolation.index_lo = static_cast<unsigned long>((contact_angle - generated_geometry.begin_ca) / dca);
    ca_interpolation.index_hi = ca_interpolation.index_lo + 1;
    ca_interpolation.weight_hi = (contact_angle - generated_geometry.ca_sets[ca_interpolation.index_lo].ca) / dca;
    ca_interpolation.weight_lo = 1.0 - ca_interpolation.weight_hi;

    // contact angle must be greater than or equal to the first point
    // and less than the last point in the generated geometry
    if (contact_angle < generated_geometry.ca_sets[ca_interpolation.index_lo].ca || contact_angle >= generated_geometry.ca_sets[ca_interpolation.index_hi].ca)
        throw InvalidContactAngleException();
}

double GeometryInterpolator::get_filling_angle(unsigned long index) const {
    if (index >= Geometry_t::n_fa)
        throw IndexOutOfBoundsException();

    return generated_geometry.ca_sets[ca_interpolation.index_lo].filling_angle[index] * ca_interpolation.weight_lo
        + generated_geometry.ca_sets[ca_interpolation.index_hi].filling_angle[index] * ca_interpolation.weight_hi;
}

double GeometryInterpolator::get_volume(unsigned long index) const {
    if (index >= Geometry_t::n_fa)
        throw IndexOutOfBoundsException();

    return generated_geometry.ca_sets[ca_interpolation.index_lo].volume[index] * ca_interpolation.weight_lo
        + generated_geometry.ca_sets[ca_interpolation.index_hi].volume[index] * ca_interpolation.weight_hi;
}

double GeometryInterpolator::get_area(unsigned long index) const {
    if (index >= Geometry_t::n_fa)
        throw IndexOutOfBoundsException();

    return generated_geometry.ca_sets[ca_interpolation.index_lo].area[index] * ca_interpolation.weight_lo
        + generated_geometry.ca_sets[ca_interpolation.index_hi].area[index] * ca_interpolation.weight_hi;
}

double GeometryInterpolator::get_kappa(unsigned long index) const {
    if (index >= Geometry_t::n_fa)
        throw IndexOutOfBoundsException();

    return generated_geometry.ca_sets[ca_interpolation.index_lo].kappa[index] * ca_interpolation.weight_lo
        + generated_geometry.ca_sets[ca_interpolation.index_hi].kappa[index] * ca_interpolation.weight_hi;
}

std::pair<double, LinearInterpolation> GeometryInterpolator::volume_to_filling_angle(double volume) const {
    unsigned long index_hi = 0;

    do {
        index_hi ++;
    } while (volume > get_volume(index_hi) && index_hi < Geometry_t::n_fa - 1);

    if (index_hi == 0 || index_hi == Geometry_t::n_fa)
        throw VolumeOutOfBoundsException();

    unsigned long index_lo = index_hi - 1;
    double dvolume = get_volume(index_hi) - get_volume(index_lo);
    double weight_hi = (volume - get_volume(index_lo)) / dvolume;
    double weight_lo = 1.0 - weight_hi;
    double filling_angle = get_filling_angle(index_lo) * weight_lo + get_filling_angle(index_hi) * weight_hi;
    return std::make_pair(
        filling_angle,
        LinearInterpolation{index_lo, index_hi, weight_lo, weight_hi}
    );
}

double GeometryInterpolator::interpolate_area(LinearInterpolation const & filling_angle_interpolation) const {
    return get_area(filling_angle_interpolation.index_lo) * filling_angle_interpolation.weight_lo
        + get_area(filling_angle_interpolation.index_hi) * filling_angle_interpolation.weight_hi;
}

double GeometryInterpolator::interpolate_kappa(LinearInterpolation const & filling_angle_interpolation) const {
    return get_kappa(filling_angle_interpolation.index_lo) * filling_angle_interpolation.weight_lo
        + get_kappa(filling_angle_interpolation.index_hi) * filling_angle_interpolation.weight_hi;
}

LinearInterpolation GeometryInterpolator::get_filling_angle_interpolation(double filling_angle) const {
    unsigned long index_hi = 0;

    do {
        index_hi ++;
    } while (filling_angle > get_filling_angle(index_hi) && index_hi < Geometry_t::n_fa - 1);

    if (index_hi == 0 || index_hi == Geometry_t::n_fa)
        throw FillingAngleOutOfBoundsException();

    unsigned long index_lo = index_hi - 1;
    double dfa = get_filling_angle(index_hi) - get_filling_angle(index_lo);
    double weight_hi = (filling_angle - get_filling_angle(index_lo)) / dfa;
    double weight_lo = 1.0 - weight_hi;
    return {index_lo, index_hi, weight_lo, weight_hi};
}

double GeometryInterpolator::interpolate_volume(double filling_angle) const {
    LinearInterpolation fa_interpolation = get_filling_angle_interpolation(filling_angle);
    return get_volume(fa_interpolation.index_lo) * fa_interpolation.weight_lo
        + get_volume(fa_interpolation.index_hi) * fa_interpolation.weight_hi;
}

double GeometryInterpolator::interpolate_area(double filling_angle) const {
    LinearInterpolation fa_interpolation = get_filling_angle_interpolation(filling_angle);
    return get_area(fa_interpolation.index_lo) * fa_interpolation.weight_lo
        + get_area(fa_interpolation.index_hi) * fa_interpolation.weight_hi;
}

double GeometryInterpolator::interpolate_kappa(double filling_angle) const {
    LinearInterpolation fa_interpolation = get_filling_angle_interpolation(filling_angle);
    return get_kappa(fa_interpolation.index_lo) * fa_interpolation.weight_lo
        + get_kappa(fa_interpolation.index_hi) * fa_interpolation.weight_hi;
}

geometry_interfaces::ConstantMeanCurvatureSurface::ConstantMeanCurvatureSurface(
    double contact_angle,
    double neck_filling_angle,
    double r_part
)
    : contact_angle{contact_angle}
    , neck_filling_angle{neck_filling_angle}
    , r_part{r_part}
    , liquid_interpolator{contact_angle}
{
    GeometryInterpolator neck_interpolator(0.0);
    neck_volume = neck_interpolator.interpolate_volume(neck_filling_angle) * r_part * r_part * r_part;
    max_liquid_volume = liquid_interpolator.interpolate_volume(generated_geometry.end_fa) * r_part * r_part * r_part - neck_volume;
}

geometry_interfaces::GeometryProps geometry_interfaces::ConstantMeanCurvatureSurface::get_liquid_props(double condensate_volume) const {
    double total_volume = condensate_volume + neck_volume;
    auto [filling_angle, interpolation] = liquid_interpolator.volume_to_filling_angle(total_volume / (r_part * r_part * r_part));
    return {
        .area = liquid_interpolator.interpolate_area(interpolation) * r_part * r_part,
        .kappa = liquid_interpolator.interpolate_kappa(interpolation) / r_part
    };
}

double geometry_interfaces::ConstantMeanCurvatureSurface::get_neck_volume() const {
    return neck_volume;
}

double geometry_interfaces::ConstantMeanCurvatureSurface::get_max_liquid_volume() const {
    return max_liquid_volume;
}


double geometry_interfaces::ConstantMeanCurvatureSurface::get_filling_angle(double condensate_volume) const {
    double total_volume = condensate_volume + neck_volume;
    auto [filling_angle, interpolation] = liquid_interpolator.volume_to_filling_angle(total_volume / (r_part * r_part * r_part));
    return filling_angle;
}

double geometry_interfaces::ConstantMeanCurvatureSurface::get_r_part() const {
    return r_part;
}

geometry_interfaces::SphericalSurface::SphericalSurface(double r_part)
    : r_part{r_part}
    , core_volume{4.0 / 3.0 * M_PI * r_part * r_part * r_part}
{}

geometry_interfaces::GeometryProps geometry_interfaces::SphericalSurface::get_liquid_props(double condensate_volume) const {
    double area = 4.0 * M_PI * r_part * r_part;
    double r_equivalent = get_equivalent_radius(condensate_volume);
    double kappa = 1.0 / r_equivalent;
    return {area, kappa};
}

double geometry_interfaces::SphericalSurface::get_core_volume() const {
    return core_volume;
}

double geometry_interfaces::SphericalSurface::get_max_liquid_volume() const { // NOLINT
    return max_liquid_volume;
}

double geometry_interfaces::SphericalSurface::get_equivalent_radius(double condensate_volume) const {
    double total_volume = condensate_volume + core_volume;
    double r_equivalent = pow(3.0 * total_volume / 4.0 / M_PI, 1.0 / 3.0);
    return r_equivalent;
}

double geometry_interfaces::SphericalSurface::get_r_part() const {
    return r_part;
}

// Unit tests below, compiled only is this is a test target
// available only in a native build
#ifdef DO_TEST

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>


TEST_CASE("GeometryInterpolator tested", "[GeometryInterpolator]") {
    // Test filling andle getters four boundary cases
    GeometryInterpolator test_interpolator(0.0);
    REQUIRE(test_interpolator.get_filling_angle(0) == generated_geometry.begin_fa);
    REQUIRE(test_interpolator.get_filling_angle(Geometry_t::n_fa - 1) == generated_geometry.end_fa);

    // Test getters
    unsigned long mid_index = Geometry_t::n_fa / 2ul;
    REQUIRE(generated_geometry.ca_sets[0].filling_angle[mid_index] == test_interpolator.get_filling_angle(mid_index));
    REQUIRE(generated_geometry.ca_sets[0].area[mid_index] == test_interpolator.get_area(mid_index));
    REQUIRE(generated_geometry.ca_sets[0].volume[mid_index] == test_interpolator.get_volume(mid_index));
    REQUIRE(generated_geometry.ca_sets[0].kappa[mid_index] == test_interpolator.get_kappa(mid_index));

    // Test interpolators
    double vol_0 = test_interpolator.get_volume(0);
    double vol_1 = test_interpolator.get_volume(1);
    double fa_0 = test_interpolator.get_filling_angle(0);
    double fa_1 = test_interpolator.get_filling_angle(1);
    double dfa = fa_1 - fa_0;
    double fa_x = fa_0 + 0.3 * dfa;
    double weight_1 = (fa_x - fa_0) / dfa;
    double weight_0 = 1.0 - weight_1;
    double vol_x = vol_0 * weight_0 + vol_1 * weight_1;
    REQUIRE(vol_x == test_interpolator.interpolate_volume(fa_x));

    double area_0 = test_interpolator.get_area(0);
    double area_1 = test_interpolator.get_area(1);
    double area_x = area_0 * weight_0 + area_1 * weight_1;
    REQUIRE(area_x == test_interpolator.interpolate_area(fa_x));

    double kappa_0 = test_interpolator.get_kappa(0);
    double kappa_1 = test_interpolator.get_kappa(1);
    double kappa_x = kappa_0 * weight_0 + kappa_1 * weight_1;
    REQUIRE(kappa_x == test_interpolator.interpolate_kappa(fa_x));

    // Test FA calculation from volume
    auto [fa_back, fa_interp] = test_interpolator.volume_to_filling_angle(vol_x);
    REQUIRE(fa_back == fa_x);
    double area_back = test_interpolator.interpolate_area(fa_interp);
    double kappa_back = test_interpolator.interpolate_kappa(fa_interp);
    REQUIRE_THAT(area_back, Catch::Matchers::WithinAbs(area_x, 0.000001));
    REQUIRE_THAT(kappa_back, Catch::Matchers::WithinAbs(kappa_x, 0.000001));

    // Test contact angle interpolation
    double dca = (generated_geometry.end_ca - generated_geometry.begin_ca) / static_cast<double>(Geometry_t::n_ca - 1);
    double ca_0 = generated_geometry.begin_ca;
    double ca_1 = ca_0 + dca;
    double ca_x = ca_0 + 0.3 * dca;
    weight_1 = (ca_x - ca_0) / dca;
    weight_0 = 1.0 - weight_1;
    double fa_mid = generated_geometry.ca_sets[ca_0].filling_angle[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].filling_angle[mid_index] * weight_1;
    double vol_mid = generated_geometry.ca_sets[ca_0].volume[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].volume[mid_index] * weight_1;
    double area_mid = generated_geometry.ca_sets[ca_0].area[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].area[mid_index] * weight_1;
    double kappa_mid = generated_geometry.ca_sets[ca_0].kappa[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].kappa[mid_index] * weight_1;
    GeometryInterpolator test_interpolator_2(ca_x);
    REQUIRE(fa_mid == test_interpolator_2.get_filling_angle(mid_index));
    REQUIRE(vol_mid == test_interpolator_2.get_volume(mid_index));
    REQUIRE(area_mid == test_interpolator_2.get_area(mid_index));
    REQUIRE(kappa_mid == test_interpolator_2.get_kappa(mid_index));

    // Test high contact and filling angles
    ca_1 = generated_geometry.end_ca;
    ca_0 = ca_1 - dca;
    ca_x = ca_0 + 0.3 * dca;
    weight_1 = (ca_x - ca_0) / dca;
    weight_0 = 1.0 - weight_1;
    fa_mid = generated_geometry.ca_sets[ca_0].filling_angle[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].filling_angle[mid_index] * weight_1;
    vol_mid = generated_geometry.ca_sets[ca_0].volume[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].volume[mid_index] * weight_1;
    area_mid = generated_geometry.ca_sets[ca_0].area[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].area[mid_index] * weight_1;
    kappa_mid = generated_geometry.ca_sets[ca_0].kappa[mid_index] * weight_0 + generated_geometry.ca_sets[ca_1].kappa[mid_index] * weight_1;
    GeometryInterpolator test_interpolator_3(ca_x);
    REQUIRE(fa_mid == test_interpolator_3.get_filling_angle(mid_index));
    REQUIRE(vol_mid == test_interpolator_3.get_volume(mid_index));
    REQUIRE(area_mid == test_interpolator_3.get_area(mid_index));
    REQUIRE(kappa_mid == test_interpolator_3.get_kappa(mid_index));

    fa_0 = test_interpolator_3.get_filling_angle(Geometry_t::n_fa - 2);
    fa_1 = test_interpolator_3.get_filling_angle(Geometry_t::n_fa - 1);
    dfa = fa_1 - fa_0;
    fa_x = fa_0 + 0.3 * dfa;
    weight_1 = (fa_x - fa_0) / dfa;
    weight_0 = 1.0 - weight_1;
    vol_0 = test_interpolator_3.get_volume(Geometry_t::n_fa - 2);
    vol_1 = test_interpolator_3.get_volume(Geometry_t::n_fa - 1);
    vol_x = vol_0 * weight_0 + vol_1 * weight_1;
    REQUIRE(vol_x == test_interpolator_3.interpolate_volume(fa_x));

    area_0 = test_interpolator_3.get_area(Geometry_t::n_fa - 2);
    area_1 = test_interpolator_3.get_area(Geometry_t::n_fa - 1);
    area_x = area_0 * weight_0 + area_1 * weight_1;
    REQUIRE(area_x == test_interpolator_3.interpolate_area(fa_x));

    kappa_0 = test_interpolator_3.get_kappa(Geometry_t::n_fa - 2);
    kappa_1 = test_interpolator_3.get_kappa(Geometry_t::n_fa - 1);
    kappa_x = kappa_0 * weight_0 + kappa_1 * weight_1;
    REQUIRE(kappa_x == test_interpolator_3.interpolate_kappa(fa_x));
}

TEST_CASE("generated_geometry tested", "[generated_geometry]") {
    // Test that area and volume are increasing with filling angle
    bool volume_increasing = true;
    bool area_increasing = true;
    for (unsigned long ca_index = 0; ca_index < Geometry_t::n_ca; ca_index ++) {
        double vol_prev = generated_geometry.ca_sets[ca_index].volume[0];
        double area_prev = generated_geometry.ca_sets[ca_index].area[0];
        for (unsigned long fa_index = 1; fa_index < Geometry_t::n_fa; fa_index ++) {
            double vol_new = generated_geometry.ca_sets[ca_index].volume[fa_index];
            double area_new = generated_geometry.ca_sets[ca_index].area[fa_index];
            double kappa_new = generated_geometry.ca_sets[ca_index].kappa[fa_index];
            if (vol_new < vol_prev)
                volume_increasing = false;
            if (area_new < area_prev)
                area_increasing = false;
            vol_prev = vol_new;
            area_prev = area_new;
        }
    }
    REQUIRE(volume_increasing);
    REQUIRE(area_increasing);

    // Test that 0 deg CA and 90 deg FA is a cylinder (if such configuration is available)
    if (generated_geometry.begin_ca <= 0 && generated_geometry.end_fa >= 90.0) {
        double target_volume = 2.0 * M_PI - 4.0 / 3.0 * M_PI;
        double target_area = 4.0 * M_PI;
        double target_kappa = 1.0 / 2.0;
        GeometryInterpolator interp(0.0);
        REQUIRE_THAT(target_volume, Catch::Matchers::WithinAbs(interp.interpolate_volume(90.0), 0.000001));
        REQUIRE_THAT(target_area, Catch::Matchers::WithinAbs(interp.interpolate_area(90.0), 0.000001));
        REQUIRE_THAT(target_kappa, Catch::Matchers::WithinAbs(interp.interpolate_kappa(90.0), 0.000001));
    }

    // TODO: add catenoid test
}

TEST_CASE("ConstantMeanCurvatureSurface tested", "[ConstantMeanCurvatureSurface]") {
    GeometryInterpolator neck_interpolator(0.0);
    double r_part = 0.1;
    double neck_volume = neck_interpolator.interpolate_volume(10.0) * r_part * r_part * r_part;
    double liquid_volume = 2.0 * r_part * M_PI * r_part * r_part - 4.0 / 3.0 * M_PI * r_part * r_part * r_part;
    double target_area = 2.0 * M_PI * r_part * 2.0 * r_part;
    double target_kappa = 1.0 / 2.0 / r_part;
    double condensate_volume = liquid_volume - neck_volume;

    geometry_interfaces::ConstantMeanCurvatureSurface surface_interface(0.0, 10.0, r_part);
    auto [area, kappa] = surface_interface.get_liquid_props(condensate_volume);
    REQUIRE_THAT(target_area, Catch::Matchers::WithinAbs(area, 0.000001));
    REQUIRE_THAT(target_kappa, Catch::Matchers::WithinAbs(kappa, 0.000001));
    REQUIRE_THAT(neck_volume, Catch::Matchers::WithinAbs(surface_interface.get_neck_volume(), 0.000001));

    // Test filling angle
    double filling_angle = surface_interface.get_filling_angle(condensate_volume);
    REQUIRE_THAT(filling_angle, Catch::Matchers::WithinAbs(90.0, 0.000001));
}

TEST_CASE("SphericalSurface tested", "[SphericalSurface]") {
    double r_part = 0.1;

    geometry_interfaces::SphericalSurface surface_interface(r_part);

    double condensate_volume = 0.0;
    double target_equivalent_radius = r_part;
    double target_area = 4.0 * M_PI * target_equivalent_radius * target_equivalent_radius;
    double target_kappa = 1.0 / target_equivalent_radius;
    auto [area, kappa] = surface_interface.get_liquid_props(condensate_volume);
    REQUIRE_THAT(target_area, Catch::Matchers::WithinAbs(area, 0.000001));
    REQUIRE_THAT(target_kappa, Catch::Matchers::WithinAbs(kappa, 0.000001));
}

#endif //DO_TEST
