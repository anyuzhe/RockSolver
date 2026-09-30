#include "rocksolver/rocksolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <charconv>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numbers>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <unordered_map>

namespace rocksolver {
namespace {

constexpr double kTiny = 1.0e-18;

Vec3 operator+(const Vec3& a, const Vec3& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 operator-(const Vec3& a, const Vec3& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

Vec3 operator*(const Vec3& value, double scale) {
    return {value.x * scale, value.y * scale, value.z * scale};
}

Vec3& operator+=(Vec3& a, const Vec3& b) {
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

double dot(const Vec3& a, const Vec3& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x,
    };
}

double norm(const Vec3& value) { return std::sqrt(dot(value, value)); }

struct TetGeometry {
    double volume{};
    std::array<Vec3, 4> gradients{};
};

double determinant(const std::array<std::array<double, 3>, 3>& m) {
    return m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
           m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
           m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
}

std::array<std::array<double, 3>, 3> inverse(
    const std::array<std::array<double, 3>, 3>& m) {
    const double det = determinant(m);
    if (std::abs(det) <= kTiny) {
        throw std::runtime_error("RockSolver encountered a zero-volume tetrahedron");
    }
    const double inv = 1.0 / det;
    return {{{
                 (m[1][1] * m[2][2] - m[1][2] * m[2][1]) * inv,
                 (m[0][2] * m[2][1] - m[0][1] * m[2][2]) * inv,
                 (m[0][1] * m[1][2] - m[0][2] * m[1][1]) * inv,
             },
             {
                 (m[1][2] * m[2][0] - m[1][0] * m[2][2]) * inv,
                 (m[0][0] * m[2][2] - m[0][2] * m[2][0]) * inv,
                 (m[0][2] * m[1][0] - m[0][0] * m[1][2]) * inv,
             },
             {
                 (m[1][0] * m[2][1] - m[1][1] * m[2][0]) * inv,
                 (m[0][1] * m[2][0] - m[0][0] * m[2][1]) * inv,
                 (m[0][0] * m[1][1] - m[0][1] * m[1][0]) * inv,
             }}};
}

TetGeometry tetra_geometry(const RockModel& model, const Zone& zone) {
    const Vec3& p0 = model.nodes[zone.nodes[0]].position;
    const Vec3& p1 = model.nodes[zone.nodes[1]].position;
    const Vec3& p2 = model.nodes[zone.nodes[2]].position;
    const Vec3& p3 = model.nodes[zone.nodes[3]].position;
    const std::array<std::array<double, 3>, 3> matrix{{
        {{p1.x - p0.x, p2.x - p0.x, p3.x - p0.x}},
        {{p1.y - p0.y, p2.y - p0.y, p3.y - p0.y}},
        {{p1.z - p0.z, p2.z - p0.z, p3.z - p0.z}},
    }};
    const double det = determinant(matrix);
    const auto inv = inverse(matrix);
    TetGeometry geometry;
    geometry.volume = std::abs(det) / 6.0;
    geometry.gradients[1] = {inv[0][0], inv[0][1], inv[0][2]};
    geometry.gradients[2] = {inv[1][0], inv[1][1], inv[1][2]};
    geometry.gradients[3] = {inv[2][0], inv[2][1], inv[2][2]};
    geometry.gradients[0] = (geometry.gradients[1] + geometry.gradients[2] +
                             geometry.gradients[3]) *
                            -1.0;
    return geometry;
}

std::array<double, 3> principal_values(const std::array<double, 6>& stress) {
    std::array<std::array<double, 3>, 3> a{{
        {{stress[0], stress[3], stress[4]}},
        {{stress[3], stress[1], stress[5]}},
        {{stress[4], stress[5], stress[2]}},
    }};
    for (int iteration = 0; iteration < 18; ++iteration) {
        int p = 0;
        int q = 1;
        double largest = std::abs(a[0][1]);
        for (const auto pair : {std::pair{0, 2}, std::pair{1, 2}}) {
            const double value = std::abs(a[pair.first][pair.second]);
            if (value > largest) {
                largest = value;
                p = pair.first;
                q = pair.second;
            }
        }
        if (largest < 1.0e-10) break;
        const double angle = 0.5 * std::atan2(2.0 * a[p][q], a[q][q] - a[p][p]);
        const double c = std::cos(angle);
        const double s = std::sin(angle);
        const double app = c * c * a[p][p] - 2.0 * s * c * a[p][q] +
                           s * s * a[q][q];
        const double aqq = s * s * a[p][p] + 2.0 * s * c * a[p][q] +
                           c * c * a[q][q];
        for (int r = 0; r < 3; ++r) {
            if (r == p || r == q) continue;
            const double arp = c * a[r][p] - s * a[r][q];
            const double arq = s * a[r][p] + c * a[r][q];
            a[r][p] = a[p][r] = arp;
            a[r][q] = a[q][r] = arq;
        }
        a[p][p] = app;
        a[q][q] = aqq;
        a[p][q] = a[q][p] = 0.0;
    }
    std::array<double, 3> values{a[0][0], a[1][1], a[2][2]};
    std::sort(values.begin(), values.end());
    return values;
}

std::array<double, 6> total_zone_strain(const RockModel& model,
                                        const Zone& zone) {
    const TetGeometry geometry = tetra_geometry(model, zone);
    std::array<double, 6> strain{};
    for (std::size_t i = 0; i < 4; ++i) {
        const Vec3& u = model.nodes[zone.nodes[i]].displacement;
        const Vec3& g = geometry.gradients[i];
        strain[0] += u.x * g.x;
        strain[1] += u.y * g.y;
        strain[2] += u.z * g.z;
        strain[3] += u.x * g.y + u.y * g.x;
        strain[4] += u.x * g.z + u.z * g.x;
        strain[5] += u.y * g.z + u.z * g.y;
    }
    return strain;
}

std::vector<std::array<double, 6>> nodal_mixed_strains(
    const RockModel& model) {
    std::vector<std::array<double, 6>> strains(model.zones.size());
    std::vector<double> nodal_mean(model.nodes.size(), 0.0);
    std::vector<double> nodal_volume(model.nodes.size(), 0.0);
    std::vector<double> volumes(model.zones.size(), 0.0);

    // FLAC3D NMD first volume-weights each surrounding Tet4 mean strain at
    // every node.  Mean strain is trace(epsilon)/3; engineering shear terms
    // are unaffected by the volumetric projection.
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const Zone& zone = model.zones[zi];
        if (!zone.active) continue;
        strains[zi] = total_zone_strain(model, zone);
        const double volume = tetra_geometry(model, zone).volume * zone.integration_weight;
        volumes[zi] = volume;
        const double mean =
            (strains[zi][0] + strains[zi][1] + strains[zi][2]) / 3.0;
        for (const std::size_t node : zone.nodes) {
            nodal_mean[node] += mean * volume;
            nodal_volume[node] += volume;
        }
    }
    for (std::size_t node = 0; node < nodal_mean.size(); ++node) {
        if (nodal_volume[node] > kTiny) nodal_mean[node] /= nodal_volume[node];
    }

    // The element mean is the arithmetic average of its four nodal values.
    // Replace only the volumetric component and preserve the raw deviator.
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const Zone& zone = model.zones[zi];
        if (!zone.active || volumes[zi] <= kTiny) continue;
        const double raw_mean =
            (strains[zi][0] + strains[zi][1] + strains[zi][2]) / 3.0;
        double mixed_mean = 0.0;
        for (const std::size_t node : zone.nodes) mixed_mean += nodal_mean[node];
        mixed_mean *= 0.25;
        const double correction = mixed_mean - raw_mean;
        strains[zi][0] += correction;
        strains[zi][1] += correction;
        strains[zi][2] += correction;
    }
    return strains;
}

const MaterialParameters& zone_material(const RockModel& model, const Zone& zone) {
    return zone.material_index == 0 ? model.config.material
        : model.config.region_materials.at(zone.material_index - 1).material;
}

void validate_water(const ModelConfig& config) {
    const auto& water=config.water;
    if(config.seepage.enabled &&
       (!std::isfinite(config.seepage.biot_coefficient) ||
        config.seepage.biot_coefficient<0.0 || config.seepage.biot_coefficient>1.0 ||
        !std::isfinite(config.seepage.biot_modulus_pa) ||
        config.seepage.biot_modulus_pa<0.0))
        throw std::invalid_argument("Biot coefficient must be in [0,1] and modulus nonnegative");
    if(water.external_enabled) {
        const auto& direction=water.external_outward_direction;
        if(!std::isfinite(water.external_surface_z_m) ||
           !std::isfinite(water.external_density_kg_m3) ||
           water.external_density_kg_m3<=0.0 ||
           !std::isfinite(config.gravity_m_s2) || config.gravity_m_s2<=0.0 ||
           !std::isfinite(direction.x) || !std::isfinite(direction.y) ||
           !std::isfinite(direction.z) || norm(direction)<=kTiny ||
           !std::isfinite(water.external_min_alignment) ||
           water.external_min_alignment<0.0 || water.external_min_alignment>1.0)
            throw std::invalid_argument("external water requires a finite level, positive density/gravity, and valid outward direction");
    }
    if (!water.enabled) return;
    if (!std::isfinite(config.water.table_z_m) ||
        !std::isfinite(config.water.density_kg_m3) ||
        config.water.density_kg_m3 <= 0.0 ||
        !std::isfinite(config.gravity_m_s2) || config.gravity_m_s2 < 0.0)
        throw std::invalid_argument("hydrostatic water settings must be finite, with positive density and nonnegative gravity");
    if(water.plane_enabled &&
       (!std::isfinite(water.plane_origin_m.x) || !std::isfinite(water.plane_origin_m.y) ||
        !std::isfinite(water.plane_origin_m.z) || !std::isfinite(water.plane_normal.x) ||
        !std::isfinite(water.plane_normal.y) || !std::isfinite(water.plane_normal.z) ||
        water.plane_normal.z<=0.0))
        throw std::invalid_argument("water-table plane requires a finite origin and upward-pointing normal");
    if(water.auto_saturated_density &&
       (water.saturated_densities_kg_m3.size()!=config.region_materials.size()+1 ||
        std::any_of(water.saturated_densities_kg_m3.begin(),water.saturated_densities_kg_m3.end(),
                    [](double density){return !std::isfinite(density)||density<=0.0;})))
        throw std::invalid_argument("saturated density is required for every material");
}

double water_surface_z(const WaterConfig& water, const Vec3& position) {
    if(!water.plane_enabled)return water.table_z_m;
    return water.plane_origin_m.z-
        (water.plane_normal.x*(position.x-water.plane_origin_m.x)+
         water.plane_normal.y*(position.y-water.plane_origin_m.y))/water.plane_normal.z;
}

double node_pore_pressure(const RockModel& model, std::size_t node_index) {
    if(model.hydraulic.valid && model.hydraulic.nodal_pore_pressure_pa.size()==model.nodes.size())
        return model.hydraulic.nodal_pore_pressure_pa.at(node_index);
    const auto& water = model.config.water;
    if (!water.enabled) return 0.0;
    const auto& position=model.nodes[node_index].position;
    return std::max(0.0, water.density_kg_m3 * model.config.gravity_m_s2 *
                             (water_surface_z(water,position)-position.z));
}

double zone_density(const RockModel& model, const Zone& zone, const Vec3& center) {
    const auto& water=model.config.water;
    if(water.auto_saturated_density && model.hydraulic.valid)
        return water.saturated_densities_kg_m3.at(zone.material_index);
    if(water.enabled && water.auto_saturated_density && center.z<=water_surface_z(water,center))
        return water.saturated_densities_kg_m3.at(zone.material_index);
    return zone_material(model,zone).density_kg_m3;
}

std::vector<double> zone_pore_pressures(const RockModel& model) {
    std::vector<double> pressure(model.zones.size(), 0.0);
    if (!model.config.water.enabled && !model.hydraulic.valid) return pressure;
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const auto& nodes = model.zones[zi].nodes;
        for (const auto node : nodes) pressure[zi] += node_pore_pressure(model, node);
        pressure[zi] /= static_cast<double>(nodes.size());
    }
    // FLAC stores pore pressure on gridpoints and exposes its physical-zone
    // average to every constitutive subzone, including both brick overlays.
    for (const auto& cell : model.physical_cells) {
        if (cell.nodes.empty()) continue;
        double average = 0.0;
        for (const auto node : cell.nodes) average += node_pore_pressure(model, node);
        average /= static_cast<double>(cell.nodes.size());
        for (const auto zi : cell.subzones) pressure[zi] = average;
    }
    return pressure;
}

double continuum_biot_coefficient(const RockModel& model) {
    return model.config.seepage.enabled ? model.config.seepage.biot_coefficient : 1.0;
}

void update_zone_state(RockModel& model, Zone& zone,
                       double pore_pressure_pa,
                       const std::array<double, 6>* prescribed_strain = nullptr) {
    const TetGeometry geometry = tetra_geometry(model, zone);
    zone.volume_m3 = geometry.volume;
    const std::array<double, 6> strain =
        prescribed_strain ? *prescribed_strain : total_zone_strain(model, zone);
    const MaterialParameters& material = zone_material(model, zone);
    std::array<double, 6> increment{};
    for (std::size_t component = 0; component < increment.size(); ++component)
        increment[component] = strain[component] - zone.strain[component];
    // Keep total stress for force equilibrium; FLAC's constitutive models see
    // effective stress (tension-positive total stress + pore pressure * I).
    std::array<double, 6> effective_stress = zone.stress_pa;
    for (int i = 0; i < 3; ++i) effective_stress[i] +=
        continuum_biot_coefficient(model)*pore_pressure_pa;
    bool updated = false;
    if (material.type == MaterialType::SMRM) {
        updated = update_smrm_increment(
            material, zone.parent_volume_m3 > 0 ? zone.parent_volume_m3 : zone.volume_m3, increment,
            effective_stress, zone.smrm_history, zone.smrm, zone.tangent_pa);
    } else if (material.type == MaterialType::MohrCoulomb) {
        updated = update_mohr_coulomb_increment(
            material, increment, effective_stress, zone.plastic_state, zone.tangent_pa);
    } else {
        updated = update_elastic_increment(
            material, increment, effective_stress, zone.tangent_pa);
    }
    if (!updated)
        throw std::runtime_error(
            "constitutive update failed or produced non-finite stress");
    for (int i = 0; i < 3; ++i) effective_stress[i] -=
        continuum_biot_coefficient(model)*pore_pressure_pa;
    zone.stress_pa = effective_stress;
    zone.strain = strain;
    if (material.type == MaterialType::SMRM) {
        zone.failure_index = zone.smrm.failure_probability;
        zone.reinforce_requirement = zone.smrm.reinforce_requirement_mpa;
        zone.plastic_state = zone.smrm.yield_state;
        if ((zone.smrm.yield_state & 2U) != 0U) {
            zone.joint_state = JointState::TensionFailure;
        } else if ((zone.smrm.yield_state & 1U) != 0U) {
            zone.joint_state = JointState::ShearFailure;
        } else {
            const bool open = std::any_of(
                zone.smrm.joint_normal_stress_mpa.begin(),
                zone.smrm.joint_normal_stress_mpa.end(),
                [](double value) { return value < 0.0; });
            const bool sliding = std::any_of(
                zone.smrm.joint_h.begin(), zone.smrm.joint_h.end(),
                [](double value) { return value > 1.0e-12; });
            zone.joint_state = open ? JointState::Open
                                     : (sliding ? JointState::Sliding
                                                : JointState::Closed);
        }
        if (zone.material_index != 0 && material.joints.empty()) zone.joint_state = JointState::None;
    } else if (material.type == MaterialType::MohrCoulomb) {
        zone.smrm = {};
        zone.failure_index = zone.plastic_state ? 1.0 : 0.0;
        zone.reinforce_requirement = 0.0;
        if ((zone.plastic_state & 0x0AU) != 0U)
            zone.joint_state = JointState::TensionFailure;
        else if ((zone.plastic_state & 0x05U) != 0U)
            zone.joint_state = JointState::ShearFailure;
        else
            zone.joint_state = JointState::None;
    } else {
        zone.smrm = {};
        zone.failure_index = 0.0;
        zone.reinforce_requirement = 0.0;
        zone.plastic_state = 0;
        zone.joint_state = JointState::None;
    }
}

void enforce_cell_mean_stress(RockModel& model, const PhysicalCell& cell) {
    double mean = 0.0;
    double volume = 0.0;
    for (const std::size_t index : cell.subzones) {
        const Zone& zone = model.zones[index];
        const double v = zone.volume_m3 * zone.integration_weight;
        mean += (zone.stress_pa[0] + zone.stress_pa[1] + zone.stress_pa[2]) * v / 3.0;
        volume += v;
    }
    if (volume <= kTiny) return;
    mean /= volume;
    for (const std::size_t index : cell.subzones) {
        Zone& zone = model.zones[index];
        const double local =
            (zone.stress_pa[0] + zone.stress_pa[1] + zone.stress_pa[2]) / 3.0;
        const double correction = mean - local;
        zone.stress_pa[0] += correction;
        zone.stress_pa[1] += correction;
        zone.stress_pa[2] += correction;
    }
}

void initialize_material_state(RockModel& model) {
    const auto pore_pressure = zone_pore_pressures(model);
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        Zone& zone = model.zones[zi];
        if (!zone.active || zone.smrm_history.step != 0) continue;
        zone.stress_pa = {model.config.initial_stress_x_pa,
                          model.config.initial_stress_y_pa,
                          model.config.initial_stress_z_pa, 0.0, 0.0, 0.0};
        const auto& material = zone_material(model, zone);
        const std::array<double, 6> zero_increment{};
        std::array<double, 6> effective_stress = zone.stress_pa;
        for (int i = 0; i < 3; ++i) effective_stress[i] +=
            continuum_biot_coefficient(model)*pore_pressure[zi];
        bool initialized = false;
        if (material.type == MaterialType::SMRM) {
            initialized = update_smrm_increment(
                material, zone.parent_volume_m3 > 0 ? zone.parent_volume_m3 : zone.volume_m3,
                zero_increment, effective_stress, zone.smrm_history, zone.smrm,
                zone.tangent_pa);
        } else if (material.type == MaterialType::MohrCoulomb) {
            initialized = update_mohr_coulomb_increment(
                material, zero_increment, effective_stress, zone.plastic_state,
                zone.tangent_pa);
        } else {
            initialized = update_elastic_increment(
                material, zero_increment, effective_stress, zone.tangent_pa);
        }
        if (!initialized)
            throw std::runtime_error(
                "failed to initialize the constitutive state");
        for (int i = 0; i < 3; ++i) effective_stress[i] -=
            continuum_biot_coefficient(model)*pore_pressure[zi];
        zone.stress_pa = effective_stress;
    }
}

void add_zone_internal_forces(RockModel& model, const Zone& zone) {
    const TetGeometry geometry = tetra_geometry(model, zone);
    for (std::size_t i = 0; i < 4; ++i) {
        const Vec3& g = geometry.gradients[i];
        const Vec3 internal{
            zone.stress_pa[0] * g.x + zone.stress_pa[3] * g.y +
                zone.stress_pa[4] * g.z,
            zone.stress_pa[3] * g.x + zone.stress_pa[1] * g.y +
                zone.stress_pa[5] * g.z,
            zone.stress_pa[4] * g.x + zone.stress_pa[5] * g.y +
                zone.stress_pa[2] * g.z,
        };
        const Vec3 contribution = internal * (-geometry.volume * zone.integration_weight);
        Node& node = model.nodes[zone.nodes[i]];
        node.force += contribution;
        node.force_reference += norm(contribution);
    }
}

void reset_mass(RockModel& model) {
    for (Node& node : model.nodes) node.mass_kg = 0.0;
    std::vector<double> densities(model.zones.size(),0.0);
    for(const auto& cell:model.physical_cells) {
        Vec3 center{};
        for(const auto node:cell.nodes)center+=model.nodes[node].position;
        center=center*(1.0/static_cast<double>(cell.nodes.size()));
        for(const auto zi:cell.subzones)
            densities[zi]=zone_density(model,model.zones[zi],center);
    }
    for (std::size_t zi=0;zi<model.zones.size();++zi) {
        Zone& zone=model.zones[zi];
        if (!zone.active) continue;
        zone.volume_m3 = tetra_geometry(model, zone).volume;
        if(densities[zi]==0.0) {
            Vec3 center{};
            for(const auto node:zone.nodes)center+=model.nodes[node].position;
            center=center*0.25;
            densities[zi]=zone_density(model,zone,center);
        }
        const double lump = densities[zi] * zone.volume_m3 * zone.integration_weight / 4.0;
        for (const std::size_t node : zone.nodes) model.nodes[node].mass_kg += lump;
    }
}

double quadratic_form(const Matrix6& matrix,
                      const std::array<double, 6>& vector) {
    double value = 0.0;
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            value += vector[row] * matrix[row][column] * vector[column];
    return value;
}

void reset_dynamic_mass(RockModel& model) {
    for (Node& node : model.nodes) node.dynamic_mass = {};
    for (const Zone& zone : model.zones) {
        if (!zone.active) continue;
        const TetGeometry geometry = tetra_geometry(model, zone);
        for (std::size_t local = 0; local < 4; ++local) {
            const Vec3& g = geometry.gradients[local];
            Node& node = model.nodes[zone.nodes[local]];
            if (model.config.solver.flac3d7_static_mass) {
                const auto& material = zone_material(model, zone);
                const double young = material.elastic_modulus_pa;
                const double poisson = material.poisson_ratio;
                const double shear = young / (2.0 * (1.0 + poisson));
                const double bulk = young / (3.0 * (1.0 - 2.0 * poisson));
                const double alpha1 = bulk + 4.0 * shear / 3.0;
                // FLAC3D 7.00 Release 126 gp.mass.inertial, checked at all
                // 3,991 nodes of the 601 slope, is 6.75 times the documented
                // alpha1 * V * max(grad N)^2 contribution. Keep this observed
                // static solver convention opt-in; it is not a material factor.
                const double scalar = 6.75 * alpha1 * geometry.volume *
                    zone.integration_weight *
                    std::max({g.x * g.x, g.y * g.y, g.z * g.z});
                node.dynamic_mass.x += scalar;
                node.dynamic_mass.y += scalar;
                node.dynamic_mass.z += scalar;
            } else {
                const std::array<double, 6> bx{g.x, 0.0, 0.0, g.y, g.z, 0.0};
                const std::array<double, 6> by{0.0, g.y, 0.0, g.x, 0.0, g.z};
                const std::array<double, 6> bz{0.0, 0.0, g.z, 0.0, g.x, g.y};
                node.dynamic_mass.x +=
                    std::abs(quadratic_form(zone.tangent_pa, bx)) * geometry.volume * zone.integration_weight;
                node.dynamic_mass.y +=
                    std::abs(quadratic_form(zone.tangent_pa, by)) * geometry.volume * zone.integration_weight;
                node.dynamic_mass.z +=
                    std::abs(quadratic_form(zone.tangent_pa, bz)) * geometry.volume * zone.integration_weight;
            }
        }
    }
    double maximum_modulus = model.config.material.elastic_modulus_pa;
    for (const auto& region : model.config.region_materials)
        maximum_modulus = std::max(maximum_modulus, region.material.elastic_modulus_pa);
    const double floor = std::max(maximum_modulus, 1.0) *
                         1.0e-12;
    for (Node& node : model.nodes) {
        node.dynamic_mass.x = std::max(node.dynamic_mass.x, floor);
        node.dynamic_mass.y = std::max(node.dynamic_mass.y, floor);
        node.dynamic_mass.z = std::max(node.dynamic_mass.z, floor);
    }
}

using FaceKey = std::array<std::size_t, 3>;

FaceKey face_key(std::size_t a, std::size_t b, std::size_t c) {
    FaceKey key{a, b, c};
    std::sort(key.begin(), key.end());
    return key;
}

const std::array<std::array<int, 3>, 4> kTetFaces{{
    {{0, 2, 1}}, {{0, 1, 3}}, {{0, 3, 2}}, {{1, 2, 3}},
}};

Vec3 zone_centroid(const RockModel& model, const Zone& zone) {
    Vec3 center{};
    for (const std::size_t node : zone.nodes) center += model.nodes[node].position;
    return center * 0.25;
}

void add_far_field_forces(RockModel& model) {
    const std::array<double, 3> normal_stress{
        model.config.initial_stress_x_pa,
        model.config.initial_stress_y_pa,
        model.config.initial_stress_z_pa,
    };
    if (std::abs(normal_stress[0]) + std::abs(normal_stress[1]) +
            std::abs(normal_stress[2]) <
        1.0) {
        return;
    }
    Vec3 low = model.nodes.front().position;
    Vec3 high = low;
    for (const Node& node : model.nodes) {
        low.x = std::min(low.x, node.position.x);
        low.y = std::min(low.y, node.position.y);
        low.z = std::min(low.z, node.position.z);
        high.x = std::max(high.x, node.position.x);
        high.y = std::max(high.y, node.position.y);
        high.z = std::max(high.z, node.position.z);
    }
    const double scale = std::max({high.x - low.x, high.y - low.y, high.z - low.z, 1.0});
    const double tolerance = scale * 1.0e-8;
    struct BoundaryTriangle {std::array<std::size_t,3> face; Vec3 center; double weight;};
    std::vector<BoundaryTriangle> boundary;
    if(model.physical_cells.empty()) {
        std::map<FaceKey,std::vector<std::pair<int,std::array<std::size_t,3>>>> faces;
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& z=model.zones[zi];if(!z.active)continue;
            for(const auto& local:kTetFaces){std::array<std::size_t,3> f{z.nodes[local[0]],z.nodes[local[1]],z.nodes[local[2]]};faces[face_key(f[0],f[1],f[2])].push_back({static_cast<int>(zi),f});}
        }
        for(const auto& [key,owners]:faces)if(owners.size()==1)boundary.push_back({owners[0].second,zone_centroid(model,model.zones[owners[0].first]),1.0});
    } else {
        const std::map<std::size_t,std::vector<std::vector<int>>> shapes{
            {4,{{0,2,1},{0,1,3},{1,2,3},{2,0,3}}},
            {5,{{0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}}},
            {6,{{0,2,1},{3,4,5},{0,1,4,3},{1,2,5,4},{2,0,3,5}}},
            {8,{{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}}}
        };
        struct Polygon {std::vector<std::size_t> nodes;Vec3 center;int count;};
        std::map<std::vector<std::size_t>,Polygon> polygons;
        for(const auto& cell:model.physical_cells){
            Vec3 center{};for(auto i:cell.nodes)center+=model.nodes[i].position;center=center*(1.0/cell.nodes.size());
            for(const auto& local:shapes.at(cell.nodes.size())){
                std::vector<std::size_t> face;for(auto i:local)face.push_back(cell.nodes[i]);
                auto key=face;std::sort(key.begin(),key.end());
                if(polygons.count(key))polygons[key].count++;else polygons[key]={face,center,1};
            }
        }
        for(const auto& [key,p]:polygons)if(p.count==1){
            const auto& f=p.nodes;
            if(f.size()==3)boundary.push_back({{f[0],f[1],f[2]},p.center,1.0});
            else for(const auto& t:std::array<std::array<int,3>,4>{{{0,1,2},{0,2,3},{0,1,3},{1,2,3}}})boundary.push_back({{f[t[0]],f[t[1]],f[t[2]]},p.center,.5});
        }
    }
    for(const auto& item:boundary){
        const auto& face=item.face;
        const Vec3& a = model.nodes[face[0]].position;
        const Vec3& b = model.nodes[face[1]].position;
        const Vec3& c = model.nodes[face[2]].position;
        Vec3 normal = cross(b - a, c - a);
        const double twice_area = norm(normal);
        if (twice_area <= kTiny) continue;
        normal = normal * (1.0 / twice_area);
        const Vec3 face_center = (a + b + c) * (1.0 / 3.0);
        if (dot(normal, item.center - face_center) >
            0.0) {
            normal = normal * -1.0;
        }
        const bool on_outer_box =
            (std::abs(a.x - low.x) < tolerance && std::abs(b.x - low.x) < tolerance &&
             std::abs(c.x - low.x) < tolerance) ||
            (std::abs(a.x - high.x) < tolerance &&
             std::abs(b.x - high.x) < tolerance &&
             std::abs(c.x - high.x) < tolerance) ||
            (std::abs(a.y - low.y) < tolerance && std::abs(b.y - low.y) < tolerance &&
             std::abs(c.y - low.y) < tolerance) ||
            (std::abs(a.y - high.y) < tolerance &&
             std::abs(b.y - high.y) < tolerance &&
             std::abs(c.y - high.y) < tolerance) ||
            (std::abs(a.z - low.z) < tolerance && std::abs(b.z - low.z) < tolerance &&
             std::abs(c.z - low.z) < tolerance) ||
            (std::abs(a.z - high.z) < tolerance &&
             std::abs(b.z - high.z) < tolerance &&
             std::abs(c.z - high.z) < tolerance);
        if (!on_outer_box) continue;
        const Vec3 traction{normal_stress[0] * normal.x,
                            normal_stress[1] * normal.y,
                            normal_stress[2] * normal.z};
        const Vec3 nodal_force = traction * (0.5 * twice_area * item.weight / 3.0);
        for (const std::size_t node : face) {
            model.nodes[node].force += nodal_force;
            model.nodes[node].force_reference += norm(nodal_force);
        }
    }
}

struct ExternalWaterFaceForce {
    std::vector<std::size_t> nodes;
    Vec3 force_per_node;
};

std::vector<ExternalWaterFaceForce> external_water_faces(const RockModel& model) {
    const auto& water=model.config.water;
    if(!water.external_enabled)return {};
    struct Face {std::vector<std::size_t> nodes;Vec3 cell_center;int owners{1};};
    std::map<std::vector<std::array<double,3>>,Face> faces;
    const auto add_face=[&](std::vector<std::size_t> nodes,const Vec3& center){
        std::vector<std::array<double,3>> key;
        key.reserve(nodes.size());
        for(const auto node:nodes){const auto& p=model.nodes[node].position;key.push_back({p.x,p.y,p.z});}
        std::sort(key.begin(),key.end());
        auto [entry,inserted]=faces.try_emplace(key,Face{std::move(nodes),center,1});
        if(!inserted)++entry->second.owners;
    };
    if(model.physical_cells.empty()) {
        for(const auto& zone:model.zones)if(zone.active) {
            const Vec3 center=zone_centroid(model,zone);
            for(const auto& local:kTetFaces)
                add_face({zone.nodes[local[0]],zone.nodes[local[1]],zone.nodes[local[2]]},center);
        }
    } else {
        const std::map<std::size_t,std::vector<std::vector<int>>> shapes{
            {4,{{0,2,1},{0,1,3},{1,2,3},{2,0,3}}},
            {5,{{0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}}},
            {6,{{0,2,1},{3,4,5},{0,1,4,3},{1,2,5,4},{2,0,3,5}}},
            {8,{{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}}},
        };
        for(const auto& cell:model.physical_cells) {
            if(std::none_of(cell.subzones.begin(),cell.subzones.end(),
                            [&](std::size_t zi){return model.zones[zi].active;}))continue;
            Vec3 center{};
            for(const auto node:cell.nodes)center+=model.nodes[node].position;
            center=center*(1.0/static_cast<double>(cell.nodes.size()));
            for(const auto& local:shapes.at(cell.nodes.size())) {
                std::vector<std::size_t> nodes;
                for(const auto i:local)nodes.push_back(cell.nodes[i]);
                add_face(std::move(nodes),center);
            }
        }
    }
    const Vec3 direction=water.external_outward_direction*
        (1.0/norm(water.external_outward_direction));
    std::vector<ExternalWaterFaceForce> selected;
    bool matched_direction=false;
    for(const auto& [key,face]:faces)if(face.owners==1) {
        Vec3 face_center{};
        for(const auto node:face.nodes)face_center+=model.nodes[node].position;
        face_center=face_center*(1.0/static_cast<double>(face.nodes.size()));
        const Vec3 p0=model.nodes[face.nodes.front()].position;
        Vec3 area_normal{};
        for(std::size_t i=1;i+1<face.nodes.size();++i)
            area_normal+=cross(model.nodes[face.nodes[i]].position-p0,
                               model.nodes[face.nodes[i+1]].position-p0);
        const double twice_area=norm(area_normal);
        if(twice_area<=kTiny)continue;
        Vec3 normal=area_normal*(1.0/twice_area);
        if(dot(normal,face.cell_center-face_center)>0.0)normal=normal*(-1.0);
        if(dot(normal,direction)<water.external_min_alignment)continue;
        matched_direction=true;
        const double depth=water.external_surface_z_m-face_center.z;
        if(depth<=0.0)continue;
        const double pressure=water.external_density_kg_m3*model.config.gravity_m_s2*depth;
        selected.push_back({face.nodes,normal*(-pressure*0.5*twice_area/face.nodes.size())});
    }
    if(!matched_direction)throw std::invalid_argument("external water direction selects no exposed mesh faces");
    return selected;
}

void add_external_water_forces(RockModel& model,
                               const std::vector<ExternalWaterFaceForce>& faces) {
    for(const auto& face:faces)for(const auto node:face.nodes) {
        model.nodes[node].force+=face.force_per_node;
        model.nodes[node].force_reference+=norm(face.force_per_node);
    }
}

void install_support(RockModel& model, double stiffness) {
    std::map<FaceKey, std::vector<std::size_t>> owners;
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const Zone& zone = model.zones[zi];
        for (const auto& local : kTetFaces) {
            owners[face_key(zone.nodes[local[0]], zone.nodes[local[1]],
                            zone.nodes[local[2]])]
                .push_back(zi);
        }
    }
    std::set<std::size_t> support_nodes;
    for (const auto& [face, zones] : owners) {
        bool has_active = false;
        bool has_inactive = false;
        for (const std::size_t zone : zones) {
            has_active = has_active || model.zones[zone].active;
            has_inactive = has_inactive || !model.zones[zone].active;
        }
        if (has_active && has_inactive) support_nodes.insert(face.begin(), face.end());
    }
    for (const std::size_t index : support_nodes) {
        Node& node = model.nodes[index];
        node.supported = true;
        node.support_reference = node.displacement;
    }
    model.config.support_stiffness_n_m = stiffness;
}

void deactivate_tunnel_segment(RockModel& model, double start, double end) {
    const double cx = 0.0;
    const double cz = model.config.height_m * 0.5;
    const double rx = std::max(model.config.tunnel_radius_x_m, 0.1);
    const double rz = std::max(model.config.tunnel_radius_z_m, 0.1);
    for (Zone& zone : model.zones) {
        if (!zone.active) continue;
        const Vec3 center = zone_centroid(model, zone);
        const double radial = ((center.x - cx) * (center.x - cx)) / (rx * rx) +
                              ((center.z - cz) * (center.z - cz)) / (rz * rz);
        if (center.y >= start - 1.0e-9 && center.y <= end + 1.0e-9 &&
            radial <= 1.0) {
            zone.active = false;
            zone.stress_pa.fill(0.0);
            zone.strain.fill(0.0);
            zone.tangent_pa = {};
            zone.smrm = {};
            zone.smrm_history = {};
        }
    }
}

#include "interface_plane.inc"
#include "cable_support.inc"

StageResult equilibrate(RockModel& model, const std::string& name,
                        int stage_index, int stage_count,
                        const ProgressCallback& progress) {
    const auto pore_pressure = zone_pore_pressures(model);
    reset_mass(model);
    const auto water_faces=external_water_faces(model);
    // A FLAC static solve uses local time stepping rather than the smallest
    // physical wave time step for every grid point.  Diagonal tangent
    // stiffness supplies the corresponding nodal inertial scaling here.
    reset_dynamic_mass(model);
    interface_dynamic_mass(model);
    cable_dynamic_mass(model);
    std::vector<std::size_t> platen_nodes;
    double platen_dynamic_mass=0.0;
    if(model.config.solver.rigid_top_platen_z){
        double z_min=model.nodes.front().position.z,z_max=z_min;
        for(const auto& node:model.nodes){
            z_min=std::min(z_min,node.position.z);
            z_max=std::max(z_max,node.position.z);
        }
        const double tolerance=std::max(1.0,z_max-z_min)*1.0e-8;
        for(std::size_t i=0;i<model.nodes.size();++i){
            const auto& node=model.nodes[i];
            if(std::abs(node.position.z-z_max)>tolerance)continue;
            if(node.fixed[2]||node.mass_kg<=0.0)
                throw std::invalid_argument("rigid top platen requires free vertical top nodes");
            platen_nodes.push_back(i);
            platen_dynamic_mass+=node.dynamic_mass.z;
        }
        if(platen_nodes.size()<3||!(platen_dynamic_mass>0.0))
            throw std::invalid_argument("rigid top platen requires a horizontal top face");
        const double initial_z=model.nodes[platen_nodes.front()].displacement.z;
        const double initial_v=model.nodes[platen_nodes.front()].velocity.z;
        for(const auto i:platen_nodes)
            if(std::abs(model.nodes[i].displacement.z-initial_z)>1.0e-9||
               std::abs(model.nodes[i].velocity.z-initial_v)>1.0e-9)
                throw std::invalid_argument("rigid top platen requires initially tied vertical state");
    }
    const double pseudo_dt = std::clamp(model.config.solver.time_step_scale,
        1.0e-4, model.config.solver.flac3d7_static_mass ? 1.0 : 0.95);
    const double local_damping =
        std::clamp(model.config.solver.damping, 0.0, 0.9999);
    StageResult result;
    result.name = name;
    double ratio = std::numeric_limits<double>::infinity();
    const int progress_interval =
        std::max(1, model.config.solver.max_cycles / 100);
    // Evaluate forces at the current displacement, including after the last advance.
    for (int cycle = 0; cycle <= model.config.solver.max_cycles; ++cycle) {
        double reference_force = 0.0;
        for (Node& node : model.nodes) {
            node.force = {0.0, 0.0, -node.mass_kg * model.config.gravity_m_s2};
            node.force_reference =
                std::abs(node.mass_kg * model.config.gravity_m_s2);
        }
        add_far_field_forces(model);
        add_external_water_forces(model,water_faces);
        std::vector<std::array<double, 6>> mixed_strain;
        if (model.config.solver.nodal_mixed_discretization) {
            mixed_strain = nodal_mixed_strains(model);
        }
        if (!model.physical_cells.empty()) {
            if (mixed_strain.empty()) {
                mixed_strain.resize(model.zones.size());
                for (std::size_t i=0;i<model.zones.size();++i)
                    mixed_strain[i] = total_zone_strain(model, model.zones[i]);
            }
            for (const auto& cell : model.physical_cells) {
                if (cell.nodes.size()==4) continue;
                double mean=0.0, volume=0.0;
                for (auto i : cell.subzones) {
                    mixed_strain[i]=total_zone_strain(model,model.zones[i]);
                    const double v=model.zones[i].volume_m3*model.zones[i].integration_weight;
                    mean+=(mixed_strain[i][0]+mixed_strain[i][1]+mixed_strain[i][2])*v/3.0;
                    volume+=v;
                }
                mean/=volume;
                // FLAC keeps each overlay tet's deviatoric strain. Only the
                // central brick pair shares shear; copying all six components
                // between overlay tets loses real nonuniform deformation.
                if (cell.nodes.size()==8 && cell.subzones.size()==10) {
                    const auto a=cell.subzones[2], b=cell.subzones[7];
                    for (int k=3;k<6;++k) {
                        const double shear=0.5*(mixed_strain[a][k]+mixed_strain[b][k]);
                        mixed_strain[a][k]=shear;
                        mixed_strain[b][k]=shear;
                    }
                }
                for (auto i : cell.subzones) {
                    const double correction=mean-(mixed_strain[i][0]+mixed_strain[i][1]+mixed_strain[i][2])/3.0;
                    for (int k=0;k<3;++k) mixed_strain[i][k]+=correction;
                }
            }
        }
        for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
            Zone& zone = model.zones[zi];
            if (!zone.active) continue;
            update_zone_state(
                model, zone, pore_pressure[zi],
                mixed_strain.empty() ? nullptr : &mixed_strain[zi]);
        }
        for (const auto& cell : model.physical_cells)
            if (cell.nodes.size()!=4) enforce_cell_mean_stress(model, cell);
        for (const Zone& zone : model.zones) {
            if (!zone.active) continue;
            add_zone_internal_forces(model, zone);
        }
        if (model.config.support_stiffness_n_m > 0.0) {
            for (Node& node : model.nodes) {
                if (!node.supported) continue;
                const Vec3 delta = node.displacement - node.support_reference;
                node.force += delta * (-model.config.support_stiffness_n_m);
                node.force_reference +=
                    norm(delta) * model.config.support_stiffness_n_m;
            }
        }
        add_interface_forces(model);
        update_cable_forces(model);
        double residual = 0.0;
        double local_ratio = 0.0;
        double maximum_component_ratio = 0.0;
        double unbalanced_maximum = 0.0;
        for (Node& node : model.nodes) {
            if (node.mass_kg <= 0.0) continue;
            if (node.fixed[0]) node.force.x = 0.0;
            if (node.fixed[1]) node.force.y = 0.0;
            if (node.fixed[2]) node.force.z = 0.0;
        }
        if(!platen_nodes.empty()){
            double total_force=0.0,total_reference=0.0;
            for(const auto i:platen_nodes){
                total_force+=model.nodes[i].force.z;
                total_reference+=model.nodes[i].force_reference;
            }
            for(const auto i:platen_nodes){
                const double share=model.nodes[i].dynamic_mass.z/platen_dynamic_mass;
                model.nodes[i].force.z=total_force*share;
                model.nodes[i].force_reference=total_reference*share;
            }
        }
        for (Node& node : model.nodes) {
            if (node.mass_kg <= 0.0) continue;
            const double node_residual = norm(node.force);
            residual += node_residual;
            reference_force += node.force_reference;
            unbalanced_maximum = std::max(unbalanced_maximum, node_residual);
            local_ratio = std::max(
                local_ratio,
                node_residual / std::max(node.force_reference, 1.0));
            maximum_component_ratio = std::max(
                maximum_component_ratio,
                std::max({std::abs(node.force.x), std::abs(node.force.y),
                          std::abs(node.force.z)}) /
                    std::max(node.force_reference, 1.0));

        }
        const auto cable_balance = cable_force_balance(model);
        residual += cable_balance.residual;
        reference_force += cable_balance.reference;
        local_ratio = std::max(local_ratio, cable_balance.local);
        maximum_component_ratio = std::max(maximum_component_ratio, cable_balance.local);
        unbalanced_maximum = std::max(unbalanced_maximum, cable_balance.maximum_force);
        ratio = residual / std::max(reference_force, 1.0);
        result.force_ratio_local = local_ratio;
        result.force_ratio_maximum = maximum_component_ratio;
        result.unbalanced_force_maximum_n = unbalanced_maximum;
        result.cycles = cycle;
        const bool converged =
            cycle > 25 && ratio <= model.config.solver.convergence_tolerance &&
            cable_balance.average() <= model.config.solver.convergence_tolerance &&
            cable_balance.scaled_local <= model.config.solver.convergence_tolerance;
        if (progress &&
            (cycle == 1 || cycle % progress_interval == 0 || converged ||
             cycle == model.config.solver.max_cycles)) {
            progress(SolveProgress{
                .stage_name = name,
                .stage_index = stage_index,
                .stage_count = stage_count,
                .cycle = cycle,
                .max_cycles = model.config.solver.max_cycles,
                .force_ratio = ratio,
                .force_ratio_local = local_ratio,
                .force_ratio_maximum = maximum_component_ratio,
                .unbalanced_force_maximum_n = unbalanced_maximum,
                .stage_converged = converged,
            });
        }
        if (converged) {
            result.converged = true;
            break;
        }
        if (cycle == model.config.solver.max_cycles) break;
        for (Node& node : model.nodes) {
            if (node.mass_kg <= 0.0) continue;
            const auto advance_component = [&](double force, double mass,
                                               double& velocity,
                                               double& displacement) {
                double damped_force = force;
                if (std::abs(velocity) > kTiny) {
                    damped_force -= local_damping * std::abs(force) *
                                    std::copysign(1.0, velocity);
                    // FLAC local damping is inactive on a velocity reversal.
                    // Detect the crossing with the tentative damped step, so
                    // the force itself advances that component through zero.
                    const double tentative =
                        velocity + damped_force / mass * pseudo_dt;
                    if (velocity * tentative < 0.0) damped_force = force;
                }
                velocity += damped_force / mass * pseudo_dt;
                displacement += velocity * pseudo_dt;
            };
            advance_component(node.force.x, node.dynamic_mass.x,
                              node.velocity.x, node.displacement.x);
            advance_component(node.force.y, node.dynamic_mass.y,
                              node.velocity.y, node.displacement.y);
            advance_component(node.force.z, node.dynamic_mass.z,
                              node.velocity.z, node.displacement.z);
            if (node.fixed[0]) node.velocity.x = 0.0;
            if (node.fixed[1]) node.velocity.y = 0.0;
            if (node.fixed[2]) node.velocity.z = 0.0;
        }
        advance_cables(model, pseudo_dt, local_damping);
    }
    result.force_ratio = ratio;
    for (const Node& node : model.nodes) {
        result.max_displacement_m = std::max(result.max_displacement_m,
                                             norm(node.displacement));
    }
    result.active_zones = static_cast<std::size_t>(
        std::count_if(model.zones.begin(), model.zones.end(),
                      [](const Zone& zone) { return zone.active; }));
    return result;
}

std::size_t node_index(int ix, int iy, int iz, int nx, int ny) {
    return static_cast<std::size_t>((iz * (ny + 1) + iy) * (nx + 1) + ix);
}

void append_cell_tetrahedra(RockModel& model, int ix, int iy, int iz, int nx,
                            int ny, int& next_id) {
    const std::size_t n000 = node_index(ix, iy, iz, nx, ny);
    const std::size_t n100 = node_index(ix + 1, iy, iz, nx, ny);
    const std::size_t n010 = node_index(ix, iy + 1, iz, nx, ny);
    const std::size_t n110 = node_index(ix + 1, iy + 1, iz, nx, ny);
    const std::size_t n001 = node_index(ix, iy, iz + 1, nx, ny);
    const std::size_t n101 = node_index(ix + 1, iy, iz + 1, nx, ny);
    const std::size_t n011 = node_index(ix, iy + 1, iz + 1, nx, ny);
    const std::size_t n111 = node_index(ix + 1, iy + 1, iz + 1, nx, ny);
    for (const auto& nodes : std::array<std::array<std::size_t, 4>, 6>{{
             {{n000, n100, n110, n111}},
             {{n000, n110, n010, n111}},
             {{n000, n010, n011, n111}},
             {{n000, n011, n001, n111}},
             {{n000, n001, n101, n111}},
             {{n000, n101, n100, n111}},
         }}) {
        Zone zone;
        zone.id = next_id++;
        zone.nodes.assign(nodes.begin(), nodes.end());
        model.zones.push_back(std::move(zone));
    }
}

std::string json_escape(const std::string& value) {
    std::ostringstream stream;
    for (const char ch : value) {
        switch (ch) {
            case '\\': stream << "\\\\"; break;
            case '"': stream << "\\\""; break;
            case '\n': stream << "\\n"; break;
            case '\r': stream << "\\r"; break;
            case '\t': stream << "\\t"; break;
            default: stream << ch;
        }
    }
    return stream.str();
}

template <typename Values>
void write_numeric_array(std::ostringstream& out, const Values& values,
                         double scale = 1.0) {
    out << '[';
    bool first = true;
    for (const auto value : values) {
        if (!first) out << ',';
        first = false;
        const double number = static_cast<double>(value) * scale;
        out << (std::isfinite(number) ? number : 0.0);
    }
    out << ']';
}

std::string trim_copy(std::string value) {
    const auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), not_space));
    value.erase(std::find_if(value.rbegin(), value.rend(), not_space).base(),
                value.end());
    return value;
}

std::string upper_copy(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(),
                   [](unsigned char ch) { return static_cast<char>(std::toupper(ch)); });
    return value;
}

std::vector<std::string> csv_fields(const std::string& line) {
    std::vector<std::string> fields;
    std::istringstream stream(line);
    std::string field;
    while (std::getline(stream, field, ',')) fields.push_back(trim_copy(field));
    return fields;
}

void append_model_stages(RockModel& model) {
    const ModelConfig& config = model.config;
    model.stages.push_back(Stage{.name = "Initial equilibrium",
                                 .type = StageType::Initial});
    if (std::any_of(config.cables.begin(),config.cables.end(),[](const auto& c){return c.after_initial;}))
        model.stages.push_back(Stage{.name="Cable installation",.type=StageType::CableInstallation});
    if(std::any_of(config.cables.begin(),config.cables.end(),[](const auto& c){return c.pretension>0;}))
        model.stages.push_back(Stage{.name="Prestress released equilibrium",.type=StageType::CableRelease});
    if (config.scenario == Scenario::Tunnel) {
        const int steps = std::clamp(config.excavation_steps, 0, 20);
        for (int step = 0; step < steps; ++step) {
            const double start = config.length_m * static_cast<double>(step) / steps;
            const double end = config.length_m * static_cast<double>(step + 1) / steps;
            model.stages.push_back(Stage{
                .name = "Excavation " + std::to_string(step + 1),
                .type = StageType::Excavation,
                .station_start_m = start,
                .station_end_m = end,
            });
            if (config.support_stiffness_n_m > 0.0) {
                model.stages.push_back(Stage{
                    .name = "Support " + std::to_string(step + 1),
                    .type = StageType::Support,
                    .station_start_m = start,
                    .station_end_m = end,
                    .support_stiffness_n_m = config.support_stiffness_n_m,
                });
            }
        }
    } else if (config.material.strength_reduction_factor > 1.0) {
        model.stages.push_back(Stage{
            .name = "Strength reduction",
            .type = StageType::StrengthReduction,
            .strength_reduction_factor = config.material.strength_reduction_factor,
        });
    }
}

}  // namespace

const char* scenario_name(Scenario scenario) {
    return scenario == Scenario::Slope ? "slope" : "tunnel";
}

const char* material_name(MaterialType material) {
    if (material == MaterialType::SMRM) return "smrm";
    if (material == MaterialType::MohrCoulomb) return "mohr-coulomb";
    return "elastic";
}

const char* joint_state_name(JointState state) {
    switch (state) {
        case JointState::Closed: return "closed";
        case JointState::Open: return "open";
        case JointState::Sliding: return "sliding";
        case JointState::TensionFailure: return "tension_failure";
        case JointState::ShearFailure: return "shear_failure";
        default: return "none";
    }
}

RockModel generate_model(const ModelConfig& input) {
    ModelConfig config = input;
    validate_water(config);
    config.nx = std::clamp(config.nx, 1, 80);
    config.ny = std::clamp(config.ny, 1, 120);
    config.nz = std::clamp(config.nz, 1, 80);
    if (config.length_m <= 0.0 || config.width_m <= 0.0 || config.height_m <= 0.0) {
        throw std::invalid_argument("model dimensions must be positive");
    }
    if (config.material.elastic_modulus_pa <= 0.0 ||
        config.material.density_kg_m3 <= 0.0) {
        throw std::invalid_argument("material density and elastic modulus must be positive");
    }
    RockModel model;
    model.name = config.scenario == Scenario::Slope ? "Parametric slope"
                                                    : "Staged tunnel";
    model.scenario = config.scenario;
    model.config = config;
    const int nx = config.nx;
    const int ny = config.ny;
    const int nz = config.nz;
    model.nodes.reserve(static_cast<std::size_t>(nx + 1) * (ny + 1) * (nz + 1));
    int next_node_id = 1;
    const double slope_drop = std::min(
        config.height_m * 0.72,
        config.length_m * std::tan(config.slope_angle_deg * std::numbers::pi / 180.0));
    for (int iz = 0; iz <= nz; ++iz) {
        for (int iy = 0; iy <= ny; ++iy) {
            for (int ix = 0; ix <= nx; ++ix) {
                const double fx = static_cast<double>(ix) / nx;
                const double fy = static_cast<double>(iy) / ny;
                const double fz = static_cast<double>(iz) / nz;
                Vec3 position;
                if (config.scenario == Scenario::Slope) {
                    position.x = -config.length_m * 0.5 + fx * config.length_m;
                    position.y = -config.width_m * 0.5 + fy * config.width_m;
                    position.z = fz * (config.height_m - slope_drop * fx);
                } else {
                    position.x = -config.width_m * 0.5 + fx * config.width_m;
                    position.y = fy * config.length_m;
                    position.z = fz * config.height_m;
                }
                Node node;
                node.id = next_node_id++;
                node.position = position;
                node.fixed = {
                    iz == 0 || ix == 0 || ix == nx,
                    iz == 0 || iy == 0 || iy == ny,
                    iz == 0,
                };
                model.nodes.push_back(node);
            }
        }
    }
    model.zones.reserve(static_cast<std::size_t>(nx) * ny * nz * 6);
    int next_zone_id = 1;
    for (int iz = 0; iz < nz; ++iz) {
        for (int iy = 0; iy < ny; ++iy) {
            for (int ix = 0; ix < nx; ++ix) {
                append_cell_tetrahedra(model, ix, iy, iz, nx, ny, next_zone_id);
            }
        }
    }
    for (Zone& zone : model.zones) zone.volume_m3 = tetra_geometry(model, zone).volume;
    append_model_stages(model);
    return model;
}

RockModel import_abaqus_c3d4_stream(std::istream& input_stream,
                                    const ModelConfig& input) {
    validate_water(input);
    if (input.material.elastic_modulus_pa <= 0.0 ||
        input.material.density_kg_m3 <= 0.0) {
        throw std::invalid_argument(
            "material density and elastic modulus must be positive");
    }
    enum class Section { None, Nodes, Tet4 };
    std::size_t node_count=4;
    Section section = Section::None;
    std::vector<std::vector<int>> pending_elements;
    std::vector<std::size_t> pending_materials;
    std::vector<std::string> pending_groups;
    std::string current_group;
    std::size_t current_material = 0;
    std::vector<bool> assigned_materials(input.region_materials.size(), false);
    std::unordered_map<int, std::size_t> node_indices;
    RockModel model;
    model.name = "Imported Abaqus C3D4 model";
    model.scenario = input.scenario;
    model.config = input;

    std::string line;
    std::size_t line_number = 0;
    while (std::getline(input_stream, line)) {
        ++line_number;
        line = trim_copy(line);
        if (line.empty() || line.rfind("**", 0) == 0) continue;
        if (line.front() == '*') {
            const std::string keyword = upper_copy(line);
            if (keyword.rfind("*NODE", 0) == 0) {
                section = Section::Nodes;
            } else if (keyword.rfind("*ELEMENT", 0) == 0) {
                section = Section::Tet4;
                node_count=0;
                for (const auto& field : csv_fields(keyword)) {
                    auto token=trim_copy(field);
                    token.erase(std::remove_if(token.begin(),token.end(),[](unsigned char c){return std::isspace(c);}),token.end());
                    if(token=="TYPE=C3D4")node_count=4;
                    if(token=="TYPE=C3D5")node_count=5;
                    if(token=="TYPE=C3D6")node_count=6;
                    if(token=="TYPE=C3D8")node_count=8;
                }
                if(!node_count)throw std::runtime_error("unsupported Abaqus element type (requires C3D4/5/6/8)");
                current_material = 0;
                current_group.clear();
                for (const auto& field : csv_fields(keyword)) {
                    const auto token = trim_copy(field);
                    if (token.rfind("ELSET=", 0) != 0) continue;
                    const auto name = trim_copy(token.substr(6));
                    current_group=name;
                    for (std::size_t i = 0; i < input.region_materials.size(); ++i) {
                        if (upper_copy(input.region_materials[i].name) == name) {
                            current_material = i + 1;
                        }
                    }
                }
            } else {
                section = Section::None;
            }
            continue;
        }
        const auto fields = csv_fields(line);
        try {
            if (section == Section::Nodes) {
                if (fields.size() < 4) {
                    throw std::runtime_error("node needs id,x,y,z");
                }
                Node node;
                node.id = std::stoi(fields[0]);
                node.position = {std::stod(fields[1]), std::stod(fields[2]),
                                 std::stod(fields[3])};
                if (!node_indices.emplace(node.id, model.nodes.size()).second) {
                    throw std::runtime_error("duplicate node id");
                }
                model.nodes.push_back(node);
            } else if (section == Section::Tet4) {
                if (fields.size()!=node_count+1)throw std::runtime_error("element node count does not match type");
                pending_materials.push_back(current_material);
                pending_groups.push_back(current_group);
                if(current_material)assigned_materials[current_material-1]=true;
                std::vector<int> row;
                for(const auto& field:fields)row.push_back(std::stoi(field));
                pending_elements.push_back(row);
            }
        } catch (const std::exception& error) {
            throw std::runtime_error("invalid Abaqus INP line " +
                                     std::to_string(line_number) + ": " +
                                     error.what());
        }
    }
    if (model.nodes.empty()) throw std::runtime_error("Abaqus INP has no nodes");
    if (pending_elements.empty()) {
        throw std::runtime_error("Abaqus INP has no C3D4 elements");
    }

    for (std::size_t i = 0; i < assigned_materials.size(); ++i) {
        if (!assigned_materials[i]) throw std::runtime_error("material group has no elements: " + input.region_materials[i].name);
    }
    const bool mixed=std::any_of(pending_elements.begin(),pending_elements.end(),[](const auto& e){return e.size()!=5;});
    if(mixed && input.scenario!=Scenario::Slope)throw std::runtime_error("mixed-element import currently supports slope calculations only");
    std::set<int> element_ids;
    const std::map<std::size_t,std::vector<std::array<int,4>>> templates{
        {4,{{0,1,2,3}}},
        {5,{{0,1,2,4},{0,2,3,4},{0,1,3,4},{1,2,3,4}}},
        {6,{{0,1,2,3},{1,2,3,4},{2,3,4,5},{1,2,0,5},{0,1,5,4},{5,0,4,3}}},
        {8,{{0,1,3,4},{1,2,3,6},{1,3,4,6},{1,4,5,6},{3,4,6,7},
            {0,1,2,5},{0,2,3,7},{2,0,5,7},{0,4,5,7},{2,5,6,7}}}
    };
    for(std::size_t ei=0;ei<pending_elements.size();++ei){
        const auto& element=pending_elements[ei];
        if(!element_ids.insert(element[0]).second)throw std::runtime_error("duplicate element id");
        PhysicalCell cell;cell.id=element[0];
        for(std::size_t k=1;k<element.size();++k){
            const auto found=node_indices.find(element[k]);
            if(found==node_indices.end())throw std::runtime_error("element references missing node");
            cell.nodes.push_back(found->second);
        }
        if(std::set<std::size_t>(cell.nodes.begin(),cell.nodes.end()).size()!=cell.nodes.size())throw std::runtime_error("element has duplicate corners");
        double parent_volume=0;
        for(const auto& tet:templates.at(cell.nodes.size())){
            Zone zone;zone.id=cell.id;zone.material_index=pending_materials[ei];zone.material_group=pending_groups[ei];
            for(int k=0;k<4;++k)zone.nodes[k]=cell.nodes[tet[k]];
            const auto a=model.nodes[zone.nodes[0]].position;
            const auto b=model.nodes[zone.nodes[1]].position;
            const auto c=model.nodes[zone.nodes[2]].position;
            const auto d=model.nodes[zone.nodes[3]].position;
            const double determinant=dot(b-a,cross(c-a,d-a));
            if(cell.nodes.size()!=4 && determinant<=kTiny)throw std::runtime_error("mixed element has inverted/degenerate integration tetrahedron");
            zone.volume_m3=tetra_geometry(model,zone).volume;
            zone.integration_weight=cell.nodes.size()==4?1.0:0.5;
            parent_volume+=zone.volume_m3*zone.integration_weight;
            cell.subzones.push_back(model.zones.size());model.zones.push_back(zone);
        }
        for(auto i:cell.subzones)model.zones[i].parent_volume_m3=parent_volume;
        if(mixed)model.physical_cells.push_back(cell);
    }

    Vec3 low = model.nodes.front().position;
    Vec3 high = low;
    for (const Node& node : model.nodes) {
        low.x = std::min(low.x, node.position.x);
        low.y = std::min(low.y, node.position.y);
        low.z = std::min(low.z, node.position.z);
        high.x = std::max(high.x, node.position.x);
        high.y = std::max(high.y, node.position.y);
        high.z = std::max(high.z, node.position.z);
    }
    const double scale = std::max({high.x - low.x, high.y - low.y,
                                   high.z - low.z, 1.0});
    const double tolerance = scale * 1.0e-7;
    for (Node& node : model.nodes) {
        const bool bottom = std::abs(node.position.z - low.z) <= tolerance;
        const bool side_x = std::abs(node.position.x - low.x) <= tolerance ||
                            std::abs(node.position.x - high.x) <= tolerance;
        const bool side_y = std::abs(node.position.y - low.y) <= tolerance ||
                            std::abs(node.position.y - high.y) <= tolerance;
        // The generated FLAC command fixes Bottom in XYZ and constrains the
        // normal component on West/East and South/North for both project types.
        // Interior tunnel/slope faces are not on the outer bounding planes and
        // therefore remain free.
        node.fixed = {bottom || side_x, bottom || side_y, bottom};
    }
    model.config.width_m = std::max(high.x - low.x, tolerance);
    model.config.length_m = std::max(high.y - low.y, tolerance);
    model.config.height_m = std::max(high.z - low.z, tolerance);
    initialize_interfaces(model);
    append_model_stages(model);
    return model;
}

RockModel import_abaqus_c3d4(const std::string& file_path,
                             const ModelConfig& config) {
    std::ifstream input_stream(file_path);
    if (!input_stream) throw std::runtime_error("cannot open Abaqus INP file");
    return import_abaqus_c3d4_stream(input_stream, config);
}

RockModel import_abaqus_c3d4_text(const std::string& inp_text,
                                  const ModelConfig& config) {
    std::istringstream input_stream(inp_text);
    return import_abaqus_c3d4_stream(input_stream, config);
}

static RockResult solve_stages(RockModel& model, const ProgressCallback& progress, bool initialize) {
    RockResult result;
    if (initialize) {
        // Incremental Interface mode follows the paper workflow: establish the
        // common equilibrium at unreduced continuum/interface strength, then
        // apply the requested reduction in its explicit stage. Legacy models
        // (scope=false, interfaces empty) retain their previous behaviour.
        if(model.config.strength_reduction_interfaces_only) {
            model.config.material.strength_reduction_factor=1.0;
            for(auto& region:model.config.region_materials)region.material.strength_reduction_factor=1.0;
            reset_interface_strength_reduction(model);
        }
        initialize_material_state(model);
        model.cables.clear();
        install_cables(model, false);
    }
    if (model.config.solver.large_strain) {
        result.warnings.emplace_back(
            "largeStrain was requested; this release preserves the package's FLAC small-strain setting");
    }
    const int stage_count = static_cast<int>(model.stages.size());
    for (std::size_t stage_position = 0; stage_position < model.stages.size();
         ++stage_position) {
        const Stage& stage = model.stages[stage_position];
        if (stage.type == StageType::CableInstallation) {
            if(result.stages.empty() || !result.stages.back().converged)
                throw std::runtime_error("初始平衡尚未收敛，未安装锚杆；请先完成初始平衡");
            install_cables(model, true);
        } else if(stage.type == StageType::CableRelease){
            if(result.stages.empty() || !result.stages.back().converged) {
                const StageResult& holding = result.stages.back();
                const auto cable_balance = cable_force_balance(model);
                std::ostringstream message;
                message << std::setprecision(6) << std::scientific
                        << "张拉阶段未收敛，未解除控制：已运行 "
                        << holding.cycles << '/' << model.config.solver.max_cycles
                        << " cycle，平均力比 " << holding.force_ratio
                        << "，锚索平均残差比 " << cable_balance.average()
                        << "，锚索缩放局部残差比 " << cable_balance.scaled_local
                        << "，目标 " << model.config.solver.convergence_tolerance
                        << "。若残差已形成平台，请放宽收敛力比或调整锚索参数；只有残差仍持续下降时才适合增加计算步数";
                throw std::runtime_error(message.str());
            }
            for(auto& c:model.cables)c.tension_control=false;
        } else if (stage.type == StageType::Excavation) {
            deactivate_tunnel_segment(model, stage.station_start_m, stage.station_end_m);
        } else if (stage.type == StageType::Support) {
            install_support(model, stage.support_stiffness_n_m);
        } else if (stage.type == StageType::StrengthReduction) {
            const double factor=std::max(stage.strength_reduction_factor, 1.0);
            if(model.config.strength_reduction_interfaces_only && model.config.strength_reduction_continuum_only) {
                // Region material bands only: keep host continuum and Interfaces at original strength.
                model.config.material.strength_reduction_factor=1.0;
                for(auto& region:model.config.region_materials)region.material.strength_reduction_factor=factor;
                reset_interface_strength_reduction(model);
            } else if(model.config.strength_reduction_interfaces_only) {
                model.config.material.strength_reduction_factor=1.0;
                for(auto& region:model.config.region_materials)region.material.strength_reduction_factor=1.0;
                set_interface_strength_reduction(model,factor);
            } else {
                model.config.material.strength_reduction_factor=factor;
                for(auto& region:model.config.region_materials)region.material.strength_reduction_factor=factor;
                if(model.config.strength_reduction_continuum_only) reset_interface_strength_reduction(model);
                else set_interface_strength_reduction(model,factor);
            }
        }
        StageResult stage_result =
            equilibrate(model, stage.name, static_cast<int>(stage_position) + 1,
                        stage_count, progress);
        result.total_cycles += stage_result.cycles;
        result.stages.push_back(stage_result);
    }
    if(!model.physical_cells.empty())for(auto& stage:result.stages)stage.active_zones=model.physical_cells.size();
    result.converged = !result.stages.empty() &&
                       std::all_of(result.stages.begin(), result.stages.end(),
                                   [](const StageResult& stage) {
                                       return stage.converged;
                                   });
    if (!result.stages.empty()) {
        result.final_force_ratio = result.stages.back().force_ratio;
        result.final_force_ratio_local =
            result.stages.back().force_ratio_local;
        result.final_force_ratio_maximum =
            result.stages.back().force_ratio_maximum;
        result.final_unbalanced_force_maximum_n =
            result.stages.back().unbalanced_force_maximum_n;
        result.max_displacement_m = result.stages.back().max_displacement_m;
        result.active_zones = result.stages.back().active_zones;
    }
    for (const Zone& zone : model.zones) {
        if (zone.active) {
            result.max_reinforce_requirement = std::max(
                result.max_reinforce_requirement,
                std::abs(zone.smrm.reinforce_requirement_mpa));
        }
    }
    if (!result.converged) {
        result.warnings.emplace_back(
            "one or more stages reached maxCycles before the requested force ratio");
    }
    return result;
}

RockResult solve(RockModel& model, const ProgressCallback& progress) {
    return solve_stages(model, progress, true);
}

RockResult solve_baseline(RockModel& model, const ProgressCallback& progress) {
    if (model.scenario != Scenario::Slope || model.config.solver.large_strain)
        throw std::invalid_argument("基准状态仅支持小应变边坡");
    model.config.cables.clear();
    model.config.material.strength_reduction_factor = 1;
    for (auto& r : model.config.region_materials) r.material.strength_reduction_factor = 1;
    reset_interface_strength_reduction(model);
    model.stages = {{"Unreinforced baseline", StageType::Initial}};
    return solve_stages(model, progress, true);
}

RockResult continue_slope(RockModel& model, const RockResult& baseline,
                         const ModelConfig& branch, bool allow_unconverged,
                         const ProgressCallback& progress) {
    if (model.scenario != Scenario::Slope || model.config.solver.large_strain || branch.solver.large_strain)
        throw std::invalid_argument("仅支持小应变边坡状态恢复");
    if (!model.cables.empty()) throw std::invalid_argument("方案必须从未加固基准恢复，不能叠加其他方案");
    if (!baseline.converged && !allow_unconverged)
        throw std::invalid_argument("基准未收敛；仅可明确选择研究模式继续");
    model.config.solver = branch.solver;
    model.config.cables = branch.cables;
    model.config.strength_reduction_interfaces_only = branch.strength_reduction_interfaces_only;
    model.config.strength_reduction_continuum_only = branch.strength_reduction_continuum_only;
    if(branch.interfaces.size()!=model.interfaces.size())
        throw std::invalid_argument("恢复分支的 Interface 配置与基准状态不一致");
    model.config.interfaces=branch.interfaces;
    for(std::size_t i=0;i<model.interfaces.size();++i)model.interfaces[i].config=branch.interfaces[i];
    for (auto& c : model.config.cables) {
        if (!c.after_initial) throw std::invalid_argument("从基准激活时，锚杆安装阶段必须为平衡后");
    }
    install_cables(model, true);
    model.stages = {{branch.cables.empty() ? "Unreinforced continuation" : "Reinforced equilibrium", StageType::Initial}};
    if(std::any_of(branch.cables.begin(),branch.cables.end(),[](const auto& c){return c.pretension>0;})){
        model.stages.front().name="Prestress holding equilibrium";
        model.stages.push_back(Stage{.name="Prestress released equilibrium",.type=StageType::CableRelease});
    }
    if (branch.material.strength_reduction_factor > 1) {
        Stage reduction; reduction.name = "Strength reduction"; reduction.type = StageType::StrengthReduction;
        reduction.strength_reduction_factor = branch.material.strength_reduction_factor;
        model.stages.push_back(reduction);
    }
    auto result = solve_stages(model, progress, false);
    if (!baseline.converged) result.warnings.push_back("研究用未收敛基准：本结果继承该基准，不能视为已收敛初始状态");
    return result;
}

RockResult continue_water_stage(RockModel& model, const RockResult& previous,
                                const WaterConfig& water,
                                const SolverSettings& solver,
                                const ProgressCallback& progress) {
    if (!previous.converged)
        throw std::invalid_argument("水位阶段要求上一阶段已收敛");
    if (model.config.solver.large_strain || solver.large_strain)
        throw std::invalid_argument("水位阶段目前仅支持小应变恢复");
    if (model.config.solver.flac3d7_static_mass != solver.flac3d7_static_mass)
        throw std::invalid_argument("水位阶段不能切换 FLAC3D 静力质量算法，请重新建立初始阶段");
    if (std::any_of(model.cables.begin(), model.cables.end(),
                    [](const auto& cable) { return cable.tension_control; }))
        throw std::invalid_argument("水位阶段前必须解除锚索恒张拉控制");

    ModelConfig next_config = model.config;
    next_config.water = water;
    next_config.solver = solver;
    validate_water(next_config);
    const auto old_pressure = zone_pore_pressures(model);
    model.config.water = water;
    model.config.solver = solver;
    const auto new_pressure = zone_pore_pressures(model);
    // FLAC3D's no-flow water-level stage: preserve the instantaneous
    // effective stress when assigning new pore pressure, then re-equilibrate
    // total stress under the new body weight and boundary water pressure.
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const double correction = continuum_biot_coefficient(model)*
            (new_pressure[zi] - old_pressure[zi]);
        for (int axis = 0; axis < 3; ++axis)
            model.zones[zi].stress_pa[axis] -= correction;
    }
    model.stages = {{"Water-level change", StageType::Initial}};
    return solve_stages(model, progress, false);
}

RockResult continue_hydraulic_stage(RockModel& model, const RockResult& previous,
                                    const HydraulicState& hydraulic,
                                    const SolverSettings& solver,
                                    const ProgressCallback& progress) {
    if (!previous.converged)
        throw std::invalid_argument("瞬态渗流力学阶段要求上一阶段已收敛");
    if (!hydraulic.valid || !hydraulic.converged ||
        hydraulic.nodal_pore_pressure_pa.size() != model.nodes.size())
        throw std::invalid_argument("瞬态渗流力学阶段需要完整且收敛的孔压场");
    if (model.config.solver.large_strain || solver.large_strain)
        throw std::invalid_argument("瞬态渗流力学阶段目前仅支持小应变恢复");
    if (model.config.solver.flac3d7_static_mass != solver.flac3d7_static_mass)
        throw std::invalid_argument("瞬态渗流力学阶段不能切换 FLAC3D 静力质量算法");
    if (std::any_of(model.cables.begin(), model.cables.end(),
                    [](const auto& cable) { return cable.tension_control; }))
        throw std::invalid_argument("瞬态渗流力学阶段前必须解除锚索恒张拉控制");

    const auto old_pressure = zone_pore_pressures(model);
    model.hydraulic = hydraulic;
    model.config.solver = solver;
    const auto new_pressure = zone_pore_pressures(model);
    // Preserve instantaneous effective stress when pore pressure changes,
    // then let the quasi-static solve equilibrate total stress and deformation.
    for (std::size_t zi = 0; zi < model.zones.size(); ++zi) {
        const double correction = continuum_biot_coefficient(model)*
            (new_pressure[zi] - old_pressure[zi]);
        for (int axis = 0; axis < 3; ++axis)
            model.zones[zi].stress_pa[axis] -= correction;
    }
    model.stages = {{"Transient seepage mechanics", StageType::Initial}};
    return solve_stages(model, progress, false);
}

RockResult reduction_trial(RockModel& model,double factor,const ProgressCallback& progress){
    if(!std::isfinite(factor)||factor<1)throw std::invalid_argument("折减系数必须不小于 1");
    if(std::any_of(model.cables.begin(),model.cables.end(),[](const auto& c){return c.tension_control;}))throw std::invalid_argument("必须解除恒定张拉控制后再分析稳定性");
    model.stages={{"Stability strength reduction",StageType::StrengthReduction}};
    model.stages[0].strength_reduction_factor=factor;
    // Stability trials advance the same restored state in bounded blocks.  A
    // model-level maxCycles value such as 30000 is appropriate for a normal
    // equilibrium solve, but using it for each of the six plateau blocks can
    // silently turn one failed RF trial into 180000 cycles.  Four thousand is
    // the calibrated block size used by the Case05 audit; smaller user limits
    // remain respected.
    const int requested_cycles=model.config.solver.max_cycles;
    model.config.solver.max_cycles=std::min(requested_cycles,4000);
    struct CycleLimitRestore {
        SolverSettings& solver;
        int value;
        ~CycleLimitRestore(){solver.max_cycles=value;}
    } restore_cycle_limit{model.config.solver,requested_cycles};
    RockResult result; double previous=0; int steady=0,total=0;
    for(int block=0;block<6;++block){
        ProgressCallback block_progress;
        if(progress)block_progress=[&,block](const SolveProgress& value){
            auto adjusted=value;
            adjusted.stage_index=block+1;
            adjusted.stage_count=6;
            progress(adjusted);
        };
        result=solve_stages(model,block_progress,false);total+=result.total_cycles;
        if(result.converged){result.total_cycles=total;result.warnings.push_back("FS_TRIAL_STABLE");return result;}
        const double ratio=result.final_force_ratio;
        if(previous>0 && std::abs(ratio/previous-1)<0.1)++steady;else steady=0;
        if(steady>=2){result.total_cycles=total;result.warnings.push_back("FS_TRIAL_UNSTABLE_PLATEAU");return result;}
        previous=ratio;
    }
    result.total_cycles=total;result.warnings.push_back("FS_TRIAL_INCONCLUSIVE");return result;
}
#include "checkpoint.inc"

ResourceEstimate estimate_resources(const ModelConfig& config) {
    ResourceEstimate estimate;
    estimate.nodes = static_cast<std::size_t>(std::max(config.nx, 1) + 1) *
                     (std::max(config.ny, 1) + 1) *
                     (std::max(config.nz, 1) + 1);
    estimate.zones = static_cast<std::size_t>(std::max(config.nx, 1)) *
                     std::max(config.ny, 1) * std::max(config.nz, 1) * 6;
    estimate.degrees_of_freedom = estimate.nodes * 3;
    const std::size_t bytes = estimate.nodes * 176 + estimate.zones * 264;
    estimate.memory_mb = static_cast<double>(bytes) / (1024.0 * 1024.0);
    estimate.recommended_backend = estimate.zones <= 200000 && estimate.memory_mb <= 512.0
                                       ? "wasm"
                                       : "native_server";
    return estimate;
}

ResourceEstimate estimate_resources(const RockModel& model) {
    ResourceEstimate estimate;
    estimate.nodes = model.nodes.size();
    estimate.zones = model.zones.size();
    estimate.degrees_of_freedom = estimate.nodes * 3;
    const std::size_t bytes = estimate.nodes * 176 + estimate.zones * 264;
    estimate.memory_mb = static_cast<double>(bytes) / (1024.0 * 1024.0);
    estimate.recommended_backend = estimate.zones <= 200000 && estimate.memory_mb <= 512.0
                                       ? "wasm"
                                       : "native_server";
    return estimate;
}

std::string result_to_json(const RockModel& model, const RockResult& result) {
    if(!model.physical_cells.empty()) {
        RockModel physical=model;
        physical.zones.clear();physical.physical_cells.clear();
        physical.hydraulic.zone_head_m.clear();
        physical.hydraulic.zone_saturation.clear();
        physical.hydraulic.zone_darcy_velocity_m_s.clear();
        physical.hydraulic.zone_conductivity_m_s.clear();
        physical.hydraulic.zone_joint_aperture_m.clear();
        physical.hydraulic.zone_joint_normal_stress_mpa.clear();
        for(auto& snapshot:physical.hydraulic.history){
            snapshot.zone_head_m.clear();
            snapshot.zone_saturation.clear();
            snapshot.zone_darcy_velocity_m_s.clear();
            snapshot.zone_conductivity_m_s.clear();
            snapshot.zone_joint_aperture_m.clear();
            snapshot.zone_joint_normal_stress_mpa.clear();
            snapshot.zone_stress_pa.clear();
            snapshot.zone_porosity.clear();
        }
        for(const auto& cell:model.physical_cells){
            auto representative=cell.subzones.front();
            for(auto i:cell.subzones)if(model.zones[i].failure_index>model.zones[representative].failure_index)representative=i;
            Zone zone=model.zones[representative];
            zone.id=cell.id;zone.nodes=cell.nodes;zone.stress_pa={};zone.strain={};zone.volume_m3=0;
            zone.plastic_state=0;
            for(auto i:cell.subzones){
                const auto& sub=model.zones[i];const double v=sub.volume_m3*sub.integration_weight;
                zone.volume_m3+=v;zone.plastic_state|=sub.plastic_state;
                for(int k=0;k<6;++k){zone.stress_pa[k]+=sub.stress_pa[k]*v;zone.strain[k]+=sub.strain[k]*v;}
            }
            for(int k=0;k<6;++k){zone.stress_pa[k]/=zone.volume_m3;zone.strain[k]/=zone.volume_m3;}
            physical.zones.push_back(zone);
            if(model.hydraulic.valid) {
                double head=0.0,saturation=0.0,total_weight=0.0;
                Vec3 velocity{};
                std::array<double,6> conductivity{};
                std::array<double,3> aperture{},normal_stress{};
                for(const auto i:cell.subzones) {
                    const double weight=model.zones[i].volume_m3*model.zones[i].integration_weight;
                    total_weight+=weight;
                    if(i<model.hydraulic.zone_head_m.size())head+=model.hydraulic.zone_head_m[i]*weight;
                    if(i<model.hydraulic.zone_saturation.size())saturation+=model.hydraulic.zone_saturation[i]*weight;
                    if(i<model.hydraulic.zone_darcy_velocity_m_s.size()) {
                        velocity.x+=model.hydraulic.zone_darcy_velocity_m_s[i].x*weight;
                        velocity.y+=model.hydraulic.zone_darcy_velocity_m_s[i].y*weight;
                        velocity.z+=model.hydraulic.zone_darcy_velocity_m_s[i].z*weight;
                    }
                    if(i<model.hydraulic.zone_conductivity_m_s.size())
                        for(std::size_t k=0;k<6;++k)conductivity[k]+=model.hydraulic.zone_conductivity_m_s[i][k]*weight;
                    if(i<model.hydraulic.zone_joint_aperture_m.size())
                        for(std::size_t k=0;k<3;++k)aperture[k]+=model.hydraulic.zone_joint_aperture_m[i][k]*weight;
                    if(i<model.hydraulic.zone_joint_normal_stress_mpa.size())
                        for(std::size_t k=0;k<3;++k)normal_stress[k]+=model.hydraulic.zone_joint_normal_stress_mpa[i][k]*weight;
                }
                if(total_weight>0.0){
                    head/=total_weight;saturation/=total_weight;
                    velocity.x/=total_weight;velocity.y/=total_weight;velocity.z/=total_weight;
                    for(auto& value:conductivity)value/=total_weight;
                    for(auto& value:aperture)value/=total_weight;
                    for(auto& value:normal_stress)value/=total_weight;
                }
                physical.hydraulic.zone_head_m.push_back(head);
                physical.hydraulic.zone_saturation.push_back(saturation);
                physical.hydraulic.zone_darcy_velocity_m_s.push_back(velocity);
                physical.hydraulic.zone_conductivity_m_s.push_back(conductivity);
                physical.hydraulic.zone_joint_aperture_m.push_back(aperture);
                physical.hydraulic.zone_joint_normal_stress_mpa.push_back(normal_stress);
            }
            for(std::size_t hi=0;hi<physical.hydraulic.history.size();++hi){
                const auto& source=model.hydraulic.history[hi];
                double head=0.0,saturation=0.0,porosity=0.0,total_weight=0.0;
                Vec3 velocity{};
                std::array<double,6> stress{};
                std::array<double,6> conductivity{};
                std::array<double,3> aperture{},normal_stress{};
                for(const auto i:cell.subzones){
                    const double weight=model.zones[i].volume_m3*model.zones[i].integration_weight;
                    total_weight+=weight;
                    if(i<source.zone_head_m.size())head+=source.zone_head_m[i]*weight;
                    if(i<source.zone_saturation.size())saturation+=source.zone_saturation[i]*weight;
                    if(i<source.zone_porosity.size())porosity+=source.zone_porosity[i]*weight;
                    if(i<source.zone_darcy_velocity_m_s.size()){
                        velocity.x+=source.zone_darcy_velocity_m_s[i].x*weight;
                        velocity.y+=source.zone_darcy_velocity_m_s[i].y*weight;
                        velocity.z+=source.zone_darcy_velocity_m_s[i].z*weight;
                    }
                    if(i<source.zone_stress_pa.size())
                        for(std::size_t k=0;k<6;++k)stress[k]+=source.zone_stress_pa[i][k]*weight;
                    if(i<source.zone_conductivity_m_s.size())
                        for(std::size_t k=0;k<6;++k)conductivity[k]+=source.zone_conductivity_m_s[i][k]*weight;
                    if(i<source.zone_joint_aperture_m.size())
                        for(std::size_t k=0;k<3;++k)aperture[k]+=source.zone_joint_aperture_m[i][k]*weight;
                    if(i<source.zone_joint_normal_stress_mpa.size())
                        for(std::size_t k=0;k<3;++k)normal_stress[k]+=source.zone_joint_normal_stress_mpa[i][k]*weight;
                }
                if(total_weight>0.0){
                    head/=total_weight;
                    saturation/=total_weight;
                    porosity/=total_weight;
                    velocity.x/=total_weight;velocity.y/=total_weight;velocity.z/=total_weight;
                    for(auto& value:stress)value/=total_weight;
                    for(auto& value:conductivity)value/=total_weight;
                    for(auto& value:aperture)value/=total_weight;
                    for(auto& value:normal_stress)value/=total_weight;
                }
                physical.hydraulic.history[hi].zone_head_m.push_back(head);
                physical.hydraulic.history[hi].zone_saturation.push_back(saturation);
                physical.hydraulic.history[hi].zone_darcy_velocity_m_s.push_back(velocity);
                if(!source.zone_stress_pa.empty())
                    physical.hydraulic.history[hi].zone_stress_pa.push_back(stress);
                if(!source.zone_porosity.empty())
                    physical.hydraulic.history[hi].zone_porosity.push_back(porosity);
                if(!source.zone_conductivity_m_s.empty())
                    physical.hydraulic.history[hi].zone_conductivity_m_s.push_back(conductivity);
                if(!source.zone_joint_aperture_m.empty())
                    physical.hydraulic.history[hi].zone_joint_aperture_m.push_back(aperture);
                if(!source.zone_joint_normal_stress_mpa.empty())
                    physical.hydraulic.history[hi].zone_joint_normal_stress_mpa.push_back(normal_stress);
            }
        }
        RockResult output=result;
        output.warnings.emplace_back("Mixed-cell stress/strain are volume averages; SMRM diagnostics represent the most critical integration subzone.");
        output.active_zones=physical.zones.size();
        for(auto& stage:output.stages)stage.active_zones=physical.zones.size();
        return result_to_json(physical,output);
    }
    std::ostringstream out;
    out << std::setprecision(std::numeric_limits<double>::max_digits10);
    const auto pore_pressure = zone_pore_pressures(model);
    out << "{\"supports\":" << cables_to_json(model) << ",\"version\":\"1.0\",\"solverVersion\":\""
        << json_escape(result.solver_version) << "\",\"scenario\":\""
        << scenario_name(model.scenario) << "\",\"material\":\""
        << material_name(model.config.material.type) << "\",\"summary\":{";
    out << "\"converged\":" << (result.converged ? "true" : "false")
        << ",\"totalCycles\":" << result.total_cycles
        << ",\"forceRatio\":" << result.final_force_ratio
        << ",\"forceRatioAverage\":" << result.final_force_ratio
        << ",\"forceRatioLocal\":" << result.final_force_ratio_local
        << ",\"forceRatioMaximum\":" << result.final_force_ratio_maximum
        << ",\"cableForceRatioAverage\":" << cable_force_balance(model).average()
        << ",\"cableForceRatioScaledLocal\":" << cable_force_balance(model).scaled_local
        << ",\"cableUnbalancedForceMaximumN\":" << cable_force_balance(model).maximum_force
        << ",\"unbalancedForceMaximumN\":"
        << result.final_unbalanced_force_maximum_n
        << ",\"maxDisplacementM\":" << result.max_displacement_m
        << ",\"maxReinforceRequirement\":" << result.max_reinforce_requirement
        << ",\"activeZones\":" << result.active_zones;
    if(model.hydraulic.valid) {
        out << ",\"hydraulicConverged\":" << (model.hydraulic.converged?"true":"false")
            << ",\"hydraulicIterations\":" << model.hydraulic.iterations
            << ",\"hydraulicRelativeResidual\":" << model.hydraulic.relative_residual
            << ",\"hydraulicInflowM3S\":" << model.hydraulic.inflow_m3_s
            << ",\"hydraulicOutflowM3S\":" << model.hydraulic.outflow_m3_s
            << ",\"hydraulicNetOutflowM3S\":" << model.hydraulic.net_outflow_m3_s
            << ",\"hydraulicMechanicalStorageChangeM3\":" << model.hydraulic.mechanical_storage_change_m3
            << ",\"hydraulicMassBalanceError\":" << model.hydraulic.relative_mass_balance_error
            << ",\"hydraulicMassBalanceResidualM3\":" << model.hydraulic.mass_balance_residual_m3
            << ",\"hydraulicMassBalanceRoundoffM3\":" << model.hydraulic.mass_balance_roundoff_m3
            << ",\"hydraulicTimeS\":" << model.hydraulic.time_s
            << ",\"hydraulicRequestedDurationS\":" << model.config.seepage.duration_s
            << ",\"hydraulicBiotCoefficient\":" << continuum_biot_coefficient(model)
            << ",\"hydraulicFailureReason\":\"" << json_escape(model.hydraulic.failure_reason) << '"'
            << ",\"rainfallInfiltrationM3S\":" << model.hydraulic.rainfall_infiltration_m3_s
            << ",\"rainfallRunoffM3S\":" << model.hydraulic.rainfall_runoff_m3_s;
    }
    out << "},\"mesh\":{";
    out << "\"nodeIds\":[";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        if (i) out << ',';
        out << model.nodes[i].id;
    }
    out << "],\"nodes\":[";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        if (i) out << ',';
        write_numeric_array(out, std::array{model.nodes[i].position.x,
                                            model.nodes[i].position.y,
                                            model.nodes[i].position.z});
    }
    out << "],\"elementIds\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << model.zones[i].id;
    }
    out << "],\"tetrahedra\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        write_numeric_array(out, model.zones[i].nodes);
    }
    out << "],\"active\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << (model.zones[i].active ? "true" : "false");
    }
    out << "],\"materialGroups\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        const auto index = model.zones[i].material_index;
        out << '\"' << json_escape(!model.zones[i].material_group.empty() ? model.zones[i].material_group : (index ? model.config.region_materials.at(index - 1).name : "")) << '\"';
    }
    out << "]},\"node\":{";
    out << "\"displacementM\":[";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        if (i) out << ',';
        write_numeric_array(out, std::array{model.nodes[i].displacement.x,
                                            model.nodes[i].displacement.y,
                                            model.nodes[i].displacement.z});
    }
    out << "],\"displacementMagnitudeM\":[";
    for (std::size_t i = 0; i < model.nodes.size(); ++i) {
        if (i) out << ',';
        out << norm(model.nodes[i].displacement);
    }
    out << "],\"hydraulicHeadM\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.hydraulic.head_m.size();++i){if(i)out<<',';out<<model.hydraulic.head_m[i];}
    out << "],\"porePressureMpa\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.hydraulic.nodal_pore_pressure_pa.size();++i){if(i)out<<',';out<<model.hydraulic.nodal_pore_pressure_pa[i]*1e-6;}
    out << "],\"matricSuctionMpa\":[";
    if(model.hydraulic.valid && model.config.seepage.unsaturated)
        for(std::size_t i=0;i<model.hydraulic.head_m.size();++i){
            if(i)out<<',';
            out<<std::max(0.0,model.nodes[i].position.z-model.hydraulic.head_m[i])*
                model.config.seepage.fluid_density_kg_m3*model.config.gravity_m_s2*1e-6;
        }
    out << "],\"saturation\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.hydraulic.nodal_saturation.size();++i){if(i)out<<',';out<<model.hydraulic.nodal_saturation[i];}
    out << "]},\"zone\":{";
    out << "\"stressMpa\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        write_numeric_array(out, model.zones[i].stress_pa, 1.0e-6);
    }
    out << "],\"porePressureMpa\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << pore_pressure[i] * 1.0e-6;
    }
    out << "],\"densityKgM3\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        Vec3 center{};
        for(const auto node:model.zones[i].nodes)center+=model.nodes[node].position;
        center=center*(1.0/static_cast<double>(model.zones[i].nodes.size()));
        out << zone_density(model,model.zones[i],center);
    }
    out << "],\"hydraulicHeadM\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';out<<(i<model.hydraulic.zone_head_m.size()?model.hydraulic.zone_head_m[i]:0.0);}
    out << "],\"saturation\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';out<<(i<model.hydraulic.zone_saturation.size()?model.hydraulic.zone_saturation[i]:0.0);}
    out << "],\"matricSuctionMpa\":[";
    if(model.hydraulic.valid && model.config.seepage.unsaturated)
        for(std::size_t i=0;i<model.zones.size();++i){
            if(i)out<<',';
            double center_z=0.0;
            for(const auto node:model.zones[i].nodes)center_z+=model.nodes[node].position.z*0.25;
            out<<std::max(0.0,center_z-model.hydraulic.zone_head_m[i])*
                model.config.seepage.fluid_density_kg_m3*model.config.gravity_m_s2*1e-6;
        }
    out << "],\"darcyVelocityMS\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';const auto v=i<model.hydraulic.zone_darcy_velocity_m_s.size()?model.hydraulic.zone_darcy_velocity_m_s[i]:Vec3{};write_numeric_array(out,std::array{v.x,v.y,v.z});}
    out << "],\"darcySpeedMS\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';const auto v=i<model.hydraulic.zone_darcy_velocity_m_s.size()?model.hydraulic.zone_darcy_velocity_m_s[i]:Vec3{};out<<norm(v);}
    out << "],\"hydraulicConductivityMS\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';write_numeric_array(out,i<model.hydraulic.zone_conductivity_m_s.size()?model.hydraulic.zone_conductivity_m_s[i]:std::array<double,6>{});}
    out << "],\"jointHydraulicApertureM\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';write_numeric_array(out,i<model.hydraulic.zone_joint_aperture_m.size()?model.hydraulic.zone_joint_aperture_m[i]:std::array<double,3>{});}
    out << "],\"jointHydraulicNormalStressMpa\":[";
    if(model.hydraulic.valid)for(std::size_t i=0;i<model.zones.size();++i){if(i)out<<',';write_numeric_array(out,i<model.hydraulic.zone_joint_normal_stress_mpa.size()?model.hydraulic.zone_joint_normal_stress_mpa[i]:std::array<double,3>{});}
    out << "],\"effectiveStressMpa\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        auto stress = model.zones[i].stress_pa;
        for (int k = 0; k < 3; ++k) stress[k] +=
            continuum_biot_coefficient(model)*pore_pressure[i];
        write_numeric_array(out, stress, 1.0e-6);
    }
    out << "],\"strain\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        write_numeric_array(out, model.zones[i].strain);
    }
    out << "],\"sigma1CompressionMpa\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << std::max(0.0, -principal_values(model.zones[i].stress_pa)[0]) * 1.0e-6;
    }
    out << "],\"plasticState\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << static_cast<int>(model.zones[i].plastic_state);
    }
    out << "],\"reinforceRequirement\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << model.zones[i].reinforce_requirement;
    }
    out << "],\"jointState\":[";
    for (std::size_t i = 0; i < model.zones.size(); ++i) {
        if (i) out << ',';
        out << '"' << joint_state_name(model.zones[i].joint_state) << '"';
    }
    out << "],\"smrm\":{";
    bool first_smrm_field = true;
    const auto scalar_field = [&](const char* name, const auto& getter) {
        if (!first_smrm_field) out << ',';
        first_smrm_field = false;
        out << '\"' << name << "\":[";
        for (std::size_t i = 0; i < model.zones.size(); ++i) {
            if (i) out << ',';
            const double value = static_cast<double>(getter(model.zones[i].smrm));
            out << (std::isfinite(value) ? value : 0.0);
        }
        out << ']';
    };
    const auto triplet_field = [&](const char* name, const auto& getter) {
        if (!first_smrm_field) out << ',';
        first_smrm_field = false;
        out << '\"' << name << "\":[";
        for (std::size_t i = 0; i < model.zones.size(); ++i) {
            if (i) out << ',';
            write_numeric_array(out, getter(model.zones[i].smrm));
        }
        out << ']';
    };
    scalar_field("rockMassEmGpa", [](const SMRMDiagnostics& d) { return d.rock_mass_em_gpa; });
    scalar_field("smrmScore", [](const SMRMDiagnostics& d) { return d.smrm_score; });
    scalar_field("failureProbability", [](const SMRMDiagnostics& d) { return d.failure_probability; });
    scalar_field("safetyRest", [](const SMRMDiagnostics& d) { return d.safety_rest; });
    scalar_field("reinforceRequirementMpa", [](const SMRMDiagnostics& d) { return d.reinforce_requirement_mpa; });
    scalar_field("stabilityCoefficient", [](const SMRMDiagnostics& d) { return d.stability_coefficient; });
    scalar_field("sigma1cminMpa", [](const SMRMDiagnostics& d) { return d.sigma1cmin_mpa; });
    scalar_field("energyDensityMjM3", [](const SMRMDiagnostics& d) { return d.energy_density_mj_m3; });
    scalar_field("u0vMjM3", [](const SMRMDiagnostics& d) { return d.u0v_mj_m3; });
    scalar_field("u0dMjM3", [](const SMRMDiagnostics& d) { return d.u0d_mj_m3; });
    scalar_field("triaxialCompressiveStrengthMpa", [](const SMRMDiagnostics& d) { return d.triaxial_compressive_strength_mpa; });
    scalar_field("limitConfiningPressureMpa", [](const SMRMDiagnostics& d) { return d.limit_confining_pressure_mpa; });
    scalar_field("rockburstIndex", [](const SMRMDiagnostics& d) { return d.rockburst_index; });
    scalar_field("rockMassFrictionDeg", [](const SMRMDiagnostics& d) { return d.rock_mass_friction_deg; });
    scalar_field("rockMassCohesionMpa", [](const SMRMDiagnostics& d) { return d.rock_mass_cohesion_mpa; });
    scalar_field("rockMassPoissonRatio", [](const SMRMDiagnostics& d) { return d.rock_mass_poisson_ratio; });
    scalar_field("rockMassUcsMpa", [](const SMRMDiagnostics& d) { return d.rock_mass_ucs_mpa; });
    scalar_field("diagSigma3JointLimitMpa", [](const SMRMDiagnostics& d) { return d.diag_sigma3_joint_limit_mpa; });
    scalar_field("diagSigma3ControllingJoint", [](const SMRMDiagnostics& d) { return d.diag_sigma3_controlling_joint; });
    triplet_field("jointK", [](const SMRMDiagnostics& d) { return d.joint_k; });
    triplet_field("jointH", [](const SMRMDiagnostics& d) { return d.joint_h; });
    triplet_field("jointN1", [](const SMRMDiagnostics& d) { return d.joint_n1; });
    triplet_field("jointN2", [](const SMRMDiagnostics& d) { return d.joint_n2; });
    triplet_field("jointN3", [](const SMRMDiagnostics& d) { return d.joint_n3; });
    triplet_field("jointUseRockSigma3", [](const SMRMDiagnostics& d) { return d.joint_use_rock_sigma3; });
    triplet_field("jointSigma3LimitMpa", [](const SMRMDiagnostics& d) { return d.joint_sigma3_limit_mpa; });
    scalar_field("bookQuasiLimitConfiningPressureMpa", [](const SMRMDiagnostics& d) { return d.book_quasi_limit_confining_pressure_mpa; });
    scalar_field("bookQuasiReinforceRequirementMpa", [](const SMRMDiagnostics& d) { return d.book_quasi_reinforce_requirement_mpa; });
    scalar_field("diagSigma3QuasiJointLimitMpa", [](const SMRMDiagnostics& d) { return d.diag_sigma3_quasi_joint_limit_mpa; });
    scalar_field("diagSigma3QuasiControllingJoint", [](const SMRMDiagnostics& d) { return d.diag_sigma3_quasi_controlling_joint; });
    triplet_field("jointQuasiSigma3LimitMpa", [](const SMRMDiagnostics& d) { return d.joint_quasi_sigma3_limit_mpa; });
    triplet_field("jointNormalStressMpa", [](const SMRMDiagnostics& d) { return d.joint_normal_stress_mpa; });
    triplet_field("jointShearStressMpa", [](const SMRMDiagnostics& d) { return d.joint_shear_stress_mpa; });
    triplet_field("jointShearStrengthMpa", [](const SMRMDiagnostics& d) { return d.joint_shear_strength_mpa; });
    scalar_field("rockMassEmXGpa", [](const SMRMDiagnostics& d) { return d.rock_mass_em_x_gpa; });
    scalar_field("rockMassEmYGpa", [](const SMRMDiagnostics& d) { return d.rock_mass_em_y_gpa; });
    scalar_field("rockMassEmZGpa", [](const SMRMDiagnostics& d) { return d.rock_mass_em_z_gpa; });
    scalar_field("confiningPressureMarginMpa", [](const SMRMDiagnostics& d) { return d.confining_pressure_margin_mpa; });
    scalar_field("bookQuasiConfiningPressureMarginMpa", [](const SMRMDiagnostics& d) { return d.book_quasi_confining_pressure_margin_mpa; });
    scalar_field("yieldState", [](const SMRMDiagnostics& d) { return d.yield_state; });
    scalar_field("hasFailed", [](const SMRMDiagnostics& d) { return d.has_failed ? 1.0 : 0.0; });
    out << "}},\"interfaces\":[";
    for(std::size_t i=0;i<model.interfaces.size();++i) {
        if(i)out<<',';
        const auto& state=model.interfaces[i];
        double max_slip=0,max_normal=0,max_shear=0;
        std::size_t shear_failed=0,shear_past=0,no_shear_failure=0,tension_failed=0;
        for(const auto& facet:state.facets) {
            const Vec3 jump=interface_average_jump(model,facet);
            const double un=dot(jump,facet.normal);
            max_slip=std::max(max_slip,norm(jump-facet.normal*un));
            max_normal=std::max(max_normal,std::abs(facet.normal_stress_pa));
            max_shear=std::max(max_shear,std::abs(facet.shear_stress_pa));
            const bool shear_now=facet.state==JointState::ShearFailure;
            const bool slipped_ever=norm(facet.plastic_shear_slip)>1.0e-14;
            shear_failed+=shear_now;
            shear_past+=!shear_now && slipped_ever;
            no_shear_failure+=!shear_now && !slipped_ever;
            tension_failed+=facet.state==JointState::TensionFailure;
        }
        out<<"{\"name\":\""<<json_escape(state.config.name)<<"\",\"sourceGroup\":\""<<json_escape(state.config.source_group)
           <<"\",\"strengthReductionFactor\":"<<state.strength_reduction_factor
           <<",\"facetCount\":"<<state.facets.size()<<",\"maxShearSlipM\":"<<max_slip
           <<",\"maxNormalStressMpa\":"<<max_normal*1e-6<<",\"maxShearStressMpa\":"<<max_shear*1e-6
           <<",\"shearFailedFacets\":"<<shear_failed
           <<",\"shearPastFacets\":"<<shear_past
           <<",\"noShearFailureFacets\":"<<no_shear_failure
           <<",\"tensionFailedFacets\":"<<tension_failed<<'}';
    }
    out << "],\"hydraulicHistory\":[";
    for(std::size_t hi=0;hi<model.hydraulic.history.size();++hi){
        if(hi)out<<',';
        const auto& snapshot=model.hydraulic.history[hi];
        out<<"{\"timeS\":"<<snapshot.time_s
           <<",\"iterations\":"<<snapshot.iterations
           <<",\"relativeResidual\":"<<snapshot.relative_residual
           <<",\"inflowM3S\":"<<snapshot.inflow_m3_s
           <<",\"outflowM3S\":"<<snapshot.outflow_m3_s
           <<",\"netOutflowM3S\":"<<snapshot.net_outflow_m3_s
           <<",\"storageChangeM3\":"<<snapshot.storage_change_m3
           <<",\"mechanicalStorageChangeM3\":"<<snapshot.mechanical_storage_change_m3
           <<",\"massBalanceError\":"<<snapshot.relative_mass_balance_error
           <<",\"massBalanceResidualM3\":"<<snapshot.mass_balance_residual_m3
           <<",\"massBalanceRoundoffM3\":"<<snapshot.mass_balance_roundoff_m3
           <<",\"rainfallInfiltrationM3S\":"<<snapshot.rainfall_infiltration_m3_s
           <<",\"rainfallRunoffM3S\":"<<snapshot.rainfall_runoff_m3_s
           <<",\"headM\":[";
        for(std::size_t i=0;i<snapshot.head_m.size();++i){if(i)out<<',';out<<snapshot.head_m[i];}
        out<<"],\"nodalBoundaryOutwardFluxM3S\":[";
        for(std::size_t i=0;i<snapshot.nodal_boundary_outward_flux_m3_s.size();++i){
            if(i)out<<',';out<<snapshot.nodal_boundary_outward_flux_m3_s[i];
        }
        out<<"],\"porePressureMpa\":[";
        for(std::size_t i=0;i<snapshot.nodal_pore_pressure_pa.size();++i){if(i)out<<',';out<<snapshot.nodal_pore_pressure_pa[i]*1.0e-6;}
        out<<"],\"saturation\":[";
        for(std::size_t i=0;i<snapshot.nodal_saturation.size();++i){if(i)out<<',';out<<snapshot.nodal_saturation[i];}
        out<<"],\"zoneHeadM\":[";
        for(std::size_t i=0;i<snapshot.zone_head_m.size();++i){if(i)out<<',';out<<snapshot.zone_head_m[i];}
        out<<"],\"zoneSaturation\":[";
        for(std::size_t i=0;i<snapshot.zone_saturation.size();++i){if(i)out<<',';out<<snapshot.zone_saturation[i];}
        out<<"],\"zonePorosity\":[";
        for(std::size_t i=0;i<snapshot.zone_porosity.size();++i){if(i)out<<',';out<<snapshot.zone_porosity[i];}
        out<<"],\"nodalBiotModulusPa\":[";
        for(std::size_t i=0;i<snapshot.nodal_biot_modulus_pa.size();++i){if(i)out<<',';out<<snapshot.nodal_biot_modulus_pa[i];}
        out<<"],\"hydraulicConductivityMS\":[";
        for(std::size_t i=0;i<snapshot.zone_conductivity_m_s.size();++i){
            if(i)out<<',';write_numeric_array(out,snapshot.zone_conductivity_m_s[i]);
        }
        out<<"],\"jointHydraulicApertureM\":[";
        for(std::size_t i=0;i<snapshot.zone_joint_aperture_m.size();++i){
            if(i)out<<',';write_numeric_array(out,snapshot.zone_joint_aperture_m[i]);
        }
        out<<"],\"jointHydraulicNormalStressMpa\":[";
        for(std::size_t i=0;i<snapshot.zone_joint_normal_stress_mpa.size();++i){
            if(i)out<<',';write_numeric_array(out,snapshot.zone_joint_normal_stress_mpa[i]);
        }
        out<<"],\"darcyVelocityMS\":[";
        for(std::size_t i=0;i<snapshot.zone_darcy_velocity_m_s.size();++i){
            if(i)out<<',';const auto& v=snapshot.zone_darcy_velocity_m_s[i];
            write_numeric_array(out,std::array{v.x,v.y,v.z});
        }
        out<<"],\"displacementM\":[";
        for(std::size_t i=0;i<snapshot.nodal_displacement_m.size();++i){
            if(i)out<<',';const auto& v=snapshot.nodal_displacement_m[i];
            write_numeric_array(out,std::array{v.x,v.y,v.z});
        }
        out<<"],\"stressMpa\":[";
        for(std::size_t i=0;i<snapshot.zone_stress_pa.size();++i){
            if(i)out<<',';write_numeric_array(out,snapshot.zone_stress_pa[i],1.0e-6);
        }
        out<<"]}";
    }
    out << "],\"stages\":[";
    for (std::size_t i = 0; i < result.stages.size(); ++i) {
        if (i) out << ',';
        const StageResult& stage = result.stages[i];
        out << "{\"name\":\"" << json_escape(stage.name) << "\",\"cycles\":"
            << stage.cycles << ",\"converged\":"
            << (stage.converged ? "true" : "false")
            << ",\"forceRatio\":" << stage.force_ratio
            << ",\"forceRatioAverage\":" << stage.force_ratio
            << ",\"forceRatioLocal\":" << stage.force_ratio_local
            << ",\"forceRatioMaximum\":" << stage.force_ratio_maximum
            << ",\"unbalancedForceMaximumN\":"
            << stage.unbalanced_force_maximum_n
            << ",\"maxDisplacementM\":" << stage.max_displacement_m
            << ",\"activeZones\":" << stage.active_zones << '}';
    }
    out << "],\"warnings\":[";
    for (std::size_t i = 0; i < result.warnings.size(); ++i) {
        if (i) out << ',';
        out << '"' << json_escape(result.warnings[i]) << '"';
    }
    out << "]}";
    return out.str();
}

}  // namespace rocksolver
