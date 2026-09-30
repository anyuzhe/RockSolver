#pragma once

#include <array>
#include <cstdint>

namespace rocksolver {

struct MaterialParameters;

using Matrix6 = std::array<std::array<double, 6>, 6>;

// Direction cosines used by the supplied reference cmodelsmrm007.dll for
// jalpha (dip direction) and jbeta (dip). Exposed to make the compatibility
// convention directly regression-testable.
std::array<double, 3> joint_direction_cosines(double dip_direction_deg,
                                              double dip_deg);

// Complete per-integration-point result set exposed by the original
// FLAC3D ModelSMRM result properties.  Stress-like values use MPa, elastic
// moduli use GPa and energy densities use MJ/m3.
struct SMRMDiagnostics {
    double rock_mass_em_gpa{};
    double smrm_score{};
    double failure_probability{};
    double safety_rest{};
    double reinforce_requirement_mpa{};
    double stability_coefficient{};
    double sigma1cmin_mpa{};
    double energy_density_mj_m3{};
    double u0v_mj_m3{};
    double u0d_mj_m3{};
    double triaxial_compressive_strength_mpa{};
    double limit_confining_pressure_mpa{};
    double rockburst_index{};
    double rock_mass_friction_deg{};
    double rock_mass_cohesion_mpa{};
    double rock_mass_poisson_ratio{};
    double rock_mass_ucs_mpa{};
    double diag_sigma3_joint_limit_mpa{};
    double diag_sigma3_controlling_joint{};
    std::array<double, 3> joint_k{};
    std::array<double, 3> joint_h{};
    std::array<double, 3> joint_n1{};
    std::array<double, 3> joint_n2{};
    std::array<double, 3> joint_n3{};
    std::array<double, 3> joint_use_rock_sigma3{};
    std::array<double, 3> joint_sigma3_limit_mpa{};
    double book_quasi_limit_confining_pressure_mpa{};
    double book_quasi_reinforce_requirement_mpa{};
    double diag_sigma3_quasi_joint_limit_mpa{};
    double diag_sigma3_quasi_controlling_joint{};
    std::array<double, 3> joint_quasi_sigma3_limit_mpa{};
    std::array<double, 3> joint_normal_stress_mpa{};
    std::array<double, 3> joint_shear_stress_mpa{};
    std::array<double, 3> joint_shear_strength_mpa{};
    double rock_mass_em_x_gpa{};
    double rock_mass_em_y_gpa{};
    double rock_mass_em_z_gpa{};
    double confining_pressure_margin_mpa{};
    double book_quasi_confining_pressure_margin_mpa{};
    std::uint8_t yield_state{};  // bit 0 shear, bit 1 tension; history is sticky
    bool has_failed{};
};

struct SMRMHistory {
    std::array<double, 3> joint_k{};
    std::array<double, 3> joint_h{};
    std::uint8_t yield_state{};
    std::uint64_t step{};
};

// FLAC-compatible material-point update. Public vectors use RockSolver order:
// [xx, yy, zz, xy, xz, yz], tension-positive stress and engineering shear
// strain. The tangent returned uses the same order and SI units (Pa).
bool update_smrm_increment(const MaterialParameters& material, double volume_m3,
                           const std::array<double, 6>& strain_increment,
                           std::array<double, 6>& stress_pa,
                           SMRMHistory& history, SMRMDiagnostics& diagnostics,
                           Matrix6& tangent_pa);

// Isotropic path with the same incremental/sign/Voigt conventions.
bool update_elastic_increment(const MaterialParameters& material,
                              const std::array<double, 6>& strain_increment,
                              std::array<double, 6>& stress_pa,
                              Matrix6& tangent_pa);

// Isotropic elastic-perfectly-plastic Mohr-Coulomb update following the
// supplied FLAC ModelMohr trial/correction sequence. Stress uses the same
// tension-positive RockSolver convention. State bits match FLAC:
// bit0 shear-now, bit1 tension-now, bit2 shear-past, bit3 tension-past.
bool update_mohr_coulomb_increment(
    const MaterialParameters& material,
    const std::array<double, 6>& strain_increment,
    std::array<double, 6>& stress_pa,
    std::uint8_t& plastic_state, Matrix6& tangent_pa);

}  // namespace rocksolver
