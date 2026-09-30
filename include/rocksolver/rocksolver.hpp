#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <string>
#include <vector>

#include "rocksolver/smrm_material.hpp"

namespace rocksolver {

constexpr const char* kVersion = "0.13.23-joint-stress-eps";

struct Vec3 {
    double x{};
    double y{};
    double z{};
};

enum class Scenario { Slope, Tunnel };
enum class MaterialType { Elastic, MohrCoulomb, SMRM };
enum class StageType { Initial, Excavation, Support, CableInstallation, StrengthReduction, CableRelease };
enum class JointState : std::uint8_t {
    None = 0,
    Closed = 1,
    Open = 2,
    Sliding = 3,
    TensionFailure = 4,
    ShearFailure = 5,
};

enum class HydraulicBoundaryKind : std::uint8_t {
    Head,
    PorePressure,
    Flux,
    Rainfall,
    Drain,
    SeepageFace,
    NoFlow,
};

enum class HydraulicBoundarySelector : std::uint8_t {
    XMin,
    XMax,
    YMin,
    YMax,
    ZMin,
    ZMax,
    Plane,
};

struct HydraulicHistoryPoint {
    double time_s{};
    double value{};
};

// Steady saturated-flow boundary. Flux is the outward Darcy flux in m/s;
// infiltration is therefore negative. A Drain imposes zero gauge pressure.
struct HydraulicBoundary {
    std::string name{"Hydraulic boundary"};
    HydraulicBoundaryKind kind{HydraulicBoundaryKind::NoFlow};
    HydraulicBoundarySelector selector{HydraulicBoundarySelector::ZMin};
    double value{};  // head (m), pressure (Pa), or outward flux (m/s)
    Vec3 plane_origin_m{};
    Vec3 plane_normal{0.0, 0.0, 1.0};
    double tolerance_m{0.0};  // 0 => mesh-scale automatic tolerance
    // Optional piecewise-linear value history used by transient flow.
    // An empty history keeps value constant.
    std::vector<HydraulicHistoryPoint> history;
};

// Hydraulic inputs are deliberately separate from JointSet's mechanical
// stiffness.  The book's aperture law uses a stress-scale parameter, while
// JointSet::normal_stiffness_pa_m has units Pa/m and is not interchangeable.
struct JointHydraulicParameters {
    bool enabled{false};
    double zero_stress_aperture_m{};
    double closure_stress_pa{};
    double connectivity{1.0};
    double minimum_aperture_m{};
    double maximum_aperture_m{};
};

struct UnsaturatedHydraulicParameters {
    double porosity{0.30};
    double residual_saturation{0.10};
    double van_genuchten_alpha_1_m{1.0};
    double van_genuchten_n{2.0};
    double mualem_l{0.5};
};

struct SteadySeepageConfig {
    bool enabled{false};
    bool solve_mechanics{false};
    double fluid_density_kg_m3{1000.0};
    // Hydraulic conductivity tensor [Kxx,Kyy,Kzz,Kxy,Kxz,Kyz], m/s.
    // Index 0 is the host material; subsequent values match region_materials.
    std::vector<std::array<double, 6>> conductivity_m_s;
    // Optional Stage 2.5 contribution based on the host SMRM joint sets:
    // K += g/(24*nu) * connectivity * lambda * aperture^3 * (I - n*n).
    bool smrm_stress_permeability{false};
    double fluid_kinematic_viscosity_m2_s{1.01e-6};
    std::vector<JointHydraulicParameters> joint_hydraulics;
    std::vector<HydraulicBoundary> boundaries;
    int max_iterations{20000};
    double relative_tolerance{1.0e-10};
    // Stage 3: saturated transient flow, integrated by backward Euler.
    bool transient{false};
    bool mechanics_each_step{false};
    double duration_s{};
    double time_step_s{};
    double initial_head_m{};
    // Specific storage Ss (1/m), indexed like conductivity_m_s.
    std::vector<double> specific_storage_1_m;
    // Stage 5 poroelastic parameters. Defaults preserve the Stage 3
    // effective-stress convention; M=0 leaves explicit Ss unchanged.
    double biot_coefficient{1.0};
    double biot_modulus_pa{0.0};
    bool bidirectional_coupling{false};
    // Opt-in small-strain, incompressible-grain (alpha=1) law:
    // n=1-(1-n0)/J and M=n0*M0/n, with J relative to t=0 equilibrium.
    bool evolving_porosity_storage{false};
    double initial_porosity{0.30};
    int coupling_max_iterations{200};
    double coupling_relative_tolerance{1.0e-7};
    // Stage 4: mixed saturated/unsaturated Richards flow. Rainfall boundary
    // values are positive precipitation intensity (m/s); infiltration is
    // capped by the book's saturated vertical capacity K33 and excess is
    // reported as runoff.
    bool unsaturated{false};
    // Opt-in FLAC3D 7 no-capillarity law: nodal p=0 below full saturation,
    // relative mobility s^2(3-2s), and explicit nodal fluid volumes.
    bool flac_no_capillarity{false};
    double initial_saturation{1.0};
    std::vector<UnsaturatedHydraulicParameters> unsaturated_materials;
    int nonlinear_max_iterations{50};
    double nonlinear_relative_tolerance{1.0e-8};
};

struct JointSet {
    std::string name{"Joint 1"};
    double dip_direction_deg{0.0};
    double dip_deg{45.0};
    double normal_density_1_m{0.0};
    double mean_radius_m{1.0};
    double cohesion_pa{0.0};
    double friction_deg{30.0};
    double tensile_strength_pa{0.0};
    double normal_stiffness_pa_m{1.0e11};
    double shear_stiffness_pa_m{5.0e10};
};

struct MaterialParameters {
    MaterialType type{MaterialType::Elastic};
    double density_kg_m3{2700.0};
    double elastic_modulus_pa{27.0e9};
    double poisson_ratio{0.25};
    double cohesion_pa{5.0e6};
    double friction_deg{45.0};
    double dilation_deg{0.0};
    double tensile_strength_pa{2.0e6};
    double ucs_pa{50.0e6};
    double fracture_toughness_mpa_sqrt_m{0.0};  // 0 => FLAC rule 0.012*UCS
    double representative_volume_m3{1.0};
    double weibull_shape{4.0};
    double strength_reduction_factor{1.0};
    // Opt-in match for the SMRM source DLL built on the 601 FLAC3D 7.0 host.
    bool flac3d7_smrm_contact{false};
    std::vector<JointSet> joints;
};

struct RegionMaterial {
    std::string name;
    MaterialParameters material;
};

struct SolverSettings {
    int max_cycles{3000};
    double convergence_tolerance{1.0e-3};
    // FLAC-style local damping coefficient (0.8 is the quasi-static default).
    double damping{0.8};
    double time_step_scale{0.26};
    // Opt-in static mass scaling observed from FLAC3D 7.00 Release 126.
    bool flac3d7_static_mass{false};
    bool large_strain{false};
    // FLAC3D compatibility option for imported/all-Tet4 grids.  It replaces
    // each element's volumetric strain with the volume-weighted nodal average
    // while retaining the element deviatoric strain (NMD/B-bar formulation).
    bool nodal_mixed_discretization{false};
    // Optional horizontal, frictionless rigid platen for Mandel benchmarks:
    // top nodes share one vertical displacement but remain free tangentially.
    bool rigid_top_platen_z{false};
};

struct Stage {
    std::string name;
    StageType type{StageType::Initial};
    double station_start_m{0.0};
    double station_end_m{0.0};
    double support_stiffness_n_m{0.0};
    double strength_reduction_factor{1.0};
};

// Straight, small-strain, fully grouted axial reinforcement. SI units.
struct CableConfig {
    int id{};
    bool after_initial{true};
    double young{}, area{}, yield_tension{}, yield_compression{};
    double grout_stiffness{}, grout_cohesion{}, grout_friction{}, perimeter{}, density{};
    double pretension{}, free_length{};
    std::vector<Vec3> points;
};
struct CableNode {
    Vec3 position, reference_displacement;
    std::array<std::size_t, 4> hosts{};
    std::array<double, 4> weights{};
    // Rebuilt from the physical cell after checkpoint restore, never serialized.
    std::array<std::size_t, 8> physical_hosts{};
    std::array<double, 8> physical_weights{};
    std::size_t physical_host_count{};
    std::size_t zone{};
    int host_id{};
    double bonded_length{}, length{}, axial{}, velocity{}, force{}, reference_force{}, stiffness{};
    double slip{}, plastic_slip{}, grout_force{}, mass{}, confinement{};
    bool sliding{false};
};
struct CableState {
    CableConfig config;
    Vec3 direction;
    std::vector<CableNode> nodes;
    std::vector<double> lengths, axial_force, plastic_extension, strain;
    std::vector<bool> yielded;
    std::vector<double> reference_extension;
    bool tension_control{false};
};

// Optional zero-thickness discontinuity. Empty `interfaces` preserves the
// pre-interface solver path exactly. source_group is only a geometric carrier:
// when replace_source_with_host is true its zones use the host continuum material.
struct InterfaceConfig {
    std::string name{"Interface"};
    bool enabled{true};
    std::string source_group;
    Vec3 moving_side_hint{-1.0, 0.0, 1.0};
    double normal_stiffness_pa_m{2.0e9};
    double shear_stiffness_pa_m{2.0e9};
    double cohesion_pa{0.0};
    double friction_deg{40.0};
    double tensile_strength_pa{0.0};
    bool replace_source_with_host{true};
};

struct InterfaceFacet {
    std::vector<std::size_t> master_nodes;
    std::vector<std::size_t> slave_nodes;
    Vec3 normal;
    double area_m2{0.0};
    Vec3 plastic_shear_slip;
    double normal_stress_pa{0.0};
    double shear_stress_pa{0.0};
    JointState state{JointState::Closed};
};

struct InterfaceState {
    InterfaceConfig config;
    double strength_reduction_factor{1.0};
    std::vector<InterfaceFacet> facets;
};

struct WaterConfig {
    bool enabled{false};
    double table_z_m{0.0};
    double density_kg_m3{1000.0};
    bool plane_enabled{false};
    Vec3 plane_origin_m{};
    Vec3 plane_normal{0.0, 0.0, 1.0};
    bool auto_saturated_density{false};
    // Index 0 is the host material; following entries match region_materials.
    std::vector<double> saturated_densities_kg_m3;
    // External water is a mechanical boundary load, independent of pore pressure.
    bool external_enabled{false};
    double external_surface_z_m{0.0};
    double external_density_kg_m3{1000.0};
    Vec3 external_outward_direction{0.0, 0.0, 1.0};
    double external_min_alignment{0.5};
};

struct ModelConfig {
    Scenario scenario{Scenario::Tunnel};
    int nx{8};
    int ny{12};
    int nz{8};
    double length_m{40.0};
    double width_m{30.0};
    double height_m{30.0};
    double slope_angle_deg{45.0};
    double tunnel_radius_x_m{4.0};
    double tunnel_radius_z_m{4.0};
    double initial_stress_x_pa{-2.0e6};
    double initial_stress_y_pa{-3.0e6};
    double initial_stress_z_pa{-4.0e6};
    double gravity_m_s2{9.81};
    int excavation_steps{3};
    double support_stiffness_n_m{2.0e8};
    MaterialParameters material;
    std::vector<RegionMaterial> region_materials;
    // Legacy/default behavior is false: strength reduction continues to act on
    // the continuum. Set true only for cases such as the paper Case 05 where
    // the reduction is defined exclusively on a prescribed discontinuity.
    bool strength_reduction_interfaces_only{false};
    bool strength_reduction_continuum_only{false};
    // Target encoding intentionally reuses the existing two persisted flags:
    // 00=all, 10=interfaces only, 01=continuum only, 11=region material bands only.
    // This preserves checkpoint compatibility while adding a fourth target.
    std::vector<InterfaceConfig> interfaces;
    SolverSettings solver;
    std::vector<CableConfig> cables;
    WaterConfig water;
    SteadySeepageConfig seepage;
};

struct HydraulicSnapshot {
    double time_s{};
    int iterations{};
    double relative_residual{};
    double inflow_m3_s{};
    double outflow_m3_s{};
    double net_outflow_m3_s{};
    double storage_change_m3{};
    double mechanical_storage_change_m3{};
    double relative_mass_balance_error{};
    double mass_balance_residual_m3{};
    double mass_balance_roundoff_m3{};
    double rainfall_infiltration_m3_s{};
    double rainfall_runoff_m3_s{};
    std::vector<double> head_m;
    std::vector<double> nodal_boundary_outward_flux_m3_s;
    std::vector<double> nodal_pore_pressure_pa;
    std::vector<double> nodal_saturation;
    std::vector<double> zone_head_m;
    std::vector<double> zone_saturation;
    std::vector<Vec3> zone_darcy_velocity_m_s;
    std::vector<std::array<double, 6>> zone_conductivity_m_s;
    std::vector<std::array<double, 3>> zone_joint_aperture_m;
    std::vector<std::array<double, 3>> zone_joint_normal_stress_mpa;
    std::vector<Vec3> nodal_displacement_m;
    std::vector<std::array<double, 6>> zone_stress_pa;
    std::vector<double> zone_porosity;
    std::vector<double> nodal_biot_modulus_pa;
};

struct HydraulicState {
    bool valid{false};
    bool converged{false};
    double time_s{};
    int iterations{};
    double relative_residual{};
    double inflow_m3_s{};
    double outflow_m3_s{};
    double net_outflow_m3_s{};
    double relative_mass_balance_error{};
    double mechanical_storage_change_m3{};
    double mass_balance_residual_m3{};
    double mass_balance_roundoff_m3{};
    std::string failure_reason; // Failed trial only; fields/history remain at last accepted time.
    double rainfall_infiltration_m3_s{};
    double rainfall_runoff_m3_s{};
    std::vector<double> head_m;                  // node values
    std::vector<double> nodal_boundary_outward_flux_m3_s; // outward-positive prescribed/Neumann reaction
    std::vector<double> nodal_pore_pressure_pa; // node values
    std::vector<double> nodal_saturation;        // node values, 0..1
    std::vector<double> zone_head_m;             // integration-zone values
    std::vector<double> zone_saturation;         // integration-zone values
    std::vector<Vec3> zone_darcy_velocity_m_s;   // integration-zone values
    // Effective tensors and apertures actually used by the hydraulic solve.
    std::vector<std::array<double, 6>> zone_conductivity_m_s;
    std::vector<std::array<double, 3>> zone_joint_aperture_m;
    std::vector<std::array<double, 3>> zone_joint_normal_stress_mpa;
    std::vector<HydraulicSnapshot> history;
};

struct Node {
    int id{};
    Vec3 position;
    Vec3 displacement;
    Vec3 velocity;
    Vec3 force;
    Vec3 dynamic_mass;
    Vec3 support_reference;
    std::array<bool, 3> fixed{false, false, false};
    double mass_kg{0.0};
    double force_reference{0.0};
    bool supported{false};
};

struct Zone {
    int id{};
    std::size_t material_index{0};  // 0 = host; otherwise region_materials[index - 1]
    std::string material_group;
    std::vector<std::size_t> nodes = std::vector<std::size_t>(4);
    double integration_weight{1.0};
    double parent_volume_m3{0.0};
    bool active{true};
    std::array<double, 6> stress_pa{};  // xx, yy, zz, xy, xz, yz
    std::array<double, 6> strain{};     // xx, yy, zz, gamma_xy, gamma_xz, gamma_yz
    Matrix6 tangent_pa{};
    SMRMHistory smrm_history;
    SMRMDiagnostics smrm;
    double volume_m3{0.0};
    double reinforce_requirement{0.0};
    double failure_index{0.0};
    std::uint8_t plastic_state{0};
    JointState joint_state{JointState::None};
};

struct PhysicalCell {
    int id{};
    std::vector<std::size_t> nodes;
    std::vector<std::size_t> subzones;
};

struct RockModel {
    std::string version{"1.0"};
    std::string name{"RockSolver model"};
    Scenario scenario{Scenario::Tunnel};
    ModelConfig config;
    std::vector<Node> nodes;
    std::vector<Zone> zones;
    std::vector<PhysicalCell> physical_cells;
    std::vector<Stage> stages;
    std::vector<CableState> cables;
    std::vector<InterfaceState> interfaces;
    HydraulicState hydraulic;
    bool initialize_interface_stresses{true};
    bool cable_physical_host_interpolation{true};
    bool cable_effective_confinement{true};
};

struct StageResult {
    std::string name;
    int cycles{};
    bool converged{false};
    double force_ratio{};
    double force_ratio_local{};
    double force_ratio_maximum{};
    double unbalanced_force_maximum_n{};
    double max_displacement_m{};
    std::size_t active_zones{};
};

struct RockResult {
    std::string solver_version{kVersion};
    bool converged{false};
    int total_cycles{};
    double final_force_ratio{};
    double final_force_ratio_local{};
    double final_force_ratio_maximum{};
    double final_unbalanced_force_maximum_n{};
    double max_displacement_m{};
    double max_reinforce_requirement{};
    std::size_t active_zones{};
    std::vector<StageResult> stages;
    std::vector<std::string> warnings;
};

struct SolveProgress {
    std::string stage_name;
    int stage_index{};
    int stage_count{};
    int cycle{};
    int max_cycles{};
    double force_ratio{};
    double force_ratio_local{};
    double force_ratio_maximum{};
    double unbalanced_force_maximum_n{};
    bool stage_converged{false};
};

using ProgressCallback = std::function<void(const SolveProgress&)>;

struct ResourceEstimate {
    std::size_t nodes{};
    std::size_t zones{};
    std::size_t degrees_of_freedom{};
    double memory_mb{};
    std::string recommended_backend;
};

RockModel generate_model(const ModelConfig& config);
RockModel import_abaqus_c3d4(const std::string& file_path,
                             const ModelConfig& config);
RockModel import_abaqus_c3d4_text(const std::string& inp_text,
                                  const ModelConfig& config);
RockResult solve(RockModel& model, const ProgressCallback& progress = {});
RockResult reduction_trial(RockModel&, double factor, const ProgressCallback& = {});
ResourceEstimate estimate_resources(const ModelConfig& config);
ResourceEstimate estimate_resources(const RockModel& model);
std::string result_to_json(const RockModel& model, const RockResult& result);
const char* scenario_name(Scenario scenario);
const char* material_name(MaterialType material);
const char* joint_state_name(JointState state);

// Versioned, lossless state storage, including integration-point history.
std::string save_checkpoint(const RockModel&, const RockResult&);
void restore_checkpoint(const std::string&, RockModel&, RockResult&);
RockResult solve_baseline(RockModel&, const ProgressCallback& = {});
RockResult continue_slope(RockModel&, const RockResult& baseline,
                         const ModelConfig& branch, bool allow_unconverged,
                         const ProgressCallback& = {});
RockResult continue_water_stage(RockModel&, const RockResult& previous,
                                const WaterConfig& water,
                                const SolverSettings& solver,
                                const ProgressCallback& = {});
RockResult continue_hydraulic_stage(RockModel&, const RockResult& previous,
                                    const HydraulicState& hydraulic,
                                    const SolverSettings& solver,
                                    const ProgressCallback& = {});
HydraulicState solve_steady_seepage(RockModel&);
RockResult solve_steady_seepage_mechanics(RockModel&, const ProgressCallback& = {});
HydraulicState solve_transient_seepage(RockModel&);
HydraulicState solve_unsaturated_seepage(RockModel&);
RockResult solve_transient_seepage_mechanics(RockModel&, const ProgressCallback& = {});

}  // namespace rocksolver
