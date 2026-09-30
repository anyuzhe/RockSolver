#include "rocksolver/rocksolver.hpp"

#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>

int main() {
    const auto remove_098_checkpoint_extension=[](std::string& checkpoint,
                                                   std::size_t region_count) {
        std::string suffix="0 0 "+std::to_string(region_count)+" ";
        for(std::size_t i=0;i<region_count;++i)suffix+="0 ";
        // Newer checkpoints append the optional rigid-platen flag after the
        // original 0.9.8 compatibility extension.
        if(checkpoint.ends_with(suffix+"0 "))checkpoint.erase(checkpoint.size()-2);
        if(!checkpoint.ends_with(suffix))
            throw std::runtime_error("unexpected 0.9.8 checkpoint extension");
        checkpoint.erase(checkpoint.size()-suffix.size());
    };
    // A frictionless rigid top platen ties only vertical displacement. An
    // uneven bottom makes the untied top settle differently across x.
    const std::string platen_mesh=R"INP(*NODE
1,0,0,0
2,1,0,0
3,2,0,0.5
4,0,1,0
5,1,1,0
6,2,1,0.5
7,0,0,2
8,1,0,2
9,2,0,2
10,0,1,2
11,1,1,2
12,2,1,2
*ELEMENT,TYPE=C3D8,ELSET=HOST
1,1,2,5,4,7,8,11,10
2,2,3,6,5,8,9,12,11
)INP";
    rocksolver::ModelConfig platen_config;
    platen_config.scenario=rocksolver::Scenario::Slope;
    platen_config.gravity_m_s2=9.81;
    platen_config.initial_stress_x_pa=0.0;
    platen_config.initial_stress_y_pa=0.0;
    platen_config.initial_stress_z_pa=0.0;
    platen_config.material.type=rocksolver::MaterialType::Elastic;
    platen_config.material.elastic_modulus_pa=1.0e8;
    platen_config.material.poisson_ratio=0.25;
    platen_config.solver.max_cycles=20000;
    platen_config.solver.convergence_tolerance=1.0e-8;
    platen_config.solver.nodal_mixed_discretization=true;
    const auto platen_spread=[&](bool rigid){
        platen_config.solver.rigid_top_platen_z=rigid;
        auto model=rocksolver::import_abaqus_c3d4_text(platen_mesh,platen_config);
        for(auto& node:model.nodes)if(node.id<=6)node.fixed[2]=true;
        const auto result=rocksolver::solve_baseline(model);
        assert(result.converged);
        double low=1.0,high=-1.0;
        for(const auto& node:model.nodes)if(node.position.z==2.0){
            low=std::min(low,node.displacement.z);
            high=std::max(high,node.displacement.z);
        }
        if(rigid){
            const auto checkpoint=rocksolver::save_checkpoint(model,result);
            rocksolver::RockModel restored;
            rocksolver::RockResult restored_result;
            rocksolver::restore_checkpoint(checkpoint,restored,restored_result);
            assert(restored.config.solver.rigid_top_platen_z);
        }
        return high-low;
    };
    assert(platen_spread(false)>8.0e-5);
    assert(platen_spread(true)<1.0e-11);
    const auto direction_zero = rocksolver::joint_direction_cosines(0.0, 90.0);
    const auto direction_ninety =
        rocksolver::joint_direction_cosines(90.0, 90.0);
    assert(std::abs(direction_zero[0] + 1.0) < 1.0e-12);
    assert(std::abs(direction_zero[1]) < 1.0e-12);
    assert(std::abs(direction_zero[2]) < 1.0e-12);
    assert(std::abs(direction_ninety[0]) < 1.0e-12);
    assert(std::abs(direction_ninety[1] - 1.0) < 1.0e-12);
    assert(std::abs(direction_ninety[2]) < 1.0e-12);

    rocksolver::MaterialParameters elastic_material;
    elastic_material.type = rocksolver::MaterialType::Elastic;
    std::array<double, 6> elastic_stress{};
    rocksolver::Matrix6 elastic_tangent{};
    const std::array<double, 6> test_increment{-1.0e-4, 2.0e-5, 0.0,
                                               1.0e-5, 0.0, 0.0};
    assert(rocksolver::update_elastic_increment(
        elastic_material, test_increment, elastic_stress, elastic_tangent));
    assert(elastic_stress[0] < 0.0);
    assert(elastic_tangent[0][0] > 0.0);
    assert(std::abs(elastic_tangent[3][4]) < 1.0e-12);

    // The parallel Mohr-Coulomb path must perform the same elastic-trial ->
    // shear/tension correction sequence as the supplied FLAC ModelMohr.
    rocksolver::MaterialParameters mc_material = elastic_material;
    mc_material.type = rocksolver::MaterialType::MohrCoulomb;
    mc_material.elastic_modulus_pa = 100.0e6;
    mc_material.poisson_ratio = 0.25;
    mc_material.cohesion_pa = 0.10e6;
    mc_material.friction_deg = 30.0;
    mc_material.dilation_deg = 0.0;
    mc_material.tensile_strength_pa = 0.02e6;
    std::array<double, 6> mc_stress{};
    rocksolver::Matrix6 mc_tangent{};
    std::uint8_t mc_state = 0;
    const std::array<double, 6> mc_compression{-0.02, 0, 0, 0, 0, 0};
    assert(rocksolver::update_mohr_coulomb_increment(
        mc_material, mc_compression, mc_stress, mc_state, mc_tangent));
    // Uniaxial STRAIN at nu=0.25 gives (-2.4,-0.8,-0.8) MPa.
    // For phi=30 deg, Nphi=3: f=3*sigma3-sigma1-2*c*sqrt(3)<0.
    // Lateral confinement keeps this trial elastic; it is not a shear test.
    assert(mc_state == 0);
    assert(std::abs(mc_stress[0] + 2.4e6) < 1e-8);
    assert(std::abs(mc_stress[1] + 0.8e6) < 1e-8);
    assert(std::abs(mc_stress[2] + 0.8e6) < 1e-8);
    // Allow Poisson expansion to obtain an unconfined elastic trial of
    // (-2,0,0) MPa, then verify return to the MC shear surface.
    mc_stress = {};
    const std::array<double, 6> mc_unconfined{-0.02, 0.005, 0.005, 0, 0, 0};
    assert(rocksolver::update_mohr_coulomb_increment(
        mc_material, mc_unconfined, mc_stress, mc_state, mc_tangent));
    assert((mc_state & 0x01U) != 0U);
    assert(std::abs(3.0*mc_stress[2]-mc_stress[0]-2.0e5*std::sqrt(3.0)) < 1e-8);
    assert(std::all_of(mc_stress.begin(), mc_stress.end(), [](double x){ return std::isfinite(x); }));
    std::array<double, 6> mc_tension_stress{};
    std::uint8_t mc_tension_state = 0;
    const std::array<double, 6> mc_tension{0.01, 0, 0, 0, 0, 0};
    assert(rocksolver::update_mohr_coulomb_increment(
        mc_material, mc_tension, mc_tension_stress, mc_tension_state, mc_tangent));
    assert((mc_tension_state & 0x02U) != 0U);

    rocksolver::MaterialParameters material_point = elastic_material;
    material_point.type = rocksolver::MaterialType::SMRM;
    material_point.cohesion_pa = 5.8e6;
    material_point.friction_deg = 60.5;
    material_point.tensile_strength_pa = 10.36e6;
    material_point.joints.push_back(rocksolver::JointSet{
        .name = "J1",
        .dip_direction_deg = 90.0,
        .dip_deg = 90.0,
        .normal_density_1_m = 0.76,
        .mean_radius_m = 1.21,
        .cohesion_pa = 10.0,
        .friction_deg = 36.55,
    });
    std::array<double, 6> smrm_stress{};
    rocksolver::SMRMHistory smrm_history;
    rocksolver::SMRMDiagnostics smrm_diagnostics;
    rocksolver::Matrix6 smrm_tangent{};
    assert(rocksolver::update_smrm_increment(
        material_point, 1.0, test_increment, smrm_stress, smrm_history,
        smrm_diagnostics, smrm_tangent));
    assert(smrm_history.step == 1);
    assert(smrm_tangent[0][0] < elastic_tangent[0][0]);
    assert(smrm_diagnostics.joint_k[0] >= 0.0 &&
           smrm_diagnostics.joint_k[0] <= 1.0);
    assert(smrm_diagnostics.rock_mass_em_x_gpa > 0.0);
    assert(std::isfinite(smrm_diagnostics.failure_probability));

    // Hydrostatic pore pressure has no tangential traction on any joint.
    auto hydrostatic_material = material_point;
    hydrostatic_material.joints[0].dip_direction_deg = 291.0;
    hydrostatic_material.joints[0].dip_deg = 37.0;
    hydrostatic_material.joints.push_back(rocksolver::JointSet{
        .name = "J2", .dip_direction_deg = 37.0, .dip_deg = 74.0,
        .normal_density_1_m = 0.75, .mean_radius_m = 0.98,
        .cohesion_pa = 1.0e-5, .friction_deg = 30.0,
    });
    hydrostatic_material.joints.push_back(rocksolver::JointSet{
        .name = "J3", .dip_direction_deg = 162.0, .dip_deg = 68.0,
        .normal_density_1_m = 0.23, .mean_radius_m = 0.49,
        .cohesion_pa = 1.0e-5, .friction_deg = 30.0,
    });
    std::array<double, 6> hydrostatic_effective_stress{
        437725.876, 437725.876, 437725.876, 0.0, 0.0, 0.0};
    rocksolver::SMRMHistory hydrostatic_history;
    rocksolver::SMRMDiagnostics hydrostatic_diagnostics;
    rocksolver::Matrix6 hydrostatic_tangent{};
    const std::array<double, 6> zero_strain{};
    assert(rocksolver::update_smrm_increment(
        hydrostatic_material, 1.0, zero_strain, hydrostatic_effective_stress,
        hydrostatic_history, hydrostatic_diagnostics, hydrostatic_tangent));
    for (std::size_t i = 0; i < hydrostatic_material.joints.size(); ++i) {
        if (hydrostatic_diagnostics.joint_shear_stress_mpa[i] >= 1.0e-12 ||
            hydrostatic_diagnostics.joint_h[i] != 0.0)
            throw std::runtime_error("hydrostatic stress created spurious joint shear");
    }

    // FLAC uses the old stress state to build each step tangent, so splitting
    // a load into cycles is intentionally path dependent.
    std::array<double, 6> split_stress{};
    rocksolver::SMRMHistory split_history;
    rocksolver::SMRMDiagnostics split_diagnostics;
    rocksolver::Matrix6 split_tangent{};
    const std::array<double, 6> compression_step{-1.0e-4, 0.0, 0.0,
                                                  0.0, 0.0, 0.0};
    assert(rocksolver::update_smrm_increment(
        material_point, 1.0, compression_step, split_stress, split_history,
        split_diagnostics, split_tangent));
    assert(rocksolver::update_smrm_increment(
        material_point, 1.0, compression_step, split_stress, split_history,
        split_diagnostics, split_tangent));
    std::array<double, 6> one_step_stress{};
    rocksolver::SMRMHistory one_step_history;
    rocksolver::SMRMDiagnostics one_step_diagnostics;
    rocksolver::Matrix6 one_step_tangent{};
    auto double_step = compression_step;
    double_step[0] *= 2.0;
    assert(rocksolver::update_smrm_increment(
        material_point, 1.0, double_step, one_step_stress, one_step_history,
        one_step_diagnostics, one_step_tangent));
    assert(std::abs(split_stress[0] - one_step_stress[0]) > 1.0);

    rocksolver::ModelConfig slope;
    slope.scenario = rocksolver::Scenario::Slope;
    slope.nx = 2;
    slope.ny = 2;
    slope.nz = 2;
    slope.initial_stress_x_pa = 0.0;
    slope.initial_stress_y_pa = 0.0;
    slope.initial_stress_z_pa = 0.0;
    auto slope_model = rocksolver::generate_model(slope);
    assert(slope_model.nodes.size() == 27);
    assert(slope_model.zones.size() == 48);
    for (const auto& zone : slope_model.zones) assert(zone.volume_m3 > 0.0);

    // FLAC does not apply local damping on a component's reversing step.
    auto reversal_config = slope;
    reversal_config.material.type = rocksolver::MaterialType::Elastic;
    reversal_config.solver.max_cycles = 1;
    reversal_config.solver.damping = 0.8;
    auto reversal_model = rocksolver::generate_model(reversal_config);
    for (auto& node : reversal_model.nodes) node.fixed = {true, true, true};
    auto& reversing_node = reversal_model.nodes[13];
    reversing_node.fixed = {true, true, false};
    constexpr double initial_upward_velocity = 1.0e-9;
    reversing_node.velocity.z = initial_upward_velocity;
    rocksolver::solve(reversal_model);
    const auto& advanced_node = reversal_model.nodes[13];
    const double undamped_step =
        -advanced_node.mass_kg * reversal_config.gravity_m_s2 /
        advanced_node.dynamic_mass.z * reversal_config.solver.time_step_scale;
    assert(initial_upward_velocity + undamped_step < 0.0);
    assert(std::abs(advanced_node.velocity.z -
                    (initial_upward_velocity + undamped_step)) < 1.0e-12);

    // FLAC feeds effective (not total) stress to SMRM. A fixed hydrostatic
    // pore pressure must weaken the joint without changing total stress.
    rocksolver::ModelConfig hydro = slope;
    hydro.nx = hydro.ny = hydro.nz = 1;
    hydro.initial_stress_x_pa = hydro.initial_stress_y_pa = hydro.initial_stress_z_pa = -2.0e6;
    hydro.material.type = rocksolver::MaterialType::SMRM;
    hydro.material.joints.push_back(rocksolver::JointSet{
        .name = "wet joint", .dip_direction_deg = 0.0, .dip_deg = 45.0,
        .normal_density_1_m = 0.1, .mean_radius_m = 1.0,
        .cohesion_pa = 1.0e5, .friction_deg = 30.0,
    });
    // solve_stages only checks convergence after cycle 25, even if fully fixed.
    hydro.solver.max_cycles = 30;
    hydro.solver.convergence_tolerance = 1.0e-10;
    auto dry_model = rocksolver::generate_model(hydro);
    for (auto& node : dry_model.nodes) node.fixed = {true, true, true};
    const auto dry_result = rocksolver::solve(dry_model);
    hydro.water = {true, 200.0, 1000.0};
    auto wet_model = rocksolver::generate_model(hydro);
    for (auto& node : wet_model.nodes) node.fixed = {true, true, true};
    const auto wet_result = rocksolver::solve(wet_model);
    assert(dry_result.converged && wet_result.converged);
    assert(std::abs(dry_model.zones[0].stress_pa[0] - wet_model.zones[0].stress_pa[0]) < 1.0e-6);
    assert(wet_model.zones[0].smrm.joint_normal_stress_mpa[0] <
           dry_model.zones[0].smrm.joint_normal_stress_mpa[0]);
    assert(rocksolver::result_to_json(wet_model, wet_result).find("\"porePressureMpa\"") != std::string::npos);
    rocksolver::RockModel restored_wet;
    rocksolver::RockResult restored_wet_result;
    rocksolver::restore_checkpoint(rocksolver::save_checkpoint(wet_model, wet_result),
                                   restored_wet, restored_wet_result);
    assert(restored_wet.config.water.enabled);
    assert(restored_wet.config.water.table_z_m == hydro.water.table_z_m);

    const auto estimate = rocksolver::estimate_resources(slope);
    assert(estimate.nodes == 27);
    assert(estimate.zones == 48);
    assert(estimate.recommended_backend == "wasm");

    const auto inp_path =
        std::filesystem::temp_directory_path() / "rocksolver-c3d4-import-test.inp";
    const std::string inp_text =
        "*HEADING\nImported test\n*NODE\n"
        "10,0,0,0\n20,1,0,0\n30,0,1,0\n40,0,0,1\n"
        "50,0.25,0.25,0\n"
        "*ELEMENT, TYPE=C3D4, ELSET=ROCK\n"
        "100,10,20,30,40\n101,50,20,30,40\n";
    {
        std::ofstream inp(inp_path);
        inp << inp_text;
    }
    rocksolver::ModelConfig imported_config;
    imported_config.scenario = rocksolver::Scenario::Tunnel;
    imported_config.excavation_steps = 0;
    auto imported = rocksolver::import_abaqus_c3d4(inp_path.string(), imported_config);
    assert(imported.nodes.size() == 5);
    assert(imported.zones.size() == 2);
    assert(imported.nodes[0].id == 10);
    assert(imported.zones[0].id == 100);
    assert(imported.zones[0].nodes[3] == 3);
    assert(imported.nodes[4].fixed[0]);
    assert(imported.nodes[4].fixed[1]);
    assert(imported.nodes[4].fixed[2]);
    assert(imported.stages.size() == 1);
    assert(rocksolver::estimate_resources(imported).zones == 2);
    const auto imported_from_text =
        rocksolver::import_abaqus_c3d4_text(inp_text, imported_config);
    assert(imported_from_text.nodes.size() == imported.nodes.size());
    assert(imported_from_text.zones.size() == imported.zones.size());

    imported_config.scenario = rocksolver::Scenario::Slope;
    auto imported_slope =
        rocksolver::import_abaqus_c3d4(inp_path.string(), imported_config);
    std::filesystem::remove(inp_path);
    assert(imported_slope.nodes[4].fixed[0]);
    assert(imported_slope.nodes[4].fixed[1]);
    assert(imported_slope.nodes[4].fixed[2]);

    // FLAC3D C3D8 uses two 5-Tet overlays with a common cell-mean stress.
    // A skew brick checks the physical-cell response; nonuniform overlay
    // stresses are covered by the mixed-cell FLAC regression in Stage 5.
    const std::string skew_brick_inp = R"INP(*NODE
1,0,0,0
2,4,0,0
3,4,1,0
4,0,1,0
5,0,0,1
6,3,0,1
7,3,1,1
8,0,1,1
*ELEMENT,TYPE=C3D8,ELSET=HOST
1,1,2,3,4,5,6,7,8
)INP";
    rocksolver::ModelConfig skew_config;
    skew_config.scenario = rocksolver::Scenario::Slope;
    skew_config.initial_stress_x_pa = 0.0;
    skew_config.initial_stress_y_pa = 0.0;
    skew_config.initial_stress_z_pa = 0.0;
    skew_config.gravity_m_s2 = 9.81;
    skew_config.material.type = rocksolver::MaterialType::Elastic;
    skew_config.material.density_kg_m3 = 2000.0;
    skew_config.material.elastic_modulus_pa = 0.1e9;
    skew_config.material.poisson_ratio = 0.3;
    skew_config.solver.max_cycles = 30000;
    skew_config.solver.convergence_tolerance = 1.0e-7;
    skew_config.solver.damping = 0.8;
    skew_config.solver.time_step_scale = 0.26;
    auto skew_model =
        rocksolver::import_abaqus_c3d4_text(skew_brick_inp, skew_config);
    const auto skew_result = rocksolver::solve(skew_model);
    assert(skew_result.converged);
    // FLAC3D 7.00 R126 reference: 9.710275454e-5 m.
    assert(std::abs(skew_result.max_displacement_m - 9.710275454e-5) < 5.0e-9);
    assert(skew_model.physical_cells.size() == 1);
    const auto& skew_subzones = skew_model.physical_cells.front().subzones;
    const auto skew_mean=[](const rocksolver::Zone& zone){
        return (zone.stress_pa[0]+zone.stress_pa[1]+zone.stress_pa[2])/3.0;
    };
    for (const auto zi : skew_subzones)
        assert(std::abs(skew_mean(skew_model.zones[zi])-
                        skew_mean(skew_model.zones[skew_subzones.front()]))<1.0e-7);

    // FLAC C3D5 subzones share the cell's mean normal stress after their
    // individual constitutive updates; this must also hold for nonlinear SMRM.
    const std::string pyramid_inp = R"INP(*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0.5,0.5,1
*ELEMENT,TYPE=C3D5,ELSET=HOST
1,1,2,3,4,5
)INP";
    auto pyramid_config = skew_config;
    pyramid_config.material = hydrostatic_material;
    pyramid_config.water.enabled = true;
    pyramid_config.water.table_z_m = 2.0;
    pyramid_config.solver.max_cycles = 2;
    pyramid_config.solver.convergence_tolerance = 1.0e-12;
    auto pyramid_model = rocksolver::import_abaqus_c3d4_text(pyramid_inp, pyramid_config);
    rocksolver::solve(pyramid_model);
    const auto& pyramid_subzones = pyramid_model.physical_cells.front().subzones;
    const auto& first_pyramid = pyramid_model.zones[pyramid_subzones.front()];
    const double first_mean = (first_pyramid.stress_pa[0] +
        first_pyramid.stress_pa[1] + first_pyramid.stress_pa[2]) / 3.0;
    for (const auto zi : pyramid_subzones) {
        const auto& stress = pyramid_model.zones[zi].stress_pa;
        if (std::abs((stress[0]+stress[1]+stress[2])/3.0 - first_mean) >= 1.0e-7)
            throw std::runtime_error("C3D5 subzone mean stresses differ");
    }

    // A Cable at the C3D8 midplane must transfer load to all eight physical
    // corners, not just the first integration tetrahedron containing it.
    auto cable_brick_config=skew_config;
    cable_brick_config.gravity_m_s2=0.0;
    rocksolver::CableConfig brick_cable;
    brick_cable.id=1;
    brick_cable.after_initial=false;
    brick_cable.young=200e9;
    brick_cable.area=0.001;
    brick_cable.yield_tension=500e3;
    brick_cable.yield_compression=500e3;
    brick_cable.grout_stiffness=2e6;
    brick_cable.grout_cohesion=100e3;
    brick_cable.grout_friction=30;
    brick_cable.perimeter=0.4;
    brick_cable.points={{0.5,0.5,0.5},{1.5,0.5,0.5}};
    cable_brick_config.cables.push_back(brick_cable);
    auto cable_brick=rocksolver::import_abaqus_c3d4_text(skew_brick_inp,cable_brick_config);
    auto cable_brick_result=rocksolver::solve(cable_brick);
    assert(cable_brick_result.converged);
    const auto& cable_host=cable_brick.cables.front().nodes.front();
    assert(cable_host.physical_host_count==8);
    for(int j=0;j<8;++j)assert(cable_host.physical_weights[j]>0);
    assert(std::abs(cable_host.physical_weights[0]-cable_host.physical_weights[3])<1e-12);
    assert(std::abs(cable_host.physical_weights[1]-cable_host.physical_weights[2])<1e-12);
    rocksolver::RockModel restored_cable_brick;
    rocksolver::RockResult restored_cable_brick_result;
    rocksolver::restore_checkpoint(rocksolver::save_checkpoint(cable_brick,cable_brick_result),
                                   restored_cable_brick,restored_cable_brick_result);
    assert(restored_cable_brick.cables.front().nodes.front().physical_host_count==8);
    assert(restored_cable_brick.cable_effective_confinement);
    assert(std::abs(restored_cable_brick.cables.front().nodes.front().physical_weights[0]-
                    cable_host.physical_weights[0])<1e-15);
    auto prior_face_host_checkpoint=rocksolver::save_checkpoint(cable_brick,cable_brick_result);
    remove_098_checkpoint_extension(prior_face_host_checkpoint,0);
    prior_face_host_checkpoint.replace(0,std::string("ROCKSTATE-1 ").size()+std::string(rocksolver::kVersion).size(),
        "ROCKSTATE-1 0.9.5-cable-effective-confinement");
    rocksolver::restore_checkpoint(prior_face_host_checkpoint,restored_cable_brick,restored_cable_brick_result);
    if(!restored_cable_brick.cable_physical_host_interpolation || !restored_cable_brick.cable_effective_confinement)
        throw std::runtime_error("0.9.5 cable checkpoint compatibility lost");
    auto old_cable_checkpoint=rocksolver::save_checkpoint(cable_brick,cable_brick_result);
    remove_098_checkpoint_extension(old_cable_checkpoint,0);
    const std::string cable_header="ROCKSTATE-1 "+std::string(rocksolver::kVersion)+" ";
    const std::string cable_flags="1 1 1 \""+std::string(rocksolver::kVersion)+"\"";
    const auto cable_flags_at=old_cable_checkpoint.rfind(cable_flags);
    assert(cable_flags_at!=std::string::npos);
    old_cable_checkpoint.erase(cable_flags_at+4,2);
    old_cable_checkpoint.replace(0,cable_header.size(),
        "ROCKSTATE-1 0.9.4-cable-cell-host-weights ");
    rocksolver::restore_checkpoint(old_cable_checkpoint,restored_cable_brick,restored_cable_brick_result);
    assert(restored_cable_brick.cable_physical_host_interpolation);
    assert(!restored_cable_brick.cable_effective_confinement);
    assert(restored_cable_brick.cables.front().nodes.front().physical_host_count==8);
    const auto old_cable_flags_at=old_cable_checkpoint.rfind(
        "1 1 \""+std::string(rocksolver::kVersion)+"\"");
    assert(old_cable_flags_at!=std::string::npos);
    old_cable_checkpoint.erase(old_cable_flags_at+2,2);
    old_cable_checkpoint.replace(0,std::string("ROCKSTATE-1 0.9.4-cable-cell-host-weights ").size(),
        "ROCKSTATE-1 0.9.3-interface-initial-stress ");
    rocksolver::restore_checkpoint(old_cable_checkpoint,restored_cable_brick,restored_cable_brick_result);
    assert(!restored_cable_brick.cable_physical_host_interpolation);
    assert(!restored_cable_brick.cable_effective_confinement);
    assert(restored_cable_brick.cables.front().nodes.front().physical_host_count==0);

    // A cable node exactly on a shared C3D8 face follows FLAC3D's
    // last-imported-zone choice, independently of zone ID.
    const std::string shared_face_nodes = R"INP(*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
9,2,0,0
10,2,1,0
11,2,0,1
12,2,1,1
)INP";
    const std::string left_cell = "*ELEMENT,TYPE=C3D8,ELSET=LEFT\n3,1,2,3,4,5,6,7,8\n";
    const std::string right_cell = "*ELEMENT,TYPE=C3D8,ELSET=RIGHT\n1,2,9,10,3,6,11,12,7\n";
    auto shared_face_config = cable_brick_config;
    shared_face_config.cables.front().points = {{0.5,0.5,0.5},{1.0,0.5,0.5},{1.5,0.5,0.5}};
    auto left_then_right = rocksolver::import_abaqus_c3d4_text(
        shared_face_nodes+left_cell+right_cell, shared_face_config);
    const auto left_then_right_result=rocksolver::solve(left_then_right);
    if(!left_then_right_result.converged || left_then_right.cables.front().nodes[1].host_id!=1)
        throw std::runtime_error("shared-face cable host: expected last imported right cell");
    auto right_interior_config=shared_face_config;
    right_interior_config.cables.front().points[1].x=1.000001;
    auto right_interior=rocksolver::import_abaqus_c3d4_text(
        shared_face_nodes+left_cell+right_cell, right_interior_config);
    if(!rocksolver::solve(right_interior).converged ||
       left_then_right.cables.front().nodes[1].zone!=right_interior.cables.front().nodes[1].zone)
        throw std::runtime_error("shared-face cable host: incorrect right-cell overlay subzone");
    auto right_then_left = rocksolver::import_abaqus_c3d4_text(
        shared_face_nodes+right_cell+left_cell, shared_face_config);
    const auto right_then_left_result=rocksolver::solve(right_then_left);
    if(!right_then_left_result.converged || right_then_left.cables.front().nodes[1].host_id!=3)
        throw std::runtime_error("shared-face cable host: expected last imported left cell");
    auto left_interior_config=shared_face_config;
    left_interior_config.cables.front().points[1].x=0.999999;
    auto left_interior=rocksolver::import_abaqus_c3d4_text(
        shared_face_nodes+right_cell+left_cell, left_interior_config);
    if(!rocksolver::solve(left_interior).converged ||
       right_then_left.cables.front().nodes[1].zone!=left_interior.cables.front().nodes[1].zone)
        throw std::runtime_error("shared-face cable host: incorrect left-cell overlay subzone");

    const std::string wedge_inp = R"INP(*NODE
1,0,0,0
2,2,0,0
3,4,0,0
4,0,0,1
5,2,0,1
6,3,0,1
7,0,0,2
8,2,0,2
9,0,1,0
10,2,1,0
11,4,1,0
12,0,1,1
13,2,1,1
14,3,1,1
15,0,1,2
16,2,1,2
*ELEMENT,TYPE=C3D6,ELSET=HOST
101,1,5,2,9,13,10
102,1,4,5,9,12,13
201,4,8,5,12,16,13
202,4,7,8,12,15,16
301,2,6,3,10,14,11
302,2,5,6,10,13,14
303,5,8,6,13,16,14
)INP";
    auto wedge_config = skew_config;
    auto wedge_model =
        rocksolver::import_abaqus_c3d4_text(wedge_inp, wedge_config);
    const auto wedge_result = rocksolver::solve(wedge_model);
    assert(wedge_result.converged);
    // FLAC3D 7.00 R126 reference for the same seven-wedge slope.
    assert(std::abs(wedge_result.max_displacement_m - 2.894907718e-4) < 5.0e-9);
    for (const auto& cell : wedge_model.physical_cells) {
        if (cell.nodes.size() != 6) continue;
        for (const auto zi : cell.subzones)
            assert(std::abs(skew_mean(wedge_model.zones[zi])-
                            skew_mean(wedge_model.zones[cell.subzones.front()]))<1.0e-7);
    }

    rocksolver::ModelConfig tunnel;
    tunnel.scenario = rocksolver::Scenario::Tunnel;
    tunnel.nx = 3;
    tunnel.ny = 3;
    tunnel.nz = 3;
    tunnel.length_m = 12.0;
    tunnel.width_m = 12.0;
    tunnel.height_m = 12.0;
    tunnel.tunnel_radius_x_m = 2.5;
    tunnel.tunnel_radius_z_m = 2.5;
    tunnel.excavation_steps = 2;
    tunnel.support_stiffness_n_m = 1.0e7;
    tunnel.solver.max_cycles = 120;
    tunnel.solver.convergence_tolerance = 2.0;
    tunnel.material.type = rocksolver::MaterialType::SMRM;
    tunnel.material.ucs_pa = 1.0e6;
    tunnel.material.joints.push_back(rocksolver::JointSet{
        .name = "J1",
        .dip_direction_deg = 30.0,
        .dip_deg = 55.0,
        .normal_density_1_m = 1.0,
        .mean_radius_m = 1.5,
        .cohesion_pa = 0.2e6,
        .friction_deg = 25.0,
        .tensile_strength_pa = 0.05e6,
    });
    auto tunnel_model = rocksolver::generate_model(tunnel);
    const auto before = tunnel_model.zones.size();
    const auto result = rocksolver::solve(tunnel_model);
    assert(result.stages.size() == 5);
    assert(result.active_zones < before);
    assert(result.total_cycles > 0);
    const std::string json = rocksolver::result_to_json(tunnel_model, result);
    assert(json.find("\"tetrahedra\"") != std::string::npos);
    assert(json.find("\"reinforceRequirement\"") != std::string::npos);
    assert(json.find("\"rockMassEmGpa\"") != std::string::npos);
    assert(json.find("\"jointShearStrengthMpa\"") != std::string::npos);
    assert(json.find("\"failureProbability\"") != std::string::npos);
    assert(json.find("\"stages\"") != std::string::npos);

    rocksolver::ModelConfig zero_tension_joint = slope;
    zero_tension_joint.material.type = rocksolver::MaterialType::SMRM;
    zero_tension_joint.material.tensile_strength_pa = 10.0e6;
    zero_tension_joint.material.joints.push_back(rocksolver::JointSet{
        .name = "zero tension",
        .dip_direction_deg = 20.0,
        .dip_deg = 60.0,
        .normal_density_1_m = 0.8,
        .mean_radius_m = 1.0,
        .cohesion_pa = 0.0,
        .friction_deg = 30.0,
        .tensile_strength_pa = 0.0,
    });
    zero_tension_joint.initial_stress_x_pa = 1.0e6;
    zero_tension_joint.solver.max_cycles = 30;
    zero_tension_joint.solver.convergence_tolerance = 10.0;
    auto zero_tension_model = rocksolver::generate_model(zero_tension_joint);
    const auto zero_tension_result = rocksolver::solve(zero_tension_model);
    assert(std::isfinite(zero_tension_result.max_reinforce_requirement));
    const std::string zero_tension_json =
        rocksolver::result_to_json(zero_tension_model, zero_tension_result);
    assert(zero_tension_json.find(":inf") == std::string::npos);

    // Optional interface is incremental: legacy import keeps its original node
    // topology; enabling it splits only the requested internal face.
    const std::string interface_inp = R"INP(*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
9,2,0,0
10,2,1,0
11,2,0,1
12,2,1,1
13,3,0,0
14,3,1,0
15,3,0,1
16,3,1,1
*ELEMENT,TYPE=C3D8,ELSET=HOST_LEFT
1,1,2,3,4,5,6,7,8
*ELEMENT,TYPE=C3D8,ELSET=SLOPE_BAND_001
2,2,9,10,3,6,11,12,7
*ELEMENT,TYPE=C3D8,ELSET=HOST_RIGHT
3,9,13,14,10,11,15,16,12
)INP";
    rocksolver::ModelConfig legacy_interface_config;
    legacy_interface_config.scenario=rocksolver::Scenario::Slope;
    legacy_interface_config.material.type=rocksolver::MaterialType::MohrCoulomb;
    legacy_interface_config.region_materials.push_back({"SLOPE_BAND_001",legacy_interface_config.material});
    auto legacy_interface_model=rocksolver::import_abaqus_c3d4_text(interface_inp,legacy_interface_config);
    assert(legacy_interface_model.nodes.size()==16 && legacy_interface_model.interfaces.empty());
    // A sloped static water plane sets gridpoint pressure, while each physical
    // zone uses its own dry/saturated density (including the material band).
    auto sloped_water_config=legacy_interface_config;
    sloped_water_config.material.type=rocksolver::MaterialType::Elastic;
    sloped_water_config.material.density_kg_m3=2000.0;
    sloped_water_config.region_materials[0].material.type=rocksolver::MaterialType::Elastic;
    sloped_water_config.region_materials[0].material.density_kg_m3=2100.0;
    sloped_water_config.gravity_m_s2=10.0;
    sloped_water_config.water.enabled=true;
    sloped_water_config.water.density_kg_m3=1000.0;
    sloped_water_config.water.plane_enabled=true;
    sloped_water_config.water.plane_origin_m={0,0,0.25};
    sloped_water_config.water.plane_normal={-0.25,0,1};
    sloped_water_config.water.auto_saturated_density=true;
    sloped_water_config.water.saturated_densities_kg_m3={2300.0,2400.0};
    sloped_water_config.solver.max_cycles=30;
    sloped_water_config.solver.convergence_tolerance=10.0;
    auto sloped_water_model=rocksolver::import_abaqus_c3d4_text(interface_inp,sloped_water_config);
    for(auto& node:sloped_water_model.nodes)node.fixed={true,true,true};
    const auto sloped_water_result=rocksolver::solve(sloped_water_model);
    double total_mass=0.0;
    for(const auto& node:sloped_water_model.nodes)total_mass+=node.mass_kg;
    assert(std::abs(total_mass-6700.0)<1.0e-5);
    const auto sloped_water_json=rocksolver::result_to_json(sloped_water_model,sloped_water_result);
    assert(sloped_water_json.find("\"densityKgM3\":[2000,2400,2300]")!=std::string::npos);
    rocksolver::RockModel restored_sloped_water_model;
    rocksolver::RockResult restored_sloped_water_result;
    rocksolver::restore_checkpoint(rocksolver::save_checkpoint(sloped_water_model,sloped_water_result),
        restored_sloped_water_model,restored_sloped_water_result);
    assert(restored_sloped_water_model.config.water.plane_enabled);
    assert(restored_sloped_water_model.config.water.auto_saturated_density);
    assert(restored_sloped_water_model.config.water.saturated_densities_kg_m3==
           sloped_water_config.water.saturated_densities_kg_m3);
    // A staged water-table change starts from the saved mechanical state:
    // total stress shifts by -delta(p), while instantaneous effective stress
    // and the prior displacement/history are retained.
    const auto old_zone_stress=restored_sloped_water_model.zones[0].stress_pa[0];
    const auto& first_cell=restored_sloped_water_model.physical_cells[0];
    double first_cell_pressure=0.0;
    for(const auto node_index:first_cell.nodes) {
        const auto& p=restored_sloped_water_model.nodes[node_index].position;
        const double surface=0.25+0.25*p.x;
        first_cell_pressure+=std::max(0.0,1000.0*10.0*(surface-p.z));
    }
    first_cell_pressure/=first_cell.nodes.size();
    auto dry_stage_water=sloped_water_config.water;
    dry_stage_water.enabled=false;
    const auto dried_result=rocksolver::continue_water_stage(
        restored_sloped_water_model,restored_sloped_water_result,
        dry_stage_water,sloped_water_config.solver);
    if(!dried_result.converged ||
       std::abs(restored_sloped_water_model.zones[first_cell.subzones[0]].stress_pa[0]-
                old_zone_stress-first_cell_pressure)>1.0e-7)
        throw std::runtime_error("water stage did not preserve effective stress");
    double dry_mass=0.0;
    for(const auto& node:restored_sloped_water_model.nodes)dry_mass+=node.mass_kg;
    if(std::abs(dry_mass-6100.0)>1.0e-5)
        throw std::runtime_error("water stage did not restore dry density");
    const auto rewetted_result=rocksolver::continue_water_stage(
        restored_sloped_water_model,dried_result,sloped_water_config.water,
        sloped_water_config.solver);
    if(!rewetted_result.converged ||
       std::abs(restored_sloped_water_model.zones[first_cell.subzones[0]].stress_pa[0]-
                old_zone_stress)>1.0e-7)
        throw std::runtime_error("water-stage reversal did not restore stress");
    auto stage_checkpoint=rocksolver::save_checkpoint(restored_sloped_water_model,rewetted_result);
    rocksolver::restore_checkpoint(stage_checkpoint,restored_sloped_water_model,restored_sloped_water_result);
    if(!restored_sloped_water_model.config.water.enabled ||
       restored_sloped_water_model.config.water.saturated_densities_kg_m3!=
           sloped_water_config.water.saturated_densities_kg_m3)
        throw std::runtime_error("water stage checkpoint lost the water field");
    auto compatible_checkpoint_model=rocksolver::import_abaqus_c3d4_text(interface_inp,sloped_water_config);
    compatible_checkpoint_model.config.solver.flac3d7_static_mass=true;
    compatible_checkpoint_model.config.material.flac3d7_smrm_contact=true;
    compatible_checkpoint_model.config.region_materials[0].material.flac3d7_smrm_contact=true;
    const auto compatible_result=rocksolver::solve(compatible_checkpoint_model);
    rocksolver::RockModel compatible_restored;
    rocksolver::RockResult compatible_restored_result;
    rocksolver::restore_checkpoint(rocksolver::save_checkpoint(compatible_checkpoint_model,compatible_result),
                                   compatible_restored,compatible_restored_result);
    if(!compatible_restored.config.solver.flac3d7_static_mass ||
       !compatible_restored.config.material.flac3d7_smrm_contact ||
       !compatible_restored.config.region_materials[0].material.flac3d7_smrm_contact)
        throw std::runtime_error("FLAC compatibility options were lost in checkpoint");
    // Boundary water is a separate mechanical load: with internal pore water
    // disabled, it must still increase downward displacement on the top face.
    auto external_config=legacy_interface_config;
    external_config.material.type=rocksolver::MaterialType::Elastic;
    external_config.region_materials[0].material.type=rocksolver::MaterialType::Elastic;
    external_config.gravity_m_s2=9.81;
    external_config.solver.max_cycles=10000;
    external_config.solver.convergence_tolerance=1.0e-6;
    auto no_external_model=rocksolver::import_abaqus_c3d4_text(interface_inp,external_config);
    const auto no_external_result=rocksolver::solve(no_external_model);
    external_config.water.external_enabled=true;
    external_config.water.external_surface_z_m=2.0;
    external_config.water.external_density_kg_m3=1000.0;
    external_config.water.external_outward_direction={0,0,1};
    auto external_model=rocksolver::import_abaqus_c3d4_text(interface_inp,external_config);
    const auto external_result=rocksolver::solve(external_model);
    assert(no_external_result.converged && external_result.converged);
    assert(!external_model.config.water.enabled);
    assert(external_model.nodes[4].displacement.z<no_external_model.nodes[4].displacement.z);
    rocksolver::RockModel restored_external_model;
    rocksolver::RockResult restored_external_result;
    rocksolver::restore_checkpoint(rocksolver::save_checkpoint(external_model,external_result),
        restored_external_model,restored_external_result);
    assert(restored_external_model.config.water.external_enabled);
    assert(restored_external_model.config.water.external_surface_z_m==2.0);
    auto interface_config=legacy_interface_config;
    interface_config.interfaces.push_back(rocksolver::InterfaceConfig{.name="test",.source_group="SLOPE_BAND_001",.moving_side_hint={-1,0,0}});
    interface_config.strength_reduction_interfaces_only=true;
    auto interface_model=rocksolver::import_abaqus_c3d4_text(interface_inp,interface_config);
    assert(interface_model.nodes.size()==20);
    assert(interface_model.interfaces.size()==1 && interface_model.interfaces[0].facets.size()==4);
    assert(interface_model.zones[interface_model.physical_cells[1].subzones.front()].material_index==0);
    // 601 FLAC3D 7.0 two-brick contact probe: at 2 MPa normal compression
    // and 1.5 MPa shear trial, target-face p=0/1 MPa gives 1.5/1.0 MPa
    // shear traction. The fixed-node test isolates the same contact law.
    auto interface_water_config=interface_config;
    interface_water_config.material.type=rocksolver::MaterialType::Elastic;
    interface_water_config.region_materials[0].material.type=rocksolver::MaterialType::Elastic;
    interface_water_config.interfaces[0].normal_stiffness_pa_m=1.0e9;
    interface_water_config.interfaces[0].shear_stiffness_pa_m=1.0e9;
    interface_water_config.interfaces[0].cohesion_pa=0.0;
    interface_water_config.interfaces[0].friction_deg=45.0;
    interface_water_config.initial_stress_x_pa=0.0;
    interface_water_config.initial_stress_y_pa=0.0;
    interface_water_config.initial_stress_z_pa=0.0;
    interface_water_config.gravity_m_s2=9.81;
    interface_water_config.solver.max_cycles=1;
    interface_water_config.solver.convergence_tolerance=10.0;
    const auto interface_slip_at_pressure=[&](bool wet) {
        auto config=interface_water_config;
        config.water.enabled=wet;
        config.water.table_z_m=1.0+1.0e6/9810.0;
        auto model=rocksolver::import_abaqus_c3d4_text(interface_inp,config);
        for(auto& node:model.nodes)node.fixed={true,true,true};
        for(const auto& facet:model.interfaces[0].facets)
            if(facet.master_nodes[0]!=facet.slave_nodes[0])
                model.nodes[facet.slave_nodes[0]].displacement={
                    facet.normal.x*(-0.002),
                    facet.normal.y*(-0.002)+0.0015,
                    facet.normal.z*(-0.002)};
        (void)rocksolver::solve(model);
        for(const auto& facet:model.interfaces[0].facets)
            if(facet.master_nodes[0]!=facet.slave_nodes[0] &&
               std::abs(model.nodes[facet.master_nodes[0]].position.z-1.0)<1.0e-9) {
                assert(std::abs(facet.normal_stress_pa+2.0e6)<1.0);
                return facet.shear_stress_pa;
            }
        return 0.0;
    };
    const double dry_interface_shear=interface_slip_at_pressure(false);
    const double wet_interface_shear=interface_slip_at_pressure(true);
    assert(std::abs(dry_interface_shear-1.5e6)<1.0);
    assert(std::abs(wet_interface_shear-1.0e6)<1.0);
    // FLAC `zone interface node initialize-stresses` preloads the host-face
    // total traction. At zero displacement a vertical split face under
    // sigma_xx=-80 kPa already carries -80 kPa, even with fixed pore water.
    auto prestress_config=interface_water_config;
    prestress_config.initial_stress_x_pa=-80000.0;
    prestress_config.water.enabled=true;
    prestress_config.water.table_z_m=6.0;
    auto prestress_model=rocksolver::import_abaqus_c3d4_text(interface_inp,prestress_config);
    for(auto& node:prestress_model.nodes)node.fixed={true,true,true};
    (void)rocksolver::solve(prestress_model);
    for(const auto& facet:prestress_model.interfaces[0].facets)
        assert(std::abs(facet.normal_stress_pa+80000.0)<1.0);
    assert(prestress_model.initialize_interface_stresses);
    // Current Interfaces are node based (one master/slave node per facet). They
    // must survive the same portable checkpoint path used by slope branches.
    interface_config.solver.max_cycles=60;
    interface_config.solver.convergence_tolerance=10.0;
    auto checkpoint_interface_model=rocksolver::import_abaqus_c3d4_text(interface_inp,interface_config);
    const auto checkpoint_interface_result=rocksolver::solve(checkpoint_interface_model);
    const auto interface_checkpoint=rocksolver::save_checkpoint(checkpoint_interface_model,checkpoint_interface_result);
    rocksolver::RockModel restored_interface_model;
    rocksolver::RockResult restored_interface_result;
    rocksolver::restore_checkpoint(interface_checkpoint,restored_interface_model,restored_interface_result);
    assert(restored_interface_model.initialize_interface_stresses);
    assert(restored_interface_model.cable_physical_host_interpolation);
    assert(restored_interface_model.cable_effective_confinement);
    assert(restored_interface_model.interfaces.size()==1);
    assert(restored_interface_model.interfaces[0].facets.size()==4);
    for(const auto& facet:restored_interface_model.interfaces[0].facets) {
        assert(facet.master_nodes.size()==1);
        assert(facet.slave_nodes.size()==1);
    }
    // Each legacy checkpoint retains the Cable and Interface mechanics that
    // were actually used in its saved equilibrium state.
    auto legacy_interface_checkpoint=interface_checkpoint;
    remove_098_checkpoint_extension(legacy_interface_checkpoint,
                                    interface_model.config.region_materials.size());
    const std::string current_header="ROCKSTATE-1 "+std::string(rocksolver::kVersion)+" ";
    assert(legacy_interface_checkpoint.starts_with(current_header));
    const std::string flags_marker="1 1 1 \""+std::string(rocksolver::kVersion)+"\"";
    const auto flags_at=legacy_interface_checkpoint.rfind(flags_marker);
    assert(flags_at!=std::string::npos);
    legacy_interface_checkpoint.erase(flags_at+4,2);
    legacy_interface_checkpoint.replace(0,current_header.size(),
        "ROCKSTATE-1 0.9.4-cable-cell-host-weights ");
    rocksolver::RockModel restored_legacy_interface_model;
    rocksolver::RockResult restored_legacy_interface_result;
    rocksolver::restore_checkpoint(legacy_interface_checkpoint,
        restored_legacy_interface_model,restored_legacy_interface_result);
    assert(restored_legacy_interface_model.initialize_interface_stresses);
    assert(restored_legacy_interface_model.cable_physical_host_interpolation);
    assert(!restored_legacy_interface_model.cable_effective_confinement);
    const auto old_flags_at=legacy_interface_checkpoint.rfind(
        "1 1 \""+std::string(rocksolver::kVersion)+"\"");
    assert(old_flags_at!=std::string::npos);
    legacy_interface_checkpoint.erase(old_flags_at+2,2);
    legacy_interface_checkpoint.replace(0,std::string("ROCKSTATE-1 0.9.4-cable-cell-host-weights ").size(),
        "ROCKSTATE-1 0.9.3-interface-initial-stress ");
    rocksolver::restore_checkpoint(legacy_interface_checkpoint,
        restored_legacy_interface_model,restored_legacy_interface_result);
    assert(restored_legacy_interface_model.initialize_interface_stresses);
    assert(!restored_legacy_interface_model.cable_physical_host_interpolation);
    assert(!restored_legacy_interface_model.cable_effective_confinement);
    const auto preload_at=legacy_interface_checkpoint.rfind(
        "1 \""+std::string(rocksolver::kVersion)+"\"");
    assert(preload_at!=std::string::npos);
    legacy_interface_checkpoint.erase(preload_at,2);
    legacy_interface_checkpoint.replace(0,std::string("ROCKSTATE-1 0.9.3-interface-initial-stress ").size(),
        "ROCKSTATE-1 0.9.2-external-water-pressure ");
    rocksolver::restore_checkpoint(legacy_interface_checkpoint,
        restored_legacy_interface_model,restored_legacy_interface_result);
    assert(!restored_legacy_interface_model.initialize_interface_stresses);
    assert(!restored_legacy_interface_model.cable_physical_host_interpolation);
    assert(!restored_legacy_interface_model.cable_effective_confinement);
    // A third, optional reduction target is additive: continuum-only keeps
    // Interface strength unchanged. The legacy/default flags remain false.
    auto continuum_interface_config=interface_config;
    continuum_interface_config.strength_reduction_interfaces_only=false;
    continuum_interface_config.strength_reduction_continuum_only=true;
    continuum_interface_config.solver.max_cycles=60;
    continuum_interface_config.solver.convergence_tolerance=10.0;
    auto continuum_interface_model=rocksolver::import_abaqus_c3d4_text(interface_inp,continuum_interface_config);
    const auto continuum_interface_result=rocksolver::solve(continuum_interface_model);
    const auto interface_json=rocksolver::result_to_json(continuum_interface_model,continuum_interface_result);
    assert(interface_json.find("\"shearPastFacets\"")!=std::string::npos);
    assert(interface_json.find("\"noShearFailureFacets\"")!=std::string::npos);
    continuum_interface_model.config.solver.max_cycles=5000;
    int reduction_block_cycles=0;
    (void)rocksolver::reduction_trial(continuum_interface_model,1.3,[&](const rocksolver::SolveProgress& progress){
        reduction_block_cycles=std::max(reduction_block_cycles,progress.max_cycles);
    });
    assert(reduction_block_cycles==4000);
    assert(continuum_interface_model.config.solver.max_cycles==5000);
    assert(continuum_interface_model.interfaces[0].strength_reduction_factor==1.0);
    assert(continuum_interface_model.config.material.strength_reduction_factor==1.3);

    // Both persisted target flags encode the fourth mode: region material bands only.
    auto band_only_config=legacy_interface_config;
    band_only_config.strength_reduction_interfaces_only=true;
    band_only_config.strength_reduction_continuum_only=true;
    band_only_config.solver.max_cycles=60;
    band_only_config.solver.convergence_tolerance=10.0;
    auto band_only_model=rocksolver::import_abaqus_c3d4_text(interface_inp,band_only_config);
    (void)rocksolver::solve(band_only_model);
    (void)rocksolver::reduction_trial(band_only_model,1.3);
    assert(band_only_model.config.material.strength_reduction_factor==1.0);
    assert(band_only_model.config.region_materials[0].material.strength_reduction_factor==1.3);

    // Stage 2 steady saturated seepage: a unit C3D8 with h=1 m at X-min and
    // h=0 at X-max has the exact linear solution h=1-x and qx=Kxx.
    const std::string seepage_cube_inp=R"INP(*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
*ELEMENT,TYPE=C3D8,ELSET=HOST
1,1,2,3,4,5,6,7,8
)INP";
    auto seepage_config=skew_config;
    seepage_config.gravity_m_s2=9.81;
    seepage_config.seepage.enabled=true;
    seepage_config.seepage.fluid_density_kg_m3=1000.0;
    seepage_config.seepage.conductivity_m_s={{{2e-6,1e-6,5e-7,0,0,0}}};
    seepage_config.seepage.relative_tolerance=1e-12;
    seepage_config.seepage.boundaries={
        {.name="upstream",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=1.0},
        {.name="downstream",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMax,.value=0.0},
        {.name="impermeable top",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::ZMax},
    };
    auto seepage_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,seepage_config);
    const auto seepage=rocksolver::solve_steady_seepage(seepage_model);
    assert(seepage.converged);
    assert(seepage.relative_residual<1e-12);
    for(std::size_t i=0;i<seepage_model.nodes.size();++i)
        assert(std::abs(seepage.head_m[i]-(1.0-seepage_model.nodes[i].position.x))<1e-10);
    for(const auto& velocity:seepage.zone_darcy_velocity_m_s){
        assert(std::abs(velocity.x-2e-6)<1e-15);
        assert(std::abs(velocity.y)<1e-15);
        assert(std::abs(velocity.z)<1e-15);
    }
    assert(std::abs(seepage.inflow_m3_s-2e-6)<1e-14);
    assert(std::abs(seepage.outflow_m3_s-2e-6)<1e-14);
    assert(seepage.relative_mass_balance_error<1e-10);

    // Pressure and drain boundaries use h=z+p/(rho*g), while outward flux is
    // positive (infiltration is negative). These two forms must be accepted
    // independently of the prescribed-head path above.
    auto pressure_config=seepage_config;
    pressure_config.seepage.boundaries={
        {.name="pressurized bottom",.kind=rocksolver::HydraulicBoundaryKind::PorePressure,
         .selector=rocksolver::HydraulicBoundarySelector::ZMin,.value=9810.0},
        {.name="drained top",.kind=rocksolver::HydraulicBoundaryKind::Drain,
         .selector=rocksolver::HydraulicBoundarySelector::ZMax},
    };
    pressure_config.seepage.conductivity_m_s={{{1e-6,1e-6,1e-6,0,0,0}}};
    auto pressure_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,pressure_config);
    const auto pressure=rocksolver::solve_steady_seepage(pressure_model);
    assert(pressure.converged);
    for(std::size_t i=0;i<pressure_model.nodes.size();++i)
        assert(std::abs(pressure.nodal_pore_pressure_pa[i]-9810.0*(1.0-pressure_model.nodes[i].position.z))<1e-6);

    auto flux_config=seepage_config;
    flux_config.seepage.conductivity_m_s={{{1e-6,1e-6,1e-6,0,0,0}}};
    flux_config.seepage.boundaries={
        {.name="infiltration",.kind=rocksolver::HydraulicBoundaryKind::Flux,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=-1e-6},
        {.name="fixed outlet head",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMax,.value=0.0},
    };
    auto flux_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,flux_config);
    const auto flux=rocksolver::solve_steady_seepage(flux_model);
    assert(flux.converged);
    for(std::size_t i=0;i<flux_model.nodes.size();++i)
        assert(std::abs(flux.head_m[i]-(1.0-flux_model.nodes[i].position.x))<1e-10);
    assert(std::abs(flux.inflow_m3_s-1e-6)<1e-14);
    assert(std::abs(flux.outflow_m3_s-1e-6)<1e-14);
    assert(flux.relative_mass_balance_error<1e-10);

    auto invalid_tensor_config=seepage_config;
    invalid_tensor_config.seepage.conductivity_m_s={{{1e-6,1e-6,-1e-6,0,0,0}}};
    auto invalid_tensor_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,invalid_tensor_config);
    bool rejected_invalid_tensor=false;
    try{(void)rocksolver::solve_steady_seepage(invalid_tensor_model);}catch(const std::invalid_argument&){rejected_invalid_tensor=true;}
    assert(rejected_invalid_tensor);

    // Stage 3 transient saturated seepage: the prescribed upstream head rises
    // and falls linearly. Interior nodes must retain a delayed response, and
    // every requested output time must be preserved.
    const std::string transient_inp=R"INP(*NODE
1,0,0,0
2,0,1,0
3,0,1,1
4,0,0,1
5,1,0,0
6,1,1,0
7,1,1,1
8,1,0,1
9,2,0,0
10,2,1,0
11,2,1,1
12,2,0,1
*ELEMENT,TYPE=C3D8,ELSET=HOST
1,1,5,6,2,4,8,7,3
2,5,9,10,6,8,12,11,7
)INP";
    auto transient_config=seepage_config;
    transient_config.seepage.transient=true;
    transient_config.seepage.duration_s=20.0;
    transient_config.seepage.time_step_s=1.0;
    transient_config.seepage.initial_head_m=0.0;
    transient_config.seepage.conductivity_m_s={{{0.1,0.1,0.1,0,0,0}}};
    transient_config.seepage.specific_storage_1_m={1.0};
    transient_config.seepage.boundaries={
        {.name="reservoir",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=0.0,
         .history={{0.0,0.0},{10.0,1.0},{20.0,0.0}}},
        {.name="closed far end",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::XMax},
    };
    auto transient_model=rocksolver::import_abaqus_c3d4_text(transient_inp,transient_config);
    const auto transient=rocksolver::solve_transient_seepage(transient_model);
    assert(transient.converged);
    assert(transient.history.size()==21);
    assert(std::abs(transient.history.front().time_s)<1e-15);
    assert(std::abs(transient.history.back().time_s-20.0)<1e-12);
    for(std::size_t i=0;i<transient_model.nodes.size();++i){
        const double x=transient_model.nodes[i].position.x;
        if(std::abs(x)<1e-12){
            assert(std::abs(transient.history[10].head_m[i]-1.0)<1e-12);
            assert(std::abs(transient.history.back().head_m[i])<1e-12);
        }else if(std::abs(x-1.0)<1e-12){
            assert(transient.history[10].head_m[i]>0.0);
            assert(transient.history[10].head_m[i]<1.0);
            assert(transient.history.back().head_m[i]>0.0);
        }
    }
    assert(transient.history[10].storage_change_m3>0.0);
    assert(transient.history.back().storage_change_m3<0.0);

    // 20 / 0.0125 is exactly 1600 steps. Repeated floating-point time
    // addition previously created a spurious near-zero 1601st step.
    auto exact_step_config=transient_config;
    exact_step_config.seepage.time_step_s=0.0125;
    auto exact_step_model=rocksolver::import_abaqus_c3d4_text(transient_inp,exact_step_config);
    const auto exact_step_result=rocksolver::solve_transient_seepage(exact_step_model);
    assert(exact_step_result.converged);
    assert(exact_step_result.history.size()==1601);
    assert(exact_step_result.history.back().time_s==20.0);

    // In a rigid saturated column Ss=rho*g/M. An explicit Biot modulus
    // must reproduce the same diffusion field as the equivalent Stage 3 Ss.
    auto biot_storage_config=transient_config;
    biot_storage_config.seepage.specific_storage_1_m={2.0};
    biot_storage_config.seepage.biot_modulus_pa=9810.0;
    auto biot_storage_model=rocksolver::import_abaqus_c3d4_text(transient_inp,biot_storage_config);
    const auto biot_storage=rocksolver::solve_transient_seepage(biot_storage_model);
    assert(biot_storage.converged);
    assert(biot_storage.history.size()==transient.history.size());
    for(std::size_t step=0;step<transient.history.size();++step)
        for(std::size_t node=0;node<transient_model.nodes.size();++node)
            assert(std::abs(biot_storage.history[step].head_m[node]-
                            transient.history[step].head_m[node])<1e-10);

    // A prescribed head change must feed the elastic volume change back into
    // fluid storage when Biot coupling is enabled. Alpha=0 is the rigid-flow
    // control; both runs must close the discrete water balance.
    auto coupled_config=transient_config;
    coupled_config.seepage.solve_mechanics=true;
    coupled_config.seepage.bidirectional_coupling=true;
    coupled_config.seepage.biot_modulus_pa=1.0e7;
    coupled_config.seepage.duration_s=2.0;
    coupled_config.seepage.time_step_s=1.0;
    coupled_config.seepage.initial_head_m=0.0;
    coupled_config.seepage.conductivity_m_s={{{1.0e-4,1.0e-4,1.0e-4,0,0,0}}};
    coupled_config.seepage.boundaries={
        {.name="pressure step",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=1.0},
        {.name="closed far end",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::XMax},
    };
    coupled_config.material.elastic_modulus_pa=1.0e6;
    coupled_config.solver.convergence_tolerance=1.0e-9;
    // The assertion below requires 1e-7 mass balance, so request a tighter
    // coupling solve rather than accepting the older 1e-5 stopping criterion.
    coupled_config.seepage.coupling_relative_tolerance=1.0e-8;
    coupled_config.seepage.coupling_max_iterations=200;
    auto uncoupled_config=coupled_config;
    uncoupled_config.seepage.biot_coefficient=0.0;
    auto uncoupled_model=rocksolver::import_abaqus_c3d4_text(transient_inp,uncoupled_config);
    const auto uncoupled_result=rocksolver::solve_transient_seepage_mechanics(uncoupled_model);
    assert(uncoupled_result.converged);
    auto coupled_model=rocksolver::import_abaqus_c3d4_text(transient_inp,coupled_config);
    const auto coupled_result=rocksolver::solve_transient_seepage_mechanics(coupled_model);
    assert(coupled_result.converged);
    assert(coupled_model.hydraulic.history.size()==3);
    for(std::size_t step=1;step<3;++step){
        assert(coupled_model.hydraulic.history[step].relative_mass_balance_error<1.0e-7);
        assert(uncoupled_model.hydraulic.history[step].relative_mass_balance_error<1.0e-7);
    }
    double head_difference=0.0;
    for(std::size_t i=0;i<coupled_model.nodes.size();++i)
        head_difference=std::max(head_difference,std::abs(coupled_model.hydraulic.head_m[i]-
                                                       uncoupled_model.hydraulic.head_m[i]));
    assert(head_difference>1.0e-5);

    // Prescribe FLAC's t=1 nodal pressure field to isolate the mechanical
    // stage response from the still-unmatched bidirectional flow update.
    auto prescribed_config=coupled_config;
    prescribed_config.material.density_kg_m3=2300.0;
    prescribed_config.solver.nodal_mixed_discretization=true;
    auto prescribed_model=rocksolver::import_abaqus_c3d4_text(transient_inp,prescribed_config);
    prescribed_model.hydraulic.valid=true;
    prescribed_model.hydraulic.nodal_pore_pressure_pa.resize(prescribed_model.nodes.size());
    for(std::size_t i=0;i<prescribed_model.nodes.size();++i)
        prescribed_model.hydraulic.nodal_pore_pressure_pa[i]=9810.0*(1.0-prescribed_model.nodes[i].position.z);
    const auto prescribed_initial=rocksolver::solve(prescribed_model);
    assert(prescribed_initial.converged);
    auto prescribed_next=prescribed_model.hydraulic;
    const std::array<double,3> flac_pressure_base={19620.0,6196.71231242926,8720.26034400602};
    for(std::size_t i=0;i<prescribed_model.nodes.size();++i)
        prescribed_next.nodal_pore_pressure_pa[i]=flac_pressure_base[static_cast<std::size_t>(prescribed_model.nodes[i].position.x)]-
            9810.0*prescribed_model.nodes[i].position.z;
    prescribed_next.converged=true;
    const auto prescribed_result=rocksolver::continue_hydraulic_stage(
        prescribed_model,prescribed_initial,prescribed_next,prescribed_config.solver);
    assert(prescribed_result.converged);
    const std::array<double,3> flac_top_uz={-0.00594167040097213,-0.00810314421588222,-0.0102646178857424};
    for(const auto& node:prescribed_model.nodes)
        if(node.position.y==0.0&&node.position.z==1.0)
            assert(std::abs(node.displacement.z-flac_top_uz[static_cast<std::size_t>(node.position.x)])<1.0e-7);

    // FLAC3D 7.00 R126 coupled two-cell column: the first step tests
    // fixed-pressure mechanical-source redistribution, while t=1/2 s test
    // the separated flow/mechanical Biot stress increments.
    auto first_config=prescribed_config;
    first_config.seepage.initial_head_m=1.0;
    first_config.seepage.boundaries[0].value=2.0;
    first_config.seepage.duration_s=0.01;
    first_config.seepage.time_step_s=0.01;
    auto first_model=rocksolver::import_abaqus_c3d4_text(transient_inp,first_config);
    const auto first_result=rocksolver::solve_transient_seepage_mechanics(first_model);
    assert(first_result.converged);
    for(std::size_t ci=0;ci<first_model.physical_cells.size();++ci){
        const auto& cell=first_model.physical_cells[ci];
        double e=0.0,v=0.0;
        for(const auto zi:cell.subzones){
            const auto& z=first_model.zones[zi];
            const double w=z.volume_m3*z.integration_weight;
            e+=(z.strain[0]+z.strain[1]+z.strain[2])*w;
            v+=w;
        }
        const std::array<double,2> flac_strain={-0.00794615648907445,-0.00838174282778547};
        assert(std::abs(e/v-flac_strain[ci])<1.0e-7);
    }
    const std::array<double,3> flac_first_pressure={19620.0,5482.13579808898,9821.99951327106};
    for(std::size_t i=0;i<first_model.nodes.size();++i)
        if(first_model.nodes[i].position.y==0.0&&first_model.nodes[i].position.z==0.0)
            assert(std::abs(first_model.hydraulic.nodal_pore_pressure_pa[i]-
                            flac_first_pressure[static_cast<std::size_t>(first_model.nodes[i].position.x)])<15.0);
    auto long_config=first_config;
    long_config.seepage.duration_s=2.0;
    // The 0.01 s increments are small: resolve inner mechanical equilibrium
    // more tightly than the outer mass-balance acceptance; do not relax it.
    long_config.solver.convergence_tolerance=1.0e-12;
    long_config.solver.max_cycles=30000;
    long_config.seepage.coupling_max_iterations=2000;
    auto long_model=rocksolver::import_abaqus_c3d4_text(transient_inp,long_config);
    const auto long_result=rocksolver::solve_transient_seepage_mechanics(long_model);
    if (!long_result.converged)
        std::cerr << "Long Biot regression stopped: " << long_model.hydraulic.failure_reason
                  << "; accepted t=" << long_model.hydraulic.time_s << '\n';
    assert(long_result.converged);
    assert(long_model.hydraulic.history.size()==201);
    const std::array<std::array<double,3>,2> flac_coupled_pressure{{
        {{19620.0,6210.21875943262,8693.44208352226}},
        {{19620.0,6798.68457721497,7832.89648695815}},
    }};
    const std::array<std::array<double,3>,2> flac_coupled_top_uz{{
        {{-0.00787047185539145,-0.00810310799800377,-0.00833574429709695}},
        {{-0.00780392472686383,-0.00804435073017895,-0.00828477656124616}},
    }};
    const std::array<double,2> flac_coupled_mid_ux={0.000232636218829794,0.00024042590794355};
    const std::array<std::array<std::array<double,3>,2>,2> flac_coupled_stress{{
        {{{{-5117.93196470914,-5207.4074333659,-11351.0919904778}},
          {{-4978.74787632454,-4889.27240766779,-11211.9079036764}}}},
        {{{{-5165.46026122431,-5257.93176406187,-11353.4223266705}},
          {{-5021.6156632354,-4929.14416039784,-11209.5777361443}}}},
    }};
    const std::array<std::array<double,2>,2> flac_coupled_qx{{
        {{0.00013669501774345,-0.0000253131837319233}},
        {{0.000130696385553998,-0.0000105424251758226}},
    }};
    for(const auto ti:{100,200}){
        const auto& snapshot=long_model.hydraulic.history[ti];
        const auto ref=static_cast<std::size_t>(ti/100-1);
        assert(snapshot.relative_mass_balance_error<1.0e-7);
        for(std::size_t i=0;i<long_model.nodes.size();++i){
            const auto& p=long_model.nodes[i].position;
            if(p.y!=0.0)continue;
            const auto x=static_cast<std::size_t>(p.x);
            if(p.z==0.0)
                assert(std::abs(snapshot.nodal_pore_pressure_pa[i]-flac_coupled_pressure[ref][x])<20.0);
            else if(p.z==1.0){
                assert(std::abs(snapshot.nodal_displacement_m[i].z-flac_coupled_top_uz[ref][x])<5.0e-7);
                if(x==1)
                    assert(std::abs(snapshot.nodal_displacement_m[i].x-flac_coupled_mid_ux[ref])<5.0e-7);
            }
        }
        for(std::size_t ci=0;ci<long_model.physical_cells.size();++ci){
            const auto& cell=long_model.physical_cells[ci];
            double sxx=0.0,syy=0.0,szz=0.0,qx=0.0,v=0.0;
            for(const auto zi:cell.subzones){
                const auto& z=long_model.zones[zi];
                const double w=z.volume_m3*z.integration_weight;
                sxx+=snapshot.zone_stress_pa[zi][0]*w;
                syy+=snapshot.zone_stress_pa[zi][1]*w;
                szz+=snapshot.zone_stress_pa[zi][2]*w;
                qx+=snapshot.zone_darcy_velocity_m_s[zi].x*w;
                v+=w;
            }
            assert(std::abs(sxx/v-flac_coupled_stress[ref][ci][0])<1.0);
            assert(std::abs(syy/v-flac_coupled_stress[ref][ci][1])<1.0);
            assert(std::abs(szz/v-flac_coupled_stress[ref][ci][2])<1.0);
            assert(std::abs(qx/v-flac_coupled_qx[ref][ci])<5.0e-7);
        }
    }

    // Refining the coupled step must give a Cauchy sequence for the pore
    // pressure and displacement, independent of the FLAC reference step.
    auto coarse_time_config=long_config;
    coarse_time_config.seepage.time_step_s=0.02;
    auto coarse_time_model=rocksolver::import_abaqus_c3d4_text(transient_inp,coarse_time_config);
    assert(rocksolver::solve_transient_seepage_mechanics(coarse_time_model).converged);
    auto fine_time_config=long_config;
    fine_time_config.seepage.time_step_s=0.005;
    auto fine_time_model=rocksolver::import_abaqus_c3d4_text(transient_inp,fine_time_config);
    assert(rocksolver::solve_transient_seepage_mechanics(fine_time_model).converged);
    double coarse_step_error=0.0,fine_step_error=0.0;
    double coarse_displacement_error=0.0,fine_displacement_error=0.0;
    for(std::size_t i=0;i<long_model.nodes.size();++i){
        coarse_step_error=std::max(coarse_step_error,std::abs(
            coarse_time_model.hydraulic.nodal_pore_pressure_pa[i]-
            long_model.hydraulic.nodal_pore_pressure_pa[i]));
        fine_step_error=std::max(fine_step_error,std::abs(
            fine_time_model.hydraulic.nodal_pore_pressure_pa[i]-
            long_model.hydraulic.nodal_pore_pressure_pa[i]));
        coarse_displacement_error=std::max(coarse_displacement_error,std::abs(
            coarse_time_model.nodes[i].displacement.z-
            long_model.nodes[i].displacement.z));
        fine_displacement_error=std::max(fine_displacement_error,std::abs(
            fine_time_model.nodes[i].displacement.z-
            long_model.nodes[i].displacement.z));
    }
    assert(fine_step_error<coarse_step_error);
    assert(fine_step_error<20.0);
    assert(fine_displacement_error<coarse_displacement_error);

    // Terzaghi excess-pressure relaxation: a 20 m laterally confined column,
    // drained at the top and sealed at its base. The small nonzero gravity
    // represents a nearly uniform initial excess pressure through the
    // existing hydraulic-head input (variation < 20 Pa over 20 m).
    std::ostringstream terzaghi_mesh;
    terzaghi_mesh<<"*NODE\n";
    for(int layer=0;layer<=20;++layer){
        const int first=4*layer+1;
        terzaghi_mesh<<first<<",0,0,"<<layer<<"\n"
                     <<first+1<<",1,0,"<<layer<<"\n"
                     <<first+2<<",1,1,"<<layer<<"\n"
                     <<first+3<<",0,1,"<<layer<<"\n";
    }
    terzaghi_mesh<<"*ELEMENT,TYPE=C3D8,ELSET=HOST\n";
    for(int layer=0;layer<20;++layer){
        const int b=4*layer+1,t=b+4;
        terzaghi_mesh<<layer+1<<","<<b<<","<<b+1<<","<<b+2<<","<<b+3
                     <<","<<t<<","<<t+1<<","<<t+2<<","<<t+3<<"\n";
    }
    auto terzaghi_config=long_config;
    terzaghi_config.gravity_m_s2=0.001;
    terzaghi_config.material.elastic_modulus_pa=9.0*5.0e8*2.0e8/(3.0*5.0e8+2.0e8);
    terzaghi_config.material.poisson_ratio=(3.0*5.0e8-2.0*2.0e8)/
        (2.0*(3.0*5.0e8+2.0e8));
    terzaghi_config.seepage.biot_modulus_pa=4.0e9;
    terzaghi_config.solver.convergence_tolerance=1.0e-10;
    // A 20:1 column needs more mechanical relaxation than the two-cell fixture.
    terzaghi_config.solver.max_cycles=200000;
    terzaghi_config.seepage.relative_tolerance=1.0e-14;
    terzaghi_config.seepage.initial_head_m=83916.0;
    terzaghi_config.seepage.conductivity_m_s={{{1.0e-10,1.0e-10,1.0e-10,0,0,0}}};
    terzaghi_config.seepage.duration_s=500.0;
    terzaghi_config.seepage.time_step_s=0.5;
    terzaghi_config.seepage.boundaries={
        {.name="drained top",.kind=rocksolver::HydraulicBoundaryKind::PorePressure,
         .selector=rocksolver::HydraulicBoundarySelector::ZMax,.value=0.0},
        {.name="sealed base",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::ZMin},
    };
    auto terzaghi_model=rocksolver::import_abaqus_c3d4_text(
        terzaghi_mesh.str(),terzaghi_config);
    const auto terzaghi_result=rocksolver::solve_transient_seepage_mechanics(terzaghi_model);
    if (!terzaghi_result.converged)
        std::cerr << "Terzaghi regression stopped: " << terzaghi_model.hydraulic.failure_reason
                  << "; accepted t=" << terzaghi_model.hydraulic.time_s
                  << "; force ratio=" << terzaghi_result.final_force_ratio
                  << "; cycles=" << terzaghi_result.total_cycles << '\n';
    assert(terzaghi_result.converged);
    std::size_t terzaghi_top_node=0;
    for(std::size_t i=0;i<terzaghi_model.nodes.size();++i)
        if(terzaghi_model.nodes[i].position.x==0.0 &&
           terzaghi_model.nodes[i].position.y==0.0 &&
           terzaghi_model.nodes[i].position.z==20.0)
            terzaghi_top_node=i;
    // With unit plan area and fixed lateral displacements, the actual
    // mechanical fluid-content increment is alpha times the top settlement
    // increment. The reported mass ledger must use the accepted mechanics,
    // not the previous fixed-point strain guess.
    for(std::size_t step=1;step<terzaghi_model.hydraulic.history.size();++step){
        const auto& current=terzaghi_model.hydraulic.history[step];
        const auto& previous=terzaghi_model.hydraulic.history[step-1];
        const double actual_mechanical=current.nodal_displacement_m[terzaghi_top_node].z-
            previous.nodal_displacement_m[terzaghi_top_node].z;
        assert(std::abs(current.mechanical_storage_change_m3-actual_mechanical)<
               1.0e-8*std::max(std::abs(actual_mechanical),1.0e-12));
    }
    const double confined=5.0e8+4.0*2.0e8/3.0;
    const double diffusivity=1.0e-10/(1.0/4.0e9+1.0/confined);
    double terzaghi_max_pressure_error=0.0;
    for(std::size_t i=0;i<terzaghi_model.nodes.size();++i){
        const double z=terzaghi_model.nodes[i].position.z;
        double series=0.0;
        for(int mode=0;mode<30;++mode){
            const double odd=2.0*mode+1.0;
            const double lambda=odd*3.14159265358979323846/40.0;
            series+=(mode%2==0?1.0:-1.0)*4.0/(odd*3.14159265358979323846)*
                std::cos(lambda*z)*std::exp(-diffusivity*lambda*lambda*500.0);
        }
        terzaghi_max_pressure_error=std::max(terzaghi_max_pressure_error,std::abs(
            terzaghi_model.hydraulic.nodal_pore_pressure_pa[i]-83916.0*series));
    }
    assert(terzaghi_max_pressure_error<0.02*83916.0);

    auto flow_only_config=prescribed_config;
    flow_only_config.seepage.biot_coefficient=0.0;
    flow_only_config.seepage.initial_head_m=1.0;
    flow_only_config.seepage.time_step_s=0.01;
    flow_only_config.seepage.boundaries[0].value=2.0;
    auto flow_only_model=rocksolver::import_abaqus_c3d4_text(transient_inp,flow_only_config);
    const auto flow_only_result=rocksolver::solve_transient_seepage_mechanics(flow_only_model);
    assert(flow_only_result.converged);
    const std::array<std::array<double,3>,2> flac_flow_only_pressure{{
        {{19620.0,10718.3856949098,9898.56454143512}},
        {{19620.0,11474.1443465575,10123.3543396867}},
    }};
    for(const auto time_index:{100,200}){
        const auto& snapshot=flow_only_model.hydraulic.history[time_index];
        for(std::size_t i=0;i<flow_only_model.nodes.size();++i)
            if(flow_only_model.nodes[i].position.y==0.0&&flow_only_model.nodes[i].position.z==0.0)
                assert(std::abs(snapshot.nodal_pore_pressure_pa[i]-
                                flac_flow_only_pressure[time_index/100-1][static_cast<std::size_t>(flow_only_model.nodes[i].position.x)])<3.0);
    }
    auto low_m_config=flow_only_config;
    low_m_config.seepage.biot_coefficient=1.0;
    low_m_config.seepage.biot_modulus_pa=1.0e5;
    low_m_config.solver.convergence_tolerance=1.0e-12;
    low_m_config.seepage.coupling_max_iterations=1000;
    low_m_config.seepage.coupling_relative_tolerance=1.0e-7;
    auto low_m_model=rocksolver::import_abaqus_c3d4_text(transient_inp,low_m_config);
    const auto low_m_result=rocksolver::solve_transient_seepage_mechanics(low_m_model);
    assert(low_m_result.converged);
    for(const auto time_index:{100,200}){
        const auto& snapshot=low_m_model.hydraulic.history[time_index];
        assert(snapshot.relative_mass_balance_error<1.0e-7);
        for(const double pressure:snapshot.nodal_pore_pressure_pa)
            assert(std::isfinite(pressure));
    }
    auto stricter_low_m_config=low_m_config;
    stricter_low_m_config.solver.convergence_tolerance=1.0e-13;
    stricter_low_m_config.seepage.coupling_relative_tolerance=1.0e-8;
    auto stricter_low_m_model=rocksolver::import_abaqus_c3d4_text(transient_inp,stricter_low_m_config);
    assert(rocksolver::solve_transient_seepage_mechanics(stricter_low_m_model).converged);
    for(std::size_t i=0;i<low_m_model.nodes.size();++i)
        assert(std::abs(low_m_model.hydraulic.nodal_pore_pressure_pa[i]-
                        stricter_low_m_model.hydraulic.nodal_pore_pressure_pa[i])<0.01);

    // Stage 4 rainfall/unsaturated flow: rainfall is positive into the model,
    // actual infiltration is capped by saturated K33, and the remainder is
    // reported as runoff.  The retained water must increase saturation while
    // preserving the fully discrete mass balance.
    auto rainfall_config=seepage_config;
    rainfall_config.seepage.transient=true;
    rainfall_config.seepage.unsaturated=true;
    rainfall_config.seepage.duration_s=0.1;
    rainfall_config.seepage.time_step_s=0.1;
    rainfall_config.seepage.initial_head_m=-1.0;
    rainfall_config.seepage.conductivity_m_s={{{1e-4,1e-4,1e-4,0,0,0}}};
    rainfall_config.seepage.specific_storage_1_m={1e-4};
    rainfall_config.seepage.unsaturated_materials={{.porosity=0.3,.residual_saturation=0.1,
        .van_genuchten_alpha_1_m=1.0,.van_genuchten_n=2.0,.mualem_l=0.5}};
    rainfall_config.seepage.nonlinear_max_iterations=100;
    rainfall_config.seepage.nonlinear_relative_tolerance=1e-9;
    rainfall_config.seepage.boundaries={{.name="rain",.kind=rocksolver::HydraulicBoundaryKind::Rainfall,
        .selector=rocksolver::HydraulicBoundarySelector::ZMax,.value=2e-4}};
    auto rainfall_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,rainfall_config);
    const auto rainfall=rocksolver::solve_unsaturated_seepage(rainfall_model);
    assert(rainfall.converged);
    assert(rainfall.history.size()==2);
    assert(std::abs(rainfall.rainfall_infiltration_m3_s-1e-4)<1e-14);
    assert(std::abs(rainfall.rainfall_runoff_m3_s-1e-4)<1e-14);
    assert(rainfall.relative_mass_balance_error<1e-8);
    assert(rainfall.zone_saturation.size()==rainfall_model.zones.size());
    assert(*std::max_element(rainfall.history.back().zone_saturation.begin(),
                            rainfall.history.back().zone_saturation.end())>
           *std::max_element(rainfall.history.front().zone_saturation.begin(),
                            rainfall.history.front().zone_saturation.end()));

    // FLAC3D 7's separate no-capillarity law keeps p=0 below saturation one.
    // This verifies its explicit nodal water volume and rain/runoff accounting
    // without changing the van Genuchten branch above.
    auto flac_unsat_config=rainfall_config;
    flac_unsat_config.seepage.flac_no_capillarity=true;
    flac_unsat_config.seepage.initial_saturation=0.4;
    auto flac_unsat_model=rocksolver::import_abaqus_c3d4_text(seepage_cube_inp,flac_unsat_config);
    const auto flac_unsat=rocksolver::solve_unsaturated_seepage(flac_unsat_model);
    assert(flac_unsat.converged);
    assert(flac_unsat.history.size()==2);
    assert(std::abs(flac_unsat.rainfall_infiltration_m3_s-1e-4)<1e-14);
    assert(std::abs(flac_unsat.rainfall_runoff_m3_s-1e-4)<1e-14);
    assert(std::abs(flac_unsat.history.back().storage_change_m3-1e-5)<1e-12);
    assert(flac_unsat.relative_mass_balance_error<1e-10);
    for(const double p:flac_unsat.nodal_pore_pressure_pa)assert(p==0.0);

    std::cout << "RockSolver native tests passed\n";
    return 0;
}
