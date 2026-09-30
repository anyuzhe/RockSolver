#include "rocksolver/rocksolver.hpp"
#include "rocksolver/smrm_material.hpp"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {

rocksolver::RockModel two_cell_model(double duration_s) {
    constexpr auto mesh=R"INP(*NODE
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
    rocksolver::ModelConfig config;
    config.scenario=rocksolver::Scenario::Slope;
    config.material.type=rocksolver::MaterialType::SMRM;
    // Match the low-normal-stress contact law of the current 601 source DLL.
    config.material.flac3d7_smrm_contact=true;
    config.material.elastic_modulus_pa=1.0e7;
    config.material.poisson_ratio=0.25;
    config.material.density_kg_m3=2300.0;
    config.material.cohesion_pa=1.0e8;
    config.material.friction_deg=40.0;
    config.material.tensile_strength_pa=1.0e8;
    config.material.joints.push_back({.name="J1",.dip_direction_deg=90.0,
        .dip_deg=90.0,.normal_density_1_m=2.0,.mean_radius_m=1.0,
        .cohesion_pa=1.0e8,.friction_deg=40.0,.normal_stiffness_pa_m=1.0e11,
        .shear_stiffness_pa_m=5.0e10});
    config.gravity_m_s2=9.81;
    config.initial_stress_x_pa=0.0;
    config.initial_stress_y_pa=0.0;
    config.initial_stress_z_pa=0.0;
    config.solver.max_cycles=30000;
    config.solver.convergence_tolerance=1.0e-10;
    config.seepage.enabled=true;
    config.seepage.transient=true;
    config.seepage.solve_mechanics=true;
    config.seepage.bidirectional_coupling=true;
    config.seepage.biot_coefficient=1.0;
    config.seepage.biot_modulus_pa=1.0e7;
    config.seepage.coupling_max_iterations=2000;
    config.seepage.coupling_relative_tolerance=1.0e-6;
    config.seepage.duration_s=duration_s;
    config.seepage.time_step_s=0.25;
    config.seepage.initial_head_m=0.0;
    config.seepage.relative_tolerance=1.0e-13;
    config.seepage.conductivity_m_s={{{1.0e-5,1.0e-5,1.0e-5,0,0,0}}};
    config.seepage.specific_storage_1_m={1.0};
    config.seepage.smrm_stress_permeability=true;
    config.seepage.fluid_kinematic_viscosity_m2_s=1.0e-6;
    config.seepage.joint_hydraulics.push_back({.enabled=true,
        .zero_stress_aperture_m=1.0e-4,.closure_stress_pa=1.0e7,
        .connectivity=1.0,.minimum_aperture_m=1.0e-8,
        .maximum_aperture_m=2.0e-4});
    config.seepage.boundaries={
        {.name="west",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=1.0},
        {.name="east",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::XMax}};
    return rocksolver::import_abaqus_c3d4_text(mesh,config);
}

} // namespace

int main() {
    // The 601 DLL compares shear stress with 1e-12 Pa. The SMRM update
    // works in MPa, so that threshold must be converted before comparison.
    {
        auto material=two_cell_model(0.25).config.material;
        for (const auto [shear_pa,expected_h]:
             {std::pair{1.0e-12,0.0},std::pair{1.0e-10,1.0},
              std::pair{1.0e-8,1.0}}) {
            std::array<double,6> stress_pa{0,0,0,shear_pa,0,0};
            rocksolver::SMRMHistory history;
            rocksolver::SMRMDiagnostics diagnostics;
            rocksolver::Matrix6 tangent{};
            assert(rocksolver::update_smrm_increment(material,1.0,{},
                stress_pa,history,diagnostics,tangent));
            assert(diagnostics.joint_h[0]==expected_h);
        }
    }
    // FLAC3D 7.00 R126, 9-cell mixed mesh under dry elastic gravity.
    // Nonuniform deformation distinguishes independent overlay tets from
    // incorrectly copying one constitutive state to every overlay pair.
    constexpr auto mixed_mesh=R"INP(*NODE
1,0,0,0
2,0,1,0
3,0,1,1
4,0,0,1
5,0.7,0,0
6,0.7,1,0
7,0.7,1,1
8,0.7,0,1
9,1.8,0,0
10,1.8,1,0
11,1.8,1,1
12,1.8,0,1
13,3,0,0
14,3,1,0
15,3,1,1
16,3,0,1
*ELEMENT,TYPE=C3D8,ELSET=HOST
1,1,5,6,2,4,8,7,3
*ELEMENT,TYPE=C3D6,ELSET=HOST
2,5,6,7,9,10,11
3,5,7,8,9,11,12
*ELEMENT,TYPE=C3D4,ELSET=HOST
4,9,13,14,15
5,9,14,10,15
6,9,10,11,15
7,9,11,12,15
8,9,12,16,15
9,9,16,13,15
)INP";
    rocksolver::ModelConfig mixed_config;
    mixed_config.scenario=rocksolver::Scenario::Slope;
    mixed_config.material.type=rocksolver::MaterialType::Elastic;
    mixed_config.material.elastic_modulus_pa=1.0e7;
    mixed_config.material.poisson_ratio=0.25;
    mixed_config.material.density_kg_m3=2300.0;
    mixed_config.initial_stress_x_pa=0.0;
    mixed_config.initial_stress_y_pa=0.0;
    mixed_config.initial_stress_z_pa=0.0;
    mixed_config.solver.max_cycles=30000;
    mixed_config.solver.convergence_tolerance=1.0e-10;
    auto mixed=rocksolver::import_abaqus_c3d4_text(mixed_mesh,mixed_config);
    const auto mixed_result=rocksolver::solve_baseline(mixed);
    assert(mixed_result.converged);
    assert(std::abs(mixed.nodes[2].displacement.z+0.000994041789305677)<1.0e-9);
    assert(std::abs(mixed.nodes[3].displacement.z+0.000889172097223665)<1.0e-9);
    assert(std::abs(mixed.zones[mixed.physical_cells[0].subzones[2]].stress_pa[0]
        +3582.39756886231)<0.05);

    // A pure C3D4 import has no physical_cells overlay. Its Biot step must
    // still use the same per-cell flow-source redistribution and mechanical
    // stress correction as the identical mesh with explicit one-tet cells.
    constexpr auto tetra_mesh=R"INP(*NODE
1,0,0,0
2,1,0,0
3,1,1,0
4,0,1,0
5,0,0,1
6,1,0,1
7,1,1,1
8,0,1,1
*ELEMENT,TYPE=C3D4,ELSET=HOST
1,1,2,3,7
2,1,3,4,7
3,1,4,8,7
4,1,8,5,7
5,1,5,6,7
6,1,6,2,7
)INP";
    auto tetra_config=mixed_config;
    tetra_config.solver.nodal_mixed_discretization=true;
    tetra_config.seepage.enabled=true;
    tetra_config.seepage.transient=true;
    tetra_config.seepage.solve_mechanics=true;
    tetra_config.seepage.bidirectional_coupling=true;
    tetra_config.seepage.biot_coefficient=1.0;
    tetra_config.seepage.biot_modulus_pa=1.0e7;
    tetra_config.seepage.coupling_max_iterations=2000;
    tetra_config.seepage.coupling_relative_tolerance=1.0e-6;
    tetra_config.seepage.duration_s=0.1;
    tetra_config.seepage.time_step_s=0.05;
    tetra_config.seepage.initial_head_m=0.0;
    tetra_config.seepage.relative_tolerance=1.0e-13;
    tetra_config.seepage.conductivity_m_s={{{1.0e-5,1.0e-5,1.0e-5,0,0,0}}};
    tetra_config.seepage.specific_storage_1_m={1.0};
    tetra_config.seepage.boundaries={
        {.name="west",.kind=rocksolver::HydraulicBoundaryKind::Head,
         .selector=rocksolver::HydraulicBoundarySelector::XMin,.value=1.0},
        {.name="east",.kind=rocksolver::HydraulicBoundaryKind::NoFlow,
         .selector=rocksolver::HydraulicBoundarySelector::XMax}};
    auto pure_tetra=rocksolver::import_abaqus_c3d4_text(tetra_mesh,tetra_config);
    assert(pure_tetra.physical_cells.empty());
    auto represented_tetra=pure_tetra;
    for(std::size_t zi=0;zi<represented_tetra.zones.size();++zi){
        rocksolver::PhysicalCell cell;
        cell.id=represented_tetra.zones[zi].id;
        cell.nodes.assign(represented_tetra.zones[zi].nodes.begin(),
            represented_tetra.zones[zi].nodes.end());
        cell.subzones={zi};
        represented_tetra.physical_cells.push_back(cell);
    }
    auto pure_tetra_evolving=pure_tetra;
    auto represented_tetra_evolving=represented_tetra;
    pure_tetra_evolving.config.seepage.evolving_porosity_storage=true;
    represented_tetra_evolving.config.seepage.evolving_porosity_storage=true;
    assert(rocksolver::solve_transient_seepage_mechanics(pure_tetra).converged);
    assert(rocksolver::solve_transient_seepage_mechanics(represented_tetra).converged);
    assert(pure_tetra.hydraulic.history.size()==represented_tetra.hydraulic.history.size());
    for(std::size_t ti=0;ti<pure_tetra.hydraulic.history.size();++ti){
        const auto& a=pure_tetra.hydraulic.history[ti];
        const auto& b=represented_tetra.hydraulic.history[ti];
        assert(a.relative_mass_balance_error<1.0e-6);
        for(std::size_t ni=0;ni<pure_tetra.nodes.size();++ni){
            assert(std::abs(a.nodal_pore_pressure_pa[ni]-b.nodal_pore_pressure_pa[ni])<1.0e-8);
            assert(std::abs(a.nodal_displacement_m[ni].z-b.nodal_displacement_m[ni].z)<1.0e-12);
        }
        for(std::size_t zi=0;zi<pure_tetra.zones.size();++zi)
            for(int axis=0;axis<3;++axis)
                assert(std::abs(a.zone_stress_pa[zi][axis]-b.zone_stress_pa[zi][axis])<1.0e-8);
    }
    assert(rocksolver::solve_transient_seepage_mechanics(pure_tetra_evolving).converged);
    assert(rocksolver::solve_transient_seepage_mechanics(represented_tetra_evolving).converged);
    assert(pure_tetra_evolving.hydraulic.history.size()==
        represented_tetra_evolving.hydraulic.history.size());
    for(std::size_t ti=0;ti<pure_tetra_evolving.hydraulic.history.size();++ti){
        const auto& a=pure_tetra_evolving.hydraulic.history[ti];
        const auto& b=represented_tetra_evolving.hydraulic.history[ti];
        assert(a.relative_mass_balance_error<1.0e-6);
        for(std::size_t ni=0;ni<pure_tetra_evolving.nodes.size();++ni){
            assert(std::abs(a.nodal_pore_pressure_pa[ni]-b.nodal_pore_pressure_pa[ni])<1.0e-8);
            assert(std::abs(a.nodal_displacement_m[ni].z-b.nodal_displacement_m[ni].z)<1.0e-12);
            assert(std::abs(a.nodal_biot_modulus_pa[ni]-b.nodal_biot_modulus_pa[ni])<1.0e-8);
        }
        for(std::size_t zi=0;zi<pure_tetra_evolving.zones.size();++zi){
            assert(std::abs(a.zone_porosity[zi]-b.zone_porosity[zi])<1.0e-12);
            for(int axis=0;axis<3;++axis)
                assert(std::abs(a.zone_stress_pa[zi][axis]-b.zone_stress_pa[zi][axis])<1.0e-8);
        }
    }

    // Dirichlet identity rows must not hide the free flow-equation residual.
    // Exact 1-D solution: h(x)=head-x, qx=K=1e-10 m/s, on two unit bricks.
    for(const double head:{1.0,100000.0}){
        auto flow=two_cell_model(0.5);
        flow.config.seepage.smrm_stress_permeability=false;
        flow.config.seepage.conductivity_m_s={{{1e-10,1e-10,1e-10,0,0,0}}};
        flow.config.seepage.relative_tolerance=1e-10;
        flow.config.seepage.boundaries[0].value=head;
        flow.config.seepage.boundaries[1].kind=rocksolver::HydraulicBoundaryKind::Flux;
        flow.config.seepage.boundaries[1].value=1e-10;
        const auto state=rocksolver::solve_steady_seepage(flow);
        assert(state.converged);
        assert(state.iterations>0);
        for(std::size_t i=0;i<flow.nodes.size();++i)
            assert(std::abs(state.head_m[i]-(head-flow.nodes[i].position.x))<1e-8);
        for(const auto& q:state.zone_darcy_velocity_m_s)
            assert(std::abs(q.x-1e-10)<1e-18);
    }
    // Fast regressions for equilibrium, saturation, input validation and
    // accepted-history preservation. No relaxed physical/relative tolerance.
    for(const auto material:{rocksolver::MaterialType::Elastic,rocksolver::MaterialType::SMRM}){
        for(const double head:{0.0,1.0}){
            auto closed=two_cell_model(0.5);
            closed.config.material.type=material;
            closed.config.seepage.smrm_stress_permeability=material==rocksolver::MaterialType::SMRM;
            closed.config.seepage.initial_head_m=head;
            closed.config.seepage.boundaries[0].kind=rocksolver::HydraulicBoundaryKind::NoFlow;
            const auto closed_result=rocksolver::solve_transient_seepage_mechanics(closed);
            assert(closed_result.converged);
            assert(closed.hydraulic.history.size()==3);
            const auto& start=closed.hydraulic.history.front();
            for(const auto& step:closed.hydraulic.history){
                assert(step.nodal_pore_pressure_pa==start.nodal_pore_pressure_pa);
                assert(step.mass_balance_residual_m3<1.0e-14);
                assert(step.nodal_saturation.size()==closed.nodes.size());
                assert(step.zone_saturation.size()==closed.zones.size());
                for(const auto saturation:step.nodal_saturation)assert(saturation==1.0);
                for(const auto saturation:step.zone_saturation)assert(saturation==1.0);
                for(std::size_t i=0;i<closed.nodes.size();++i){
                    assert(step.nodal_displacement_m[i].x==start.nodal_displacement_m[i].x);
                    assert(step.nodal_displacement_m[i].z==start.nodal_displacement_m[i].z);
                }
            }
        }
    }
    auto incomplete=two_cell_model(0.5);
    incomplete.config.seepage.coupling_max_iterations=1;
    incomplete.config.seepage.boundaries[0].history={{0.0,0.0},{0.25,0.0},{0.5,1.0}};
    const auto incomplete_result=rocksolver::solve_transient_seepage_mechanics(incomplete);
    assert(!incomplete_result.converged);
    assert(!incomplete.hydraulic.converged);
    assert(incomplete.hydraulic.history.size()==2);
    assert(incomplete.hydraulic.time_s==0.25);
    assert(incomplete.hydraulic.history.back().time_s==0.25);
    assert(incomplete.hydraulic.failure_reason.find("t=0.5")!=std::string::npos);
    assert(incomplete.hydraulic.nodal_pore_pressure_pa==
        incomplete.hydraulic.history.back().nodal_pore_pressure_pa);
    for(const double invalid:{0.0,-1.0,std::numeric_limits<double>::quiet_NaN()}){
        auto bad=two_cell_model(invalid);
        bool rejected=false;
        try{(void)rocksolver::solve_transient_seepage_mechanics(bad);}
        catch(const std::invalid_argument&){rejected=true;}
        assert(rejected);
        assert(!bad.hydraulic.valid);
    }
    auto first=two_cell_model(0.25);
    const auto first_result=rocksolver::solve_transient_seepage_mechanics(first);
    assert(first_result.converged);
    for(const auto& zone:first.zones)if(zone.active){
        const auto& joint=first.config.material.joints.front();
        const double normal=std::max(0.0,zone.smrm.joint_normal_stress_mpa[0]);
        const double closed=joint.cohesion_pa*1.0e-6+
            normal*std::tan(joint.friction_deg*std::acos(-1.0)/180.0);
        const double effective=(1.0-zone.smrm.joint_k[0])*closed;
        assert(std::abs(zone.smrm.joint_shear_strength_mpa[0]-effective)<1.0e-9);
    }
    for(const auto* version:{"0.13.0-rainfall-unsaturated","0.13.13-biot-equilibrium-history",
                            "0.13.14-flow-free-residual","0.13.20-json-roundtrip",
                            "0.13.21-node-boundary-flux",
                            "0.13.22-contact-shear-strength"}){
        auto legacy_checkpoint=rocksolver::save_checkpoint(first,first_result);
        legacy_checkpoint.replace(std::string("ROCKSTATE-1 ").size(),
            std::string(rocksolver::kVersion).size(),version);
        rocksolver::RockModel restored;
        rocksolver::RockResult restored_result;
        rocksolver::restore_checkpoint(legacy_checkpoint,restored,restored_result);
        assert(restored.nodes.size()==first.nodes.size());
        assert(restored.zones.front().stress_pa==first.zones.front().stress_pa);
        assert(restored.initialize_interface_stresses==first.initialize_interface_stresses);
        assert(restored.cable_physical_host_interpolation==first.cable_physical_host_interpolation);
        assert(restored.cable_effective_confinement==first.cable_effective_confinement);
    }

    auto model=two_cell_model(2.0);
    const auto result=rocksolver::solve_transient_seepage_mechanics(model);
    assert(result.converged);
    assert(model.hydraulic.history.size()==9);
    assert(model.hydraulic.history.front().zone_conductivity_m_s.size()==model.zones.size());
    assert(model.hydraulic.history.front().zone_joint_aperture_m.size()==model.zones.size());
    assert(model.hydraulic.history.front().zone_joint_normal_stress_mpa.size()==model.zones.size());
    assert(model.hydraulic.history.back().zone_conductivity_m_s==model.hydraulic.zone_conductivity_m_s);
    assert(model.hydraulic.history.front().zone_conductivity_m_s!=
           model.hydraulic.history.back().zone_conductivity_m_s);
    for(std::size_t step=1;step<9;++step)
        assert(model.hydraulic.history[step].relative_mass_balance_error<1.0e-6);
    // FLAC3D 7.00 Release 126, same two bricks and dt=0.25 s at t=2 s.
    // Node ordering here follows the Abaqus mesh above, not FLAC gridpoint IDs.
    constexpr std::array<std::array<double,3>,12> flac_nodes{{
        {{9810,0,0}},{{9810,0,0}},
        {{0,0,-0.000789631366491662}},{{0,0,-0.000789657113441742}},
        {{-2176.94680097695,0,0}},{{-2176.94680097695,0,0}},
        {{-11986.946800977,0.000124925673372421,-0.000914566082214327}},
        {{-11986.946800977,0.000124927822480526,-0.000914575893170984}},
        {{-76.6187923828857,0,0}},{{-76.6187923828857,0,0}},
        {{-9886.61879238286,0,-0.00103949139888825}},
        {{-9886.61879238286,0,-0.00103950407278818}},
    }};
    assert(model.nodes.size()==flac_nodes.size());
    for(std::size_t i=0;i<model.nodes.size();++i){
        assert(std::abs(model.hydraulic.nodal_pore_pressure_pa[i]-flac_nodes[i][0])<20.0);
        assert(std::abs(model.nodes[i].displacement.x-flac_nodes[i][1])<1.0e-7);
        assert(std::abs(model.nodes[i].displacement.z-flac_nodes[i][2])<1.0e-7);
    }
    constexpr std::array<std::array<double,3>,2> flac_stress{{
        {{-4332.58313819254,-3039.47706629339,-11649.1510409579}},
        {{-3597.28105915728,-737.492114091272,-10913.8489623966}},
    }};
    assert(model.physical_cells.size()==flac_stress.size());
    for(std::size_t ci=0;ci<model.physical_cells.size();++ci){
        const auto& cell=model.physical_cells[ci];
        std::array<double,3> stress{};
        double volume=0.0;
        for(const auto zi:cell.subzones){
            const auto& zone=model.zones[zi];
            const double weight=zone.volume_m3*zone.integration_weight;
            for(int axis=0;axis<3;++axis)stress[axis]+=zone.stress_pa[axis]*weight;
            volume+=weight;
        }
        for(int axis=0;axis<3;++axis)
            assert(std::abs(stress[axis]/volume-flac_stress[ci][axis])<5.0);
    }
    assert(std::abs(model.hydraulic.zone_conductivity_m_s.front()[0]-
                    first.hydraulic.zone_conductivity_m_s.front()[0])>1.0e-14);
    assert(std::find(result.warnings.begin(),result.warnings.end(),
        "BIDIRECTIONAL_BIOT_SMRM_LAGGED_PERMEABILITY_FIXED_STORAGE")!=
        result.warnings.end());

    auto evolving=two_cell_model(2.0);
    evolving.config.seepage.evolving_porosity_storage=true;
    evolving.config.seepage.initial_porosity=0.30;
    const auto evolving_result=rocksolver::solve_transient_seepage_mechanics(evolving);
    assert(evolving_result.converged);
    if(std::getenv("ROCKSOLVER_STAGE5_DUMP_JSON")){
        std::cout<<rocksolver::result_to_json(evolving,evolving_result)<<std::endl;
        return 0;
    }
    assert(evolving.hydraulic.history.size()==9);
    const auto& initial=evolving.hydraulic.history.front();
    const auto& final=evolving.hydraulic.history.back();
    assert(initial.nodal_boundary_outward_flux_m3_s.size()==evolving.nodes.size());
    for(const double flux:initial.nodal_boundary_outward_flux_m3_s)assert(flux==0.0);
    bool observed_boundary_flux=false;
    for(std::size_t step=1;step<evolving.hydraulic.history.size();++step){
        const auto& accepted=evolving.hydraulic.history[step];
        assert(accepted.nodal_boundary_outward_flux_m3_s.size()==evolving.nodes.size());
        double net=0.0;
        for(std::size_t ni=0;ni<evolving.nodes.size();++ni){
            const double flux=accepted.nodal_boundary_outward_flux_m3_s[ni];
            net+=flux;
            observed_boundary_flux=observed_boundary_flux||flux!=0.0;
            if(evolving.nodes[ni].position.x!=0.0)assert(flux==0.0);
        }
        const double scale=std::max(accepted.inflow_m3_s+accepted.outflow_m3_s,1.0e-15);
        assert(std::abs(net-accepted.net_outflow_m3_s)<1.0e-12*scale);
    }
    assert(observed_boundary_flux);
    assert(rocksolver::result_to_json(evolving,evolving_result).find(
        "\"nodalBoundaryOutwardFluxM3S\"")!=std::string::npos);
    assert(initial.zone_porosity.size()==evolving.zones.size());
    assert(final.nodal_biot_modulus_pa.size()==evolving.nodes.size());
    for(const double porosity:initial.zone_porosity)assert(std::abs(porosity-0.30)<1.0e-12);
    for(const double modulus:initial.nodal_biot_modulus_pa)
        assert(std::abs(modulus-1.0e7)<1.0e-6);
    for(std::size_t step=1;step<9;++step){
        const auto& previous=evolving.hydraulic.history[step-1];
        const auto& current=evolving.hydraulic.history[step];
        double physical_storage_change=0.0;
        for(std::size_t ni=0;ni<evolving.nodes.size();++ni){
            const double nodal_volume=evolving.nodes[ni].position.x==1.0?0.25:0.125;
            physical_storage_change+=nodal_volume*(
                current.nodal_pore_pressure_pa[ni]/current.nodal_biot_modulus_pa[ni]-
                previous.nodal_pore_pressure_pa[ni]/previous.nodal_biot_modulus_pa[ni]);
        }
        const double boundary_volume=current.net_outflow_m3_s*
            (current.time_s-previous.time_s);
        const double residual=std::abs(physical_storage_change+
            current.mechanical_storage_change_m3+boundary_volume);
        const double scale=std::max({std::abs(physical_storage_change),
            std::abs(current.mechanical_storage_change_m3),
            std::abs(boundary_volume)});
        assert(std::abs(physical_storage_change-current.storage_change_m3)<1.0e-9);
        assert(residual/scale<1.0e-6);
        assert(current.relative_mass_balance_error<1.0e-6);
    }
    // FLAC3D 7.00 R126, conservative n/M FISH reference at t=2 s:
    // D2_storage_conservative_stress_v2.dat (SHA-256 56b5643d9058989227b76aba8bdc5129aa18e041ecb874e6dbfa7fef4151ea7d).
    // The parameter-update substep preserves p/M on free nodes and applies
    // the corresponding total-stress correction before mechanical rebalance.
    // Values below are mapped by coordinates, not FLAC gridpoint IDs.
    constexpr std::array<std::array<double,3>,12> evolving_flac_nodes{{
        {{9810,0,0}},{{9810,0,0}},
        {{0,0,-0.000789615011906799}},{{0,0,-0.00078957661717553}},
        {{-2176.96865276826,0,0}},{{-2176.96865276826,0,0}},
        {{-11984.4028185892,0.000124922962472359,-0.000914524778299803}},
        {{-11984.4028185892,0.000124930249046239,-0.000914520058816692}},
        {{-77.0732741988501,0,0}},{{-77.0732741988501,0,0}},
        {{-9887.02624870282,0,-0.00103945420068757}},
        {{-9887.02624870282,0,-0.00103944384972466}},
    }};
    for(std::size_t i=0;i<evolving.nodes.size();++i){
        assert(std::abs(evolving.hydraulic.nodal_pore_pressure_pa[i]-
            evolving_flac_nodes[i][0])<20.0);
        // Both references conserve fluid content; the staggered parameter
        // update and the coupled Rock step still differ at finite dt.
        assert(std::abs(evolving.nodes[i].displacement.x-
            evolving_flac_nodes[i][1])<2.0e-7);
        assert(std::abs(evolving.nodes[i].displacement.z-
            evolving_flac_nodes[i][2])<2.0e-7);
    }
    constexpr std::array<double,12> evolving_flac_modulus{{
        9994195.58973624,9994195.58973624,9994195.58973624,9994195.58973624,
        9997108.98149048,9997108.98149048,9997108.98149048,9997108.98149048,
        10000024.0722962,10000024.0722962,10000024.0722962,10000024.0722962,
    }};
    for(std::size_t i=0;i<evolving.nodes.size();++i)
        assert(std::abs(final.nodal_biot_modulus_pa[i]-evolving_flac_modulus[i])<100.0);
    assert(std::abs(final.zone_porosity[evolving.physical_cells[0].subzones[0]]-
        0.30017423344015)<1.0e-6);
    assert(std::abs(final.zone_porosity[evolving.physical_cells[1].subzones[0]]-
        0.299999277832852)<1.0e-6);
    constexpr std::array<std::array<double,3>,2> evolving_flac_stress{{
        {{-4332.97125759755,-3039.97456444657,-11649.1506138521}},
        {{-3597.67002865459,-737.992245508128,-10913.8493815776}},
    }};
    for(std::size_t ci=0;ci<evolving.physical_cells.size();++ci){
        const auto& cell=evolving.physical_cells[ci];
        std::array<double,3> stress{};
        double volume=0.0;
        for(const auto zi:cell.subzones){
            const auto& zone=evolving.zones[zi];
            const double weight=zone.volume_m3*zone.integration_weight;
            for(int axis=0;axis<3;++axis)stress[axis]+=zone.stress_pa[axis]*weight;
            volume+=weight;
        }
        for(int axis=0;axis<3;++axis)
            assert(std::abs(stress[axis]/volume-evolving_flac_stress[ci][axis])<5.0);
    }
    assert(std::find(evolving_result.warnings.begin(),evolving_result.warnings.end(),
        "BIDIRECTIONAL_BIOT_INCOMPRESSIBLE_GRAINS_EVOLVING_STORAGE")!=
        evolving_result.warnings.end());
}
