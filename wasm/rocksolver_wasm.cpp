#include "rocksolver/rocksolver.hpp"

#include <sstream>
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <emscripten/emscripten.h>

EM_JS(void, rocksolver_emit_progress,
      (int task_id, const char* stage_name, int stage_index, int stage_count,
       int cycle, int max_cycles, double force_ratio, double force_ratio_local,
       double force_ratio_maximum, double unbalanced_force_maximum_n,
       int stage_converged), {
    if (typeof postMessage !== 'function') return;
    postMessage({
        type: 'rocksolver-progress',
        id: task_id,
        progress: {
            phase: 'solving',
            stageName: UTF8ToString(stage_name),
            stageIndex: stage_index,
            stageCount: stage_count,
            cycle: cycle,
            maxCycles: max_cycles,
            forceRatio: force_ratio,
            forceRatioAverage: force_ratio,
            forceRatioLocal: force_ratio_local,
            forceRatioMaximum: force_ratio_maximum,
            unbalancedForceMaximumN: unbalanced_force_maximum_n,
            converged: Boolean(stage_converged),
        },
    });
});

EM_JS(void, rocksolver_emit_phase,
      (int task_id, const char* phase, const char* label, int percent), {
    if (typeof postMessage !== 'function') return;
    postMessage({
        type: 'rocksolver-progress',
        id: task_id,
        progress: {
            phase: UTF8ToString(phase),
            label: UTF8ToString(label),
            percent: percent,
        },
    });
});

namespace {

std::string last_result;
std::string last_checkpoint;

std::string json_escape(const std::string& value) {
    std::ostringstream output;
    for (const char ch : value) {
        switch (ch) {
            case '\\': output << "\\\\"; break;
            case '"': output << "\\\""; break;
            case '\n': output << "\\n"; break;
            case '\r': output << "\\r"; break;
            case '\t': output << "\\t"; break;
            default: output << ch;
        }
    }
    return output.str();
}

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> parts;
    std::istringstream stream(value);
    std::string part;
    while (std::getline(stream, part, delimiter)) parts.push_back(part);
    return parts;
}

struct WebConfig {
    std::unordered_map<std::string, std::string> values;
    std::vector<std::string> joints;
    std::vector<std::string> regions;
    std::vector<std::string> cables;
    std::vector<std::string> interfaces;
    std::vector<std::string> hydraulic_materials;
    std::vector<std::string> hydraulic_joints;
    std::vector<std::string> hydraulic_boundaries;
    std::vector<std::string> hydraulic_boundary_history;
    std::vector<std::string> hydraulic_storage;
    std::vector<std::string> hydraulic_unsaturated_materials;
};

WebConfig parse_web_config(const std::string& text) {
    WebConfig parsed;
    std::istringstream stream(text);
    std::string line;
    while (std::getline(stream, line)) {
        if (line.empty()) continue;
        const auto separator = line.find('=');
        if (separator == std::string::npos) continue;
        const std::string key = line.substr(0, separator);
        const std::string value = line.substr(separator + 1);
        if (key == "joint") parsed.joints.push_back(value);
        else if (key == "cable") parsed.cables.push_back(value);
        else if (key == "regionMaterial") parsed.regions.push_back(value);
        else if (key == "interface") parsed.interfaces.push_back(value);
        else if (key == "hydraulicMaterial") parsed.hydraulic_materials.push_back(value);
        else if (key == "hydraulicJoint") parsed.hydraulic_joints.push_back(value);
        else if (key == "hydraulicBoundary") parsed.hydraulic_boundaries.push_back(value);
        else if (key == "hydraulicBoundaryHistory") parsed.hydraulic_boundary_history.push_back(value);
        else if (key == "hydraulicStorage") parsed.hydraulic_storage.push_back(value);
        else if (key == "hydraulicUnsaturatedMaterial") parsed.hydraulic_unsaturated_materials.push_back(value);
        else parsed.values[key] = value;
    }
    return parsed;
}

double number(const WebConfig& config, const std::string& key, double fallback) {
    const auto found = config.values.find(key);
    if (found == config.values.end()) return fallback;
    std::size_t end = 0;
    const double value = std::stod(found->second, &end);
    if (end != found->second.size()) {
        throw std::invalid_argument("invalid numeric Web config: " + key);
    }
    return value;
}

int integer(const WebConfig& config, const std::string& key, int fallback) {
    const auto found = config.values.find(key);
    if (found == config.values.end()) return fallback;
    std::size_t end = 0;
    const int value = std::stoi(found->second, &end);
    if (end != found->second.size()) {
        throw std::invalid_argument("invalid integer Web config: " + key);
    }
    return value;
}

rocksolver::JointSet parse_joint(const std::string& value) {
    const auto parts = split(value, '|');
    if (parts.size() != 10) {
        throw std::invalid_argument("Web joint requires 10 pipe-separated fields");
    }
    rocksolver::JointSet joint;
    joint.name = parts[0];
    joint.dip_direction_deg = std::stod(parts[1]);
    joint.dip_deg = std::stod(parts[2]);
    joint.normal_density_1_m = std::stod(parts[3]);
    joint.mean_radius_m = std::stod(parts[4]);
    joint.cohesion_pa = std::stod(parts[5]) * 1.0e6;
    joint.friction_deg = std::stod(parts[6]);
    joint.tensile_strength_pa = std::stod(parts[7]) * 1.0e6;
    joint.normal_stiffness_pa_m = std::stod(parts[8]) * 1.0e9;
    joint.shear_stiffness_pa_m = std::stod(parts[9]) * 1.0e9;
    return joint;
}

rocksolver::HydraulicBoundaryKind hydraulic_boundary_kind(const std::string& value) {
    if(value=="head")return rocksolver::HydraulicBoundaryKind::Head;
    if(value=="pressure")return rocksolver::HydraulicBoundaryKind::PorePressure;
    if(value=="flux")return rocksolver::HydraulicBoundaryKind::Flux;
    if(value=="rainfall")return rocksolver::HydraulicBoundaryKind::Rainfall;
    if(value=="drain")return rocksolver::HydraulicBoundaryKind::Drain;
    if(value=="seepage-face")return rocksolver::HydraulicBoundaryKind::SeepageFace;
    if(value=="no-flow")return rocksolver::HydraulicBoundaryKind::NoFlow;
    throw std::invalid_argument("invalid hydraulic boundary kind");
}

rocksolver::HydraulicBoundarySelector hydraulic_boundary_selector(const std::string& value) {
    if(value=="x-min")return rocksolver::HydraulicBoundarySelector::XMin;
    if(value=="x-max")return rocksolver::HydraulicBoundarySelector::XMax;
    if(value=="y-min")return rocksolver::HydraulicBoundarySelector::YMin;
    if(value=="y-max")return rocksolver::HydraulicBoundarySelector::YMax;
    if(value=="z-min")return rocksolver::HydraulicBoundarySelector::ZMin;
    if(value=="z-max")return rocksolver::HydraulicBoundarySelector::ZMax;
    if(value=="plane")return rocksolver::HydraulicBoundarySelector::Plane;
    throw std::invalid_argument("invalid hydraulic boundary selector");
}

rocksolver::ModelConfig model_config(const std::string& text) {
    const WebConfig web = parse_web_config(text);
    rocksolver::ModelConfig config;
    const auto scenario = web.values.find("scenario");
    config.scenario = scenario != web.values.end() && scenario->second == "slope"
                          ? rocksolver::Scenario::Slope
                          : rocksolver::Scenario::Tunnel;
    const auto material = web.values.find("material");
    const std::string material_name = material == web.values.end() ? "smrm" : material->second;
    if (material_name == "elastic") config.material.type = rocksolver::MaterialType::Elastic;
    else if (material_name == "mohr-coulomb" || material_name == "mohr" || material_name == "mc")
        config.material.type = rocksolver::MaterialType::MohrCoulomb;
    else if (material_name == "smrm") config.material.type = rocksolver::MaterialType::SMRM;
    else throw std::invalid_argument("invalid constitutive material");
    config.material.density_kg_m3 =
        number(web, "densityKgM3", config.material.density_kg_m3);
    config.material.elastic_modulus_pa =
        number(web, "elasticModulusGpa", 27.0) * 1.0e9;
    config.material.poisson_ratio =
        number(web, "poissonRatio", config.material.poisson_ratio);
    config.material.cohesion_pa = number(web, "cohesionMpa", 5.0) * 1.0e6;
    config.material.friction_deg =
        number(web, "frictionDeg", config.material.friction_deg);
    config.material.dilation_deg =
        number(web, "dilationDeg", config.material.dilation_deg);
    config.material.tensile_strength_pa =
        number(web, "tensileStrengthMpa", 2.0) * 1.0e6;
    config.material.ucs_pa = number(web, "ucsMpa", 50.0) * 1.0e6;
    config.material.fracture_toughness_mpa_sqrt_m =
        number(web, "fractureToughnessMpaSqrtM", 0.0);
    config.material.representative_volume_m3 =
        number(web, "representativeVolumeM3", 1.0);
    config.material.weibull_shape = number(web, "weibullShape", 4.0);
    config.material.strength_reduction_factor =
        number(web, "strengthReductionFactor", 1.0);
    const bool flac3d700_compatibility =
        integer(web, "flac3d700Compatibility", 0) != 0;
    config.material.flac3d7_smrm_contact = flac3d700_compatibility;
    if (config.material.type == rocksolver::MaterialType::MohrCoulomb &&
        config.material.dilation_deg > config.material.friction_deg) {
        throw std::invalid_argument("Mohr-Coulomb dilation angle cannot exceed friction angle");
    }
    config.material.joints.reserve(web.joints.size());
    for (const std::string& joint : web.joints) {
        config.material.joints.push_back(parse_joint(joint));
    }
    std::set<std::string> region_names;
    for (const auto& row : web.regions) {
        const auto parts = split(row, '|');
        if ((parts.size() != 12 && parts.size() != 13) || parts[0].rfind("SLOPE_BAND_", 0) != 0 || !region_names.insert(parts[0]).second)
            throw std::invalid_argument("invalid or duplicate region material");
        rocksolver::RegionMaterial region;
        region.name = parts[0];
        auto& m = region.material;
        if (parts[1] == "elastic") m.type = rocksolver::MaterialType::Elastic;
        else if (parts[1] == "mohr-coulomb" || parts[1] == "mohr" || parts[1] == "mc") m.type = rocksolver::MaterialType::MohrCoulomb;
        else if (parts[1] == "smrm") m.type = rocksolver::MaterialType::SMRM;
        else throw std::invalid_argument("invalid region constitutive model");
        std::vector<double> values;
        for (std::size_t i = 2; i < parts.size(); ++i) {
            std::size_t end = 0;
            double value = std::stod(parts[i], &end);
            if (end != parts[i].size() || !std::isfinite(value)) throw std::invalid_argument("invalid region material number");
            values.push_back(value);
        }
        const bool with_dilation = values.size() == 11;
        m.density_kg_m3 = values[0];
        m.elastic_modulus_pa = values[1] * 1e9;
        m.poisson_ratio = values[2];
        m.cohesion_pa = values[3] * 1e6;
        m.friction_deg = values[4];
        m.dilation_deg = with_dilation ? values[5] : 0.0;
        const std::size_t offset = with_dilation ? 1 : 0;
        m.tensile_strength_pa = values[5 + offset] * 1e6;
        m.ucs_pa = values[6 + offset] * 1e6;
        m.fracture_toughness_mpa_sqrt_m = values[7 + offset];
        m.representative_volume_m3 = values[8 + offset];
        m.weibull_shape = values[9 + offset];
        m.strength_reduction_factor = config.material.strength_reduction_factor;
        m.flac3d7_smrm_contact = flac3d700_compatibility;
        if (m.density_kg_m3 <= 0 || m.elastic_modulus_pa <= 0 || m.poisson_ratio <= -1 || m.poisson_ratio >= 0.5 ||
            m.cohesion_pa < 0 || m.friction_deg < 0 || m.friction_deg >= 90 || m.dilation_deg < 0 || m.dilation_deg >= 90 ||
            m.tensile_strength_pa < 0 || m.ucs_pa <= 0 || m.fracture_toughness_mpa_sqrt_m < 0 ||
            m.representative_volume_m3 <= 0 || m.weibull_shape <= 0 ||
            (m.type == rocksolver::MaterialType::MohrCoulomb && m.dilation_deg > m.friction_deg))
            throw std::invalid_argument("region material out of range");
        config.region_materials.push_back(region);
    }
    config.gravity_m_s2 = number(web, "gravityMS2", config.gravity_m_s2);
    config.water.enabled = integer(web, "waterEnabled", 0) != 0;
    config.water.table_z_m = number(web, "waterTableZM", 0.0);
    config.water.density_kg_m3 = number(web, "waterDensityKgM3", 1000.0);
    config.water.plane_enabled = integer(web, "waterPlaneEnabled", 0) != 0;
    config.water.plane_origin_m = {number(web, "waterPlaneOriginX", 0.0),
                                   number(web, "waterPlaneOriginY", 0.0),
                                   number(web, "waterPlaneOriginZ", 0.0)};
    config.water.plane_normal = {number(web, "waterPlaneNormalX", 0.0),
                                 number(web, "waterPlaneNormalY", 0.0),
                                 number(web, "waterPlaneNormalZ", 1.0)};
    config.water.auto_saturated_density = integer(web, "waterAutoSaturatedDensity", 0) != 0;
    if(config.water.auto_saturated_density) {
        const auto found=web.values.find("waterSaturatedDensitiesKgM3");
        if(found==web.values.end())throw std::invalid_argument("missing saturated densities");
        for(const auto& part:split(found->second,',')) {
            std::size_t end=0;
            const double density=std::stod(part,&end);
            if(end!=part.size())throw std::invalid_argument("invalid saturated density");
            config.water.saturated_densities_kg_m3.push_back(density);
        }
    }
    config.water.external_enabled = integer(web, "externalWaterEnabled", 0) != 0;
    config.water.external_surface_z_m = number(web, "externalWaterSurfaceZM", 0.0);
    config.water.external_density_kg_m3 = number(web, "externalWaterDensityKgM3", 1000.0);
    config.water.external_outward_direction = {
        number(web, "externalWaterOutwardX", 0.0),
        number(web, "externalWaterOutwardY", 0.0),
        number(web, "externalWaterOutwardZ", 1.0)};
    config.water.external_min_alignment = number(web, "externalWaterMinAlignment", 0.5);
    config.seepage.enabled=integer(web,"seepageEnabled",0)!=0;
    config.seepage.solve_mechanics=integer(web,"seepageSolveMechanics",0)!=0;
    config.seepage.fluid_density_kg_m3=number(web,"seepageFluidDensityKgM3",1000.0);
    config.seepage.smrm_stress_permeability=integer(web,"seepageSmrmStressPermeability",0)!=0;
    config.seepage.fluid_kinematic_viscosity_m2_s=number(web,"seepageKinematicViscosityM2S",1.01e-6);
    config.seepage.max_iterations=integer(web,"seepageMaxIterations",20000);
    config.seepage.relative_tolerance=number(web,"seepageRelativeTolerance",1e-10);
    config.seepage.transient=integer(web,"seepageTransient",0)!=0;
    config.seepage.mechanics_each_step=integer(web,"seepageMechanicsEachStep",0)!=0;
    config.seepage.duration_s=number(web,"seepageDurationS",0.0);
    config.seepage.time_step_s=number(web,"seepageTimeStepS",0.0);
    config.seepage.initial_head_m=number(web,"seepageInitialHeadM",0.0);
    config.seepage.biot_coefficient=number(web,"seepageBiotCoefficient",1.0);
    config.seepage.biot_modulus_pa=number(web,"seepageBiotModulusPa",0.0);
    config.seepage.bidirectional_coupling=integer(web,"seepageBidirectionalCoupling",0)!=0;
    config.seepage.evolving_porosity_storage=integer(web,"seepageEvolvingPorosityStorage",0)!=0;
    config.seepage.initial_porosity=number(web,"seepageInitialPorosity",0.30);
    config.seepage.coupling_max_iterations=integer(web,"seepageCouplingMaxIterations",200);
    config.seepage.coupling_relative_tolerance=number(web,"seepageCouplingRelativeTolerance",1e-7);
    config.seepage.unsaturated=integer(web,"seepageUnsaturated",0)!=0;
    config.seepage.flac_no_capillarity=integer(web,"seepageFlacNoCapillarity",0)!=0;
    config.seepage.initial_saturation=number(web,"seepageInitialSaturation",1.0);
    config.seepage.nonlinear_max_iterations=integer(web,"seepageNonlinearMaxIterations",50);
    config.seepage.nonlinear_relative_tolerance=number(web,"seepageNonlinearRelativeTolerance",1e-8);
    for(const auto& row:web.hydraulic_materials) {
        const auto fields=split(row,'|');
        if(fields.size()!=6)throw std::invalid_argument("hydraulic material requires six conductivity components");
        std::array<double,6> tensor{};
        for(std::size_t i=0;i<tensor.size();++i) {
            std::size_t end=0;tensor[i]=std::stod(fields[i],&end);
            if(end!=fields[i].size()||!std::isfinite(tensor[i]))throw std::invalid_argument("invalid hydraulic conductivity component");
        }
        config.seepage.conductivity_m_s.push_back(tensor);
    }
    for(const auto& row:web.hydraulic_joints) {
        const auto fields=split(row,'|');
        if(fields.size()!=6)throw std::invalid_argument("hydraulic joint requires six fields");
        rocksolver::JointHydraulicParameters hydraulic;
        hydraulic.enabled=std::stoi(fields[0])!=0;
        hydraulic.zero_stress_aperture_m=std::stod(fields[1]);
        hydraulic.closure_stress_pa=std::stod(fields[2]);
        hydraulic.connectivity=std::stod(fields[3]);
        hydraulic.minimum_aperture_m=std::stod(fields[4]);
        hydraulic.maximum_aperture_m=std::stod(fields[5]);
        config.seepage.joint_hydraulics.push_back(hydraulic);
    }
    for(const auto& row:web.hydraulic_boundaries) {
        const auto fields=split(row,'|');
        if(fields.size()!=11)throw std::invalid_argument("hydraulic boundary requires eleven fields");
        rocksolver::HydraulicBoundary boundary;
        boundary.name=fields[0];
        boundary.kind=hydraulic_boundary_kind(fields[1]);
        boundary.selector=hydraulic_boundary_selector(fields[2]);
        boundary.value=std::stod(fields[3]);
        boundary.plane_origin_m={std::stod(fields[4]),std::stod(fields[5]),std::stod(fields[6])};
        boundary.plane_normal={std::stod(fields[7]),std::stod(fields[8]),std::stod(fields[9])};
        boundary.tolerance_m=std::stod(fields[10]);
        config.seepage.boundaries.push_back(std::move(boundary));
    }
    for(const auto& row:web.hydraulic_boundary_history) {
        const auto fields=split(row,'|');
        if(fields.size()!=3)throw std::invalid_argument("hydraulic boundary history requires boundary index, time, and value");
        const auto boundary_index=static_cast<std::size_t>(std::stoul(fields[0]));
        if(boundary_index>=config.seepage.boundaries.size())
            throw std::invalid_argument("hydraulic boundary history index is out of range");
        config.seepage.boundaries[boundary_index].history.push_back(
            {std::stod(fields[1]),std::stod(fields[2])});
    }
    for(const auto& row:web.hydraulic_storage)
        config.seepage.specific_storage_1_m.push_back(std::stod(row));
    for(const auto& row:web.hydraulic_unsaturated_materials) {
        const auto fields=split(row,'|');
        if(fields.size()!=5)throw std::invalid_argument("unsaturated hydraulic material requires five fields");
        rocksolver::UnsaturatedHydraulicParameters material;
        material.porosity=std::stod(fields[0]);
        material.residual_saturation=std::stod(fields[1]);
        material.van_genuchten_alpha_1_m=std::stod(fields[2]);
        material.van_genuchten_n=std::stod(fields[3]);
        material.mualem_l=std::stod(fields[4]);
        config.seepage.unsaturated_materials.push_back(material);
    }
    config.nx = integer(web, "nx", config.nx);
    config.ny = integer(web, "ny", config.ny);
    config.nz = integer(web, "nz", config.nz);
    config.length_m = number(web, "lengthM", config.length_m);
    config.width_m = number(web, "widthM", config.width_m);
    config.height_m = number(web, "heightM", config.height_m);
    config.slope_angle_deg = number(web, "slopeAngleDeg", config.slope_angle_deg);
    config.initial_stress_x_pa = number(web, "initialXMpa", -2.0) * 1.0e6;
    config.initial_stress_y_pa = number(web, "initialYMpa", -3.0) * 1.0e6;
    config.initial_stress_z_pa = number(web, "initialZMpa", -4.0) * 1.0e6;
    config.tunnel_radius_x_m = number(web, "tunnelRadiusXM", 4.0);
    config.tunnel_radius_z_m = number(web, "tunnelRadiusZM", 4.0);
    config.excavation_steps = integer(web, "excavationSteps", 0);
    config.support_stiffness_n_m = number(web, "supportStiffnessMnM", 0.0) * 1.0e6;
    config.solver.max_cycles = integer(web, "maxCycles", 3000);
    config.solver.convergence_tolerance =
        number(web, "convergenceTolerance", 1.0e-3);
    config.solver.damping = number(web, "damping", 0.8);
    config.solver.time_step_scale = flac3d700_compatibility
        ? 1.0 : number(web, "timeStepScale", 0.26);
    config.solver.flac3d7_static_mass = flac3d700_compatibility;
    config.solver.large_strain = integer(web, "largeStrain", 0) != 0;
    config.solver.nodal_mixed_discretization =
        integer(web, "nodalMixedDiscretization", 0) != 0;
    config.solver.rigid_top_platen_z =
        integer(web, "rigidTopPlatenZ", 0) != 0;
    config.strength_reduction_interfaces_only =
        integer(web, "strengthReductionInterfacesOnly", 0) != 0;
    config.strength_reduction_continuum_only =
        integer(web, "strengthReductionContinuumOnly", 0) != 0;
    // Both flags set is the persisted encoding for region-material-bands-only reduction.
    for (const auto& line : web.interfaces) {
        const auto fields=split(line,'|');
        if(fields.size()!=12) throw std::invalid_argument("invalid interface record");
        rocksolver::InterfaceConfig value;
        value.name=fields[0];
        value.enabled=std::stoi(fields[1])!=0;
        value.source_group=fields[2];
        value.moving_side_hint={std::stod(fields[3]),std::stod(fields[4]),std::stod(fields[5])};
        value.normal_stiffness_pa_m=std::stod(fields[6]);
        value.shear_stiffness_pa_m=std::stod(fields[7]);
        value.cohesion_pa=std::stod(fields[8]);
        value.friction_deg=std::stod(fields[9]);
        value.tensile_strength_pa=std::stod(fields[10]);
        value.replace_source_with_host=std::stoi(fields[11])!=0;
        config.interfaces.push_back(std::move(value));
    }
    for (const auto& line : web.cables) {
        std::vector<std::string> fields;
        std::istringstream input(line); std::string field;
        while (std::getline(input, field, '|')) fields.push_back(field);
        if (fields.size() != 12 && fields.size()!=14) throw std::invalid_argument("invalid cable record");
        rocksolver::CableConfig c;
        c.id = std::stoi(fields[0]);
        if (fields[1] != "before" && fields[1] != "after") throw std::invalid_argument("invalid cable stage");
        c.after_initial = fields[1] == "after";
        c.young=std::stod(fields[2]); c.area=std::stod(fields[3]);
        c.yield_tension=std::stod(fields[4]); c.yield_compression=std::stod(fields[5]);
        c.grout_stiffness=std::stod(fields[6]); c.grout_cohesion=std::stod(fields[7]);
        c.grout_friction=std::stod(fields[8]); c.perimeter=std::stod(fields[9]); c.density=std::stod(fields[10]);
        if(fields.size()==14){c.pretension=std::stod(fields[12]);c.free_length=std::stod(fields[13]);}
        std::istringstream points(fields[11]);
        while(std::getline(points,field,';')) {
            std::replace(field.begin(),field.end(),',',' ');
            std::istringstream xyz(field); rocksolver::Vec3 point;
            if(!(xyz>>point.x>>point.y>>point.z)) throw std::invalid_argument("invalid cable point");
            c.points.push_back(point);
        }
        config.cables.push_back(c);
    }
    return config;
}

void set_error(const std::exception& error) {
    last_checkpoint.clear();
    last_result = "{\"error\":\"" + json_escape(error.what()) + "\"}";
}

void store_state(const rocksolver::RockModel& model,
                 const rocksolver::RockResult& result, int task_id) {
    rocksolver_emit_phase(task_id,"serializing","正在生成恢复状态",96);
    last_checkpoint=rocksolver::save_checkpoint(model,result);
    rocksolver_emit_phase(task_id,"serializing","正在生成计算结果",98);
    last_result=rocksolver::result_to_json(model,result);
}

rocksolver::ProgressCallback web_progress(int task_id) {
    return [task_id](const rocksolver::SolveProgress& progress) {
        rocksolver_emit_progress(
            task_id, progress.stage_name.c_str(), progress.stage_index,
            progress.stage_count, progress.cycle, progress.max_cycles,
            progress.force_ratio, progress.force_ratio_local,
            progress.force_ratio_maximum,
            progress.unbalanced_force_maximum_n,
            progress.stage_converged ? 1 : 0);
    };
}

}  // namespace

extern "C" {

EMSCRIPTEN_KEEPALIVE const char* rocksolver_version() { return rocksolver::kVersion; }

EMSCRIPTEN_KEEPALIVE int rocksolver_run_preset(int scenario, int material, int nx,
                                               int ny, int nz, int task_id) {
    try {
        last_checkpoint.clear();
        rocksolver::ModelConfig config;
        config.scenario = scenario == 0 ? rocksolver::Scenario::Slope
                                        : rocksolver::Scenario::Tunnel;
        config.material.type = material == 0 ? rocksolver::MaterialType::Elastic
                             : material == 2 ? rocksolver::MaterialType::MohrCoulomb
                                             : rocksolver::MaterialType::SMRM;
        config.nx = nx;
        config.ny = ny;
        config.nz = nz;
        auto model = rocksolver::generate_model(config);
        const auto result = rocksolver::solve(model, web_progress(task_id));
        last_result = rocksolver::result_to_json(model, result);
        return result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        set_error(error);
        return -1;
    }
}

EMSCRIPTEN_KEEPALIVE int rocksolver_run_config(const char* config_text,
                                               int task_id) {
    try {
        last_checkpoint.clear();
        if (config_text == nullptr)
            throw std::invalid_argument("WASM model config is missing");
        const auto config = model_config(config_text);
        auto model = rocksolver::generate_model(config);
        const auto result = rocksolver::solve(model, web_progress(task_id));
        last_result = rocksolver::result_to_json(model, result);
        return result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        set_error(error);
        return -1;
    }
}

EMSCRIPTEN_KEEPALIVE int rocksolver_run_imported(const char* inp_text,
                                                 const char* config_text,
                                                 int task_id) {
    try {
        last_checkpoint.clear();
        if (inp_text == nullptr || config_text == nullptr) {
            throw std::invalid_argument("WASM imported model input is missing");
        }
        const auto config = model_config(config_text);
        auto model = rocksolver::import_abaqus_c3d4_text(inp_text, config);
        const auto result = rocksolver::solve(model, web_progress(task_id));
        last_result = rocksolver::result_to_json(model, result);
        return result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        set_error(error);
        return -1;
    }
}

EMSCRIPTEN_KEEPALIVE int rocksolver_run_steady_seepage(const char* inp_text,
                                                       const char* config_text,
                                                       int task_id) {
    try {
        last_checkpoint.clear();
        if(!inp_text||!config_text)throw std::invalid_argument("steady seepage input is missing");
        auto model=rocksolver::import_abaqus_c3d4_text(inp_text,model_config(config_text));
        const bool transient=model.config.seepage.transient;
        rocksolver_emit_phase(task_id,"seepage",transient?"正在求解瞬态饱和渗流":"正在求解稳态饱和渗流",10);
        const auto result=transient
            ?rocksolver::solve_transient_seepage_mechanics(model,web_progress(task_id))
            :rocksolver::solve_steady_seepage_mechanics(model,web_progress(task_id));
        rocksolver_emit_phase(task_id,"serializing","正在生成渗流结果",98);
        last_result=rocksolver::result_to_json(model,result);
        return result.converged?0:1;
    }catch(const std::exception& error){set_error(error);return -1;}
}

EMSCRIPTEN_KEEPALIVE int rocksolver_run_baseline(const char* inp, const char* config_text, int task_id) {
    try {
        if(!inp || !config_text) throw std::invalid_argument("缺少基准模型");
        auto model=rocksolver::import_abaqus_c3d4_text(inp,model_config(config_text));
        auto result=rocksolver::solve_baseline(model,web_progress(task_id));
        store_state(model,result,task_id);
        return result.converged ? 0 : 1;
    } catch(const std::exception& e) {set_error(e);return -1;}
}
EMSCRIPTEN_KEEPALIVE int rocksolver_run_initial_state(const char* inp, const char* config_text, int task_id) {
    try {
        if(!inp || !config_text) throw std::invalid_argument("缺少基准模型");
        auto model=rocksolver::import_abaqus_c3d4_text(inp,model_config(config_text));
        auto result=rocksolver::solve(model,web_progress(task_id));
        store_state(model,result,task_id);
        return result.converged ? 0 : 1;
    } catch(const std::exception& e) {set_error(e);return -1;}
}
EMSCRIPTEN_KEEPALIVE int rocksolver_reduction_trial(const char* checkpoint,double factor,int task_id){
    try{
        last_checkpoint.clear();
        rocksolver::RockModel model;rocksolver::RockResult saved;
        rocksolver::restore_checkpoint(checkpoint,model,saved);
        if(!saved.converged)throw std::invalid_argument("稳定性分析要求已收敛的基准");
        auto result=rocksolver::reduction_trial(model,factor,web_progress(task_id));
        last_result=rocksolver::result_to_json(model,result);
        return result.converged?0:1;
    }catch(const std::exception& e){set_error(e);return -1;}
}
EMSCRIPTEN_KEEPALIVE int rocksolver_resume(const char* checkpoint, const char* config_text, int allow_unconverged, int task_id) {
    try {
        if(!checkpoint || !config_text) throw std::invalid_argument("缺少恢复状态");
        rocksolver::RockModel model;rocksolver::RockResult baseline;
        rocksolver::restore_checkpoint(checkpoint,model,baseline);
        auto result=rocksolver::continue_slope(model,baseline,model_config(config_text),allow_unconverged!=0,web_progress(task_id));
        store_state(model,result,task_id);
        return result.converged ? 0 : 1;
    } catch(const std::exception& e) {set_error(e);return -1;}
}

EMSCRIPTEN_KEEPALIVE int rocksolver_water_stage(const char* checkpoint,
                                                const char* config_text,
                                                int task_id) {
    try {
        if (!checkpoint || !config_text)
            throw std::invalid_argument("缺少水位阶段恢复状态或配置");
        rocksolver::RockModel model;
        rocksolver::RockResult previous;
        rocksolver::restore_checkpoint(checkpoint, model, previous);
        const auto config = model_config(config_text);
        auto result = rocksolver::continue_water_stage(
            model, previous, config.water, config.solver, web_progress(task_id));
        store_state(model, result, task_id);
        return result.converged ? 0 : 1;
    } catch (const std::exception& error) {
        set_error(error);
        return -1;
    }
}

EMSCRIPTEN_KEEPALIVE const char* rocksolver_result_ptr() { return last_result.c_str(); }

EMSCRIPTEN_KEEPALIVE int rocksolver_result_size() {
    return static_cast<int>(last_result.size());
}

EMSCRIPTEN_KEEPALIVE const char* rocksolver_checkpoint_ptr() {
    return last_checkpoint.c_str();
}

EMSCRIPTEN_KEEPALIVE int rocksolver_checkpoint_size() {
    return static_cast<int>(last_checkpoint.size());
}

}
