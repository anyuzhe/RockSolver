#include "rocksolver/rocksolver.hpp"

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

std::string require_value(int& index, int argc, char** argv) {
    if (index + 1 >= argc) throw std::invalid_argument("missing value for " + std::string(argv[index]));
    return argv[++index];
}

int as_int(const std::string& value, const std::string& name) {
    try {
        std::size_t end = 0;
        const int parsed = std::stoi(value, &end);
        if (end != value.size()) throw std::invalid_argument("trailing");
        return parsed;
    } catch (...) {
        throw std::invalid_argument(name + " must be an integer");
    }
}

double as_double(const std::string& value, const std::string& name) {
    try {
        std::size_t end = 0;
        const double parsed = std::stod(value, &end);
        if (end != value.size()) throw std::invalid_argument("trailing");
        return parsed;
    } catch (...) {
        throw std::invalid_argument(name + " must be numeric");
    }
}

std::vector<std::string> split(const std::string& value, char delimiter) {
    std::vector<std::string> parts;
    std::istringstream stream(value);
    std::string part;
    while (std::getline(stream, part, delimiter)) parts.push_back(part);
    return parts;
}

rocksolver::JointSet parse_joint(const std::string& value) {
    const auto parts = split(value, ',');
    if (parts.size() != 10) {
        throw std::invalid_argument(
            "--joint needs name,dipDirection,dip,density,radius,cohesionMPa,friction,tensileMPa,knGPaPerM,ksGPaPerM");
    }
    rocksolver::JointSet joint;
    joint.name = parts[0];
    joint.dip_direction_deg = as_double(parts[1], "joint dip direction");
    joint.dip_deg = as_double(parts[2], "joint dip");
    joint.normal_density_1_m = as_double(parts[3], "joint density");
    joint.mean_radius_m = as_double(parts[4], "joint radius");
    joint.cohesion_pa = as_double(parts[5], "joint cohesion") * 1.0e6;
    joint.friction_deg = as_double(parts[6], "joint friction");
    joint.tensile_strength_pa = as_double(parts[7], "joint tensile strength") * 1.0e6;
    joint.normal_stiffness_pa_m = as_double(parts[8], "joint normal stiffness") * 1.0e9;
    joint.shear_stiffness_pa_m = as_double(parts[9], "joint shear stiffness") * 1.0e9;
    return joint;
}

void print_help() {
    std::cout
        << "RockSolver native — FLAC-compatible SMRM core\n"
        << "  --scenario slope|tunnel --material elastic|mohr-coulomb|smrm\n"
        << "  --nx N --ny N --nz N --length M --width M --height M\n"
        << "  --e-gpa GPa --density kg/m3 --nu value --cohesion-mpa MPa\n"
        << "  --friction-deg deg --dilation-deg deg --tensile-mpa MPa --ucs-mpa MPa\n"
        << "  --kic-mpa-sqrt-m value --weibull value --representative-volume-m3 value\n"
        << "  --joint name,dipdir,dip,density,radius,c,phi,tension,kn,ks\n"
        << "  --inp FILE (use an existing Abaqus C3D4 mesh)\n"
        << "  --output FILE (write result JSON instead of stdout) --progress\n"
        << "  --nodal-mixed-discretization (FLAC-compatible Tet4 volumetric averaging)\n"
        << "  --estimate (return resource estimate without solving)\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        rocksolver::ModelConfig config;
        bool estimate_only = false;
        bool show_progress = false;
        std::string inp_path;
        std::string output_path;
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                print_help();
                return 0;
            } else if (arg == "--estimate") {
                estimate_only = true;
            } else if (arg == "--progress") {
                show_progress = true;
            } else if (arg == "--nodal-mixed-discretization") {
                config.solver.nodal_mixed_discretization = true;
            } else if (arg == "--output") {
                output_path = require_value(i, argc, argv);
            } else if (arg == "--inp") {
                inp_path = require_value(i, argc, argv);
            } else if (arg == "--scenario") {
                const std::string value = require_value(i, argc, argv);
                if (value == "slope") config.scenario = rocksolver::Scenario::Slope;
                else if (value == "tunnel") config.scenario = rocksolver::Scenario::Tunnel;
                else throw std::invalid_argument("scenario must be slope or tunnel");
            } else if (arg == "--material") {
                const std::string value = require_value(i, argc, argv);
                if (value == "elastic") config.material.type = rocksolver::MaterialType::Elastic;
                else if (value == "mohr-coulomb" || value == "mohr" || value == "mc") config.material.type = rocksolver::MaterialType::MohrCoulomb;
                else if (value == "smrm") config.material.type = rocksolver::MaterialType::SMRM;
                else throw std::invalid_argument("material must be elastic, mohr-coulomb or smrm");
            } else if (arg == "--nx") config.nx = as_int(require_value(i, argc, argv), arg);
            else if (arg == "--ny") config.ny = as_int(require_value(i, argc, argv), arg);
            else if (arg == "--nz") config.nz = as_int(require_value(i, argc, argv), arg);
            else if (arg == "--length") config.length_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--width") config.width_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--height") config.height_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--slope-angle") config.slope_angle_deg = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--tunnel-rx") config.tunnel_radius_x_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--tunnel-rz") config.tunnel_radius_z_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--initial-x-mpa") config.initial_stress_x_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--initial-y-mpa") config.initial_stress_y_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--initial-z-mpa") config.initial_stress_z_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--gravity") config.gravity_m_s2 = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--density") config.material.density_kg_m3 = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--e-gpa") config.material.elastic_modulus_pa = as_double(require_value(i, argc, argv), arg) * 1.0e9;
            else if (arg == "--nu") config.material.poisson_ratio = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--cohesion-mpa") config.material.cohesion_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--friction-deg") config.material.friction_deg = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--dilation-deg") config.material.dilation_deg = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--tensile-mpa") config.material.tensile_strength_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--ucs-mpa") config.material.ucs_pa = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else if (arg == "--kic-mpa-sqrt-m") config.material.fracture_toughness_mpa_sqrt_m = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--weibull") config.material.weibull_shape = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--representative-volume-m3") config.material.representative_volume_m3 = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--strength-reduction") config.material.strength_reduction_factor = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--joint") config.material.joints.push_back(parse_joint(require_value(i, argc, argv)));
            else if (arg == "--cycles") config.solver.max_cycles = as_int(require_value(i, argc, argv), arg);
            else if (arg == "--tolerance") config.solver.convergence_tolerance = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--damping") config.solver.damping = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--time-step-scale") config.solver.time_step_scale = as_double(require_value(i, argc, argv), arg);
            else if (arg == "--large-strain") config.solver.large_strain = true;
            else if (arg == "--excavation-steps") config.excavation_steps = as_int(require_value(i, argc, argv), arg);
            else if (arg == "--support-stiffness-mn-m") config.support_stiffness_n_m = as_double(require_value(i, argc, argv), arg) * 1.0e6;
            else throw std::invalid_argument("unknown argument: " + arg);
        }
        if (estimate_only) {
            const auto estimate = inp_path.empty()
                                      ? rocksolver::estimate_resources(config)
                                      : rocksolver::estimate_resources(
                                            rocksolver::import_abaqus_c3d4(inp_path, config));
            std::cout << "{\"version\":\"" << rocksolver::kVersion
                      << "\",\"nodes\":" << estimate.nodes
                      << ",\"zones\":" << estimate.zones
                      << ",\"degreesOfFreedom\":" << estimate.degrees_of_freedom
                      << ",\"memoryMb\":" << estimate.memory_mb
                      << ",\"recommendedBackend\":\"" << estimate.recommended_backend
                      << "\"}\n";
            return 0;
        }
        auto model = inp_path.empty()
                         ? rocksolver::generate_model(config)
                         : rocksolver::import_abaqus_c3d4(inp_path, config);
        const rocksolver::ProgressCallback progress = show_progress
            ? rocksolver::ProgressCallback([](const rocksolver::SolveProgress& value) {
                  std::cerr << std::setprecision(10)
                            << "{\"stage\":\"" << value.stage_name
                            << "\",\"cycle\":" << value.cycle
                            << ",\"maxCycles\":" << value.max_cycles
                            << ",\"forceRatioAverage\":" << value.force_ratio
                            << ",\"forceRatioLocal\":" << value.force_ratio_local
                            << ",\"forceRatioMaximum\":"
                            << value.force_ratio_maximum
                            << ",\"unbalancedForceMaximumN\":"
                            << value.unbalanced_force_maximum_n << "}\n";
              })
            : rocksolver::ProgressCallback{};
        const auto result = rocksolver::solve(model, progress);
        const std::string json = rocksolver::result_to_json(model, result);
        if (output_path.empty()) {
            std::cout << json << '\n';
        } else {
            std::ofstream output(output_path, std::ios::binary | std::ios::trunc);
            if (!output) throw std::runtime_error("cannot open result output file");
            output << json << '\n';
            if (!output) throw std::runtime_error("failed to write result output file");
        }
        return result.converged ? 0 : 3;
    } catch (const std::exception& error) {
        std::cerr << "{\"error\":\"" << error.what() << "\"}\n";
        return 2;
    }
}
