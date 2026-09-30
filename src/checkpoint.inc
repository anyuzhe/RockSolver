// Text format is portable across native/WASM. Explicit fields, full double precision.
namespace {
template<class A> void checkpoint_fields(A& a, Vec3& v) { a(v.x,v.y,v.z); }
template<class A> void checkpoint_fields(A& a, JointSet& v) { a(v.name,v.dip_direction_deg,v.dip_deg,v.normal_density_1_m,v.mean_radius_m,v.cohesion_pa,v.friction_deg,v.tensile_strength_pa,v.normal_stiffness_pa_m,v.shear_stiffness_pa_m); }
template<class A> void checkpoint_fields(A& a, MaterialParameters& v) { a(v.type,v.density_kg_m3,v.elastic_modulus_pa,v.poisson_ratio,v.cohesion_pa,v.friction_deg,v.dilation_deg,v.tensile_strength_pa,v.ucs_pa,v.fracture_toughness_mpa_sqrt_m,v.representative_volume_m3,v.weibull_shape,v.strength_reduction_factor,v.joints); }
template<class A> void checkpoint_fields(A& a, RegionMaterial& v) { a(v.name,v.material); }
template<class A> void checkpoint_fields(A& a, SolverSettings& v) { a(v.max_cycles,v.convergence_tolerance,v.damping,v.time_step_scale,v.large_strain,v.nodal_mixed_discretization); }
template<class A> void checkpoint_fields(A& a, Stage& v) { a(v.name,v.type,v.station_start_m,v.station_end_m,v.support_stiffness_n_m,v.strength_reduction_factor); }
template<class A> void checkpoint_fields(A& a, CableConfig& v) { a(v.id,v.after_initial,v.young,v.area,v.yield_tension,v.yield_compression,v.grout_stiffness,v.grout_cohesion,v.grout_friction,v.perimeter,v.density,v.pretension,v.free_length,v.points); }
template<class A> void checkpoint_fields(A& a, CableNode& v) { a(v.position,v.reference_displacement,v.hosts,v.weights,v.zone,v.host_id,v.bonded_length,v.length,v.axial,v.velocity,v.force,v.reference_force,v.stiffness,v.slip,v.plastic_slip,v.grout_force,v.mass,v.confinement,v.sliding); }
template<class A> void checkpoint_fields(A& a, CableState& v) { a(v.config,v.direction,v.nodes,v.lengths,v.axial_force,v.plastic_extension,v.strain,v.yielded,v.reference_extension,v.tension_control); }
template<class A> void checkpoint_fields(A& a, InterfaceConfig& v) { a(v.name,v.enabled,v.source_group,v.moving_side_hint,v.normal_stiffness_pa_m,v.shear_stiffness_pa_m,v.cohesion_pa,v.friction_deg,v.tensile_strength_pa,v.replace_source_with_host); }
template<class A> void checkpoint_fields(A& a, InterfaceFacet& v) { a(v.master_nodes,v.slave_nodes,v.normal,v.area_m2,v.plastic_shear_slip,v.normal_stress_pa,v.shear_stress_pa,v.state); }
template<class A> void checkpoint_fields(A& a, InterfaceState& v) { a(v.config,v.strength_reduction_factor,v.facets); }
template<class A> void checkpoint_fields(A& a, WaterConfig& v) { a(v.enabled,v.table_z_m,v.density_kg_m3,v.plane_enabled,v.plane_origin_m,v.plane_normal,v.auto_saturated_density,v.saturated_densities_kg_m3,v.external_enabled,v.external_surface_z_m,v.external_density_kg_m3,v.external_outward_direction,v.external_min_alignment); }
template<class A> void checkpoint_fields(A& a, ModelConfig& v) { a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,v.material,v.region_materials,v.strength_reduction_interfaces_only,v.strength_reduction_continuum_only,v.interfaces,v.solver,v.cables,v.water); }
template<class A> void checkpoint_fields(A& a, Node& v) { a(v.id,v.position,v.displacement,v.velocity,v.force,v.dynamic_mass,v.support_reference,v.fixed,v.mass_kg,v.force_reference,v.supported); }
template<class A> void checkpoint_fields(A& a, Zone& v) { a(v.id,v.material_index,v.material_group,v.nodes,v.integration_weight,v.parent_volume_m3,v.active,v.stress_pa,v.strain,v.tangent_pa,v.smrm_history,v.smrm,v.volume_m3,v.reinforce_requirement,v.failure_index,v.plastic_state,v.joint_state); }
template<class A> void checkpoint_fields(A& a, PhysicalCell& v) { a(v.id,v.nodes,v.subzones); }
template<class A> void checkpoint_fields(A& a, RockModel& v) { a(v.version,v.name,v.scenario,v.config,v.nodes,v.zones,v.physical_cells,v.stages,v.cables,v.interfaces,v.initialize_interface_stresses,v.cable_physical_host_interpolation,v.cable_effective_confinement); }
template<class A> void checkpoint_fields(A& a, StageResult& v) { a(v.name,v.cycles,v.converged,v.force_ratio,v.force_ratio_local,v.force_ratio_maximum,v.unbalanced_force_maximum_n,v.max_displacement_m,v.active_zones); }
template<class A> void checkpoint_fields(A& a, RockResult& v) { a(v.solver_version,v.converged,v.total_cycles,v.final_force_ratio,v.final_force_ratio_local,v.final_force_ratio_maximum,v.final_unbalanced_force_maximum_n,v.max_displacement_m,v.max_reinforce_requirement,v.active_zones,v.stages,v.warnings); }
template<class A> void checkpoint_fields(A& a, SMRMHistory& v) { a(v.joint_k,v.joint_h,v.yield_state,v.step); }
template<class A> void checkpoint_fields(A& a, SMRMDiagnostics& v) { a(v.rock_mass_em_gpa,v.smrm_score,v.failure_probability,v.safety_rest,v.reinforce_requirement_mpa,v.stability_coefficient,v.sigma1cmin_mpa,v.energy_density_mj_m3,v.u0v_mj_m3,v.u0d_mj_m3,v.triaxial_compressive_strength_mpa,v.limit_confining_pressure_mpa,v.rockburst_index,v.rock_mass_friction_deg,v.rock_mass_cohesion_mpa,v.rock_mass_poisson_ratio,v.rock_mass_ucs_mpa,v.diag_sigma3_joint_limit_mpa,v.diag_sigma3_controlling_joint,v.joint_k,v.joint_h,v.joint_n1,v.joint_n2,v.joint_n3,v.joint_use_rock_sigma3,v.joint_sigma3_limit_mpa,v.book_quasi_limit_confining_pressure_mpa,v.book_quasi_reinforce_requirement_mpa,v.diag_sigma3_quasi_joint_limit_mpa,v.diag_sigma3_quasi_controlling_joint,v.joint_quasi_sigma3_limit_mpa,v.joint_normal_stress_mpa,v.joint_shear_stress_mpa,v.joint_shear_strength_mpa,v.rock_mass_em_x_gpa,v.rock_mass_em_y_gpa,v.rock_mass_em_z_gpa,v.confining_pressure_margin_mpa,v.book_quasi_confining_pressure_margin_mpa,v.yield_state,v.has_failed); }
struct CheckpointWriter {
    std::string& stream;
    template<class... T> void operator()(T&... v) { (one(v), ...); }
    void one(std::string& v) {
        stream.push_back('"');
        for (const char ch : v) {
            if (ch == '\\' || ch == '"') stream.push_back('\\');
            stream.push_back(ch);
        }
        stream.append("\" ");
    }
    template<class T> void number(T v) {
        char buffer[128];
        std::to_chars_result converted;
        if constexpr(std::is_floating_point_v<T>) {
            converted=std::to_chars(buffer,buffer+sizeof(buffer),v,
                                    std::chars_format::general,
                                    std::numeric_limits<T>::max_digits10);
        } else {
            converted=std::to_chars(buffer,buffer+sizeof(buffer),v);
        }
        if(converted.ec!=std::errc{}) throw std::runtime_error("恢复状态序列化失败");
        stream.append(buffer,converted.ptr);
        stream.push_back(' ');
    }
    template<class T> void one(T& v) {
        if constexpr(std::is_enum_v<T>) number(static_cast<std::underlying_type_t<T>>(v));
        else if constexpr(std::is_same_v<T,bool>) number(v ? 1 : 0);
        else if constexpr(std::is_same_v<T,std::uint8_t>) number(static_cast<unsigned>(v));
        else if constexpr(std::is_arithmetic_v<T>) number(v);
        else checkpoint_fields(*this, v);
    }
    template<class T, std::size_t N> void one(std::array<T,N>& v) { for(auto& x:v) one(x); }
    template<class T> void one(std::vector<T>& v) { number(v.size()); for(auto& x:v) one(x); }
    void one(std::vector<bool>& v) { number(v.size()); for(bool x:v) one(x); }
};
struct CheckpointReader {
    std::istream& stream;
    std::size_t remaining;
    template<class... T> void operator()(T&... v) { (one(v), ...); }
    void one(std::string& v) { stream >> std::quoted(v); check(); }
    void check() { if(!stream) throw std::invalid_argument("恢复文件不完整或字段无效"); }
    template<class T> void one(T& v) {
        if constexpr(std::is_enum_v<T>) { int x; stream >> x; v=static_cast<T>(x); }
        else if constexpr(std::is_same_v<T,std::uint8_t>) { unsigned x; stream >> x; if(x>255) throw std::invalid_argument("状态标记无效"); v=x; }
        else if constexpr(std::is_arithmetic_v<T>) {
            stream >> v;
            if constexpr(std::is_floating_point_v<T>) if(!std::isfinite(v)) throw std::invalid_argument("恢复文件含非有限值");
        } else checkpoint_fields(*this,v);
        check();
    }
    template<class T,std::size_t N> void one(std::array<T,N>& v) { for(auto& x:v) one(x); }
    std::size_t count() { std::size_t n; stream >> n; check(); if(n>remaining || n>5000000) throw std::invalid_argument("恢复文件对象数量无效"); remaining-=n; return n; }
    template<class T> void one(std::vector<T>& v) { v.resize(count()); for(auto& x:v) one(x); }
    void one(std::vector<bool>& v) { v.resize(count()); for(std::size_t i=0;i<v.size();++i){bool x;one(x);v[i]=x;} }
};

// 0.9.4 distributed Cable load through the physical C3D8 cell but did not
// subtract fixed pore pressure from the Cable grout's effective confinement.
void read_checkpoint_094(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario,m.config,m.nodes,m.zones,m.physical_cells,
      m.stages,m.cables,m.interfaces,m.initialize_interface_stresses,
      m.cable_physical_host_interpolation);
    a(r);
}

// 0.9.3 preloaded Interface stress but still distributed Cable load through
// one four-node integration tetrahedron instead of the full C3D8 cell.
void read_checkpoint_093(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario,m.config,m.nodes,m.zones,m.physical_cells,
      m.stages,m.cables,m.interfaces,m.initialize_interface_stresses);
    a(r);
}

// 0.9.2 stored external water pressure but did not preload initialized
// continuum stress onto Interfaces. Keep that old equilibrium state unchanged
// when a legacy checkpoint is resumed; new comparisons need a fresh baseline.
void read_checkpoint_092(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario,m.config,m.nodes,m.zones,m.physical_cells,
      m.stages,m.cables,m.interfaces);
    a(r);
}

// 0.9.1 stored a sloped plane and density switch, but no external traction.
void read_checkpoint_091(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,
      v.material,v.region_materials,v.strength_reduction_interfaces_only,
      v.strength_reduction_continuum_only,v.interfaces,v.solver,v.cables,
      v.water.enabled,v.water.table_z_m,v.water.density_kg_m3,
      v.water.plane_enabled,v.water.plane_origin_m,v.water.plane_normal,
      v.water.auto_saturated_density,v.water.saturated_densities_kg_m3);
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables,m.interfaces);
    a(r);
}

// 0.9.0 stored only a horizontal static water table.
void read_checkpoint_090(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,
      v.material,v.region_materials,v.strength_reduction_interfaces_only,
      v.strength_reduction_continuum_only,v.interfaces,v.solver,v.cables,
      v.water.enabled,v.water.table_z_m,v.water.density_kg_m3);
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables,m.interfaces);
    a(r);
}

// 0.8.2 had the same model layout except for the new hydrostatic-water fields.
void read_checkpoint_082(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,
      v.material,v.region_materials,v.strength_reduction_interfaces_only,
      v.strength_reduction_continuum_only,v.interfaces,v.solver,v.cables);
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables,m.interfaces);
    a(r);
}

// 0.8.0 serialized Interfaces but did not yet have the continuum-only target.
void read_checkpoint_080(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,
      v.material,v.region_materials,v.strength_reduction_interfaces_only,v.interfaces,
      v.solver,v.cables);
    v.strength_reduction_continuum_only=false;
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables,m.interfaces);
    a(r);
}

// 0.7.0 did not serialize optional interfaces. Read that exact layout and
// supply the additive 0.8 defaults so saved pre-Interface studies remain usable.
void read_checkpoint_070(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m,
      v.material,v.region_materials,v.solver,v.cables);
    v.strength_reduction_interfaces_only=false;
    v.strength_reduction_continuum_only=false;
    v.interfaces.clear();
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables);
    m.interfaces.clear();
    a(r);
}

void read_material_061(CheckpointReader& a, MaterialParameters& v) {
    a(v.type,v.density_kg_m3,v.elastic_modulus_pa,v.poisson_ratio,v.cohesion_pa,
      v.friction_deg,v.tensile_strength_pa,v.ucs_pa,v.fracture_toughness_mpa_sqrt_m,
      v.representative_volume_m3,v.weibull_shape,v.strength_reduction_factor,v.joints);
    v.dilation_deg=0.0;
}

void read_checkpoint_061(CheckpointReader& a, RockModel& m, RockResult& r) {
    a(m.version,m.name,m.scenario);
    auto& v=m.config;
    a(v.scenario,v.nx,v.ny,v.nz,v.length_m,v.width_m,v.height_m,v.slope_angle_deg,
      v.tunnel_radius_x_m,v.tunnel_radius_z_m,v.initial_stress_x_pa,v.initial_stress_y_pa,
      v.initial_stress_z_pa,v.gravity_m_s2,v.excavation_steps,v.support_stiffness_n_m);
    read_material_061(a,v.material);
    v.region_materials.resize(a.count());
    for(auto& region:v.region_materials){a(region.name);read_material_061(a,region.material);}
    a(v.solver,v.cables);
    v.strength_reduction_interfaces_only=false;
    v.strength_reduction_continuum_only=false;
    v.interfaces.clear();
    a(m.nodes,m.zones,m.physical_cells,m.stages,m.cables);
    m.interfaces.clear();
    a(r);
}
}
std::string save_checkpoint(const RockModel& model, const RockResult& result) {
    std::size_t cable_nodes=0;
    for(const auto& cable:model.cables) cable_nodes+=cable.nodes.size();
    std::string out;
    out.reserve(1024 + model.nodes.size()*800 + model.zones.size()*1900 +
                model.physical_cells.size()*128 + cable_nodes*512);
    out.append("ROCKSTATE-1 "); out.append(kVersion); out.push_back(' ');
    // The writer never modifies values; this avoids copying the entire mesh.
    CheckpointWriter a{out}; a(const_cast<RockModel&>(model),const_cast<RockResult&>(result));
    // 0.9.8 extension: these opt-in solver/material settings existed in
    // 0.9.7 but were inadvertently omitted from its serialized layout.
    bool static_mass=model.config.solver.flac3d7_static_mass;
    bool host_contact=model.config.material.flac3d7_smrm_contact;
    std::vector<bool> region_contact;
    for(const auto& region:model.config.region_materials)
        region_contact.push_back(region.material.flac3d7_smrm_contact);
    a(static_mass,host_contact,region_contact);
    bool rigid_top_platen=model.config.solver.rigid_top_platen_z;
    a(rigid_top_platen);
    return out;
}
void restore_checkpoint(const std::string& text, RockModel& model, RockResult& result) {
    std::istringstream in(text); std::string format, version; in >> format >> version;
    if(format!="ROCKSTATE-1") throw std::invalid_argument("恢复文件格式无效");
    RockModel m; RockResult r; CheckpointReader a{in,text.size()/2};
    // Hydraulic-only revisions preserve the mechanical checkpoint layout.
    const bool current_layout=version==kVersion||version=="0.13.22-contact-shear-strength"||
        version=="0.13.21-node-boundary-flux"||
        version=="0.13.20-json-roundtrip"||
        version=="0.13.15-mixed-overlay-subzones"||
        version=="0.13.14-flow-free-residual"||
        version=="0.13.13-biot-equilibrium-history"||
        version=="0.13.0-rainfall-unsaturated";
    if(current_layout) {
        a(m,r);
        bool static_mass=false,host_contact=false;
        std::vector<bool> region_contact;
        a(static_mass,host_contact,region_contact);
        if(region_contact.size()!=m.config.region_materials.size())
            throw std::invalid_argument("恢复文件材料兼容标志数量无效");
        m.config.solver.flac3d7_static_mass=static_mass;
        m.config.material.flac3d7_smrm_contact=host_contact;
        for(std::size_t i=0;i<region_contact.size();++i)
            m.config.region_materials[i].material.flac3d7_smrm_contact=region_contact[i];
        in >> std::ws;
        if(!in.eof()) a(m.config.solver.rigid_top_platen_z);
    }
    else if(version=="0.9.7-smrm-flac3d700" ||
            version=="0.9.5-cable-effective-confinement") a(m,r);
    else if(version=="0.9.4-cable-cell-host-weights") read_checkpoint_094(a,m,r);
    else if(version=="0.9.3-interface-initial-stress") read_checkpoint_093(a,m,r);
    else if(version=="0.9.2-external-water-pressure") read_checkpoint_092(a,m,r);
    else if(version=="0.9.1-static-water-plane-density") read_checkpoint_091(a,m,r);
    else if(version=="0.9.0-hydrostatic-water") read_checkpoint_090(a,m,r);
    else if(version=="0.8.2-material-band-target") read_checkpoint_082(a,m,r);
    else if(version=="0.8.0-interface") read_checkpoint_080(a,m,r);
    else if(version=="0.7.0-mohr-coulomb") read_checkpoint_070(a,m,r);
    else if(version=="0.6.1-support-equilibrium") read_checkpoint_061(a,m,r);
    else throw std::invalid_argument("恢复文件版本与求解器不兼容，请使用原版本或重新计算基准");
    if(version!="0.9.3-interface-initial-stress" && version!="0.9.4-cable-cell-host-weights" && version!="0.9.5-cable-effective-confinement" && version!="0.9.7-smrm-flac3d700" && !current_layout)m.initialize_interface_stresses=false;
    if(version!="0.9.4-cable-cell-host-weights" && version!="0.9.5-cable-effective-confinement" && version!="0.9.7-smrm-flac3d700" && !current_layout)m.cable_physical_host_interpolation=false;
    if(version!="0.9.5-cable-effective-confinement" && version!="0.9.7-smrm-flac3d700" && !current_layout)m.cable_effective_confinement=false;
    in >> std::ws; if(!in.eof()) throw std::invalid_argument("恢复文件包含多余数据");
    if(m.nodes.empty() || m.zones.empty()) throw std::invalid_argument("恢复文件缺少网格");
    validate_water(m.config);
    for(const auto& z:m.zones) {
        if(z.nodes.size()!=4 || z.material_index>m.config.region_materials.size() || z.volume_m3<=0) throw std::invalid_argument("恢复单元无效");
        for(auto n:z.nodes) if(n>=m.nodes.size()) throw std::invalid_argument("恢复单元节点引用无效");
    }
    for(const auto& c:m.physical_cells) {
        for(auto n:c.nodes) if(n>=m.nodes.size()) throw std::invalid_argument("物理单元节点引用无效");
        for(auto z:c.subzones) if(z>=m.zones.size()) throw std::invalid_argument("积分子单元引用无效");
    }
    for(const auto& c:m.cables) for(const auto& n:c.nodes) {
        if(n.zone>=m.zones.size()) throw std::invalid_argument("锚杆单元引用无效");
        for(auto h:n.hosts) if(h>=m.nodes.size()) throw std::invalid_argument("锚杆节点引用无效");
    }
    for(auto& c:m.cables)for(auto& n:c.nodes)initialize_cable_physical_hosts(m,n);
    if(m.interfaces.size()!=m.config.interfaces.size()) throw std::invalid_argument("Interface 恢复状态数量无效");
    for(std::size_t si=0;si<m.interfaces.size();++si) {
        const auto& state=m.interfaces[si];
        for(std::size_t fi=0;fi<state.facets.size();++fi) {
            const auto& facet=state.facets[fi];
            const std::size_t node_count=facet.master_nodes.size();
            const bool valid_shape=node_count==1 || node_count==3 || node_count==4;
            if(!valid_shape || node_count!=facet.slave_nodes.size() || !(facet.area_m2>0)) {
                std::ostringstream message;
                message << "Interface 恢复面无效: interface=" << si
                        << " facet=" << fi
                        << " masterNodes=" << node_count
                        << " slaveNodes=" << facet.slave_nodes.size()
                        << " area=" << std::setprecision(17) << facet.area_m2;
                throw std::invalid_argument(message.str());
            }
            for(auto n:facet.master_nodes) if(n>=m.nodes.size()) throw std::invalid_argument("Interface 主面节点引用无效");
            for(auto n:facet.slave_nodes) if(n>=m.nodes.size()) throw std::invalid_argument("Interface 从面节点引用无效");
        }
    }
    model=std::move(m);result=std::move(r);
}
