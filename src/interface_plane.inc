// Optional zero-thickness interface support. Included inside rocksolver.cpp's
// anonymous namespace so it can share Vec3 helpers and the explicit solver state.
namespace {
using CellFaceKey = std::vector<std::size_t>;

const std::vector<std::vector<int>>& interface_cell_faces(std::size_t count) {
    static const std::map<std::size_t,std::vector<std::vector<int>>> shapes{
        {4,{{0,2,1},{0,1,3},{1,2,3},{2,0,3}}},
        {5,{{0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}}},
        {6,{{0,2,1},{3,4,5},{0,1,4,3},{1,2,5,4},{2,0,3,5}}},
        {8,{{0,3,2,1},{4,5,6,7},{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7}}}
    };
    return shapes.at(count);
}

CellFaceKey interface_face_key(const std::vector<std::size_t>& face) {
    auto key=face; std::sort(key.begin(),key.end()); return key;
}

Vec3 interface_cell_center(const RockModel& model,const PhysicalCell& cell) {
    Vec3 c{}; for(auto n:cell.nodes)c+=model.nodes[n].position;
    return c*(1.0/static_cast<double>(cell.nodes.size()));
}
std::pair<Vec3,double> interface_face_geometry(
    const RockModel& model,const std::vector<std::size_t>& face) {
    if(face.size()<3 || face.size()>4)throw std::runtime_error("interface face must have 3 or 4 nodes");
    const Vec3 p0=model.nodes[face[0]].position;
    Vec3 area_vector{};
    for(std::size_t i=1;i+1<face.size();++i) {
        area_vector+=cross(model.nodes[face[i]].position-p0,
                           model.nodes[face[i+1]].position-p0)*0.5;
    }
    const double area=norm(area_vector);
    if(!(area>kTiny))throw std::runtime_error("interface face has zero area");
    return {area_vector*(1.0/area),area};
}

void validate_interface_config(const InterfaceConfig& c) {
    if(!c.enabled)return;
    if(c.source_group.empty())throw std::runtime_error("interface source group is empty");
    if(!(c.normal_stiffness_pa_m>0) || !(c.shear_stiffness_pa_m>0) ||
       c.cohesion_pa<0 || c.tensile_strength_pa<0 ||
       c.friction_deg<0 || c.friction_deg>=90 || norm(c.moving_side_hint)<=kTiny)
        throw std::runtime_error("interface properties are invalid");
}

struct InterfaceFaceOwner {
    std::size_t cell{};
    std::vector<std::size_t> face;
};
std::string interface_upper(std::string value) {
    std::transform(value.begin(),value.end(),value.begin(),
                   [](unsigned char c){return static_cast<char>(std::toupper(c));});
    return value;
}

std::string interface_cell_group(const RockModel& model,const PhysicalCell& cell) {
    if(cell.subzones.empty())return {};
    return model.zones.at(cell.subzones.front()).material_group;
}

void initialize_one_interface(RockModel& model,const InterfaceConfig& config) {
    validate_interface_config(config);
    if(!config.enabled)return;
    if(model.physical_cells.empty())
        throw std::runtime_error("zero-thickness interface currently requires imported mixed C3D5/6/8 cells");
    std::map<CellFaceKey,std::vector<InterfaceFaceOwner>> face_owners;
    for(std::size_t ci=0;ci<model.physical_cells.size();++ci) {
        const auto& cell=model.physical_cells[ci];
        for(const auto& local:interface_cell_faces(cell.nodes.size())) {
            std::vector<std::size_t> face; face.reserve(local.size());
            for(int i:local)face.push_back(cell.nodes[i]);
            face_owners[interface_face_key(face)].push_back({ci,std::move(face)});
        }
    }
    struct SelectedFace { CellFaceKey key; std::size_t source{}, moving{}; std::vector<std::size_t> face; };
    std::vector<SelectedFace> selected;
    const Vec3 hint=config.moving_side_hint*(1.0/norm(config.moving_side_hint));
    const auto source_name=interface_upper(config.source_group);
    for(const auto& [key,owners]:face_owners) {
        if(owners.size()!=2)continue;
        const auto ga=interface_upper(interface_cell_group(model,model.physical_cells[owners[0].cell]));
        const auto gb=interface_upper(interface_cell_group(model,model.physical_cells[owners[1].cell]));
        const bool a_source=ga==source_name,b_source=gb==source_name;
        if(a_source==b_source)continue;
        const auto& source=a_source?owners[0]:owners[1];
        const auto& other=a_source?owners[1]:owners[0];
        const Vec3 delta=interface_cell_center(model,model.physical_cells[other.cell])-
                         interface_cell_center(model,model.physical_cells[source.cell]);
        auto [selection_normal,selection_area]=interface_face_geometry(model,source.face);
        (void)selection_area;
        if(dot(selection_normal,delta)<0)selection_normal=selection_normal*(-1.0);
        // moving_side_hint selects the actual source-group face normal, not only
        // the neighbour centroid direction. This rejects strip end-caps (toe/base)
        // while retaining curved/stair-stepped faces that genuinely point to the
        // requested moving side.
        if(dot(selection_normal,hint)>1.0e-8)
            selected.push_back({key,source.cell,other.cell,source.face});
    }
    if(selected.empty())throw std::runtime_error("interface source group has no boundary on the requested moving side");
    std::set<CellFaceKey> cuts; for(const auto& item:selected)cuts.insert(item.key);

    // Split each interface node locally. This permits an interface to terminate
    // inside the mesh: crack-front nodes stay shared instead of requiring the
    // discontinuity to cut the entire domain into two disconnected components.
    std::set<std::size_t> interface_nodes;
    for(const auto& item:selected)interface_nodes.insert(item.face.begin(),item.face.end());
    std::map<std::size_t,std::set<std::size_t>> moving_cells_by_node;
    for(auto node:interface_nodes) {
        std::map<std::size_t,std::vector<std::size_t>> adjacency;
        for(const auto& [key,owners]:face_owners) {
            if(owners.size()!=2 || cuts.count(key) ||
               std::find(key.begin(),key.end(),node)==key.end())continue;
            adjacency[owners[0].cell].push_back(owners[1].cell);
            adjacency[owners[1].cell].push_back(owners[0].cell);
        }
        std::set<std::size_t> visited,forbidden;
        std::vector<std::size_t> queue;
        for(const auto& item:selected)if(std::find(item.face.begin(),item.face.end(),node)!=item.face.end()) {
            if(visited.insert(item.moving).second)queue.push_back(item.moving);
            forbidden.insert(item.source);
        }
        for(std::size_t q=0;q<queue.size();++q)for(auto next:adjacency[queue[q]])
            if(visited.insert(next).second)queue.push_back(next);
        bool crack_front=false;for(auto source:forbidden)if(visited.count(source)){crack_front=true;break;}
        if(!crack_front)moving_cells_by_node[node]=std::move(visited);
    }

    int next_id=0; for(const auto& node:model.nodes)next_id=std::max(next_id,node.id);
    std::map<std::size_t,std::size_t> duplicate;
    for(const auto& [original,cells]:moving_cells_by_node) {
        if(cells.empty())continue;
        Node node=model.nodes[original]; node.id=++next_id;
        duplicate[original]=model.nodes.size(); model.nodes.push_back(node);
    }

    for(std::size_t ci=0;ci<model.physical_cells.size();++ci) {
        auto& cell=model.physical_cells[ci];
        for(auto& n:cell.nodes)if(auto it=duplicate.find(n);it!=duplicate.end() && moving_cells_by_node[n].count(ci))n=it->second;
        for(auto zi:cell.subzones)for(auto& n:model.zones[zi].nodes)
            if(auto it=duplicate.find(n);it!=duplicate.end() && moving_cells_by_node[n].count(ci))n=it->second;
    }
    if(config.replace_source_with_host) {
        for(auto& zone:model.zones)
            if(interface_upper(zone.material_group)==source_name)zone.material_index=0;
    }

    // FLAC3D interfaces are node based: each interface node owns its own
    // normal/shear spring and plastic slip history.  Accumulate the tributary
    // area and an area-weighted normal from all selected faces that meet at the
    // same split node, rather than averaging the jump over an entire face.
    struct NodalInterfaceContribution {
        std::size_t master{};
        std::size_t slave{};
        Vec3 weighted_normal{};
        double area_m2{0.0};
    };
    std::map<std::pair<std::size_t,std::size_t>,NodalInterfaceContribution> nodal;
    for(const auto& item:selected) {
        auto [normal,area]=interface_face_geometry(model,item.face);
        const Vec3 direction=interface_cell_center(model,model.physical_cells[item.moving])-
                             interface_cell_center(model,model.physical_cells[item.source]);
        if(dot(normal,direction)<0)normal=normal*(-1.0);
        // FLAC3D defines the interface-node normal as the area-weighted average
        // of the triangular interface elements connected to that node.  A
        // quadrilateral zone face is represented by two interface triangles;
        // using one face normal with equal node weights is only equivalent for
        // a single planar face, not at a kink shared by adjacent faces.
        std::map<std::size_t,Vec3> node_normal_weight;
        const Vec3 p0=model.nodes[item.face[0]].position;
        for(std::size_t ti=1;ti+1<item.face.size();++ti) {
            const auto a=item.face[0],b=item.face[ti],c=item.face[ti+1];
            Vec3 area_vector=cross(model.nodes[b].position-p0,
                                   model.nodes[c].position-p0)*0.5;
            if(dot(area_vector,direction)<0)area_vector=area_vector*(-1.0);
            node_normal_weight[a]+=area_vector;
            node_normal_weight[b]+=area_vector;
            node_normal_weight[c]+=area_vector;
        }
        const double share=area/static_cast<double>(item.face.size());
        for(const auto master:item.face) {
            const auto found=duplicate.find(master);
            const std::size_t slave=found==duplicate.end()?master:found->second;
            const auto key=std::pair{master,slave};
            auto& contribution=nodal[key];
            contribution.master=master;
            contribution.slave=slave;
            contribution.weighted_normal+=node_normal_weight.at(master);
            contribution.area_m2+=share;
        }
    }

    InterfaceState state; state.config=config;
    state.facets.reserve(nodal.size());
    for(const auto& [key,contribution]:nodal) {
        (void)key;
        if(!(contribution.area_m2>kTiny))continue;
        const double normal_length=norm(contribution.weighted_normal);
        if(!(normal_length>kTiny))throw std::runtime_error("interface node has indeterminate normal");
        InterfaceFacet facet;
        facet.master_nodes={contribution.master};
        facet.slave_nodes={contribution.slave};
        facet.normal=contribution.weighted_normal*(1.0/normal_length);
        facet.area_m2=contribution.area_m2;
        state.facets.push_back(std::move(facet));
    }
    model.interfaces.push_back(std::move(state));
}

void initialize_interfaces(RockModel& model) {
    model.interfaces.clear();
    for(const auto& config:model.config.interfaces)initialize_one_interface(model,config);
}
void interface_dynamic_mass(RockModel& model) {
    for(const auto& state:model.interfaces)for(const auto& facet:state.facets) {
        const double share=facet.area_m2/std::max<std::size_t>(facet.master_nodes.size(),1);
        const Vec3 n=facet.normal;
        const double kn=state.config.normal_stiffness_pa_m*share;
        const double ks=state.config.shear_stiffness_pa_m*share;
        const Vec3 diagonal{ks+(kn-ks)*n.x*n.x,
                            ks+(kn-ks)*n.y*n.y,
                            ks+(kn-ks)*n.z*n.z};
        for(std::size_t i=0;i<facet.master_nodes.size();++i) {
            if(facet.master_nodes[i]==facet.slave_nodes[i])continue;
            model.nodes[facet.master_nodes[i]].dynamic_mass+=diagonal;
            model.nodes[facet.slave_nodes[i]].dynamic_mass+=diagonal;
        }
    }
}

Vec3 interface_average_jump(const RockModel& model,const InterfaceFacet& facet) {
    Vec3 jump{};
    for(std::size_t i=0;i<facet.master_nodes.size();++i)
        jump+=model.nodes[facet.slave_nodes[i]].displacement-
              model.nodes[facet.master_nodes[i]].displacement;
    return jump*(1.0/static_cast<double>(facet.master_nodes.size()));
}

void add_interface_forces(RockModel& model) {
    for(auto& state:model.interfaces)for(auto& facet:state.facets) {
        const double reduction=std::max(state.strength_reduction_factor,1.0);
        const double phi=std::atan(std::tan(state.config.friction_deg*std::numbers::pi/180.0)/reduction);
        const Vec3 jump=interface_average_jump(model,facet);
        const double un=dot(jump,facet.normal);
        // FLAC's `zone interface node initialize-stresses` transfers the
        // initialized host-zone traction to the interface.  The Web exporter
        // issues that command after `zone initialize stress`; without this
        // offset an initially compressed split face must penetrate by
        // |sigma_n|/kn before carrying its prescribed in-situ stress.
        const Vec3 initial_traction=model.initialize_interface_stresses?Vec3{
            model.config.initial_stress_x_pa*facet.normal.x,
            model.config.initial_stress_y_pa*facet.normal.y,
            model.config.initial_stress_z_pa*facet.normal.z}:Vec3{};
        const double initial_normal=dot(initial_traction,facet.normal);
        const Vec3 initial_shear=initial_traction-facet.normal*initial_normal;
        double normal_stress=initial_normal+state.config.normal_stiffness_pa_m*un;
        const double tension=state.config.tensile_strength_pa/reduction;
        Vec3 shear_stress{};
        if(normal_stress>tension) {
            normal_stress=tension;
            facet.state=JointState::TensionFailure;
        } else {
            const Vec3 shear_displacement=jump-facet.normal*un;
            const Vec3 elastic_shear=shear_displacement-facet.plastic_shear_slip;
            const Vec3 trial=initial_shear+elastic_shear*state.config.shear_stiffness_pa_m;
            const double trial_magnitude=norm(trial);
            // FLAC3D's interface slip condition uses the target-face pore
            // pressure to reduce compressive normal force, not the contact
            // force itself. The source-side master node is the target face;
            // split master/slave nodes coincide geometrically for the current
            // single static water field, so their present pressures agree.
            const double pore_pressure=node_pore_pressure(model,facet.master_nodes.front());
            const double strength=state.config.cohesion_pa/reduction+
                std::max(0.0,-normal_stress-pore_pressure)*std::tan(phi);
            if(trial_magnitude>strength+1.0e-12) {
                shear_stress=trial_magnitude>0?trial*(strength/trial_magnitude):Vec3{};
                facet.plastic_shear_slip=shear_displacement+
                    (initial_shear-shear_stress)*(1.0/state.config.shear_stiffness_pa_m);
                facet.state=JointState::ShearFailure;
            } else {
                shear_stress=trial;
                facet.state=JointState::Closed;
            }
        }
        if(facet.state==JointState::TensionFailure && tension<=0.0)shear_stress={};
        facet.normal_stress_pa=normal_stress;
        facet.shear_stress_pa=norm(shear_stress);
        const Vec3 traction=facet.normal*normal_stress+shear_stress;
        const double share=facet.area_m2/static_cast<double>(facet.master_nodes.size());
        const Vec3 slave_force=traction*(-share);
        const Vec3 master_force=slave_force*(-1.0);
        for(std::size_t i=0;i<facet.master_nodes.size();++i) {
            if(facet.master_nodes[i]==facet.slave_nodes[i])continue;
            model.nodes[facet.slave_nodes[i]].force+=slave_force;
            model.nodes[facet.slave_nodes[i]].force_reference+=norm(slave_force);
            model.nodes[facet.master_nodes[i]].force+=master_force;
            model.nodes[facet.master_nodes[i]].force_reference+=norm(master_force);
        }
    }
}

void set_interface_strength_reduction(RockModel& model,double factor) {
    for(auto& state:model.interfaces)state.strength_reduction_factor=std::max(factor,1.0);
}

void reset_interface_strength_reduction(RockModel& model) {
    set_interface_strength_reduction(model,1.0);
}
} // namespace
