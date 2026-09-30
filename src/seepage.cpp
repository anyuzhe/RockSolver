#include "rocksolver/rocksolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace rocksolver {
namespace {

constexpr double kSmall = 1.0e-30;

template <typename Callback>
void for_each_coupled_cell(const RockModel& model,Callback&& callback) {
    if(!model.physical_cells.empty()) {
        for(const auto& cell:model.physical_cells)
            callback(cell.nodes,cell.subzones);
    } else {
        // Pure C3D4 imports store each physical element directly as one Zone.
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];
            if(!zone.active)continue;
            callback(zone.nodes,std::array<std::size_t,1>{zi});
        }
    }
}

Vec3 add(const Vec3& a, const Vec3& b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
Vec3 sub(const Vec3& a, const Vec3& b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
Vec3 scale(const Vec3& a, double s) { return {a.x*s,a.y*s,a.z*s}; }
double dot3(const Vec3& a, const Vec3& b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
Vec3 cross3(const Vec3& a, const Vec3& b) {
    return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};
}
double length(const Vec3& value) { return std::sqrt(dot3(value,value)); }

struct TetData {
    double volume{};
    std::array<Vec3,4> gradient{};
};

TetData tet_data(const RockModel& model, const Zone& zone) {
    if(zone.nodes.size()!=4)throw std::invalid_argument("steady seepage requires Tet4 integration zones");
    const Vec3 p0=model.nodes.at(zone.nodes[0]).position;
    const Vec3 a=sub(model.nodes.at(zone.nodes[1]).position,p0);
    const Vec3 b=sub(model.nodes.at(zone.nodes[2]).position,p0);
    const Vec3 c=sub(model.nodes.at(zone.nodes[3]).position,p0);
    const double det=dot3(a,cross3(b,c));
    if(!std::isfinite(det)||std::abs(det)<=kSmall)throw std::invalid_argument("steady seepage found a degenerate Tet4");
    // Rows of J^-1 are the physical gradients of xi, eta and zeta.
    const Vec3 g1=scale(cross3(b,c),1.0/det);
    const Vec3 g2=scale(cross3(c,a),1.0/det);
    const Vec3 g3=scale(cross3(a,b),1.0/det);
    return {std::abs(det)/6.0,{scale(add(add(g1,g2),g3),-1.0),g1,g2,g3}};
}

double deformed_tet_signed_volume(const RockModel& model,const Zone& zone) {
    const auto displaced=[&](std::size_t index){
        const auto& node=model.nodes.at(zone.nodes.at(index));
        return add(node.position,node.displacement);
    };
    const Vec3 p0=displaced(0);
    return dot3(sub(displaced(1),p0),cross3(sub(displaced(2),p0),
        sub(displaced(3),p0)))/6.0;
}

using Tensor=std::array<double,6>;

Vec3 multiply(const Tensor& k,const Vec3& v) {
    return {k[0]*v.x+k[3]*v.y+k[4]*v.z,
            k[3]*v.x+k[1]*v.y+k[5]*v.z,
            k[4]*v.x+k[5]*v.y+k[2]*v.z};
}

double determinant(const Tensor& k) {
    return k[0]*(k[1]*k[2]-k[5]*k[5])
          -k[3]*(k[3]*k[2]-k[4]*k[5])
          +k[4]*(k[3]*k[5]-k[4]*k[1]);
}

void validate_tensor(const Tensor& k) {
    if(std::any_of(k.begin(),k.end(),[](double v){return !std::isfinite(v);}) ||
       k[0]<=0.0 || k[0]*k[1]-k[3]*k[3]<=0.0 || determinant(k)<=0.0)
        throw std::invalid_argument("hydraulic conductivity tensor must be finite, symmetric positive definite, and use m/s");
}

const Tensor& zone_tensor(const std::vector<Tensor>& tensors,std::size_t zone_index) {
    return tensors.at(zone_index);
}

double fixed_water_surface_z(const WaterConfig& water,const Vec3& position) {
    if(!water.plane_enabled)return water.table_z_m;
    return water.plane_origin_m.z-
        (water.plane_normal.x*(position.x-water.plane_origin_m.x)+
         water.plane_normal.y*(position.y-water.plane_origin_m.y))/water.plane_normal.z;
}

double existing_zone_pore_pressure(const RockModel& model,const Zone& zone) {
    double value=0.0;
    if(model.hydraulic.valid && model.hydraulic.nodal_pore_pressure_pa.size()==model.nodes.size()){
        for(const auto node:zone.nodes)value+=model.hydraulic.nodal_pore_pressure_pa[node];
        return value/static_cast<double>(zone.nodes.size());
    }
    if(!model.config.water.enabled)return 0.0;
    for(const auto node:zone.nodes){
        const auto& position=model.nodes[node].position;
        value+=std::max(0.0,model.config.water.density_kg_m3*model.config.gravity_m_s2*
            (fixed_water_surface_z(model.config.water,position)-position.z));
    }
    return value/static_cast<double>(zone.nodes.size());
}

struct DerivedHydraulics {
    std::vector<Tensor> conductivity;
    std::vector<std::array<double,3>> aperture;
    std::vector<std::array<double,3>> normal_stress_mpa;
};

DerivedHydraulics derive_zone_hydraulics(const RockModel& model) {
    const auto& config=model.config.seepage;
    DerivedHydraulics result;
    result.conductivity.resize(model.zones.size());
    result.aperture.resize(model.zones.size());
    result.normal_stress_mpa.resize(model.zones.size());
    if(config.smrm_stress_permeability){
        if(model.config.material.type!=MaterialType::SMRM)
            throw std::invalid_argument("SMRM stress-dependent permeability requires an SMRM host material");
        if(config.joint_hydraulics.size()!=model.config.material.joints.size() ||
           config.joint_hydraulics.size()>3)
            throw std::invalid_argument("SMRM hydraulic parameters are required for every host joint set");
        if(!std::isfinite(config.fluid_kinematic_viscosity_m2_s)||config.fluid_kinematic_viscosity_m2_s<=0.0)
            throw std::invalid_argument("fluid kinematic viscosity must be positive");
        for(const auto& hydraulic:config.joint_hydraulics)if(hydraulic.enabled){
            if(!std::isfinite(hydraulic.zero_stress_aperture_m)||hydraulic.zero_stress_aperture_m<=0.0 ||
               !std::isfinite(hydraulic.closure_stress_pa)||hydraulic.closure_stress_pa<=0.0 ||
               !std::isfinite(hydraulic.connectivity)||hydraulic.connectivity<0.0||hydraulic.connectivity>1.0 ||
               !std::isfinite(hydraulic.minimum_aperture_m)||hydraulic.minimum_aperture_m<=0.0 ||
               !std::isfinite(hydraulic.maximum_aperture_m)||hydraulic.maximum_aperture_m<hydraulic.minimum_aperture_m)
                throw std::invalid_argument("invalid SMRM joint hydraulic parameters");
        }
    }
    for(std::size_t zi=0;zi<model.zones.size();++zi){
        const auto& zone=model.zones[zi];
        Tensor tensor=config.conductivity_m_s.at(zone.material_index);
        if(config.smrm_stress_permeability && zone.material_index==0){
            std::array<double,6> stress=zone.stress_pa;
            if(zone.smrm_history.step==0)
                stress={model.config.initial_stress_x_pa,model.config.initial_stress_y_pa,
                        model.config.initial_stress_z_pa,0.0,0.0,0.0};
            const double pore_pressure=existing_zone_pore_pressure(model,zone);
            for(int i=0;i<3;++i)stress[i]+=pore_pressure;
            for(std::size_t ji=0;ji<config.joint_hydraulics.size();++ji){
                const auto& hydraulic=config.joint_hydraulics[ji];
                if(!hydraulic.enabled)continue;
                const auto& joint=model.config.material.joints[ji];
                const auto n=joint_direction_cosines(joint.dip_direction_deg,joint.dip_deg);
                const double normal_tension_positive=
                    n[0]*(stress[0]*n[0]+stress[3]*n[1]+stress[4]*n[2])+
                    n[1]*(stress[3]*n[0]+stress[1]*n[1]+stress[5]*n[2])+
                    n[2]*(stress[4]*n[0]+stress[5]*n[1]+stress[2]*n[2]);
                const double normal_compression=-normal_tension_positive;
                const double exponent=std::clamp(-normal_compression/hydraulic.closure_stress_pa,-700.0,700.0);
                const double aperture=std::clamp(hydraulic.zero_stress_aperture_m*std::exp(exponent),
                                                 hydraulic.minimum_aperture_m,
                                                 hydraulic.maximum_aperture_m);
                result.aperture[zi][ji]=aperture;
                result.normal_stress_mpa[zi][ji]=normal_compression*1.0e-6;
                const double contribution=model.config.gravity_m_s2/(24.0*config.fluid_kinematic_viscosity_m2_s)*
                    hydraulic.connectivity*std::max(0.0,joint.normal_density_1_m)*aperture*aperture*aperture;
                tensor[0]+=contribution*(1.0-n[0]*n[0]);
                tensor[1]+=contribution*(1.0-n[1]*n[1]);
                tensor[2]+=contribution*(1.0-n[2]*n[2]);
                tensor[3]-=contribution*n[0]*n[1];
                tensor[4]-=contribution*n[0]*n[2];
                tensor[5]-=contribution*n[1]*n[2];
            }
        }
        validate_tensor(tensor);
        result.conductivity[zi]=tensor;
    }
    return result;
}

struct Bounds {Vec3 low,high;double tolerance;};
Bounds bounds(const RockModel& model) {
    if(model.nodes.empty())throw std::invalid_argument("steady seepage model has no nodes");
    Vec3 low=model.nodes.front().position,high=low;
    for(const auto& node:model.nodes){
        low.x=std::min(low.x,node.position.x);low.y=std::min(low.y,node.position.y);low.z=std::min(low.z,node.position.z);
        high.x=std::max(high.x,node.position.x);high.y=std::max(high.y,node.position.y);high.z=std::max(high.z,node.position.z);
    }
    const double span=std::max({high.x-low.x,high.y-low.y,high.z-low.z,1.0});
    return {low,high,span*1.0e-8};
}

bool selected(const Vec3& point,const HydraulicBoundary& boundary,const Bounds& box) {
    const double tolerance=boundary.tolerance_m>0.0?boundary.tolerance_m:box.tolerance;
    switch(boundary.selector){
        case HydraulicBoundarySelector::XMin:return std::abs(point.x-box.low.x)<=tolerance;
        case HydraulicBoundarySelector::XMax:return std::abs(point.x-box.high.x)<=tolerance;
        case HydraulicBoundarySelector::YMin:return std::abs(point.y-box.low.y)<=tolerance;
        case HydraulicBoundarySelector::YMax:return std::abs(point.y-box.high.y)<=tolerance;
        case HydraulicBoundarySelector::ZMin:return std::abs(point.z-box.low.z)<=tolerance;
        case HydraulicBoundarySelector::ZMax:return std::abs(point.z-box.high.z)<=tolerance;
        case HydraulicBoundarySelector::Plane:{
            const double normal_length=length(boundary.plane_normal);
            if(!std::isfinite(normal_length)||normal_length<=kSmall)
                throw std::invalid_argument("hydraulic plane boundary requires a nonzero normal");
            return std::abs(dot3(sub(point,boundary.plane_origin_m),boundary.plane_normal)/normal_length)<=tolerance;
        }
    }
    return false;
}

struct BoundaryFace {std::vector<std::size_t> nodes;std::size_t owner_zone{};Vec3 normal;double area{};};

std::vector<BoundaryFace> boundary_faces(const RockModel& model) {
    static const std::array<std::array<int,3>,4> tet_faces{{{0,2,1},{0,1,3},{1,2,3},{2,0,3}}};
    struct Candidate {std::array<std::size_t,3> nodes;std::size_t owner{};int count{1};};
    std::map<std::array<std::size_t,3>,Candidate> faces;
    for(std::size_t zi=0;zi<model.zones.size();++zi){
        const auto& zone=model.zones[zi];if(!zone.active)continue;
        for(const auto& local:tet_faces){
            std::array<std::size_t,3> oriented{zone.nodes[local[0]],zone.nodes[local[1]],zone.nodes[local[2]]};
            auto key=oriented;std::sort(key.begin(),key.end());
            auto [it,inserted]=faces.emplace(key,Candidate{oriented,zi,1});
            if(!inserted)++it->second.count;
        }
    }
    std::vector<BoundaryFace> result;
    for(const auto& [key,candidate]:faces)if(candidate.count==1){
        const Vec3 a=model.nodes[candidate.nodes[0]].position;
        const Vec3 b=model.nodes[candidate.nodes[1]].position;
        const Vec3 c=model.nodes[candidate.nodes[2]].position;
        Vec3 area_vector=cross3(sub(b,a),sub(c,a));
        const double twice_area=length(area_vector);if(twice_area<=kSmall)continue;
        Vec3 normal=scale(area_vector,1.0/twice_area);
        Vec3 center{};for(const auto ni:model.zones[candidate.owner].nodes)center=add(center,model.nodes[ni].position);
        center=scale(center,0.25);
        const Vec3 face_center=scale(add(add(a,b),c),1.0/3.0);
        if(dot3(normal,sub(center,face_center))>0.0)normal=scale(normal,-1.0);
        // Mixed physical cells (C3D6/C3D8) use two overlapping tetrahedral
        // overlays. Their volume terms are averaged with integration_weight;
        // exterior face fluxes must use the same weight or a prescribed
        // discharge is applied once per overlay (twice for C3D8).
        const double area=0.5*twice_area*model.zones[candidate.owner].integration_weight;
        result.push_back({{candidate.nodes[0],candidate.nodes[1],candidate.nodes[2]},candidate.owner,normal,area});
    }
    return result;
}

using SparseRows=std::vector<std::unordered_map<std::size_t,double>>;

std::vector<double> multiply(const SparseRows& rows,const std::vector<double>& x) {
    std::vector<double> value(rows.size(),0.0);
    for(std::size_t i=0;i<rows.size();++i)
        for(const auto& [j,a]:rows[i])value[i]+=a*x[j];
    return value;
}

double dot_vector(const std::vector<double>& a,const std::vector<double>& b) {
    double value=0.0;for(std::size_t i=0;i<a.size();++i)value+=a[i]*b[i];return value;
}

struct CgResult {bool converged{};int iterations{};double relative_residual{};};
CgResult conjugate_gradient(const SparseRows& rows,const std::vector<double>& rhs,
                            std::vector<double>& x,int max_iterations,double tolerance,
                            const std::map<std::size_t,double>& prescribed) {
    auto ax=multiply(rows,x);
    std::vector<double> r(rhs.size()),z(rhs.size()),p(rhs.size()),diagonal(rhs.size());
    for(std::size_t i=0;i<rhs.size();++i){
        r[i]=rhs[i]-ax[i];
        const auto found=rows[i].find(i);
        diagonal[i]=found==rows[i].end()?0.0:found->second;
        if(!(diagonal[i]>0.0)||!std::isfinite(diagonal[i]))throw std::runtime_error("steady seepage matrix has a nonpositive diagonal");
        z[i]=r[i]/diagonal[i];p[i]=z[i];
    }
    // Dirichlet identity rows contain heads (m), whereas free rows contain
    // water rates (m3/s). Including prescribed heads in the residual scale can
    // declare an unsolved low-conductivity field converged at iteration zero.
    double free_rhs_squared=0.0;
    for(std::size_t i=0;i<rhs.size();++i)
        if(prescribed.find(i)==prescribed.end())free_rhs_squared+=rhs[i]*rhs[i];
    const double rhs_norm=std::sqrt(std::max(free_rhs_squared,kSmall));
    double rz=dot_vector(r,z);
    double residual=std::sqrt(std::max(dot_vector(r,r),0.0))/rhs_norm;
    if(residual<=tolerance)return {true,0,residual};
    for(int iteration=1;iteration<=max_iterations;++iteration){
        const auto ap=multiply(rows,p);
        const double denominator=dot_vector(p,ap);
        if(!(denominator>0.0)||!std::isfinite(denominator))throw std::runtime_error("steady seepage matrix is singular or not positive definite");
        const double alpha=rz/denominator;
        for(std::size_t i=0;i<x.size();++i){x[i]+=alpha*p[i];r[i]-=alpha*ap[i];}
        residual=std::sqrt(std::max(dot_vector(r,r),0.0))/rhs_norm;
        if(residual<=tolerance)return {true,iteration,residual};
        for(std::size_t i=0;i<z.size();++i)z[i]=r[i]/diagonal[i];
        const double next_rz=dot_vector(r,z);
        const double beta=next_rz/rz;
        for(std::size_t i=0;i<p.size();++i)p[i]=z[i]+beta*p[i];
        rz=next_rz;
    }
    return {false,max_iterations,residual};
}

double boundary_value_at(const HydraulicBoundary& boundary,double time_s) {
    if(boundary.history.empty())return boundary.value;
    for(std::size_t i=0;i<boundary.history.size();++i){
        const auto& point=boundary.history[i];
        if(!std::isfinite(point.time_s)||point.time_s<0.0||!std::isfinite(point.value)||
           (i>0&&point.time_s<=boundary.history[i-1].time_s))
            throw std::invalid_argument("hydraulic boundary history requires finite, strictly increasing nonnegative times");
    }
    if(time_s<=boundary.history.front().time_s)return boundary.history.front().value;
    if(time_s>=boundary.history.back().time_s)return boundary.history.back().value;
    const auto upper=std::upper_bound(boundary.history.begin(),boundary.history.end(),time_s,
        [](double time,const HydraulicHistoryPoint& point){return time<point.time_s;});
    const auto& right=*upper;
    const auto& left=*(upper-1);
    const double fraction=(time_s-left.time_s)/(right.time_s-left.time_s);
    return left.value+(right.value-left.value)*fraction;
}

struct UnsaturatedPoint {
    double saturation{};
    double water_content{};
    double relative_conductivity{};
    double capacity_1_m{};
};

UnsaturatedPoint unsaturated_point(const UnsaturatedHydraulicParameters& material,
                                   double specific_storage_1_m,
                                   double pressure_head_m) {
    if(!(material.porosity>0.0&&material.porosity<=1.0) ||
       !(material.residual_saturation>=0.0&&material.residual_saturation<1.0) ||
       !(material.van_genuchten_alpha_1_m>0.0) || !(material.van_genuchten_n>1.0) ||
       !std::isfinite(material.mualem_l))
        throw std::invalid_argument("invalid van Genuchten-Mualem unsaturated parameters");
    if(pressure_head_m>=0.0)
        return {1.0,material.porosity+specific_storage_1_m*pressure_head_m,1.0,specific_storage_1_m};
    const double m=1.0-1.0/material.van_genuchten_n;
    const double a=material.van_genuchten_alpha_1_m;
    const double x=std::max(a*(-pressure_head_m),1.0e-30);
    const double xn=std::pow(x,material.van_genuchten_n);
    const double effective=std::pow(1.0+xn,-m);
    const double saturation=material.residual_saturation+
        (1.0-material.residual_saturation)*effective;
    const double derivative=m*material.van_genuchten_n*a*
        std::pow(x,material.van_genuchten_n-1.0)*std::pow(1.0+xn,-m-1.0);
    const double capacity=material.porosity*(1.0-material.residual_saturation)*derivative;
    const double inner=1.0-std::pow(1.0-std::pow(effective,1.0/m),m);
    const double relative=std::clamp(std::pow(effective,material.mualem_l)*inner*inner,1.0e-12,1.0);
    return {saturation,material.porosity*saturation,relative,std::max(capacity,1.0e-16)};
}

Tensor scaled_tensor(const Tensor& tensor,double scale_value) {
    Tensor result=tensor;
    for(auto& value:result)value*=scale_value;
    return result;
}

} // namespace

HydraulicState solve_steady_seepage(RockModel& model) {
    const auto& config=model.config.seepage;
    if(!config.enabled)throw std::invalid_argument("steady seepage is not enabled");
    for(const auto& boundary:config.boundaries)
        if(boundary.kind==HydraulicBoundaryKind::Rainfall||
           boundary.kind==HydraulicBoundaryKind::SeepageFace)
            throw std::invalid_argument("rainfall and seepage-face boundaries require transient unsaturated seepage");
    if(!std::isfinite(config.fluid_density_kg_m3)||config.fluid_density_kg_m3<=0.0 ||
       !std::isfinite(model.config.gravity_m_s2)||model.config.gravity_m_s2<=0.0)
        throw std::invalid_argument("steady seepage requires positive fluid density and gravity");
    if(config.conductivity_m_s.size()!=model.config.region_materials.size()+1)
        throw std::invalid_argument("steady seepage requires one conductivity tensor for every material");
    for(const auto& tensor:config.conductivity_m_s)validate_tensor(tensor);
    if(config.max_iterations<=0||!std::isfinite(config.relative_tolerance)||config.relative_tolerance<=0.0)
        throw std::invalid_argument("steady seepage iteration controls are invalid");
    if(config.solve_mechanics && model.config.water.auto_saturated_density &&
       model.config.water.saturated_densities_kg_m3.size()!=model.config.region_materials.size()+1)
        throw std::invalid_argument("steady seepage mechanics requires saturated density for every material");

    const auto derived=derive_zone_hydraulics(model);
    SparseRows matrix(model.nodes.size());
    std::vector<double> rhs(model.nodes.size(),0.0);
    for(std::size_t zi=0;zi<model.zones.size();++zi){
        const auto& zone=model.zones[zi];if(!zone.active)continue;
        const auto geometry=tet_data(model,zone);
        const auto& tensor=zone_tensor(derived.conductivity,zi);
        const double weight=geometry.volume*zone.integration_weight;
        for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j)
            matrix[zone.nodes[i]][zone.nodes[j]]+=weight*dot3(geometry.gradient[i],multiply(tensor,geometry.gradient[j]));
    }

    const auto box=bounds(model);
    const auto faces=boundary_faces(model);
    std::map<std::size_t,double> prescribed;
    const double rho_g=config.fluid_density_kg_m3*model.config.gravity_m_s2;
    for(const auto& boundary:config.boundaries){
        if(!std::isfinite(boundary.value))throw std::invalid_argument("hydraulic boundary value must be finite");
        std::size_t matches=0;
        if(boundary.kind==HydraulicBoundaryKind::Head ||
           boundary.kind==HydraulicBoundaryKind::PorePressure ||
           boundary.kind==HydraulicBoundaryKind::Drain){
            for(std::size_t i=0;i<model.nodes.size();++i)if(selected(model.nodes[i].position,boundary,box)){
                ++matches;
                const double head=boundary.kind==HydraulicBoundaryKind::Head?boundary.value:
                    model.nodes[i].position.z+(boundary.kind==HydraulicBoundaryKind::Drain?0.0:boundary.value)/rho_g;
                auto [entry,inserted]=prescribed.emplace(i,head);
                if(!inserted&&std::abs(entry->second-head)>std::max(1.0,std::abs(head))*1.0e-10)
                    throw std::invalid_argument("conflicting prescribed hydraulic heads on one node");
            }
        }else{
            for(const auto& face:faces){
                bool face_selected=true;for(const auto node:face.nodes)face_selected&=selected(model.nodes[node].position,boundary,box);
                if(!face_selected)continue;++matches;
                if(boundary.kind==HydraulicBoundaryKind::Flux){
                    // Positive value is outward flow: A h = - integral(N q_n).
                    const double nodal=-boundary.value*face.area/3.0;
                    for(const auto node:face.nodes)rhs[node]+=nodal;
                }
            }
        }
        if(matches==0)throw std::invalid_argument("hydraulic boundary selects no mesh nodes or exterior faces: "+boundary.name);
    }
    if(prescribed.empty())throw std::invalid_argument("steady seepage requires at least one head, pore-pressure, or drain boundary");

    // Eliminate prescribed heads while preserving symmetry and positive definiteness.
    for(std::size_t row=0;row<matrix.size();++row){
        if(const auto fixed=prescribed.find(row);fixed!=prescribed.end())continue;
        for(auto entry=matrix[row].begin();entry!=matrix[row].end();){
            const auto fixed=prescribed.find(entry->first);
            if(fixed!=prescribed.end()){
                rhs[row]-=entry->second*fixed->second;
                entry=matrix[row].erase(entry);
            }else ++entry;
        }
    }
    for(const auto& [node,head]:prescribed){matrix[node].clear();matrix[node][node]=1.0;rhs[node]=head;}

    const double initial=std::accumulate(prescribed.begin(),prescribed.end(),0.0,
        [](double sum,const auto& item){return sum+item.second;})/prescribed.size();
    std::vector<double> head(model.nodes.size(),initial);
    for(const auto& [node,value]:prescribed)head[node]=value;
    const auto cg=conjugate_gradient(matrix,rhs,head,config.max_iterations,config.relative_tolerance,prescribed);

    HydraulicState state;
    state.valid=true;state.converged=cg.converged;state.iterations=cg.iterations;state.relative_residual=cg.relative_residual;
    state.head_m=head;
    state.nodal_pore_pressure_pa.resize(model.nodes.size());
    for(std::size_t i=0;i<model.nodes.size();++i)
        state.nodal_pore_pressure_pa[i]=rho_g*(head[i]-model.nodes[i].position.z);
    state.zone_head_m.resize(model.zones.size());
    state.zone_darcy_velocity_m_s.resize(model.zones.size());
    for(std::size_t zi=0;zi<model.zones.size();++zi){
        const auto& zone=model.zones[zi];if(!zone.active)continue;
        const auto geometry=tet_data(model,zone);Vec3 gradient{};double average=0.0;
        for(std::size_t local=0;local<4;++local){
            average+=head[zone.nodes[local]]*0.25;
            gradient=add(gradient,scale(geometry.gradient[local],head[zone.nodes[local]]));
        }
        state.zone_head_m[zi]=average;
        state.zone_darcy_velocity_m_s[zi]=scale(multiply(zone_tensor(derived.conductivity,zi),gradient),-1.0);
    }
    for(const auto& face:faces){
        const double outward=dot3(state.zone_darcy_velocity_m_s[face.owner_zone],face.normal)*face.area;
        if(outward>=0.0)state.outflow_m3_s+=outward;else state.inflow_m3_s-=outward;
    }
    state.net_outflow_m3_s=state.outflow_m3_s-state.inflow_m3_s;
    state.relative_mass_balance_error=std::abs(state.net_outflow_m3_s)/
        std::max({state.inflow_m3_s,state.outflow_m3_s,kSmall});
    state.zone_conductivity_m_s=derived.conductivity;
    state.zone_joint_aperture_m=derived.aperture;
    state.zone_joint_normal_stress_mpa=derived.normal_stress_mpa;
    model.hydraulic=state;
    return state;
}

HydraulicState solve_flac_no_capillarity_seepage(RockModel& model) {
    const auto& config=model.config.seepage;
    if(!config.enabled||!config.transient||!config.unsaturated||!config.flac_no_capillarity)
        throw std::invalid_argument("FLAC no-capillarity seepage is not enabled");
    if(!std::isfinite(model.config.gravity_m_s2)||model.config.gravity_m_s2<=0.0||
       !std::isfinite(config.fluid_density_kg_m3)||config.fluid_density_kg_m3<=0.0||
       !std::isfinite(config.initial_saturation)||config.initial_saturation<0.0||config.initial_saturation>1.0||
       !std::isfinite(config.duration_s)||config.duration_s<=0.0||
       !std::isfinite(config.time_step_s)||config.time_step_s<=0.0)
        throw std::invalid_argument("invalid FLAC no-capillarity initial state or time controls");
    const std::size_t material_count=model.config.region_materials.size()+1;
    if(config.conductivity_m_s.size()!=material_count||
       config.unsaturated_materials.size()!=material_count||
       config.specific_storage_1_m.size()!=material_count)
        throw std::invalid_argument("FLAC no-capillarity seepage requires conductivity, porosity, and storage for every material");
    for(std::size_t i=0;i<material_count;++i){
        validate_tensor(config.conductivity_m_s[i]);
        const double porosity=config.unsaturated_materials[i].porosity;
        if(!std::isfinite(porosity)||porosity<=0.0||porosity>1.0||
           !std::isfinite(config.specific_storage_1_m[i])||config.specific_storage_1_m[i]<=0.0)
            throw std::invalid_argument("FLAC no-capillarity porosity and specific storage must be positive");
    }
    for(const auto& boundary:config.boundaries)
        if(boundary.kind!=HydraulicBoundaryKind::Flux&&boundary.kind!=HydraulicBoundaryKind::Rainfall&&
           boundary.kind!=HydraulicBoundaryKind::NoFlow)
            throw std::invalid_argument("FLAC no-capillarity branch currently supports flux, rainfall, and no-flow boundaries only");

    const auto derived=derive_zone_hydraulics(model);
    const auto faces=boundary_faces(model);
    const auto box=bounds(model);
    const double rho_g=config.fluid_density_kg_m3*model.config.gravity_m_s2;
    const std::size_t count=model.nodes.size();
    std::vector<double> volume(count,0.0),porosity_sum(count,0.0),storage_sum(count,0.0);
    for(const auto& zone:model.zones)if(zone.active){
        const double weight=tet_data(model,zone).volume*zone.integration_weight*0.25;
        const double porosity=config.unsaturated_materials[zone.material_index].porosity;
        const double storage=config.specific_storage_1_m[zone.material_index];
        for(const auto node:zone.nodes){
            volume[node]+=weight;
            porosity_sum[node]+=porosity*weight;
            storage_sum[node]+=storage*weight;
        }
    }
    std::vector<double> porosity(count),storage(count),saturation(count,config.initial_saturation),pressure(count,0.0);
    for(std::size_t i=0;i<count;++i){
        if(!(volume[i]>0.0))throw std::invalid_argument("FLAC no-capillarity mesh has a node without fluid volume");
        porosity[i]=porosity_sum[i]/volume[i];
        storage[i]=storage_sum[i]/volume[i];
        if(config.initial_saturation==1.0)
            pressure[i]=rho_g*std::max(0.0,config.initial_head_m-model.nodes[i].position.z);
    }
    HydraulicState state;
    state.valid=true;state.converged=true;
    state.zone_conductivity_m_s=derived.conductivity;
    state.zone_joint_aperture_m=derived.aperture;
    state.zone_joint_normal_stress_mpa=derived.normal_stress_mpa;
    const auto populate=[&](){
        state.nodal_saturation=saturation;
        state.nodal_pore_pressure_pa=pressure;
        state.head_m.resize(count);
        for(std::size_t i=0;i<count;++i)state.head_m[i]=model.nodes[i].position.z+pressure[i]/rho_g;
        state.zone_head_m.assign(model.zones.size(),0.0);
        state.zone_saturation.assign(model.zones.size(),0.0);
        state.zone_darcy_velocity_m_s.assign(model.zones.size(),Vec3{});
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];if(!zone.active)continue;
            const auto geometry=tet_data(model,zone);
            double mean_s=0.0;
            for(const auto node:zone.nodes)mean_s+=saturation[node]*0.25;
            Vec3 gradient{};
            for(std::size_t j=0;j<4;++j){
                const auto node=zone.nodes[j];
                gradient=add(gradient,scale(geometry.gradient[j],
                    pressure[node]/rho_g+mean_s*model.nodes[node].position.z));
                state.zone_head_m[zi]+=state.head_m[node]*0.25;
            }
            double upstream_s=0.0;int upstream_count=0;
            for(std::size_t j=0;j<4;++j)
                if(dot3(geometry.gradient[j],multiply(derived.conductivity[zi],gradient))>0.0){
                    upstream_s+=saturation[zone.nodes[j]];++upstream_count;
                }
            const double upstream=upstream_count?upstream_s/upstream_count:mean_s;
            const double kr=upstream*upstream*(3.0-2.0*upstream);
            state.zone_saturation[zi]=mean_s;
            // FLAC3D 7 zone.flow is a zone diagnostic distinct from its
            // conservative, upstream-weighted nodal Q.  Direct 601 runs at
            // s=0.2 and a varying s≈0.4 field give q=-K*sbar^2*e_z exactly
            // when every node has p=0.  Keep the cubic law for nodal updates.
            const bool zero_pressure=std::all_of(zone.nodes.begin(),zone.nodes.end(),
                [&](std::size_t node){return pressure[node]==0.0;});
            state.zone_darcy_velocity_m_s[zi]=zero_pressure?
                scale(multiply(derived.conductivity[zi],Vec3{0.0,0.0,1.0}),-mean_s*mean_s):
                scale(multiply(derived.conductivity[zi],gradient),-kr);
        }
        // FLAC reports one zone.flow per physical brick/prism, using its
        // physical-node mean before the nonlinear square.  Averaging squared
        // tetrahedral values would otherwise introduce a Jensen bias.
        for(const auto& cell:model.physical_cells){
            if(cell.nodes.empty()||cell.subzones.empty())continue;
            const bool zero_pressure=std::all_of(cell.nodes.begin(),cell.nodes.end(),
                [&](std::size_t node){return pressure[node]==0.0;});
            if(!zero_pressure)continue;
            double mean_s=0.0,tensor_weight=0.0;
            Tensor tensor{};
            for(const auto node:cell.nodes)mean_s+=saturation[node]/cell.nodes.size();
            for(const auto zi:cell.subzones){
                const auto& zone=model.zones[zi];
                if(!zone.active)continue;
                const double weight=tet_data(model,zone).volume*zone.integration_weight;
                tensor_weight+=weight;
                for(std::size_t k=0;k<tensor.size();++k)tensor[k]+=derived.conductivity[zi][k]*weight;
            }
            if(!(tensor_weight>0.0))continue;
            for(auto& value:tensor)value/=tensor_weight;
            const Vec3 diagnostic=scale(multiply(tensor,Vec3{0.0,0.0,1.0}),-mean_s*mean_s);
            for(const auto zi:cell.subzones)state.zone_darcy_velocity_m_s[zi]=diagnostic;
        }
    };
    const auto snapshot=[&](double storage_change){
        HydraulicSnapshot out;
        out.time_s=state.time_s;out.iterations=state.iterations;
        out.relative_residual=state.relative_residual;
        out.inflow_m3_s=state.inflow_m3_s;out.outflow_m3_s=state.outflow_m3_s;
        out.net_outflow_m3_s=state.net_outflow_m3_s;
        out.storage_change_m3=storage_change;
        out.relative_mass_balance_error=state.relative_mass_balance_error;
        out.rainfall_infiltration_m3_s=state.rainfall_infiltration_m3_s;
        out.rainfall_runoff_m3_s=state.rainfall_runoff_m3_s;
        out.head_m=state.head_m;out.nodal_pore_pressure_pa=state.nodal_pore_pressure_pa;
        out.nodal_saturation=state.nodal_saturation;out.zone_head_m=state.zone_head_m;
        out.zone_saturation=state.zone_saturation;out.zone_darcy_velocity_m_s=state.zone_darcy_velocity_m_s;
        return out;
    };
    populate();state.history.push_back(snapshot(0.0));
    double time=0.0;
    while(time+config.duration_s*1.0e-12<config.duration_s){
        const double dt=std::min(config.time_step_s,config.duration_s-time);
        const double next_time=time+dt;
        std::vector<double> nodal_outflow(count,0.0),load(count,0.0);
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];if(!zone.active)continue;
            const auto geometry=tet_data(model,zone);
            double mean_s=0.0;
            for(const auto node:zone.nodes)mean_s+=saturation[node]*0.25;
            Vec3 gradient{};
            for(std::size_t j=0;j<4;++j){
                const auto node=zone.nodes[j];
                gradient=add(gradient,scale(geometry.gradient[j],
                    pressure[node]/rho_g+mean_s*model.nodes[node].position.z));
            }
            const auto tensor_flux=multiply(derived.conductivity[zi],gradient);
            std::array<double,4> local_outflow{};
            double inflow_s=0.0;int inflow_count=0;
            for(std::size_t j=0;j<4;++j){
                local_outflow[j]=geometry.volume*zone.integration_weight*
                    dot3(geometry.gradient[j],tensor_flux);
                if(local_outflow[j]>0.0){inflow_s+=saturation[zone.nodes[j]];++inflow_count;}
            }
            const double upstream=inflow_count?inflow_s/inflow_count:mean_s;
            const double kr=upstream*upstream*(3.0-2.0*upstream);
            for(std::size_t j=0;j<4;++j)nodal_outflow[zone.nodes[j]]+=kr*local_outflow[j];
        }
        state.rainfall_infiltration_m3_s=0.0;state.rainfall_runoff_m3_s=0.0;
        for(const auto& boundary:config.boundaries){
            if(boundary.kind==HydraulicBoundaryKind::NoFlow)continue;
            const double value=boundary_value_at(boundary,next_time);
            std::size_t matches=0;
            for(const auto& face:faces){
                bool face_selected=true;
                for(const auto node:face.nodes)
                    face_selected&=selected(model.nodes[node].position,boundary,box);
                if(!face_selected||(boundary.kind==HydraulicBoundaryKind::Rainfall&&face.normal.z<=0.0))continue;
                ++matches;
                double inward=boundary.kind==HydraulicBoundaryKind::Flux?-value*face.area:
                    std::min(std::max(0.0,value),derived.conductivity[face.owner_zone][2])*face.normal.z*face.area;
                for(const auto node:face.nodes)load[node]+=inward/3.0;
                if(boundary.kind==HydraulicBoundaryKind::Rainfall){
                    state.rainfall_infiltration_m3_s+=inward;
                    state.rainfall_runoff_m3_s+=std::max(0.0,value)*face.normal.z*face.area-inward;
                }
            }
            if(matches==0)throw std::invalid_argument("FLAC no-capillarity boundary selects no exterior faces: "+boundary.name);
        }
        double storage_change=0.0;
        for(std::size_t i=0;i<count;++i){
            const double old_water=volume[i]*(porosity[i]*saturation[i]+storage[i]*pressure[i]/rho_g);
            const double new_water=old_water+dt*(load[i]-nodal_outflow[i]);
            if(new_water< -1.0e-12)throw std::runtime_error("FLAC no-capillarity explicit time step depleted a node below zero saturation");
            const double content=std::max(0.0,new_water)/volume[i];
            if(content<porosity[i]){
                saturation[i]=content/porosity[i];pressure[i]=0.0;
            }else{
                saturation[i]=1.0;pressure[i]=rho_g*(content-porosity[i])/storage[i];
            }
            storage_change+=new_water-old_water;
        }
        state.inflow_m3_s=0.0;state.outflow_m3_s=0.0;
        for(const double value:load){if(value>=0.0)state.inflow_m3_s+=value;else state.outflow_m3_s-=value;}
        state.net_outflow_m3_s=state.outflow_m3_s-state.inflow_m3_s;
        state.relative_mass_balance_error=std::abs(storage_change+state.net_outflow_m3_s*dt)/
            std::max({std::abs(storage_change),std::abs(state.net_outflow_m3_s*dt),1.0e-20});
        state.time_s=next_time;state.iterations=1;state.relative_residual=0.0;
        populate();state.history.push_back(snapshot(storage_change));
        time=next_time;
    }
    model.hydraulic=state;
    return state;
}

HydraulicState solve_unsaturated_seepage(RockModel& model) {
    const auto& config=model.config.seepage;
    if(config.flac_no_capillarity)return solve_flac_no_capillarity_seepage(model);
    if(!config.enabled||!config.transient||!config.unsaturated)
        throw std::invalid_argument("unsaturated transient seepage is not enabled");
    if(!std::isfinite(config.fluid_density_kg_m3)||config.fluid_density_kg_m3<=0.0||
       !std::isfinite(model.config.gravity_m_s2)||model.config.gravity_m_s2<=0.0)
        throw std::invalid_argument("unsaturated seepage requires positive fluid density and gravity");
    const std::size_t material_count=model.config.region_materials.size()+1;
    if(config.conductivity_m_s.size()!=material_count||
       config.specific_storage_1_m.size()!=material_count||
       config.unsaturated_materials.size()!=material_count)
        throw std::invalid_argument("unsaturated seepage requires conductivity, storage, and retention parameters for every material");
    for(const auto& tensor:config.conductivity_m_s)validate_tensor(tensor);
    if(std::any_of(config.specific_storage_1_m.begin(),config.specific_storage_1_m.end(),
                   [](double value){return !std::isfinite(value)||value<=0.0;}))
        throw std::invalid_argument("unsaturated seepage requires positive specific storage for every material");
    for(std::size_t i=0;i<material_count;++i)
        (void)unsaturated_point(config.unsaturated_materials[i],config.specific_storage_1_m[i],-1.0);
    if(!std::isfinite(config.duration_s)||config.duration_s<=0.0||
       !std::isfinite(config.time_step_s)||config.time_step_s<=0.0||
       !std::isfinite(config.initial_head_m)||config.nonlinear_max_iterations<=0||
       !std::isfinite(config.nonlinear_relative_tolerance)||config.nonlinear_relative_tolerance<=0.0)
        throw std::invalid_argument("unsaturated seepage time and nonlinear controls are invalid");

    const auto derived=derive_zone_hydraulics(model);
    const auto box=bounds(model);
    const auto faces=boundary_faces(model);
    const double rho_g=config.fluid_density_kg_m3*model.config.gravity_m_s2;
    std::vector<std::size_t> seepage_face_nodes;
    for(const auto& boundary:config.boundaries)
        if(boundary.kind==HydraulicBoundaryKind::SeepageFace)
            for(std::size_t i=0;i<model.nodes.size();++i)
                if(selected(model.nodes[i].position,boundary,box))seepage_face_nodes.push_back(i);
    std::sort(seepage_face_nodes.begin(),seepage_face_nodes.end());
    seepage_face_nodes.erase(std::unique(seepage_face_nodes.begin(),seepage_face_nodes.end()),seepage_face_nodes.end());
    std::vector<bool> seepage_face_active(model.nodes.size(),false);
    std::vector<std::size_t> rainfall_nodes;
    for(const auto& boundary:config.boundaries){
        if(boundary.kind!=HydraulicBoundaryKind::Rainfall)continue;
        for(const auto& face:faces){
            bool matches=true;
            for(const auto node:face.nodes)
                matches&=selected(model.nodes[node].position,boundary,box);
            if(matches&&face.normal.z>0.0)
                rainfall_nodes.insert(rainfall_nodes.end(),face.nodes.begin(),face.nodes.end());
        }
    }
    std::sort(rainfall_nodes.begin(),rainfall_nodes.end());
    rainfall_nodes.erase(std::unique(rainfall_nodes.begin(),rainfall_nodes.end()),rainfall_nodes.end());
    std::vector<bool> rainfall_ponded(model.nodes.size(),false);

    const auto boundary_at=[&](double time_s,std::vector<double>* load,
                               std::vector<double>* rainfall_capacity,
                               std::vector<double>* rainfall_precipitation){
        std::map<std::size_t,double> prescribed;
        for(const auto& boundary:config.boundaries){
            const double value=boundary_value_at(boundary,time_s);
            std::size_t matches=0;
            if(boundary.kind==HydraulicBoundaryKind::Head||
               boundary.kind==HydraulicBoundaryKind::PorePressure||
               boundary.kind==HydraulicBoundaryKind::Drain||
               boundary.kind==HydraulicBoundaryKind::SeepageFace){
                for(std::size_t i=0;i<model.nodes.size();++i)if(selected(model.nodes[i].position,boundary,box)){
                    ++matches;
                    if(boundary.kind==HydraulicBoundaryKind::SeepageFace&&!seepage_face_active[i])continue;
                    const double head=boundary.kind==HydraulicBoundaryKind::Head?value:
                        model.nodes[i].position.z+
                        (boundary.kind==HydraulicBoundaryKind::Drain||
                         boundary.kind==HydraulicBoundaryKind::SeepageFace?0.0:value)/rho_g;
                    auto [entry,inserted]=prescribed.emplace(i,head);
                    if(!inserted&&std::abs(entry->second-head)>std::max(1.0,std::abs(head))*1.0e-10)
                        throw std::invalid_argument("conflicting unsaturated hydraulic heads on one node");
                }
            }else{
                for(const auto& face:faces){
                    bool face_selected=true;
                    for(const auto node:face.nodes)
                        face_selected&=selected(model.nodes[node].position,boundary,box);
                    if(!face_selected)continue;
                    if(boundary.kind==HydraulicBoundaryKind::Rainfall&&face.normal.z<=0.0)continue;
                    ++matches;
                    double outward=0.0;
                    if(boundary.kind==HydraulicBoundaryKind::Flux)outward=value;
                    if(boundary.kind==HydraulicBoundaryKind::Rainfall){
                        const double precipitation=std::max(0.0,value);
                        // Book Eq. (13.11)-(13.12): W=min(F,K33).  K33 is
                        // the saturated vertical conductivity, including the
                        // Stage 2.5 SMRM fracture contribution.
                        const double capacity=std::max(0.0,derived.conductivity[face.owner_zone][2]);
                        const double infiltration=std::min(precipitation,capacity);
                        const double projection=face.normal.z;
                        outward=-infiltration*projection;
                        for(const auto node:face.nodes){
                            if(rainfall_capacity)(*rainfall_capacity)[node]+=infiltration*projection*face.area/3.0;
                            if(rainfall_precipitation)(*rainfall_precipitation)[node]+=precipitation*projection*face.area/3.0;
                            if(rainfall_ponded[node]){
                                const double head=model.nodes[node].position.z;
                                auto [entry,inserted]=prescribed.emplace(node,head);
                                if(!inserted&&std::abs(entry->second-head)>1.0e-10)
                                    throw std::invalid_argument("conflicting ponded rainfall head on one node");
                            }
                        }
                    }
                    if(load&&(boundary.kind==HydraulicBoundaryKind::Flux||
                              boundary.kind==HydraulicBoundaryKind::Rainfall)){
                        const double nodal=-outward*face.area/3.0;
                        for(const auto node:face.nodes)
                            if(boundary.kind!=HydraulicBoundaryKind::Rainfall||!rainfall_ponded[node])
                                (*load)[node]+=nodal;
                    }
                }
            }
            if(matches==0)
                throw std::invalid_argument("unsaturated hydraulic boundary selects no mesh nodes or exterior faces: "+boundary.name);
        }
        return prescribed;
    };

    const auto zone_pressure_head=[&](const Zone& zone,const std::vector<double>& head){
        double value=0.0;
        for(const auto node:zone.nodes)value+=head[node]-model.nodes[node].position.z;
        return value/4.0;
    };
    const auto make_effective_tensors=[&](const std::vector<double>& head){
        std::vector<Tensor> tensors(model.zones.size());
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];
            const auto point=unsaturated_point(config.unsaturated_materials.at(zone.material_index),
                config.specific_storage_1_m.at(zone.material_index),zone_pressure_head(zone,head));
            tensors[zi]=scaled_tensor(derived.conductivity[zi],point.relative_conductivity);
        }
        return tensors;
    };
    const auto populate_state=[&](HydraulicState& state,const std::vector<double>& head,
                                  const std::vector<Tensor>& tensors){
        state.head_m=head;
        state.nodal_pore_pressure_pa.resize(model.nodes.size());
        state.nodal_saturation.assign(model.nodes.size(),0.0);
        std::vector<double> nodal_weight(model.nodes.size(),0.0);
        for(std::size_t i=0;i<model.nodes.size();++i)
            state.nodal_pore_pressure_pa[i]=rho_g*std::max(0.0,head[i]-model.nodes[i].position.z);
        state.zone_head_m.assign(model.zones.size(),0.0);
        state.zone_saturation.assign(model.zones.size(),0.0);
        state.zone_darcy_velocity_m_s.assign(model.zones.size(),Vec3{});
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];if(!zone.active)continue;
            const auto geometry=tet_data(model,zone);
            Vec3 gradient{};double average=0.0;
            for(std::size_t local=0;local<4;++local){
                average+=head[zone.nodes[local]]*0.25;
                gradient=add(gradient,scale(geometry.gradient[local],head[zone.nodes[local]]));
            }
            const auto point=unsaturated_point(config.unsaturated_materials.at(zone.material_index),
                config.specific_storage_1_m.at(zone.material_index),zone_pressure_head(zone,head));
            state.zone_head_m[zi]=average;
            state.zone_saturation[zi]=point.saturation;
            state.zone_darcy_velocity_m_s[zi]=scale(multiply(tensors[zi],gradient),-1.0);
            const double weight=geometry.volume*zone.integration_weight;
            for(const auto node:zone.nodes){
                state.nodal_saturation[node]+=point.saturation*weight;
                nodal_weight[node]+=weight;
            }
        }
        for(std::size_t i=0;i<state.nodal_saturation.size();++i)
            if(nodal_weight[i]>0.0)state.nodal_saturation[i]/=nodal_weight[i];
        state.zone_conductivity_m_s=tensors;
        state.zone_joint_aperture_m=derived.aperture;
        state.zone_joint_normal_stress_mpa=derived.normal_stress_mpa;
    };
    const auto make_snapshot=[](const HydraulicState& state,double storage_change){
        HydraulicSnapshot snapshot;
        snapshot.time_s=state.time_s;snapshot.iterations=state.iterations;
        snapshot.relative_residual=state.relative_residual;
        snapshot.inflow_m3_s=state.inflow_m3_s;snapshot.outflow_m3_s=state.outflow_m3_s;
        snapshot.net_outflow_m3_s=state.net_outflow_m3_s;
        snapshot.storage_change_m3=storage_change;
        snapshot.relative_mass_balance_error=state.relative_mass_balance_error;
        snapshot.rainfall_infiltration_m3_s=state.rainfall_infiltration_m3_s;
        snapshot.rainfall_runoff_m3_s=state.rainfall_runoff_m3_s;
        snapshot.head_m=state.head_m;snapshot.nodal_pore_pressure_pa=state.nodal_pore_pressure_pa;
        snapshot.nodal_saturation=state.nodal_saturation;snapshot.zone_head_m=state.zone_head_m;
        snapshot.zone_saturation=state.zone_saturation;
        snapshot.zone_darcy_velocity_m_s=state.zone_darcy_velocity_m_s;
        return snapshot;
    };

    std::vector<double> head(model.nodes.size(),config.initial_head_m);
    for(const auto& [node,value]:boundary_at(0.0,nullptr,nullptr,nullptr))head[node]=value;
    HydraulicState state;state.valid=true;state.converged=true;state.time_s=0.0;
    auto effective_tensors=make_effective_tensors(head);
    populate_state(state,head,effective_tensors);
    state.history.push_back(make_snapshot(state,0.0));

    double time=0.0;
    std::size_t step_index=0;
    while(time<config.duration_s-std::max(1.0,config.duration_s)*1.0e-14){
        const double next_time=std::min(config.duration_s,
            static_cast<double>(++step_index)*config.time_step_s);
        const double dt=next_time-time;
        const auto old_head=head;
        std::vector<double> guess=head;
        CgResult final_cg{};
        SparseRows final_conductivity,final_capacity;
        std::vector<double> final_boundary_load;
        std::vector<double> final_rainfall_capacity,final_rainfall_precipitation;
        std::map<std::size_t,double> final_prescribed;
        double nonlinear_residual=std::numeric_limits<double>::infinity();
        int linear_iterations=0;
        bool nonlinear_converged=false;
        for(int nonlinear=0;nonlinear<config.nonlinear_max_iterations;++nonlinear){
            SparseRows conductivity(model.nodes.size()),capacity(model.nodes.size());
            auto tensors=make_effective_tensors(guess);
            for(std::size_t zi=0;zi<model.zones.size();++zi){
                const auto& zone=model.zones[zi];if(!zone.active)continue;
                const auto geometry=tet_data(model,zone);
                const double volume=geometry.volume*zone.integration_weight;
                const double old_psi=zone_pressure_head(zone,old_head);
                const double new_psi=zone_pressure_head(zone,guess);
                const auto& parameters=config.unsaturated_materials.at(zone.material_index);
                const double storage=config.specific_storage_1_m.at(zone.material_index);
                const auto old_point=unsaturated_point(parameters,storage,old_psi);
                const auto new_point=unsaturated_point(parameters,storage,new_psi);
                const double delta=new_psi-old_psi;
                const double secant=std::abs(delta)>1.0e-12?
                    (new_point.water_content-old_point.water_content)/delta:new_point.capacity_1_m;
                for(std::size_t i=0;i<4;++i)for(std::size_t j=0;j<4;++j){
                    conductivity[zone.nodes[i]][zone.nodes[j]]+=volume*
                        dot3(geometry.gradient[i],multiply(tensors[zi],geometry.gradient[j]));
                    capacity[zone.nodes[i]][zone.nodes[j]]+=std::max(secant,1.0e-16)*volume*
                        (i==j?2.0:1.0)/20.0;
                }
            }
            SparseRows matrix=conductivity;
            std::vector<double> rhs(model.nodes.size(),0.0);
            std::vector<double> rainfall_capacity(model.nodes.size(),0.0);
            std::vector<double> rainfall_precipitation(model.nodes.size(),0.0);
            const auto prescribed=boundary_at(next_time,&rhs,&rainfall_capacity,&rainfall_precipitation);
            const auto boundary_load=rhs;
            for(std::size_t i=0;i<matrix.size();++i)
                for(const auto& [j,value]:capacity[i]){
                    matrix[i][j]+=value/dt;
                    rhs[i]+=value/dt*old_head[j];
                }
            for(std::size_t row=0;row<matrix.size();++row){
                if(prescribed.find(row)!=prescribed.end())continue;
                for(auto entry=matrix[row].begin();entry!=matrix[row].end();){
                    const auto fixed=prescribed.find(entry->first);
                    if(fixed!=prescribed.end()){
                        rhs[row]-=entry->second*fixed->second;
                        entry=matrix[row].erase(entry);
                    }else ++entry;
                }
            }
            for(const auto& [node,value]:prescribed){matrix[node].clear();matrix[node][node]=1.0;rhs[node]=value;}
            auto next=guess;
            for(const auto& [node,value]:prescribed)next[node]=value;
            final_cg=conjugate_gradient(matrix,rhs,next,config.max_iterations,config.relative_tolerance,prescribed);
            linear_iterations+=final_cg.iterations;
            bool boundary_sets_stable=true;
            if(!seepage_face_nodes.empty()||!rainfall_nodes.empty()){
                const auto conductivity_next=multiply(conductivity,next);
                std::vector<double> head_delta(next.size());
                for(std::size_t i=0;i<next.size();++i)head_delta[i]=next[i]-old_head[i];
                const auto capacity_delta=multiply(capacity,head_delta);
                for(const auto node:seepage_face_nodes){
                    const double pressure_head=next[node]-model.nodes[node].position.z;
                    const double outward=-(conductivity_next[node]+capacity_delta[node]/dt-boundary_load[node]);
                    const bool active=seepage_face_active[node]?
                        outward>=-1.0e-12:pressure_head>1.0e-10;
                    boundary_sets_stable&=active==seepage_face_active[node];
                    seepage_face_active[node]=active;
                }
                for(const auto node:rainfall_nodes){
                    const double pressure_head=next[node]-model.nodes[node].position.z;
                    const double inward=conductivity_next[node]+capacity_delta[node]/dt-boundary_load[node];
                    const double tolerance=std::max(1.0e-15,rainfall_capacity[node]*1.0e-6);
                    const bool active=rainfall_ponded[node]?
                        inward<=rainfall_capacity[node]+tolerance:pressure_head>1.0e-10;
                    boundary_sets_stable&=active==rainfall_ponded[node];
                    rainfall_ponded[node]=active;
                }
            }
            double difference=0.0,norm_value=0.0;
            for(std::size_t i=0;i<next.size();++i){
                difference+=std::pow(next[i]-guess[i],2);
                norm_value+=std::pow(next[i],2);
            }
            nonlinear_residual=std::sqrt(difference/std::max(norm_value,kSmall));
            guess=std::move(next);
            final_conductivity=std::move(conductivity);final_capacity=std::move(capacity);
            final_boundary_load=boundary_load;final_prescribed=prescribed;
            final_rainfall_capacity=std::move(rainfall_capacity);
            final_rainfall_precipitation=std::move(rainfall_precipitation);
            if(final_cg.converged&&boundary_sets_stable&&
               nonlinear_residual<=config.nonlinear_relative_tolerance){
                nonlinear_converged=true;break;
            }
        }
        head=guess;
        effective_tensors=make_effective_tensors(head);
        state.converged=state.converged&&nonlinear_converged;
        state.iterations=linear_iterations;
        state.relative_residual=std::max(final_cg.relative_residual,nonlinear_residual);
        state.time_s=next_time;
        populate_state(state,head,effective_tensors);
        const auto conductivity_head=multiply(final_conductivity,head);
        std::vector<double> head_change(head.size(),0.0);
        for(std::size_t i=0;i<head.size();++i)head_change[i]=head[i]-old_head[i];
        const auto stored_head=multiply(final_capacity,head_change);
        state.rainfall_infiltration_m3_s=0.0;
        state.rainfall_runoff_m3_s=0.0;
        for(const auto node:rainfall_nodes){
            const double infiltration=rainfall_ponded[node]?
                std::max(0.0,conductivity_head[node]+stored_head[node]/dt-final_boundary_load[node]):
                final_rainfall_capacity[node];
            state.rainfall_infiltration_m3_s+=infiltration;
            state.rainfall_runoff_m3_s+=std::max(0.0,final_rainfall_precipitation[node]-infiltration);
        }
        state.inflow_m3_s=0.0;state.outflow_m3_s=0.0;
        const auto add_outward=[&](double outward){if(outward>=0.0)state.outflow_m3_s+=outward;else state.inflow_m3_s-=outward;};
        for(const double load:final_boundary_load)add_outward(-load);
        for(const auto& [node,value]:final_prescribed){
            (void)value;
            add_outward(-(conductivity_head[node]+stored_head[node]/dt-final_boundary_load[node]));
        }
        state.net_outflow_m3_s=state.outflow_m3_s-state.inflow_m3_s;
        double storage_change=0.0,domain_volume=0.0;
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];if(!zone.active)continue;
            const double volume=tet_data(model,zone).volume*zone.integration_weight;
            domain_volume+=volume;
            const auto& parameters=config.unsaturated_materials.at(zone.material_index);
            const double storage=config.specific_storage_1_m.at(zone.material_index);
            const auto old_point=unsaturated_point(parameters,storage,zone_pressure_head(zone,old_head));
            const auto new_point=unsaturated_point(parameters,storage,zone_pressure_head(zone,head));
            storage_change+=volume*(new_point.water_content-old_point.water_content);
        }
        const double boundary_volume=state.net_outflow_m3_s*dt;
        const double rainfall_volume=std::accumulate(final_rainfall_precipitation.begin(),
            final_rainfall_precipitation.end(),0.0)*dt;
        state.relative_mass_balance_error=std::abs(storage_change+boundary_volume)/
            std::max({std::abs(storage_change),std::abs(boundary_volume),
                      rainfall_volume,domain_volume*1.0e-14,kSmall});
        state.history.push_back(make_snapshot(state,storage_change));
        time=next_time;
        if(!nonlinear_converged)break;
    }
    model.hydraulic=state;
    return state;
}

HydraulicState solve_transient_saturated_impl(
    RockModel& model,const std::vector<double>* initial_head,
    const std::vector<double>* zone_volumetric_strain_increment,double time_offset_s,
    const std::vector<double>* zone_biot_modulus_pa=nullptr,
    const std::vector<double>* previous_zone_biot_modulus_pa=nullptr) {
    const auto& config=model.config.seepage;
    if(config.unsaturated)return solve_unsaturated_seepage(model);
    if(!config.enabled||!config.transient)
        throw std::invalid_argument("transient saturated seepage is not enabled");
    for(const auto& boundary:config.boundaries)
        if(boundary.kind==HydraulicBoundaryKind::Rainfall||
           boundary.kind==HydraulicBoundaryKind::SeepageFace)
            throw std::invalid_argument("rainfall and seepage-face boundaries require transient unsaturated seepage");
    if(!std::isfinite(config.fluid_density_kg_m3)||config.fluid_density_kg_m3<=0.0||
       !std::isfinite(model.config.gravity_m_s2)||model.config.gravity_m_s2<=0.0)
        throw std::invalid_argument("transient seepage requires positive fluid density and gravity");
    if(config.conductivity_m_s.size()!=model.config.region_materials.size()+1)
        throw std::invalid_argument("transient seepage requires one conductivity tensor for every material");
    for(const auto& tensor:config.conductivity_m_s)validate_tensor(tensor);
    if(config.specific_storage_1_m.size()!=model.config.region_materials.size()+1||
       std::any_of(config.specific_storage_1_m.begin(),config.specific_storage_1_m.end(),
                   [](double value){return !std::isfinite(value)||value<=0.0;}))
        throw std::invalid_argument("transient seepage requires positive specific storage for every material");
    if(!std::isfinite(config.biot_modulus_pa)||config.biot_modulus_pa<0.0)
        throw std::invalid_argument("Biot modulus must be nonnegative");
    if(!std::isfinite(config.duration_s)||config.duration_s<=0.0||
       !std::isfinite(config.time_step_s)||config.time_step_s<=0.0||
       !std::isfinite(config.initial_head_m))
        throw std::invalid_argument("transient seepage duration, time step, and initial head are invalid");
    if(config.max_iterations<=0||!std::isfinite(config.relative_tolerance)||config.relative_tolerance<=0.0)
        throw std::invalid_argument("transient seepage iteration controls are invalid");
    if(initial_head&&initial_head->size()!=model.nodes.size())
        throw std::invalid_argument("coupled seepage initial head size does not match mesh");
    if(zone_volumetric_strain_increment&&
       (zone_volumetric_strain_increment->size()!=model.zones.size()||
        std::abs(config.duration_s-config.time_step_s)>config.duration_s*1.0e-12))
        throw std::invalid_argument("coupled seepage strain source requires one increment per zone and one time step");
    if(zone_biot_modulus_pa&&
       (zone_biot_modulus_pa->size()!=model.zones.size()||
        std::any_of(zone_biot_modulus_pa->begin(),zone_biot_modulus_pa->end(),
            [](double value){return !std::isfinite(value)||value<=0.0;})))
        throw std::invalid_argument("evolving Biot modulus requires one positive value per zone");
    if(previous_zone_biot_modulus_pa&&
       (!zone_biot_modulus_pa||!initial_head||
        previous_zone_biot_modulus_pa->size()!=model.zones.size()||
        std::any_of(previous_zone_biot_modulus_pa->begin(),
            previous_zone_biot_modulus_pa->end(),
            [](double value){return !std::isfinite(value)||value<=0.0;})))
        throw std::invalid_argument("evolving storage source requires previous positive Biot moduli and head");

    const auto derived=derive_zone_hydraulics(model);
    const double rho_g=config.fluid_density_kg_m3*model.config.gravity_m_s2;
    SparseRows conductivity(model.nodes.size());
    SparseRows capacity(model.nodes.size());
    std::vector<double> lumped_capacity(model.nodes.size(),0.0);
    std::vector<double> mechanical_source_m3_s(model.nodes.size(),0.0);
    std::vector<double> storage_parameter_source_m3_s(model.nodes.size(),0.0);
    double mechanical_storage_change_m3=0.0;
    double storage_parameter_change_m3=0.0;
    for(std::size_t zi=0;zi<model.zones.size();++zi){
        const auto& zone=model.zones[zi];if(!zone.active)continue;
        const auto geometry=tet_data(model,zone);
        const auto& tensor=zone_tensor(derived.conductivity,zi);
        const double volume=geometry.volume*zone.integration_weight;
        const double modulus=zone_biot_modulus_pa?zone_biot_modulus_pa->at(zi):
            config.biot_modulus_pa;
        const double storage=modulus>0.0?
            rho_g/modulus:
            config.specific_storage_1_m.at(zone.material_index);
        if(previous_zone_biot_modulus_pa){
            const double change=1.0/modulus-1.0/previous_zone_biot_modulus_pa->at(zi);
            for(const auto node:zone.nodes){
                const double old_pressure=rho_g*(initial_head->at(node)-
                    model.nodes[node].position.z);
                const double content=volume*old_pressure*change/4.0;
                storage_parameter_change_m3+=content;
                storage_parameter_source_m3_s[node]-=content/config.duration_s;
            }
        }
        if(zone_volumetric_strain_increment){
            const double content_change=config.biot_coefficient*volume*
                zone_volumetric_strain_increment->at(zi);
            mechanical_storage_change_m3+=content_change;
            for(const auto node:zone.nodes)
                mechanical_source_m3_s[node]-=content_change/(4.0*config.duration_s);
        }
        for(std::size_t i=0;i<4;++i){
            lumped_capacity[zone.nodes[i]]+=storage*volume/4.0;
            for(std::size_t j=0;j<4;++j){
                conductivity[zone.nodes[i]][zone.nodes[j]]+=
                    volume*dot3(geometry.gradient[i],multiply(tensor,geometry.gradient[j]));
                capacity[zone.nodes[i]][zone.nodes[j]]+=storage*volume*(i==j?2.0:1.0)/20.0;
            }
        }
    }
    if(zone_volumetric_strain_increment){
        // FLAC3D's coupled continuity equation is nodal finite volume:
        // V_n=sum(V_t/4), storage=(rho*g/M)*V_n. Keep the existing
        // consistent-mass Stage 3 discretization for uncoupled runs.
        for(std::size_t i=0;i<capacity.size();++i){
            capacity[i].clear();
            capacity[i][i]=lumped_capacity[i];
        }
    }

    const auto box=bounds(model);
    const auto faces=boundary_faces(model);
    const auto prescribed_at=[&](double time_s,std::vector<double>* load){
        std::map<std::size_t,double> prescribed;
        for(const auto& boundary:config.boundaries){
            const double value=boundary_value_at(boundary,time_s);
            std::size_t matches=0;
            if(boundary.kind==HydraulicBoundaryKind::Head||
               boundary.kind==HydraulicBoundaryKind::PorePressure||
               boundary.kind==HydraulicBoundaryKind::Drain){
                for(std::size_t i=0;i<model.nodes.size();++i)if(selected(model.nodes[i].position,boundary,box)){
                    ++matches;
                    const double head=boundary.kind==HydraulicBoundaryKind::Head?value:
                        model.nodes[i].position.z+
                        (boundary.kind==HydraulicBoundaryKind::Drain?0.0:value)/rho_g;
                    auto [entry,inserted]=prescribed.emplace(i,head);
                    if(!inserted&&std::abs(entry->second-head)>std::max(1.0,std::abs(head))*1.0e-10)
                        throw std::invalid_argument("conflicting transient hydraulic heads on one node");
                }
            }else{
                for(const auto& face:faces){
                    bool face_selected=true;
                    for(const auto node:face.nodes)
                        face_selected&=selected(model.nodes[node].position,boundary,box);
                    if(!face_selected)continue;
                    ++matches;
                    if(load&&boundary.kind==HydraulicBoundaryKind::Flux){
                        const double nodal=-value*face.area/3.0;
                        for(const auto node:face.nodes)(*load)[node]+=nodal;
                    }
                }
            }
            if(matches==0)
                throw std::invalid_argument("transient hydraulic boundary selects no mesh nodes or exterior faces: "+boundary.name);
        }
        return prescribed;
    };

    auto populate_state=[&](HydraulicState& state,const std::vector<double>& head){
        state.head_m=head;
        state.nodal_saturation.assign(model.nodes.size(),1.0);
        state.zone_saturation.assign(model.zones.size(),1.0);
        state.nodal_pore_pressure_pa.resize(model.nodes.size());
        for(std::size_t i=0;i<model.nodes.size();++i)
            state.nodal_pore_pressure_pa[i]=rho_g*(head[i]-model.nodes[i].position.z);
        state.zone_head_m.assign(model.zones.size(),0.0);
        state.zone_darcy_velocity_m_s.assign(model.zones.size(),Vec3{});
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            const auto& zone=model.zones[zi];if(!zone.active)continue;
            const auto geometry=tet_data(model,zone);
            Vec3 gradient{};double average=0.0;
            for(std::size_t local=0;local<4;++local){
                average+=head[zone.nodes[local]]*0.25;
                gradient=add(gradient,scale(geometry.gradient[local],head[zone.nodes[local]]));
            }
            state.zone_head_m[zi]=average;
            state.zone_darcy_velocity_m_s[zi]=
                scale(multiply(zone_tensor(derived.conductivity,zi),gradient),-1.0);
        }
        state.inflow_m3_s=0.0;state.outflow_m3_s=0.0;
        for(const auto& face:faces){
            const double outward=dot3(state.zone_darcy_velocity_m_s[face.owner_zone],face.normal)*face.area;
            if(outward>=0.0)state.outflow_m3_s+=outward;else state.inflow_m3_s-=outward;
        }
        state.net_outflow_m3_s=state.outflow_m3_s-state.inflow_m3_s;
        state.zone_conductivity_m_s=derived.conductivity;
        state.zone_joint_aperture_m=derived.aperture;
        state.zone_joint_normal_stress_mpa=derived.normal_stress_mpa;
    };
    const auto make_snapshot=[](const HydraulicState& state,double storage_change){
        HydraulicSnapshot snapshot;
        snapshot.time_s=state.time_s;
        snapshot.iterations=state.iterations;
        snapshot.relative_residual=state.relative_residual;
        snapshot.inflow_m3_s=state.inflow_m3_s;
        snapshot.outflow_m3_s=state.outflow_m3_s;
        snapshot.net_outflow_m3_s=state.net_outflow_m3_s;
        snapshot.storage_change_m3=storage_change;
        snapshot.mechanical_storage_change_m3=state.mechanical_storage_change_m3;
        snapshot.relative_mass_balance_error=state.relative_mass_balance_error;
        snapshot.mass_balance_residual_m3=state.mass_balance_residual_m3;
        snapshot.mass_balance_roundoff_m3=state.mass_balance_roundoff_m3;
        snapshot.rainfall_infiltration_m3_s=state.rainfall_infiltration_m3_s;
        snapshot.rainfall_runoff_m3_s=state.rainfall_runoff_m3_s;
        snapshot.head_m=state.head_m;
        snapshot.nodal_boundary_outward_flux_m3_s=state.nodal_boundary_outward_flux_m3_s;
        snapshot.nodal_pore_pressure_pa=state.nodal_pore_pressure_pa;
        snapshot.nodal_saturation=state.nodal_saturation;
        snapshot.zone_head_m=state.zone_head_m;
        snapshot.zone_saturation=state.zone_saturation;
        snapshot.zone_darcy_velocity_m_s=state.zone_darcy_velocity_m_s;
        snapshot.zone_conductivity_m_s=state.zone_conductivity_m_s;
        snapshot.zone_joint_aperture_m=state.zone_joint_aperture_m;
        snapshot.zone_joint_normal_stress_mpa=state.zone_joint_normal_stress_mpa;
        return snapshot;
    };

    std::vector<double> head=initial_head?*initial_head:
        std::vector<double>(model.nodes.size(),config.initial_head_m);
    // A coupled step receives the actual previous-time head. Do not replace
    // it with a newly prescribed boundary value before assembling C*dh/dt:
    // that would erase the storage/reaction of a boundary pressure jump.
    if(!initial_head)
        for(const auto& [node,value]:prescribed_at(time_offset_s,nullptr))head[node]=value;

    HydraulicState state;
    state.valid=true;state.converged=true;state.time_s=time_offset_s;
    populate_state(state,head);
    state.nodal_boundary_outward_flux_m3_s.assign(head.size(),0.0);
    state.relative_mass_balance_error=0.0;
    state.history.push_back(make_snapshot(state,0.0));

    double time=0.0;
    std::size_t step_index=0;
    while(time<config.duration_s-std::max(1.0,config.duration_s)*1.0e-14){
        const double next_time=std::min(config.duration_s,
            static_cast<double>(++step_index)*config.time_step_s);
        const double dt=next_time-time;
        SparseRows matrix=conductivity;
        std::vector<double> rhs(model.nodes.size(),0.0);
        const auto prescribed=prescribed_at(time_offset_s+next_time,&rhs);
        const auto boundary_load=rhs;
        auto effective_mechanical_source=mechanical_source_m3_s;
        if(zone_volumetric_strain_increment){
            // FLAC3D 7 redistributes mechanical pore-volume change that would
            // land on fixed-pressure nodes to the free nodes of that physical
            // cell. Keep the total source unchanged for the boundary balance.
            for_each_coupled_cell(model,[&](const auto& nodes,const auto& subzones){
                std::size_t free_count=0;
                for(const auto node:nodes)
                    if(prescribed.find(node)==prescribed.end())++free_count;
                if(free_count==0)return;
                double redirected=0.0;
                for(const auto zi:subzones){
                    const auto& zone=model.zones[zi];
                    const double content=config.biot_coefficient*zone.volume_m3*
                        zone.integration_weight*zone_volumetric_strain_increment->at(zi);
                    const double contribution=-content/(4.0*dt);
                    for(const auto node:zone.nodes)
                        if(prescribed.find(node)!=prescribed.end()){
                            effective_mechanical_source[node]-=contribution;
                            redirected+=contribution;
                        }
                }
                for(const auto node:nodes)
                    if(prescribed.find(node)==prescribed.end())
                        effective_mechanical_source[node]+=redirected/static_cast<double>(free_count);
            });
        }
        for(std::size_t i=0;i<rhs.size();++i)
            rhs[i]+=effective_mechanical_source[i]+storage_parameter_source_m3_s[i];
        for(std::size_t i=0;i<matrix.size();++i)
            for(const auto& [j,value]:capacity[i]){
                matrix[i][j]+=value/dt;
                rhs[i]+=value/dt*head[j];
            }

        for(std::size_t row=0;row<matrix.size();++row){
            if(prescribed.find(row)!=prescribed.end())continue;
            for(auto entry=matrix[row].begin();entry!=matrix[row].end();){
                const auto fixed=prescribed.find(entry->first);
                if(fixed!=prescribed.end()){
                    rhs[row]-=entry->second*fixed->second;
                    entry=matrix[row].erase(entry);
                }else ++entry;
            }
        }
        for(const auto& [node,value]:prescribed){
            matrix[node].clear();matrix[node][node]=1.0;rhs[node]=value;
        }

        const auto old_head=head;
        for(const auto& [node,value]:prescribed)head[node]=value;
        const auto cg=conjugate_gradient(matrix,rhs,head,config.max_iterations,config.relative_tolerance,prescribed);
        state.converged=state.converged&&cg.converged;
        state.iterations=cg.iterations;
        state.relative_residual=cg.relative_residual;
        state.time_s=time_offset_s+next_time;
        populate_state(state,head);
        // Recover boundary discharge from the unmodified finite-element
        // equations.  A Tet4 pressure gradient is only elementwise constant,
        // so integrating it over exterior faces does not give a conservative
        // Dirichlet-boundary reaction on a transient mesh.
        state.inflow_m3_s=0.0;
        state.outflow_m3_s=0.0;
        state.nodal_boundary_outward_flux_m3_s.assign(head.size(),0.0);
        const auto conductivity_head=multiply(conductivity,head);
        std::vector<double> head_change(head.size(),0.0);
        for(std::size_t i=0;i<head.size();++i)head_change[i]=head[i]-old_head[i];
        const auto stored_head=multiply(capacity,head_change);
        const auto add_outward=[&](double outward){
            if(outward>=0.0)state.outflow_m3_s+=outward;
            else state.inflow_m3_s-=outward;
        };
        // Neumann load uses rhs = -q_out.  At prescribed-head nodes the
        // reaction is K h + C dh/dt - rhs, hence q_out = -reaction.
        for(std::size_t i=0;i<boundary_load.size();++i){
            state.nodal_boundary_outward_flux_m3_s[i]=-boundary_load[i];
            add_outward(-boundary_load[i]);
        }
        for(const auto& [node,value]:prescribed){
            (void)value;
            const double reaction=conductivity_head[node]+stored_head[node]/dt-
                boundary_load[node]-effective_mechanical_source[node]-
                storage_parameter_source_m3_s[node];
            state.nodal_boundary_outward_flux_m3_s[node]-=reaction;
            add_outward(-reaction);
        }
        state.net_outflow_m3_s=state.outflow_m3_s-state.inflow_m3_s;
        double storage_change=0.0;
        for(std::size_t i=0;i<head.size();++i)
            storage_change+=lumped_capacity[i]*(head[i]-old_head[i]);
        storage_change+=storage_parameter_change_m3;
        const double boundary_volume=state.net_outflow_m3_s*dt;
        state.mechanical_storage_change_m3=mechanical_storage_change_m3;
        state.mass_balance_residual_m3=std::abs(storage_change+
            mechanical_storage_change_m3+boundary_volume);
        // Floating-point cancellation floor, in m3, not a solver/physical
        // tolerance. Bound the row products before cancellation (C*h and
        // dt*K*h), using gamma_n = n*eps/(1-n*eps). Keep the raw relative
        // residual visible even when all physical increments approach zero.
        long double roundoff=0.0;
        for(std::size_t i=0;i<head.size();++i){
            long double terms=dt*(std::abs(boundary_load[i])+
                std::abs(effective_mechanical_source[i])+
                std::abs(storage_parameter_source_m3_s[i]));
            for(const auto& [j,value]:capacity[i])
                terms+=std::abs(value)*(std::abs(head[j])+std::abs(old_head[j]));
            for(const auto& [j,value]:conductivity[i])
                terms+=dt*std::abs(value*head[j]);
            const double neps=(16.0+2.0*(capacity[i].size()+conductivity[i].size()))*
                std::numeric_limits<double>::epsilon();
            roundoff+=neps/(1.0-neps)*terms;
        }
        state.mass_balance_roundoff_m3=static_cast<double>(roundoff);
        state.relative_mass_balance_error=state.mass_balance_residual_m3/
            std::max({std::abs(storage_change),std::abs(mechanical_storage_change_m3),
                      std::abs(boundary_volume),kSmall});
        state.history.push_back(make_snapshot(state,storage_change));
        time=next_time;
        if(!cg.converged)break;
    }
    model.hydraulic=state;
    return state;
}

HydraulicState solve_transient_seepage(RockModel& model) {
    return solve_transient_saturated_impl(model,nullptr,nullptr,0.0);
}

RockResult solve_coupled_saturated_seepage_mechanics(
    RockModel& model,const ProgressCallback& progress) {
    const auto config=model.config.seepage;
    if(!config.enabled||!config.transient||config.unsaturated||
       !config.solve_mechanics||!config.bidirectional_coupling||
       !(config.biot_modulus_pa>0.0)||!std::isfinite(config.biot_modulus_pa)||
       !std::isfinite(config.biot_coefficient)||config.biot_coefficient<0.0||
       config.biot_coefficient>1.0||config.coupling_max_iterations<1||
       !std::isfinite(config.coupling_relative_tolerance)||
       !(config.coupling_relative_tolerance>0.0))
        throw std::invalid_argument("bidirectional coupling requires saturated transient mechanics and valid Biot/iteration parameters");
    // Validate before the time loop: invalid/zero duration must not bypass
    // the per-step fluid validator and return a successful initial state.
    if(!std::isfinite(config.duration_s)||config.duration_s<=0.0||
       !std::isfinite(config.time_step_s)||config.time_step_s<=0.0||
       !std::isfinite(config.initial_head_m)||
       !std::isfinite(config.fluid_density_kg_m3)||config.fluid_density_kg_m3<=0.0||
       !std::isfinite(model.config.gravity_m_s2)||model.config.gravity_m_s2<=0.0)
        throw std::invalid_argument("transient seepage duration, time step, initial head, density or gravity is invalid");
    if((model.config.material.type!=MaterialType::Elastic&&
        model.config.material.type!=MaterialType::SMRM)||
       !model.config.region_materials.empty()||!model.cables.empty()||
       !model.interfaces.empty()||model.config.solver.large_strain)
        throw std::invalid_argument("bidirectional solver supports small-strain homogeneous elastic or SMRM continuum only");
    if(config.evolving_porosity_storage&&
       (std::abs(config.biot_coefficient-1.0)>1.0e-12||
        !std::isfinite(config.initial_porosity)||
        !(config.initial_porosity>0.0&&config.initial_porosity<1.0)))
        throw std::invalid_argument("evolving porosity/storage requires alpha=1 and 0<n0<1");

    const double rho_g=config.fluid_density_kg_m3*model.config.gravity_m_s2;
    HydraulicState initial;
    initial.valid=true;initial.converged=true;
    initial.head_m.assign(model.nodes.size(),config.initial_head_m);
    initial.nodal_pore_pressure_pa.resize(model.nodes.size());
    initial.nodal_saturation.assign(model.nodes.size(),1.0);
    for(std::size_t i=0;i<model.nodes.size();++i)
        initial.nodal_pore_pressure_pa[i]=rho_g*(config.initial_head_m-model.nodes[i].position.z);
    initial.zone_head_m.assign(model.zones.size(),config.initial_head_m);
    initial.zone_saturation.assign(model.zones.size(),1.0);
    initial.zone_darcy_velocity_m_s.assign(model.zones.size(),Vec3{});
    model.hydraulic=initial;
    model.stages={{"Coupled initial equilibrium",StageType::Initial}};
    auto previous_result=solve(model,progress);
    if(!previous_result.converged)
        throw std::runtime_error("bidirectional initial mechanical equilibrium did not converge");
    // The initial solve is quasi-static: its explicit relaxation velocity is
    // not part of the accepted physical state for the transient Biot steps.
    for(auto& node:model.nodes)node.velocity=Vec3{};

    std::vector<double> reference_signed_volume;
    if(config.evolving_porosity_storage){
        reference_signed_volume.resize(model.zones.size());
        for(std::size_t zi=0;zi<model.zones.size();++zi){
            reference_signed_volume[zi]=deformed_tet_signed_volume(model,model.zones[zi]);
            if(!std::isfinite(reference_signed_volume[zi])||
               std::abs(reference_signed_volume[zi])<=kSmall)
                throw std::invalid_argument("evolving porosity/storage found a degenerate reference cell");
        }
    }
    struct StorageState {
        std::vector<double> zone_porosity;
        std::vector<double> zone_modulus_pa;
        std::vector<double> nodal_modulus_pa;
    };
    const auto storage_at=[&](const RockModel& current){
        StorageState state;
        state.zone_modulus_pa.assign(current.zones.size(),config.biot_modulus_pa);
        if(!config.evolving_porosity_storage)return state;
        state.zone_porosity.assign(current.zones.size(),config.initial_porosity);
        const double fluid_bulk=config.initial_porosity*config.biot_modulus_pa;
        for_each_coupled_cell(current,[&](const auto&,const auto& subzones){
            double reference=0.0,deformed=0.0;
            for(const auto zi:subzones){
                const double now=deformed_tet_signed_volume(current,current.zones[zi]);
                const double before=reference_signed_volume.at(zi);
                if(!std::isfinite(now)||now*before<=0.0)
                    throw std::runtime_error("evolving porosity/storage found an inverted cell");
                const double weight=current.zones[zi].integration_weight;
                reference+=std::abs(before)*weight;
                deformed+=std::abs(now)*weight;
            }
            const double j=deformed/reference;
            const double porosity=1.0-(1.0-config.initial_porosity)/j;
            if(!std::isfinite(porosity)||!(porosity>0.0&&porosity<1.0))
                throw std::runtime_error("evolving porosity/storage left the physical porosity range");
            for(const auto zi:subzones){
                state.zone_porosity[zi]=porosity;
                state.zone_modulus_pa[zi]=fluid_bulk/porosity;
            }
        });
        std::vector<double> nodal_volume(current.nodes.size(),0.0);
        std::vector<double> nodal_porosity_volume(current.nodes.size(),0.0);
        state.nodal_modulus_pa.assign(current.nodes.size(),config.biot_modulus_pa);
        for(std::size_t zi=0;zi<current.zones.size();++zi){
            const auto& zone=current.zones[zi];if(!zone.active)continue;
            const double share=zone.volume_m3*zone.integration_weight/4.0;
            for(const auto ni:zone.nodes){
                nodal_volume[ni]+=share;
                nodal_porosity_volume[ni]+=share*state.zone_porosity[zi];
            }
        }
        for(std::size_t ni=0;ni<current.nodes.size();++ni)
            if(nodal_volume[ni]>0.0)
                state.nodal_modulus_pa[ni]=fluid_bulk*nodal_volume[ni]/
                    nodal_porosity_volume[ni];
        return state;
    };

    auto capture=[&](HydraulicSnapshot& snapshot,const RockModel& current){
        snapshot.nodal_displacement_m.resize(current.nodes.size());
        for(std::size_t i=0;i<current.nodes.size();++i)
            snapshot.nodal_displacement_m[i]=current.nodes[i].displacement;
        snapshot.zone_stress_pa.resize(current.zones.size());
        for(std::size_t i=0;i<current.zones.size();++i)
            snapshot.zone_stress_pa[i]=current.zones[i].stress_pa;
    };
    HydraulicSnapshot start;
    start.head_m=initial.head_m;
    start.nodal_boundary_outward_flux_m3_s.assign(model.nodes.size(),0.0);
    start.nodal_pore_pressure_pa=initial.nodal_pore_pressure_pa;
    start.nodal_saturation=initial.nodal_saturation;
    start.zone_head_m=initial.zone_head_m;
    start.zone_saturation=initial.zone_saturation;
    start.zone_darcy_velocity_m_s=initial.zone_darcy_velocity_m_s;
    const auto initial_hydraulics=derive_zone_hydraulics(model);
    start.zone_conductivity_m_s=initial_hydraulics.conductivity;
    start.zone_joint_aperture_m=initial_hydraulics.aperture;
    start.zone_joint_normal_stress_mpa=initial_hydraulics.normal_stress_mpa;
    capture(start,model);
    if(config.evolving_porosity_storage){
        const auto initial_storage=storage_at(model);
        start.zone_porosity=initial_storage.zone_porosity;
        start.nodal_biot_modulus_pa=initial_storage.nodal_modulus_pa;
    }
    std::vector<HydraulicSnapshot> history{std::move(start)};

    const auto& material=model.config.material;
    const double confined_modulus=material.elastic_modulus_pa*(1.0-material.poisson_ratio)/
        ((1.0+material.poisson_ratio)*(1.0-2.0*material.poisson_ratio));
    if(!(confined_modulus>0.0)||!std::isfinite(confined_modulus))
        throw std::invalid_argument("bidirectional coupling requires positive elastic confined modulus");
    std::vector<double> nodal_fluid_volume(model.nodes.size(),0.0);
    if(config.evolving_porosity_storage)
        for(const auto& zone:model.zones){
            if(!zone.active)continue;
            const double share=zone.volume_m3*zone.integration_weight/4.0;
            for(const auto ni:zone.nodes)nodal_fluid_volume[ni]+=share;
        }
    double time=0.0;
    std::size_t step_index=0;
    try {
    while(time<config.duration_s-std::max(1.0,config.duration_s)*1.0e-14){
        const double next_time=std::min(config.duration_s,
            static_cast<double>(++step_index)*config.time_step_s);
        const double dt=next_time-time;
        const RockModel step_start=model;
        const auto step_storage=storage_at(step_start);
        auto storage_guess=step_storage;
        const double largest_modulus=*std::max_element(step_storage.zone_modulus_pa.begin(),
            step_storage.zone_modulus_pa.end());
        const double stiffness_ratio=config.biot_coefficient*config.biot_coefficient*
            largest_modulus/confined_modulus;
        const double relaxation=1.0/(1.0+stiffness_ratio);
        std::vector<double> strain_guess(model.zones.size(),0.0);
        const std::vector<double> zero_strain(model.zones.size(),0.0);
        const auto flow_predictor=[&](const StorageState& storage){
            RockModel flow_only_model=step_start;
            flow_only_model.config.seepage.duration_s=dt;
            flow_only_model.config.seepage.time_step_s=dt;
            return solve_transient_saturated_impl(
                flow_only_model,&step_start.hydraulic.head_m,&zero_strain,time,
                config.evolving_porosity_storage?&storage.zone_modulus_pa:nullptr,
                config.evolving_porosity_storage?&step_storage.zone_modulus_pa:nullptr);
        };
        auto fluid_only=flow_predictor(storage_guess);
        if(!fluid_only.converged)
            throw std::runtime_error("bidirectional flow-only predictor did not converge");
        bool accepted=false;
        double last_difference=0.0,last_magnitude=0.0,last_mass_error=0.0;
        double last_mass_residual=0.0,last_mass_scale=0.0,last_fluid_mass_error=0.0;
        for(int iteration=0;iteration<config.coupling_max_iterations;++iteration){
            if(iteration>0&&config.evolving_porosity_storage){
                fluid_only=flow_predictor(storage_guess);
                if(!fluid_only.converged)
                    throw std::runtime_error("bidirectional flow-only predictor did not converge");
            }
            RockModel flow_model=step_start;
            flow_model.config.seepage.duration_s=dt;
            flow_model.config.seepage.time_step_s=dt;
            const auto fluid=solve_transient_saturated_impl(
                flow_model,&step_start.hydraulic.head_m,&strain_guess,time,
                config.evolving_porosity_storage?&storage_guess.zone_modulus_pa:nullptr,
                config.evolving_porosity_storage?&step_storage.zone_modulus_pa:nullptr);
            if(!fluid.converged)
                throw std::runtime_error("bidirectional fluid step did not converge");
            RockModel mechanical_model=step_start;
            // FLAC's flow loop applies -alpha*dp_flow to total stress. Its
            // mechanical loop applies +alpha^2*M*depsilon_v; the latter is
            // not generally -alpha*average(dp_mechanical) when a fixed-pore-
            // pressure boundary redistributes generated pore volume.
            for_each_coupled_cell(mechanical_model,[&](const auto& nodes,const auto& subzones){
                if(nodes.empty()||subzones.empty())return;
                double mechanical_pressure=0.0;
                for(const auto node:nodes)
                    mechanical_pressure+=fluid.nodal_pore_pressure_pa[node]-
                        fluid_only.nodal_pore_pressure_pa[node];
                mechanical_pressure/=static_cast<double>(nodes.size());
                double cell_strain=0.0,cell_volume=0.0,cell_modulus=0.0;
                for(const auto zi:subzones){
                    const double volume=mechanical_model.zones[zi].volume_m3*
                        mechanical_model.zones[zi].integration_weight;
                    cell_strain+=strain_guess[zi]*volume;
                    cell_modulus+=storage_guess.zone_modulus_pa[zi]*volume;
                    cell_volume+=volume;
                }
                cell_strain/=cell_volume;
                cell_modulus/=cell_volume;
                const double offset=config.biot_coefficient*mechanical_pressure+
                    config.biot_coefficient*config.biot_coefficient*
                        cell_modulus*cell_strain;
                for(const auto zi:subzones)
                    for(int axis=0;axis<3;++axis)
                        mechanical_model.zones[zi].stress_pa[axis]+=offset;
            });
            // An unchanged hydraulic load on an already accepted equilibrium
            // requires no new mechanical relaxation. Re-solving it injects
            // force-tolerance noise into a physically zero strain increment.
            const bool unchanged=fluid.nodal_pore_pressure_pa==
                step_start.hydraulic.nodal_pore_pressure_pa &&
                std::all_of(strain_guess.begin(),strain_guess.end(),
                    [](double value){return value==0.0;}) &&
                storage_guess.zone_modulus_pa==step_storage.zone_modulus_pa;
            auto result=previous_result;
            if(unchanged)mechanical_model.hydraulic=fluid;
            else result=continue_hydraulic_stage(
                mechanical_model,previous_result,fluid,model.config.solver,progress);
            if(!result.converged)
                throw std::runtime_error("bidirectional mechanical step did not converge");
            double difference=0.0,magnitude=0.0;
            std::vector<double> new_strain(strain_guess.size());
            for(std::size_t zi=0;zi<new_strain.size();++zi){
                for(int axis=0;axis<3;++axis)
                    new_strain[zi]+=mechanical_model.zones[zi].strain[axis]-
                        step_start.zones[zi].strain[axis];
                difference=std::max(difference,std::abs(new_strain[zi]-strain_guess[zi]));
                magnitude=std::max(magnitude,std::abs(new_strain[zi]));
            }
            last_difference=difference;last_magnitude=magnitude;
            double pressure_scale=1.0;
            for(std::size_t node=0;node<fluid.nodal_pore_pressure_pa.size();++node)
                pressure_scale=std::max({pressure_scale,
                    std::abs(fluid.nodal_pore_pressure_pa[node]),
                    std::abs(step_start.hydraulic.nodal_pore_pressure_pa[node]),
                    std::abs(fluid.nodal_pore_pressure_pa[node]-
                        step_start.hydraulic.nodal_pore_pressure_pa[node])});
            // A strain increment can approach zero while the mechanical solve
            // retains finite roundoff/force-tolerance noise. Test the induced
            // Biot pressure error against the pressure state, not only this
            // step's pressure increment, which can itself approach zero.
            const double current_largest_modulus=*std::max_element(
                storage_guess.zone_modulus_pa.begin(),storage_guess.zone_modulus_pa.end());
            const double pressure_mismatch=current_largest_modulus*
                config.biot_coefficient*difference;
            double actual_mechanical_storage=0.0;
            double mechanical_roundoff=0.0;
            for(std::size_t zi=0;zi<mechanical_model.zones.size();++zi){
                const auto& zone=mechanical_model.zones[zi];
                if(!zone.active)continue;
                actual_mechanical_storage+=config.biot_coefficient*
                    tet_data(mechanical_model,zone).volume*zone.integration_weight*
                    new_strain[zi];
                for(int axis=0;axis<3;++axis)
                    mechanical_roundoff+=64.0*std::numeric_limits<double>::epsilon()*
                        config.biot_coefficient*zone.volume_m3*zone.integration_weight*
                        (std::abs(zone.strain[axis])+
                         std::abs(step_start.zones[zi].strain[axis]));
            }
            const auto& fluid_snapshot=fluid.history.back();
            const double boundary_volume=fluid_snapshot.net_outflow_m3_s*dt;
            const auto updated_storage=config.evolving_porosity_storage?
                storage_at(mechanical_model):StorageState{};
            double actual_storage_change=fluid_snapshot.storage_change_m3;
            double mass_roundoff=fluid_snapshot.mass_balance_roundoff_m3+mechanical_roundoff;
            double modulus_mismatch=0.0;
            if(config.evolving_porosity_storage){
                actual_storage_change=0.0;
                for(std::size_t ni=0;ni<model.nodes.size();++ni){
                    actual_storage_change+=nodal_fluid_volume[ni]*(
                        fluid.nodal_pore_pressure_pa[ni]/
                            updated_storage.nodal_modulus_pa[ni]-
                        step_start.hydraulic.nodal_pore_pressure_pa[ni]/
                            step_storage.nodal_modulus_pa[ni]);
                    mass_roundoff+=64.0*std::numeric_limits<double>::epsilon()*
                        nodal_fluid_volume[ni]*(
                            std::abs(fluid.nodal_pore_pressure_pa[ni]/updated_storage.nodal_modulus_pa[ni])+
                            std::abs(step_start.hydraulic.nodal_pore_pressure_pa[ni]/step_storage.nodal_modulus_pa[ni]));
                }
                for(std::size_t zi=0;zi<model.zones.size();++zi)
                    modulus_mismatch=std::max(modulus_mismatch,std::abs(
                        updated_storage.zone_modulus_pa[zi]-
                        storage_guess.zone_modulus_pa[zi])/
                        updated_storage.zone_modulus_pa[zi]);
            }
            const double actual_mass_residual=std::abs(actual_storage_change+
                actual_mechanical_storage+boundary_volume);
            const double actual_mass_scale=std::max({std::abs(actual_storage_change),
                std::abs(actual_mechanical_storage),std::abs(boundary_volume),kSmall});
            const double actual_mass_error=actual_mass_residual/actual_mass_scale;
            last_mass_error=actual_mass_error;
            last_mass_residual=actual_mass_residual;
            last_mass_scale=actual_mass_scale;
            last_fluid_mass_error=fluid_snapshot.relative_mass_balance_error;
            // Near steady state the pressure scale can remain large while the
            // current step's flow and strain increments are tiny. Require the
            // accepted mechanical state to satisfy the same mass tolerance.
            if(pressure_mismatch<=config.coupling_relative_tolerance*pressure_scale &&
               actual_mass_residual<=config.coupling_relative_tolerance*actual_mass_scale+mass_roundoff &&
               modulus_mismatch<=config.coupling_relative_tolerance){
                model=std::move(mechanical_model);
                previous_result=std::move(result);
                auto snapshot=fluid.history.back();
                snapshot.storage_change_m3=actual_storage_change;
                snapshot.mechanical_storage_change_m3=actual_mechanical_storage;
                snapshot.relative_mass_balance_error=actual_mass_error;
                snapshot.mass_balance_residual_m3=actual_mass_residual;
                snapshot.mass_balance_roundoff_m3=mass_roundoff;
                model.hydraulic.mechanical_storage_change_m3=actual_mechanical_storage;
                model.hydraulic.relative_mass_balance_error=snapshot.relative_mass_balance_error;
                model.hydraulic.mass_balance_residual_m3=actual_mass_residual;
                model.hydraulic.mass_balance_roundoff_m3=mass_roundoff;
                capture(snapshot,model);
                if(config.evolving_porosity_storage){
                    snapshot.zone_porosity=updated_storage.zone_porosity;
                    snapshot.nodal_biot_modulus_pa=updated_storage.nodal_modulus_pa;
                }
                history.push_back(std::move(snapshot));
                accepted=true;
                break;
            }
            for(std::size_t zi=0;zi<strain_guess.size();++zi)
                strain_guess[zi]+=relaxation*(new_strain[zi]-strain_guess[zi]);
            if(config.evolving_porosity_storage)storage_guess=updated_storage;
        }
        if(!accepted){
            std::ostringstream message;
            message<<std::setprecision(12)<<"bidirectional Biot iteration did not converge within couplingMaxIterations at t="
                <<time+dt<<" s: strain difference="<<last_difference<<", magnitude="<<last_magnitude
                <<", mass error="<<last_mass_error<<", mass residual="<<last_mass_residual
                <<" m3, mass scale="<<last_mass_scale<<" m3, fluid mass error="
                <<last_fluid_mass_error;
            throw std::runtime_error(message.str());
        }
        time=next_time;
    }
    } catch(const std::runtime_error& error) {
        // Trial models never overwrite the accepted model. Preserve its
        // complete history for native/WASM inspection, but never report the
        // requested run as converged or include the rejected trial field.
        model.hydraulic.history=std::move(history);
        model.hydraulic.converged=false;
        model.hydraulic.failure_reason=error.what();
        previous_result.converged=false;
        previous_result.warnings.push_back("BIDIRECTIONAL_BIOT_INCOMPLETE: "+std::string(error.what()));
        return previous_result;
    }
    model.hydraulic.history=std::move(history);
    previous_result.warnings.push_back(config.evolving_porosity_storage?
        "BIDIRECTIONAL_BIOT_INCOMPRESSIBLE_GRAINS_EVOLVING_STORAGE":
        model.config.material.type==MaterialType::SMRM?
            "BIDIRECTIONAL_BIOT_SMRM_LAGGED_PERMEABILITY_FIXED_STORAGE":
            "BIDIRECTIONAL_BIOT_ELASTIC_EACH_STEP");
    return previous_result;
}

RockResult solve_steady_seepage_mechanics(RockModel& model,const ProgressCallback& progress) {
    const auto hydraulic=solve_steady_seepage(model);
    if(!hydraulic.converged)throw std::runtime_error("steady seepage linear solver did not converge");
    if(!model.config.seepage.solve_mechanics){
        RockResult result;
        result.solver_version=kVersion;result.converged=true;result.total_cycles=hydraulic.iterations;
        result.warnings.push_back("STEADY_SEEPAGE_ONLY");
        return result;
    }
    return solve(model,progress);
}

RockResult solve_transient_seepage_mechanics(RockModel& model,const ProgressCallback& progress) {
    if(model.config.seepage.bidirectional_coupling)
        return solve_coupled_saturated_seepage_mechanics(model,progress);
    auto hydraulic=solve_transient_seepage(model);
    if(!hydraulic.converged)
        throw std::runtime_error("transient seepage linear solver did not converge");
    if(!model.config.seepage.solve_mechanics){
        RockResult result;
        result.solver_version=kVersion;
        result.converged=true;
        result.total_cycles=std::accumulate(hydraulic.history.begin(),hydraulic.history.end(),0,
            [](int sum,const HydraulicSnapshot& snapshot){return sum+snapshot.iterations;});
        result.warnings.push_back("TRANSIENT_SEEPAGE_ONLY");
        return result;
    }

    auto state_from_snapshot=[&](const HydraulicSnapshot& snapshot){
        HydraulicState state;
        state.valid=true;state.converged=true;state.time_s=snapshot.time_s;
        state.iterations=snapshot.iterations;state.relative_residual=snapshot.relative_residual;
        state.inflow_m3_s=snapshot.inflow_m3_s;state.outflow_m3_s=snapshot.outflow_m3_s;
        state.net_outflow_m3_s=snapshot.net_outflow_m3_s;
        state.mechanical_storage_change_m3=snapshot.mechanical_storage_change_m3;
        state.relative_mass_balance_error=snapshot.relative_mass_balance_error;
        state.rainfall_infiltration_m3_s=snapshot.rainfall_infiltration_m3_s;
        state.rainfall_runoff_m3_s=snapshot.rainfall_runoff_m3_s;
        state.head_m=snapshot.head_m;
        state.nodal_boundary_outward_flux_m3_s=snapshot.nodal_boundary_outward_flux_m3_s;
        state.nodal_pore_pressure_pa=snapshot.nodal_pore_pressure_pa;
        state.nodal_saturation=snapshot.nodal_saturation;
        state.zone_head_m=snapshot.zone_head_m;
        state.zone_saturation=snapshot.zone_saturation;
        state.zone_darcy_velocity_m_s=snapshot.zone_darcy_velocity_m_s;
        state.zone_conductivity_m_s=hydraulic.zone_conductivity_m_s;
        state.zone_joint_aperture_m=hydraulic.zone_joint_aperture_m;
        state.zone_joint_normal_stress_mpa=hydraulic.zone_joint_normal_stress_mpa;
        return state;
    };
    const auto capture_mechanics=[&](HydraulicSnapshot& snapshot){
        snapshot.nodal_displacement_m.resize(model.nodes.size());
        for(std::size_t i=0;i<model.nodes.size();++i)
            snapshot.nodal_displacement_m[i]=model.nodes[i].displacement;
        snapshot.zone_stress_pa.resize(model.zones.size());
        for(std::size_t i=0;i<model.zones.size();++i)
            snapshot.zone_stress_pa[i]=model.zones[i].stress_pa;
    };

    auto history=hydraulic.history;
    if(!model.config.seepage.mechanics_each_step){
        model.hydraulic=hydraulic;
        auto result=solve(model,progress);
        capture_mechanics(history.back());
        model.hydraulic.history=std::move(history);
        result.warnings.push_back("TRANSIENT_MECHANICS_FINAL_TIME_ONLY");
        return result;
    }

    model.hydraulic=state_from_snapshot(history.front());
    auto result=solve(model,progress);
    capture_mechanics(history.front());
    for(std::size_t i=1;i<history.size();++i){
        const auto next=state_from_snapshot(history[i]);
        result=continue_hydraulic_stage(model,result,next,model.config.solver,progress);
        capture_mechanics(history[i]);
    }
    model.hydraulic.history=std::move(history);
    result.warnings.push_back("TRANSIENT_ONE_WAY_COUPLING_EACH_STEP");
    return result;
}

} // namespace rocksolver
