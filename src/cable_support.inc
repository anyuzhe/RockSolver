// Included inside rocksolver's anonymous namespace; no independent build target.
std::array<double,4> cable_weights(const RockModel& model, const Zone& z, const Vec3& p) {
    const auto g=tetra_geometry(model,z);
    const auto d=p-model.nodes[z.nodes[0]].position;
    std::array<double,4> w{0,dot(g.gradients[1],d),dot(g.gradients[2],d),dot(g.gradients[3],d)};
    w[0]=1-w[1]-w[2]-w[3]; return w;
}
template<class F> void for_each_cable_host(const CableNode& n, F&& action) {
    if(n.physical_host_count==8) {
        for(std::size_t j=0;j<8;++j)action(n.physical_hosts[j],n.physical_weights[j]);
    } else {
        for(std::size_t j=0;j<4;++j)action(n.hosts[j],n.weights[j]);
    }
}
void initialize_cable_physical_hosts(const RockModel& model, CableNode& n) {
    n.physical_host_count=0;
    if(!model.cable_physical_host_interpolation)return;
    const PhysicalCell* cell=nullptr;
    for(const auto& candidate:model.physical_cells)
        if(std::find(candidate.subzones.begin(),candidate.subzones.end(),n.zone)!=candidate.subzones.end()) {
            cell=&candidate;break;
        }
    if(!cell || cell->nodes.size()!=8)return;
    // Abaqus C3D8 corner order: bottom 0..3, top 4..7.
    constexpr std::array<std::array<double,3>,8> signs{{
        {{-1,-1,-1}},{{1,-1,-1}},{{1,1,-1}},{{-1,1,-1}},
        {{-1,-1,1}},{{1,-1,1}},{{1,1,1}},{{-1,1,1}}
    }};
    std::array<double,3> local{};
    double scale=0;
    for(auto index:cell->nodes)scale=std::max(scale,norm(model.nodes[index].position-n.position));
    for(int iteration=0;iteration<20;++iteration) {
        Vec3 mapped{};
        std::array<std::array<double,3>,3> jacobian{};
        for(std::size_t j=0;j<8;++j) {
            const auto& s=signs[j];const auto& p=model.nodes[cell->nodes[j]].position;
            const double f0=1+s[0]*local[0],f1=1+s[1]*local[1],f2=1+s[2]*local[2];
            const double w=f0*f1*f2/8;
            mapped+=p*w;
            const std::array<double,3> dw{{s[0]*f1*f2/8,s[1]*f0*f2/8,s[2]*f0*f1/8}};
            for(int k=0;k<3;++k){jacobian[0][k]+=p.x*dw[k];jacobian[1][k]+=p.y*dw[k];jacobian[2][k]+=p.z*dw[k];}
        }
        const Vec3 residual=mapped-n.position;
        if(norm(residual)<=1e-12*std::max(scale,1.0))break;
        const auto inv=inverse(jacobian);
        for(int k=0;k<3;++k)local[k]-=inv[k][0]*residual.x+inv[k][1]*residual.y+inv[k][2]*residual.z;
        if(iteration==19)throw std::runtime_error("锚杆 C3D8 宿主单元插值未收敛");
    }
    for(double coordinate:local)if(coordinate < -1-1e-8 || coordinate > 1+1e-8)
        throw std::runtime_error("锚杆节点超出 C3D8 宿主单元");
    for(std::size_t j=0;j<8;++j) {
        const auto& s=signs[j];
        n.physical_hosts[j]=cell->nodes[j];
        n.physical_weights[j]=(1+s[0]*local[0])*(1+s[1]*local[1])*(1+s[2]*local[2])/8;
    }
    n.physical_host_count=8;
}
Vec3 cable_host_displacement(const RockModel& model,const CableNode& n) {
    Vec3 u{};for_each_cable_host(n,[&](std::size_t index,double weight){u+=model.nodes[index].displacement*weight;});return u;
}
double cable_host_pore_pressure(const RockModel& model,const CableNode& n) {
    double pressure=0;
    for_each_cable_host(n,[&](std::size_t index,double weight){pressure+=node_pore_pressure(model,index)*weight;});
    return pressure;
}
void install_cables(RockModel& model,bool after) {
    if(model.config.cables.empty())return;
    if(model.scenario!=Scenario::Slope || model.config.solver.large_strain)
        throw std::runtime_error("锚杆仅支持小应变边坡计算");
    // Cable count is bounded by available memory, not a fixed 200-rod cap.
    std::set<int> ids;
    for(const auto& c:model.config.cables) {
        if(c.id<1 || !ids.insert(c.id).second)throw std::runtime_error("锚杆 ID 无效或重复");
        for(double value:{c.young,c.area,c.yield_tension,c.yield_compression,c.grout_stiffness,c.grout_cohesion,c.grout_friction,c.perimeter,c.density,c.pretension,c.free_length})
            if(!std::isfinite(value)||value<0)throw std::runtime_error("锚杆参数必须为有限非负数");
        if(c.young<=0||c.area<=0||c.grout_stiffness<=0||c.perimeter<=0||c.grout_friction>=90||c.points.size()<2||c.points.size()>201)
            throw std::runtime_error("锚杆参数或分段数无效");
        if(c.pretension>c.yield_tension || (c.pretension>0 && !c.after_initial))throw std::runtime_error("预应力必须在裸坡平衡后施加，且不得超过拉伸屈服力");
        if(c.after_initial!=after)continue;
        CableState state;state.config=c;
        const auto delta=c.points.back()-c.points.front();const double total=norm(delta);
        if(!(total>1e-6))throw std::runtime_error("锚杆长度无效");
        if(c.free_length>=total)throw std::runtime_error("自由段必须短于总长度");
        state.tension_control=c.pretension>0;
        state.direction=delta*(1/total);
        for(std::size_t i=0;i<c.points.size();++i) {
            auto p=c.points[i];
            if(!std::isfinite(p.x)||!std::isfinite(p.y)||!std::isfinite(p.z)||norm(cross(p-c.points[0],state.direction))>1e-7*std::max(1.0,total))
                throw std::runtime_error("第一版仅支持直线锚杆");
            CableNode n;n.position=p;bool found=false;int containing_cells=0;
            for(std::size_t zi=0;zi<model.zones.size();++zi) {
                const auto& z=model.zones[zi];if(!z.active)continue;
                auto w=cable_weights(model,z,p);
                if(*std::min_element(w.begin(),w.end()) < -1e-8)continue;
                double sum=0;for(auto& v:w){v=std::max(v,0.0);sum+=v;}for(auto& v:w)v/=sum;
                // FLAC3D selects the last imported physical cell when a cable
                // node lies on a shared face.
                if(found && z.id==n.host_id)continue;
                ++containing_cells;
                n.zone=zi;n.host_id=z.id;n.weights=w;std::copy_n(z.nodes.begin(),4,n.hosts.begin());found=true;
            }
            if(!found)throw std::runtime_error("锚杆 "+std::to_string(c.id)+" 有节点在坡体外，请调整方向或长度");
            if(containing_cells>1 && !model.physical_cells.empty()) {
                // Resolve C3D8 overlay-subzone ties just inside the selected
                // cell. At the exact face several subzones can contain p but
                // have different stresses; FLAC's host uses the interior side.
                for(const auto& cell:model.physical_cells)if(cell.id==n.host_id) {
                    Vec3 center{};
                    for(auto index:cell.nodes)center+=model.nodes[index].position;
                    center=center*(1.0/cell.nodes.size());
                    const Vec3 probe=p+(center-p)*1e-6;
                    for(auto zi:cell.subzones) {
                        auto probe_weights=cable_weights(model,model.zones[zi],probe);
                        if(*std::min_element(probe_weights.begin(),probe_weights.end()) < -1e-8)continue;
                        auto w=cable_weights(model,model.zones[zi],p);
                        double sum=0;for(auto& v:w){v=std::max(v,0.0);sum+=v;}for(auto& v:w)v/=sum;
                        n.zone=zi;n.weights=w;std::copy_n(model.zones[zi].nodes.begin(),4,n.hosts.begin());
                        break;
                    }
                    break;
                }
            }
            if(!model.physical_cells.empty())for(const auto& cell:model.physical_cells)
                if(std::find(cell.subzones.begin(),cell.subzones.end(),n.zone)!=cell.subzones.end()){n.host_id=cell.id;break;}
            initialize_cable_physical_hosts(model,n);
            n.reference_displacement=cable_host_displacement(model,n);
            state.nodes.push_back(n);
            if(i) {
                const double length=dot(p-c.points[i-1],state.direction);
                if(length<1e-7)throw std::runtime_error("锚杆存在倒序或重复节点");
                state.lengths.push_back(length);state.nodes[i-1].length+=length/2;state.nodes[i].length+=length/2;
            }
        }
        // Exact interval coverage rejects a cable traversing a void, even if endpoints are inside.
        std::vector<std::pair<double,double>> intervals;
        for(const auto& z:model.zones) {
            if(!z.active)continue;
            auto a=cable_weights(model,z,c.points.front()),b=cable_weights(model,z,c.points.back());
            double lo=0,hi=1;
            for(int k=0;k<4;++k){double d=b[k]-a[k];if(std::abs(d)<1e-14){if(a[k]<-1e-9){hi=-1;break;}}else if(d>0)lo=std::max(lo,-a[k]/d);else hi=std::min(hi,-a[k]/d);}
            if(hi>=lo)intervals.emplace_back(lo,hi);
        }
        std::sort(intervals.begin(),intervals.end());double end=0;
        for(auto [lo,hi]:intervals){if(lo>end+1e-8)break;end=std::max(end,hi);}
        if(end<1-1e-8)throw std::runtime_error("锚杆穿出坡体或穿过空洞，未建立支护");
        for(std::size_t i=0;i<state.nodes.size();++i){
            double x=dot(state.nodes[i].position-c.points.front(),state.direction);
            double lo=i?x-state.lengths[i-1]/2:0,hi=i<state.lengths.size()?x+state.lengths[i]/2:total;
            state.nodes[i].bonded_length=std::max(0.0,hi-std::max(lo,c.free_length));
        }
        state.reference_extension.resize(state.lengths.size());
        state.axial_force.resize(state.lengths.size());state.plastic_extension.resize(state.lengths.size());
        state.strain.resize(state.lengths.size());state.yielded.resize(state.lengths.size());
        for(auto& n:state.nodes){n.mass=c.density*c.area*n.length;n.stiffness=c.grout_stiffness*n.bonded_length;}
        for(std::size_t i=0;i<state.lengths.size();++i){double k=c.young*c.area/state.lengths[i];state.nodes[i].stiffness+=k;state.nodes[i+1].stiffness+=k;}
        model.cables.push_back(std::move(state));
    }
}
void cable_dynamic_mass(RockModel& model) {
    for(const auto& c:model.cables)for(const auto& n:c.nodes)for_each_cable_host(n,[&](std::size_t index,double weight) {
        // A conservative row-sum bound includes coupling to the independent cable DOF.
        double k=2*c.config.grout_stiffness*n.bonded_length*weight;
        if(c.config.free_length>0 && &n==&c.nodes.front())k+=2*c.config.young*c.config.area/c.lengths.front()*weight;
        model.nodes[index].dynamic_mass+=Vec3{k,k,k};
    });
}
void update_cable_forces(RockModel& model) {
    for(auto& c:model.cables) {
        auto& p=c.config;const auto d=c.direction;
        if(p.free_length>0)c.nodes.front().axial=dot(cable_host_displacement(model,c.nodes.front())-c.nodes.front().reference_displacement,d);
        for(auto& n:c.nodes) {
            const Vec3 gravity{0,0,-n.mass*model.config.gravity_m_s2};
            n.force=dot(gravity,d);n.reference_force=std::abs(n.force);
            const auto& s=model.zones[n.zone].stress_pa;
            const double axial=s[0]*d.x*d.x+s[1]*d.y*d.y+s[2]*d.z*d.z+2*(s[3]*d.x*d.y+s[4]*d.x*d.z+s[5]*d.y*d.z);
            const double total_confinement=-0.5*(s[0]+s[1]+s[2]-axial);
            n.confinement=std::max(0.0,total_confinement-
                (model.cable_effective_confinement?cable_host_pore_pressure(model,n):0.0));
            const double host=dot(cable_host_displacement(model,n)-n.reference_displacement,d);
            n.slip=host-n.axial;
            const double k=p.grout_stiffness*n.bonded_length;
            const double limit=(p.grout_cohesion+n.confinement*p.perimeter*std::tan(p.grout_friction*std::numbers::pi/180))*n.bonded_length;
            const double trial=k*(n.slip-n.plastic_slip);
            n.grout_force=std::clamp(trial,-limit,limit);n.sliding=std::abs(trial)>limit;
            n.plastic_slip=k>0?n.slip-n.grout_force/k:0;
            n.force+=n.grout_force;n.reference_force+=std::abs(n.grout_force);
            const Vec3 reaction=d*(-n.grout_force)+(gravity-d*dot(gravity,d));
            for_each_cable_host(n,[&](std::size_t index,double weight){auto& host_node=model.nodes[index];host_node.force+=reaction*weight;host_node.force_reference+=norm(reaction)*weight;});
        }
        for(std::size_t i=0;i<c.lengths.size();++i) {
            double extension=c.nodes[i+1].axial-c.nodes[i].axial;
            double k=p.young*p.area/c.lengths[i];
            double trial=k*(extension-c.plastic_extension[i]-c.reference_extension[i]);
            // Match FLAC3D's pretension boundary condition: the maintained
            // tension is applied at the cable head/end element only.  Locking
            // every segment in the free length to the same pretension makes
            // adjacent axial forces cancel exactly and leaves the free cable
            // nodes unable to equilibrate their own body/interaction forces.
            if(c.tension_control && i==0){
                trial=p.pretension;
                c.reference_extension[i]=extension-c.plastic_extension[i]-trial/k;
            }
            double force=std::clamp(trial,-p.yield_compression,p.yield_tension);
            c.axial_force[i]=force;c.yielded[i]=trial>p.yield_tension||trial< -p.yield_compression;
            c.plastic_extension[i]=extension-c.reference_extension[i]-force/k;c.strain[i]=extension/c.lengths[i];
            c.nodes[i].force+=force;c.nodes[i+1].force-=force;
            c.nodes[i].reference_force+=std::abs(force);c.nodes[i+1].reference_force+=std::abs(force);
        }
        if(p.free_length>0){
            auto& head=c.nodes.front();
            for_each_cable_host(head,[&](std::size_t index,double weight){auto& host=model.nodes[index];host.force+=d*(head.force*weight);host.force_reference+=std::abs(head.force)*weight;});
        }
    }
}
// Independent cable equilibrium: do not let the much larger rock-force sum
// hide a support residual. A mean-reference floor avoids tiny unloaded end
// forces dominating the local ratio. Rigidly attached heads have no cable DOF.
struct CableForceBalance {
    double residual = 0, reference = 0, local = 0, scaled_local = 0, maximum_force = 0;
    std::size_t count = 0;
    double average() const { return residual / std::max(reference, 1.0); }
};
CableForceBalance cable_force_balance(const RockModel& model) {
    CableForceBalance balance;
    for (const auto& c : model.cables) for (std::size_t i = c.config.free_length > 0 ? 1 : 0; i < c.nodes.size(); ++i) {
        const auto& n = c.nodes[i];
        const double f = std::abs(n.force);
        balance.residual += f;
        balance.reference += n.reference_force;
        balance.local = std::max(balance.local, f / std::max(n.reference_force, 1.0));
        balance.maximum_force = std::max(balance.maximum_force, f);
        ++balance.count;
    }
    const double floor = std::max(balance.reference / std::max<std::size_t>(balance.count, 1), 1.0);
    for (const auto& c : model.cables) for (std::size_t i = c.config.free_length > 0 ? 1 : 0; i < c.nodes.size(); ++i) {
        const auto& n = c.nodes[i];
        balance.scaled_local = std::max(balance.scaled_local, std::abs(n.force) / std::max(n.reference_force, floor));
    }
    return balance;
}
void advance_cables(RockModel& model,double dt,double damping) {
    for(auto& c:model.cables)for(auto& n:c.nodes) {
        if(c.config.free_length>0 && &n==&c.nodes.front())continue;
        double f=std::abs(n.force);
        double damped=n.force;if(std::abs(n.velocity)>1e-18)damped-=damping*f*std::copysign(1.0,n.velocity);
        n.velocity+=damped/(2*n.stiffness)*dt;n.axial+=n.velocity*dt;
        if(!std::isfinite(n.axial))throw std::runtime_error("锚杆迭代发散，请检查参数与时间步");
    }
}
std::string cables_to_json(const RockModel& model) {
    std::ostringstream out;out<<std::setprecision(16)<<'[';
    for(std::size_t bi=0;bi<model.cables.size();++bi) {
        if(bi)out<<',';const auto& c=model.cables[bi];out<<"{\"id\":"<<c.config.id<<",\"targetPretensionN\":"<<c.config.pretension<<",\"tensionControlActive\":"<<(c.tension_control?"true":"false")<<",\"freeLengthM\":"<<c.config.free_length<<",\"nodes\":[";
        for(std::size_t i=0;i<c.nodes.size();++i){if(i)out<<',';const auto& n=c.nodes[i];
            auto host_delta=cable_host_displacement(model,n)-n.reference_displacement;
            auto displacement=n.reference_displacement+host_delta+c.direction*(n.axial-dot(host_delta,c.direction));
            out<<"{\"index\":"<<i<<",\"position\":["<<n.position.x<<','<<n.position.y<<','<<n.position.z<<"],\"displacementM\":["<<displacement.x<<','<<displacement.y<<','<<displacement.z<<"],\"referenceDisplacementM\":["<<n.reference_displacement.x<<','<<n.reference_displacement.y<<','<<n.reference_displacement.z<<"],\"tributaryLengthM\":"<<n.length<<",\"hostZoneId\":"<<n.host_id<<",\"slipM\":"<<n.slip<<",\"plasticSlipM\":"<<n.plastic_slip<<",\"groutForceN\":"<<n.grout_force<<",\"confinementPa\":"<<n.confinement<<",\"sliding\":"<<(n.sliding?"true":"false")<<'}';}
        out<<"],\"segments\":[";
        for(std::size_t i=0;i<c.lengths.size();++i){if(i)out<<',';out<<"{\"index\":"<<i<<",\"lengthM\":"<<c.lengths[i]<<",\"axialForceN\":"<<c.axial_force[i]<<",\"strain\":"<<c.strain[i]<<",\"plasticExtensionM\":"<<c.plastic_extension[i]<<",\"yielded\":"<<(c.yielded[i]?"true":"false")<<'}';}
        out<<"]}";
    }
    out<<']';return out.str();
}
