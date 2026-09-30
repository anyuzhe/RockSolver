#include "rocksolver/rocksolver.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numbers>

namespace rocksolver {

std::array<double, 3> joint_direction_cosines(double dip_direction_deg,
                                              double dip_deg) {
    const double alpha = dip_direction_deg * std::numbers::pi / 180.0;
    const double beta = dip_deg * std::numbers::pi / 180.0;
    return {-std::cos(alpha) * std::sin(beta),
            std::sin(alpha) * std::sin(beta), std::cos(beta)};
}

namespace {

constexpr double kPi = 3.1415926535897932384626433832795;
constexpr double kEps = 1.0e-12;
constexpr double kFlacStressEpsMpa = 1.0e-18; // 601 DLL uses 1e-12 Pa.
constexpr double kInf = std::numeric_limits<double>::infinity();

struct JointWork {
    double nx{};
    double ny{};
    double nz{};
    double density{};
    double radius{};
    double average_radius{};
    double cohesion_mpa{};
    double friction_deg{};
    double k{};
    double h{};
    double normal_mpa{};
    double shear_mpa{};
    double shear_strength_mpa{};
    double n1{};
    double n2{};
    double n3{};
};

struct StrengthWork {
    double friction_deg{};
    double cohesion_mpa{};
    double tension_mpa{};
    double ucs_mpa{};
    double kic_mpa_sqrt_m{};
    double tan2_theta{};
};

struct PrincipalFrame {
    std::array<double, 3> compression{};
    std::array<std::array<double, 3>, 3> axes{};  // column per principal value
};

double square(double value) { return value * value; }

double clamp01(double value) { return std::clamp(value, 0.0, 1.0); }

double radians(double degrees) { return degrees * kPi / 180.0; }

std::array<double, 6> to_matrix_order(const std::array<double, 6>& value) {
    return {value[0], value[1], value[2], value[3], value[5], value[4]};
}

Matrix6 to_rocksolver_order(const Matrix6& value) {
    constexpr std::array<int, 6> p{0, 1, 2, 3, 5, 4};
    Matrix6 result{};
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            result[row][column] = value[p[row]][p[column]];
    return result;
}

StrengthWork strength_parameters(const MaterialParameters& material) {
    const double reduction = std::max(material.strength_reduction_factor, 1.0);
    StrengthWork value;
    value.friction_deg = std::atan(std::tan(radians(material.friction_deg)) /
                                   reduction) * 180.0 / kPi;
    value.cohesion_mpa = material.cohesion_pa * 1.0e-6 / reduction;
    value.tension_mpa = material.tensile_strength_pa * 1.0e-6 / reduction;
    const double tan_theta =
        std::tan(radians(45.0 + value.friction_deg * 0.5));
    value.tan2_theta = std::max(square(tan_theta), kEps);
    // ModelSMRM derives UCS and K1c from c/phi; the separate UI UCS is only
    // retained as an import/display field.
    value.ucs_mpa = std::max(0.0, 2.0 * value.cohesion_mpa * tan_theta);
    value.kic_mpa_sqrt_m = material.fracture_toughness_mpa_sqrt_m > 0.0
                                ? material.fracture_toughness_mpa_sqrt_m /
                                      reduction
                                : 0.012 * value.ucs_mpa;
    return value;
}

std::array<JointWork, 3> make_joints(const MaterialParameters& material) {
    std::array<JointWork, 3> result{};
    const double reduction = std::max(material.strength_reduction_factor, 1.0);
    const double volume = std::max(material.representative_volume_m3, kEps);
    const std::size_t count = std::min<std::size_t>(3, material.joints.size());
    for (std::size_t index = 0; index < count; ++index) {
        const JointSet& source = material.joints[index];
        JointWork& joint = result[index];
        // Coordinate convention compiled into the supplied
        // cmodelsmrm007.dll (2026-07-08). Direct disassembly of its
        // BasicParameters constructor gives this exact direction cosine
        // order. A later source-tree revision rotated the horizontal
        // projection by 90 degrees and is not compatible with that DLL.
        const auto normal = joint_direction_cosines(
            source.dip_direction_deg, source.dip_deg);
        joint.nx = normal[0];
        joint.ny = normal[1];
        joint.nz = normal[2];
        joint.density = std::max(0.0, source.normal_density_1_m);
        joint.radius = std::max(0.0, source.mean_radius_m);
        joint.average_radius = joint.radius;
        if (joint.radius > kEps) {
            const double lambda_v =
                joint.density / (2.0 * kPi * square(joint.radius));
            if (lambda_v > std::numbers::e)
                joint.average_radius = joint.radius * std::log(lambda_v * volume);
        }
        joint.cohesion_mpa = source.cohesion_pa * 1.0e-6 / reduction;
        joint.friction_deg =
            std::atan(std::tan(radians(source.friction_deg)) / reduction) *
            180.0 / kPi;
    }
    return result;
}

int active_joint_count(const std::array<JointWork, 3>& joints) {
    return static_cast<int>(std::count_if(
        joints.begin(), joints.end(), [](const JointWork& joint) {
            return joint.density > kEps && joint.average_radius > kEps;
        }));
}

void update_joint_state(JointWork& joint,
                        const std::array<double, 6>& stress_pa,
                        const StrengthWork& rock,
                        bool flac3d7_smrm_contact) {
    if (joint.density <= kEps || joint.average_radius <= kEps) return;
    // Convert FLAC/RockSolver tension-positive stress to the
    // compression-positive tensor used by the SMRM book equations.
    const double sxx = -stress_pa[0] * 1.0e-6;
    const double syy = -stress_pa[1] * 1.0e-6;
    const double szz = -stress_pa[2] * 1.0e-6;
    const double sxy = -stress_pa[3] * 1.0e-6;
    const double sxz = -stress_pa[4] * 1.0e-6;
    const double syz = -stress_pa[5] * 1.0e-6;
    const double tx = sxx * joint.nx + sxy * joint.ny + sxz * joint.nz;
    const double ty = sxy * joint.nx + syy * joint.ny + syz * joint.nz;
    const double tz = sxz * joint.nx + syz * joint.ny + szz * joint.nz;
    joint.normal_mpa = joint.nx * tx + joint.ny * ty + joint.nz * tz;
    // Project the traction onto the plane before taking its norm.  The
    // squared-magnitude difference creates spurious shear under isotropic
    // pore pressure through catastrophic cancellation.
    const double shear_x = tx - joint.normal_mpa * joint.nx;
    const double shear_y = ty - joint.normal_mpa * joint.ny;
    const double shear_z = tz - joint.normal_mpa * joint.nz;
    const double shear = std::sqrt(square(shear_x) + square(shear_y) +
                                   square(shear_z));
    const double scale = std::max({std::abs(sxx), std::abs(syy), std::abs(szz),
                                   std::abs(sxy), std::abs(sxz), std::abs(syz)});
    const double roundoff = 64.0 * std::numeric_limits<double>::epsilon() * scale;
    joint.shear_mpa = shear <= roundoff ? 0.0 : shear;

    const double compression = std::max(joint.normal_mpa, 0.0);
    const double raw_k = clamp01(std::exp(
        -10.0 * compression / std::max(rock.ucs_mpa, kEps)));
    if (flac3d7_smrm_contact) {
        // Current 601 source DLL smoothly suppresses contact at very low
        // confinement. The older exported Huzhou DLL uses raw_k directly.
        const double x = clamp01(compression /
            (0.005 * std::max(rock.ucs_mpa, kEps)));
        const double transition = x * x * (3.0 - 2.0 * x);
        joint.k = clamp01(1.0 - transition * (1.0 - raw_k));
    } else {
        joint.k = compression > kEps ? raw_k : 1.0;
    }
    const double closed_strength =
        joint.cohesion_mpa + compression * std::tan(radians(joint.friction_deg));
    // The current 601 DLL reports contact-weighted strength and uses that
    // same value for h evolution.
    const double effective_shear_strength =
        (1.0 - joint.k) * std::max(0.0, closed_strength);
    joint.shear_strength_mpa = effective_shear_strength;
    const double shear_eps_mpa =
        flac3d7_smrm_contact ? kFlacStressEpsMpa : kEps;
    if (joint.shear_mpa <= shear_eps_mpa) {
        joint.h = 0.0;
    } else if (effective_shear_strength <= shear_eps_mpa) {
        joint.h = 1.0;
    } else if (joint.shear_mpa <= effective_shear_strength) {
        joint.h = 0.0;
    } else {
        joint.h = clamp01(1.0 - effective_shear_strength / joint.shear_mpa);
    }
}

Matrix6 compliance_matrix(const MaterialParameters& material,
                          const std::array<JointWork, 3>& joints) {
    Matrix6 compliance{};
    const double e = material.elastic_modulus_pa;
    const double nu = material.poisson_ratio;
    const double normal = 1.0 / e;
    const double poisson = -nu / e;
    const double shear = 2.0 * (1.0 + nu) / e;
    for (int row = 0; row < 3; ++row)
        for (int column = 0; column < 3; ++column)
            compliance[row][column] = row == column ? normal : poisson;
    compliance[3][3] = compliance[4][4] = compliance[5][5] = shear;

    const double alpha = 8.0 * (1.0 - nu * nu) / kPi;
    const double beta = 2.0 / (2.0 - nu);
    for (const JointWork& joint : joints) {
        if (joint.density <= kEps || joint.average_radius <= kEps) continue;
        const double tmp1 = beta * square(joint.h);
        const double tmp2 = square(joint.k) - tmp1;
        const double coeff =
            alpha * joint.density * joint.average_radius / e;
        const double nx = joint.nx, ny = joint.ny, nz = joint.nz;
        const double nx2 = nx * nx, ny2 = ny * ny, nz2 = nz * nz;
        const double nxny = nx * ny, nxnz = nx * nz, nynz = ny * nz;
        compliance[0][0] += coeff * (tmp2 * nx2 + tmp1) * nx2;
        compliance[0][1] += coeff * tmp2 * nxny * nxny;
        compliance[0][2] += coeff * tmp2 * nxnz * nxnz;
        compliance[0][3] += coeff * (2.0 * tmp2 * nx2 + tmp1) * nxny;
        compliance[0][4] += coeff * 2.0 * tmp2 * nxny * nxnz;
        compliance[0][5] += coeff * (2.0 * tmp2 * nx2 + tmp1) * nxnz;
        compliance[1][1] += coeff * (tmp2 * ny2 + tmp1) * ny2;
        compliance[1][2] += coeff * tmp2 * nynz * nynz;
        compliance[1][3] += coeff * (2.0 * tmp2 * ny2 + tmp1) * nxny;
        compliance[1][4] += coeff * (2.0 * tmp2 * ny2 + tmp1) * nynz;
        compliance[1][5] += coeff * 2.0 * tmp2 * ny2 * nxnz;
        compliance[2][2] += coeff * (tmp2 * nz2 + tmp1) * nz2;
        compliance[2][3] += coeff * 2.0 * tmp2 * nz2 * nxny;
        compliance[2][4] += coeff * (2.0 * tmp2 * nz2 + tmp1) * nynz;
        compliance[2][5] += coeff * (2.0 * tmp2 * nz2 + tmp1) * nxnz;
        compliance[3][3] +=
            coeff * (4.0 * tmp2 * nx2 * ny2 + tmp1 * (nx2 + ny2));
        compliance[3][4] += coeff * (4.0 * tmp2 * ny2 + tmp1) * nxnz;
        compliance[3][5] += coeff * (4.0 * tmp2 * nx2 + tmp1) * nynz;
        compliance[4][4] +=
            coeff * (4.0 * tmp2 * ny2 * nz2 + tmp1 * (ny2 + nz2));
        compliance[4][5] += coeff * (4.0 * tmp2 * nz2 + tmp1) * nxny;
        compliance[5][5] +=
            coeff * (4.0 * tmp2 * nx2 * nz2 + tmp1 * (nx2 + nz2));
    }
    for (int row = 1; row < 6; ++row)
        for (int column = 0; column < row; ++column)
            compliance[row][column] = compliance[column][row];
    return compliance;
}

bool invert6(Matrix6& matrix) {
    double augmented[6][12]{};
    for (int row = 0; row < 6; ++row) {
        for (int column = 0; column < 6; ++column)
            augmented[row][column] = matrix[row][column];
        augmented[row][6 + row] = 1.0;
    }
    for (int column = 0; column < 6; ++column) {
        int pivot_row = column;
        double pivot_size = std::abs(augmented[column][column]);
        for (int row = column + 1; row < 6; ++row) {
            if (std::abs(augmented[row][column]) > pivot_size) {
                pivot_size = std::abs(augmented[row][column]);
                pivot_row = row;
            }
        }
        if (!std::isfinite(pivot_size) || pivot_size < 1.0e-30) return false;
        if (pivot_row != column)
            for (int entry = 0; entry < 12; ++entry)
                std::swap(augmented[pivot_row][entry],
                          augmented[column][entry]);
        const double pivot = augmented[column][column];
        for (double& entry : augmented[column]) entry /= pivot;
        for (int row = 0; row < 6; ++row) {
            if (row == column) continue;
            const double factor = augmented[row][column];
            for (int entry = 0; entry < 12; ++entry)
                augmented[row][entry] -= factor * augmented[column][entry];
        }
    }
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            matrix[row][column] = augmented[row][6 + column];
    return true;
}

PrincipalFrame principal_frame(const std::array<double, 6>& stress_pa) {
    // Matrix order is xx,yy,zz,xy,yz,xz. Work in MPa, tension positive.
    const auto stress = to_matrix_order(stress_pa);
    double a[3][3]{{stress[0] * 1.0e-6, stress[3] * 1.0e-6,
                    stress[5] * 1.0e-6},
                   {stress[3] * 1.0e-6, stress[1] * 1.0e-6,
                    stress[4] * 1.0e-6},
                   {stress[5] * 1.0e-6, stress[4] * 1.0e-6,
                    stress[2] * 1.0e-6}};
    double vectors[3][3]{{1.0, 0.0, 0.0},
                         {0.0, 1.0, 0.0},
                         {0.0, 0.0, 1.0}};
    for (int sweep = 0; sweep < 18; ++sweep) {
        int p = 0, q = 1;
        double largest = std::abs(a[0][1]);
        for (const auto pair : {std::pair{0, 2}, std::pair{1, 2}}) {
            if (std::abs(a[pair.first][pair.second]) > largest) {
                largest = std::abs(a[pair.first][pair.second]);
                p = pair.first;
                q = pair.second;
            }
        }
        if (largest < 1.0e-14) break;
        const double angle =
            0.5 * std::atan2(2.0 * a[p][q], a[p][p] - a[q][q]);
        const double c = std::cos(angle), s = std::sin(angle);
        const double app = a[p][p] * c * c + a[q][q] * s * s +
                           2.0 * a[p][q] * s * c;
        const double aqq = a[p][p] * s * s + a[q][q] * c * c -
                           2.0 * a[p][q] * s * c;
        a[p][p] = app;
        a[q][q] = aqq;
        a[p][q] = a[q][p] = 0.0;
        for (int row = 0; row < 3; ++row) {
            if (row != p && row != q) {
                const double arp = a[row][p] * c + a[row][q] * s;
                const double arq = -a[row][p] * s + a[row][q] * c;
                a[row][p] = a[p][row] = arp;
                a[row][q] = a[q][row] = arq;
            }
            const double vrp = vectors[row][p] * c + vectors[row][q] * s;
            const double vrq = -vectors[row][p] * s + vectors[row][q] * c;
            vectors[row][p] = vrp;
            vectors[row][q] = vrq;
        }
    }
    std::array<int, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return a[lhs][lhs] < a[rhs][rhs];
    });
    PrincipalFrame frame;
    for (int principal = 0; principal < 3; ++principal) {
        frame.compression[principal] = -a[order[principal]][order[principal]];
        for (int component = 0; component < 3; ++component)
            frame.axes[principal][component] = vectors[component][order[principal]];
    }
    return frame;
}

bool scaled_quadratic(double A, double B, double C, double& a, double& b,
                      double& c, double& root) {
    if (!std::isfinite(A) || !std::isfinite(B) || !std::isfinite(C)) return false;
    const double scale = std::max({std::abs(A), std::abs(B), std::abs(C)});
    if (scale <= 0.0) return false;
    a = A / scale;
    b = B / scale;
    c = C / scale;
    const double discriminant = std::fma(-4.0 * a, c, b * b);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                             std::max({1.0, std::abs(b * b),
                                       std::abs(4.0 * a * c)});
    if (!std::isfinite(discriminant) || discriminant < -tolerance) return false;
    root = std::sqrt(std::max(0.0, discriminant));
    return true;
}

double quadratic_upper(double A, double B, double C) {
    if (std::abs(A) <= kEps) return kInf;
    double a, b, c, root;
    if (!scaled_quadratic(A, B, C, a, b, c, root)) return kInf;
    if (b >= 0.0 && b + root > 0.0) return -2.0 * c / (b + root);
    return (-b + root) / (2.0 * a);
}

double quadratic_negative_linear_minus(double A, double B, double C) {
    if (std::abs(A) <= kEps) return kInf;
    double a, b, c, root;
    if (!scaled_quadratic(A, B, C, a, b, c, root)) return kInf;
    if (b >= 0.0 && b + root > 0.0) return 2.0 * c / (b + root);
    return (b - root) / (2.0 * a);
}

double sigma1_limit(const JointWork& j, double s2, double s3, double beta,
                    double kic) {
    const double n1s = square(j.n1), n2s = square(j.n2), n3s = square(j.n3);
    const double t1 = beta * square(j.h), t2 = square(j.k) - t1;
    const double lateral = s2 * n2s + s3 * n3s;
    const double value = quadratic_upper(
        n1s * (t2 * n1s + t1), 2.0 * n1s * t2 * lateral,
        t2 * lateral * lateral +
            t1 * (square(s2 * j.n2) + square(s3 * j.n3)) -
            kPi * square(kic) / (4.0 * j.average_radius));
    return std::isfinite(value) && value > kEps ? value : kInf;
}

double sigma3_limit(const JointWork& j, double s1, double s2, double beta,
                    double kic) {
    const double n1s = square(j.n1), n2s = square(j.n2), n3s = square(j.n3);
    const double t1 = beta * square(j.h), t2 = square(j.k) - t1;
    const double lateral = s1 * n1s + s2 * n2s;
    const double value = quadratic_negative_linear_minus(
        n3s * (t2 * n3s + t1), 2.0 * n3s * t2 * lateral,
        t2 * lateral * lateral +
            t1 * (square(s1 * j.n1) + square(s2 * j.n2)) -
            kPi * square(kic) / (4.0 * j.average_radius));
    return std::isfinite(value) && value > kEps && value <= s1 + kEps
               ? value
               : kInf;
}

double quasi_sigma3_limit(const JointWork& j, double s1, double beta,
                          double kic) {
    const double n1 = std::abs(j.n1), n3 = std::abs(j.n3);
    const double friction = std::tan(radians(j.friction_deg));
    const double denominator = n1 * n3 + friction * (1.0 - n1 * n1);
    if (std::abs(denominator) <= kEps || beta <= kEps) return kInf;
    const double strength =
        0.5 * kic * std::sqrt(kPi / (beta * j.average_radius)) +
        j.cohesion_mpa;
    const double value = (n1 * n3 - friction * n1 * n1) / denominator * s1 -
                         strength / denominator;
    return std::isfinite(value) && value <= s1 + kEps ? value : kInf;
}

double uniaxial_limit(const JointWork& j, double beta, double kic) {
    const double n1s = square(j.n1);
    const double radical = j.average_radius *
                           ((square(j.k) - beta * square(j.h)) * n1s +
                            beta * square(j.h));
    if (radical <= 0.0 || std::abs(j.n1) <= kEps) return kInf;
    const double denominator = 2.0 * std::abs(j.n1) * std::sqrt(radical);
    return denominator > kEps ? std::sqrt(kPi) * kic / denominator : kInf;
}

double tau31_limit(const JointWork& j, double s11, double s22, double s33,
                   double beta, double kic) {
    const double n1s = square(j.nx), n2s = square(j.ny), n3s = square(j.nz);
    const double sigma0 = s11 * n1s + s22 * n2s + s33 * n3s;
    const double t1 = beta * square(j.h), t2 = square(j.k) - t1;
    const double A = square(2.0 * j.nx * j.nz) * t2 + t1 * (n1s + n3s);
    const double B = 2.0 * j.nx * j.nz *
                     (2.0 * t2 * sigma0 + t1 * (s11 + s33));
    const double C = t2 * sigma0 * sigma0 +
                     t1 * (square(s11 * j.nx) + square(s22 * j.ny) +
                           square(s33 * j.nz)) -
                     kPi * square(kic) / (4.0 * j.average_radius);
    return quadratic_upper(A, B, C);
}

double tau31_slope(const JointWork& j, double s11, double s22, double s33,
                   double beta, double kic) {
    const double scale = std::max({1.0, std::abs(s11), std::abs(s22),
                                   std::abs(s33)});
    const double step = std::max(1.0, 1.0e-5 * scale);
    const double plus = tau31_limit(j, s11, s22, s33 + step, beta, kic);
    const double minus = tau31_limit(j, s11, s22, s33 - step, beta, kic);
    if (std::isfinite(plus) && std::isfinite(minus))
        return (plus - minus) / (2.0 * step);
    const double base = tau31_limit(j, s11, s22, s33, beta, kic);
    if (std::isfinite(plus) && std::isfinite(base)) return (plus - base) / step;
    if (std::isfinite(minus) && std::isfinite(base)) return (base - minus) / step;
    return kInf;
}

double joint_failure_probability(const JointWork& j, double beta, double kic,
                                 double volume_m3) {
    if (j.density <= kEps || j.radius <= kEps || volume_m3 <= kEps) return 0.0;
    const double drive = square(j.k) * square(j.normal_mpa) +
                         beta * square(j.h) * square(j.shear_mpa);
    if (drive <= kEps) return 0.0;
    const double critical_radius = kPi * square(kic) / (4.0 * drive);
    const double lambda_v = j.density / (2.0 * kPi * square(j.radius));
    const double survival = std::pow(
        std::max(0.0, 1.0 - std::exp(-critical_radius / j.radius)),
        lambda_v * volume_m3);
    return clamp01(1.0 - survival);
}

void compute_diagnostics(const MaterialParameters& material, double volume_m3,
                         const std::array<double, 6>& stress_pa,
                         std::array<JointWork, 3>& joints,
                         SMRMHistory& history, SMRMDiagnostics& out) {
    out = {};
    const StrengthWork rock = strength_parameters(material);
    const PrincipalFrame frame = principal_frame(stress_pa);
    const double s1 = frame.compression[0];
    const double s2 = frame.compression[1];
    const double s3 = frame.compression[2];
    const double alpha = 8.0 * (1.0 - square(material.poisson_ratio)) / kPi;
    const double beta = 2.0 / (2.0 - material.poisson_ratio);
    const Matrix6 compliance = compliance_matrix(material, joints);
    const int active = active_joint_count(joints);
    if (active > 0) {
        out.rock_mass_em_x_gpa = 1.0 / compliance[0][0] * 1.0e-9;
        out.rock_mass_em_y_gpa = 1.0 / compliance[1][1] * 1.0e-9;
        out.rock_mass_em_z_gpa = 1.0 / compliance[2][2] * 1.0e-9;
    } else {
        out.rock_mass_em_x_gpa = out.rock_mass_em_y_gpa =
            out.rock_mass_em_z_gpa = material.elastic_modulus_pa * 1.0e-9;
    }

    double em_sum = 0.0, s1111_sum = 0.0, s1133_sum = 0.0;
    double s1_joint = kInf, ucs_joint = kInf, s3_joint = -kInf;
    double s3_quasi = -kInf, tau_limit = kInf, tau_zero = kInf;
    double tau_slope = kInf, survival = 1.0;
    bool joint_shear = false, joint_tension = false;
    for (std::size_t index = 0; index < joints.size(); ++index) {
        JointWork& joint = joints[index];
        if (joint.density <= kEps || joint.average_radius <= kEps) continue;
        joint.n1 = joint.nx * frame.axes[0][0] +
                   joint.ny * frame.axes[0][1] +
                   joint.nz * frame.axes[0][2];
        joint.n2 = joint.nx * frame.axes[1][0] +
                   joint.ny * frame.axes[1][1] +
                   joint.nz * frame.axes[1][2];
        joint.n3 = joint.nx * frame.axes[2][0] +
                   joint.ny * frame.axes[2][1] +
                   joint.nz * frame.axes[2][2];
        out.joint_k[index] = joint.k;
        out.joint_h[index] = joint.h;
        out.joint_n1[index] = joint.n1;
        out.joint_n2[index] = joint.n2;
        out.joint_n3[index] = joint.n3;
        out.joint_normal_stress_mpa[index] = joint.normal_mpa;
        out.joint_shear_stress_mpa[index] = joint.shear_mpa;
        out.joint_shear_strength_mpa[index] = joint.shear_strength_mpa;

        const double n1s = square(joint.n1), n3s = square(joint.n3);
        const double t1 = beta * square(joint.h), t2 = square(joint.k) - t1;
        const double density_radius = joint.density * joint.average_radius;
        em_sum += density_radius * (t2 * n1s * n1s + t1 * n1s);
        s1111_sum += density_radius * (t2 * n1s * n1s + t1 * n1s);
        s1133_sum += density_radius * t2 * n1s * n3s;

        const bool locked = joint.normal_mpa >= 0.0 &&
                            (1.0 - joint.k) >= 0.95 && joint.h <= kEps;
        const double phi = radians(joint.friction_deg);
        const bool use_rock = locked || std::abs(joint.n1) > std::cos(phi) ||
                              std::abs(joint.n3) < std::sin(phi);
        out.joint_use_rock_sigma3[index] = use_rock ? 1.0 : 0.0;
        if (!locked) {
            s1_joint = std::min(s1_joint,
                                sigma1_limit(joint, s2, s3, beta,
                                             rock.kic_mpa_sqrt_m));
            ucs_joint = std::min(ucs_joint,
                                 uniaxial_limit(joint, beta,
                                                rock.kic_mpa_sqrt_m));
            const double t31 = tau31_limit(
                joint, -stress_pa[0] * 1.0e-6, -stress_pa[1] * 1.0e-6,
                -stress_pa[2] * 1.0e-6, beta, rock.kic_mpa_sqrt_m);
            const double t31_zero = tau31_limit(
                joint, -stress_pa[0] * 1.0e-6, -stress_pa[1] * 1.0e-6, 0.0,
                beta, rock.kic_mpa_sqrt_m);
            const double t31_slope = tau31_slope(
                joint, -stress_pa[0] * 1.0e-6, -stress_pa[1] * 1.0e-6,
                -stress_pa[2] * 1.0e-6, beta, rock.kic_mpa_sqrt_m);
            tau_limit = std::min(tau_limit, t31);
            tau_zero = std::min(tau_zero, t31_zero);
            tau_slope = std::min(tau_slope, t31_slope);
        }
        if (!use_rock) {
            const double limit = sigma3_limit(
                joint, s1, s2, beta, rock.kic_mpa_sqrt_m);
            if (std::isfinite(limit)) {
                out.joint_sigma3_limit_mpa[index] = limit;
                if (!std::isfinite(s3_joint) || limit >= s3_joint) {
                    s3_joint = limit;
                    out.diag_sigma3_joint_limit_mpa = limit;
                    out.diag_sigma3_controlling_joint = index + 1.0;
                }
            }
            const double quasi = quasi_sigma3_limit(
                joint, s1, beta, rock.kic_mpa_sqrt_m);
            if (std::isfinite(quasi)) {
                out.joint_quasi_sigma3_limit_mpa[index] = quasi;
                if (!std::isfinite(s3_quasi) || quasi >= s3_quasi) {
                    s3_quasi = quasi;
                    out.diag_sigma3_quasi_joint_limit_mpa = quasi;
                    out.diag_sigma3_quasi_controlling_joint = index + 1.0;
                }
            }
        }
        survival *= 1.0 - joint_failure_probability(
                              joint, beta, rock.kic_mpa_sqrt_m,
                              std::max(volume_m3, kEps));
        joint_tension = joint_tension || joint.normal_mpa < 0.0;
        joint_shear = joint_shear ||
                      (joint.normal_mpa >= 0.0 && joint.h > kEps);
    }

    const double denominator = 1.0 + alpha * em_sum;
    out.rock_mass_em_gpa = denominator > kEps
                               ? material.elastic_modulus_pa / denominator *
                                     1.0e-9
                               : 0.0;
    const double s1111 = 1.0 + alpha * s1111_sum;
    const double s1133 = -material.poisson_ratio + alpha * s1133_sum;
    out.rock_mass_poisson_ratio =
        std::abs(s1111) > kEps ? std::abs(s1133 / s1111) : 0.0;
    out.triaxial_compressive_strength_mpa =
        active == 0
            ? rock.ucs_mpa
            : std::max(0.0,
                       std::min(rock.ucs_mpa + s3 * rock.tan2_theta, s1_joint));
    out.rock_mass_ucs_mpa =
        std::max(0.0, std::min(rock.ucs_mpa, ucs_joint));
    out.sigma1cmin_mpa = out.rock_mass_ucs_mpa;
    const double rock_s3 = (s1 - rock.ucs_mpa) / rock.tan2_theta;
    out.limit_confining_pressure_mpa = rock_s3;
    if (std::isfinite(s3_joint))
        out.limit_confining_pressure_mpa =
            std::max(out.limit_confining_pressure_mpa, s3_joint);
    out.book_quasi_limit_confining_pressure_mpa = rock_s3;
    if (std::isfinite(s3_quasi))
        out.book_quasi_limit_confining_pressure_mpa =
            std::max(out.book_quasi_limit_confining_pressure_mpa, s3_quasi);
    const double rock_tan_phi = std::tan(radians(rock.friction_deg));
    out.rock_mass_friction_deg =
        std::atan(std::max(0.0, std::min(rock_tan_phi, tau_slope))) *
        180.0 / kPi;
    out.rock_mass_cohesion_mpa =
        std::max(0.0, std::min(rock.cohesion_mpa, tau_zero));
    (void)tau_limit;

    const double span = s1 - s3;
    out.stability_coefficient =
        std::abs(span) > kEps
            ? std::max(0.0,
                       (s1 - out.limit_confining_pressure_mpa) / span)
            : 0.0;
    out.safety_rest = std::max(0.0, out.stability_coefficient - 1.0);
    // The reference DLL stores required compressive confinement with FLAC's
    // negative stress sign.  Its signed requirement is therefore current
    // compression-positive sigma3 plus that signed limit (not sigma3-limit).
    out.confining_pressure_margin_mpa =
        s3 + out.limit_confining_pressure_mpa;
    // In the reference DLL the two *Reinforce_Requirement properties are
    // signed margins (negative: pressure required, positive: reserve).  The
    // later source revision added separate margin fields and clamped these to
    // zero; those fields are absent from the DLL and its exported values prove
    // the original signed convention.
    out.reinforce_requirement_mpa = out.confining_pressure_margin_mpa;
    out.book_quasi_confining_pressure_margin_mpa =
        s3 + out.book_quasi_limit_confining_pressure_mpa;
    out.book_quasi_reinforce_requirement_mpa =
        out.book_quasi_confining_pressure_margin_mpa;

    const double triaxial = rock.ucs_mpa + s3 * rock.tan2_theta;
    double rock_failure = 0.0;
    if (triaxial > kEps) {
        const double utilization = std::max(0.0, s1 / triaxial);
        rock_failure =
            1.0 - std::exp(-std::min(std::pow(utilization,
                                               std::max(material.weibull_shape,
                                                        kEps)),
                                      700.0));
    } else if (s1 > 0.0) {
        rock_failure = 1.0;
    }
    out.failure_probability =
        clamp01(1.0 - (1.0 - clamp01(rock_failure)) * survival);
    const double radical = out.rock_mass_ucs_mpa *
                           (out.rock_mass_ucs_mpa + 8.0 * s3);
    const double griffith =
        out.rock_mass_ucs_mpa > kEps && radical >= 0.0
            ? std::max(0.0, s3 + 0.5 * (out.rock_mass_ucs_mpa +
                                        std::sqrt(radical)))
            : 0.0;
    out.rockburst_index = griffith > kEps ? s1 / griffith : 0.0;
    if (out.rock_mass_em_gpa > kEps) {
        out.smrm_score = out.rock_mass_em_gpa > 10.0
                             ? (out.rock_mass_em_gpa + 100.0) * 0.5
                             : 40.0 * std::log10(out.rock_mass_em_gpa) + 10.0;
        out.smrm_score = std::clamp(out.smrm_score, 0.0, 100.0);
        const double nu = std::clamp(out.rock_mass_poisson_ratio, 0.0, 0.499999);
        const double sum = s1 + s2 + s3;
        out.u0v_mj_m3 = (1.0 - 2.0 * nu) * square(sum) /
                        (6.0 * out.rock_mass_em_gpa * 1000.0);
        out.u0d_mj_m3 =
            (1.0 + nu) *
            (square(s1 - s2) + square(s2 - s3) + square(s3 - s1)) /
            (6.0 * out.rock_mass_em_gpa * 1000.0);
        out.energy_density_mj_m3 = out.u0v_mj_m3 + out.u0d_mj_m3;
    }

    const bool rock_shear =
        s1 > out.triaxial_compressive_strength_mpa + kEps;
    const bool rock_tension = s3 < -rock.tension_mpa - kEps;
    std::uint8_t now = 0;
    if (rock_shear || joint_shear) now |= 1;
    if (rock_tension || joint_tension) now |= 2;
    history.yield_state = static_cast<std::uint8_t>(history.yield_state | now);
    out.yield_state = history.yield_state;
    out.has_failed = out.failure_probability >= 0.99 || out.yield_state != 0;
}

}  // namespace

bool update_elastic_increment(const MaterialParameters& material,
                              const std::array<double, 6>& strain_increment,
                              std::array<double, 6>& stress_pa,
                              Matrix6& tangent_pa) {
    if (!(material.elastic_modulus_pa > 0.0) ||
        !(material.poisson_ratio > -1.0 && material.poisson_ratio < 0.5))
        return false;
    Matrix6 compliance = compliance_matrix(material, {});
    if (!invert6(compliance)) return false;
    tangent_pa = to_rocksolver_order(compliance);
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            stress_pa[row] += tangent_pa[row][column] * strain_increment[column];
    return std::all_of(stress_pa.begin(), stress_pa.end(),
                       [](double value) { return std::isfinite(value); });
}

bool update_smrm_increment(const MaterialParameters& material, double volume_m3,
                           const std::array<double, 6>& strain_increment,
                           std::array<double, 6>& stress_pa,
                           SMRMHistory& history, SMRMDiagnostics& diagnostics,
                           Matrix6& tangent_pa) {
    if (!(material.elastic_modulus_pa > 0.0) ||
        !(material.poisson_ratio >= 0.0 && material.poisson_ratio < 0.5))
        return false;
    (void)volume_m3;  // Geometry volume is not ModelSMRM's statistical volume.
    std::array<JointWork, 3> joints = make_joints(material);
    const StrengthWork rock = strength_parameters(material);
    for (JointWork& joint : joints)
        update_joint_state(joint, stress_pa, rock,
                           material.flac3d7_smrm_contact);
    Matrix6 compliance = compliance_matrix(material, joints);
    if (!invert6(compliance)) return false;
    tangent_pa = to_rocksolver_order(compliance);
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            stress_pa[row] += tangent_pa[row][column] * strain_increment[column];
    if (!std::all_of(stress_pa.begin(), stress_pa.end(),
                     [](double value) { return std::isfinite(value); }))
        return false;
    for (JointWork& joint : joints)
        update_joint_state(joint, stress_pa, rock,
                           material.flac3d7_smrm_contact);
    // ModelSMRM's statistical joint probability uses its fixed representative
    // rock volume (1 m3 in this package), not each FLAC zone's geometric
    // volume.  Passing the often-thousands-of-m3 Tet4 volume grossly inflated
    // Web failure probabilities.
    compute_diagnostics(material, material.representative_volume_m3,
                        stress_pa, joints, history,
                        diagnostics);
    for (std::size_t index = 0; index < joints.size(); ++index) {
        history.joint_k[index] = joints[index].k;
        history.joint_h[index] = joints[index].h;
    }
    ++history.step;
    return true;
}

}  // namespace rocksolver
