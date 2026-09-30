#include "rocksolver/rocksolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <numbers>
#include <utility>

namespace rocksolver {
namespace {
constexpr double kEps = 1.0e-18;
constexpr std::uint8_t kShearNow = 0x01;
constexpr std::uint8_t kTensionNow = 0x02;
constexpr std::uint8_t kShearPast = 0x04;
constexpr std::uint8_t kTensionPast = 0x08;

struct PrincipalFrame {
    std::array<double, 3> values{}; // ascending, tension positive
    std::array<std::array<double, 3>, 3> axes{};
};

double radians(double degrees) {
    return degrees * std::numbers::pi / 180.0;
}

Matrix6 isotropic_tangent(const MaterialParameters& material) {
    const double e = material.elastic_modulus_pa;
    const double nu = material.poisson_ratio;
    const double g = e / (2.0 * (1.0 + nu));
    const double lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    Matrix6 d{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            d[i][j] = (i == j ? lambda + 2.0 * g : lambda);
    d[3][3] = d[4][4] = d[5][5] = g;
    return d;
}

PrincipalFrame principal_frame(const std::array<double, 6>& stress) {
    double a[3][3]{{stress[0], stress[3], stress[4]},
                   {stress[3], stress[1], stress[5]},
                   {stress[4], stress[5], stress[2]}};
    double vectors[3][3]{{1.0, 0.0, 0.0},
                         {0.0, 1.0, 0.0},
                         {0.0, 0.0, 1.0}};
    for (int sweep = 0; sweep < 24; ++sweep) {
        int p = 0, q = 1;
        double largest = std::abs(a[0][1]);
        for (const auto pair : {std::pair{0, 2}, std::pair{1, 2}}) {
            if (std::abs(a[pair.first][pair.second]) > largest) {
                largest = std::abs(a[pair.first][pair.second]);
                p = pair.first; q = pair.second;
            }
        }
        if (largest < 1.0e-9) break;
        const double angle = 0.5 * std::atan2(2.0 * a[p][q], a[p][p] - a[q][q]);
        const double c = std::cos(angle), s = std::sin(angle);
        const double app = a[p][p] * c * c + a[q][q] * s * s + 2.0 * a[p][q] * s * c;
        const double aqq = a[p][p] * s * s + a[q][q] * c * c - 2.0 * a[p][q] * s * c;
        a[p][p] = app; a[q][q] = aqq; a[p][q] = a[q][p] = 0.0;
        for (int row = 0; row < 3; ++row) {
            if (row != p && row != q) {
                const double arp = a[row][p] * c + a[row][q] * s;
                const double arq = -a[row][p] * s + a[row][q] * c;
                a[row][p] = a[p][row] = arp; a[row][q] = a[q][row] = arq;
            }
            const double vrp = vectors[row][p] * c + vectors[row][q] * s;
            const double vrq = -vectors[row][p] * s + vectors[row][q] * c;
            vectors[row][p] = vrp; vectors[row][q] = vrq;
        }
    }
    std::array<int, 3> order{0, 1, 2};
    std::sort(order.begin(), order.end(), [&](int lhs, int rhs) {
        return a[lhs][lhs] < a[rhs][rhs];
    });
    PrincipalFrame frame;
    for (int i = 0; i < 3; ++i) {
        frame.values[i] = a[order[i]][order[i]];
        for (int k = 0; k < 3; ++k)
            frame.axes[i][k] = vectors[k][order[i]];
    }
    return frame;
}

void reconstruct(const PrincipalFrame& frame, std::array<double, 6>& stress) {
    double s[3][3]{};
    for (int i = 0; i < 3; ++i)
        for (int r = 0; r < 3; ++r)
            for (int c = 0; c < 3; ++c)
                s[r][c] += frame.values[i] * frame.axes[i][r] * frame.axes[i][c];
    stress = {s[0][0], s[1][1], s[2][2], s[0][1], s[0][2], s[1][2]};
}

}  // namespace

bool update_mohr_coulomb_increment(
    const MaterialParameters& material,
    const std::array<double, 6>& strain_increment,
    std::array<double, 6>& stress_pa,
    std::uint8_t& plastic_state, Matrix6& tangent_pa) {
    if (!(material.elastic_modulus_pa > 0.0) ||
        !(material.poisson_ratio > -1.0 && material.poisson_ratio < 0.5) ||
        material.cohesion_pa < 0.0 || material.tensile_strength_pa < 0.0 ||
        !(material.friction_deg >= 0.0 && material.friction_deg < 90.0) ||
        !(material.dilation_deg >= 0.0 && material.dilation_deg < 90.0) ||
        material.dilation_deg > material.friction_deg)
        return false;

    tangent_pa = isotropic_tangent(material);
    for (int row = 0; row < 6; ++row)
        for (int column = 0; column < 6; ++column)
            stress_pa[row] += tangent_pa[row][column] * strain_increment[column];

    if (plastic_state & kShearNow) plastic_state |= kShearPast;
    if (plastic_state & kTensionNow) plastic_state |= kTensionPast;
    plastic_state &= static_cast<std::uint8_t>(~(kShearNow | kTensionNow));

    const double reduction = std::max(material.strength_reduction_factor, 1.0);
    const double cohesion = material.cohesion_pa / reduction;
    const double phi = std::atan(std::tan(radians(material.friction_deg)) / reduction);
    const double psi = std::atan(std::tan(radians(material.dilation_deg)) / reduction);
    double tension = material.tensile_strength_pa / reduction;

    const double e = material.elastic_modulus_pa;
    const double nu = material.poisson_ratio;
    const double g = e / (2.0 * (1.0 + nu));
    const double bulk = e / (3.0 * (1.0 - 2.0 * nu));
    const double e1 = bulk + 4.0 * g / 3.0;
    const double e2 = bulk - 2.0 * g / 3.0;
    const double sin_phi = std::sin(phi);
    const double nph = (1.0 + sin_phi) / std::max(1.0 - sin_phi, kEps);
    const double csn = 2.0 * cohesion * std::sqrt(nph);
    const double sin_psi = std::sin(psi);
    const double nps = (1.0 + sin_psi) / std::max(1.0 - sin_psi, kEps);
    const double rc = std::sqrt(1.0 + nph * nph);
    const double ra = e1 - nps * e2;
    const double rb = e2 - nps * e1;
    const double rd = ra - rb * nph;
    if (!std::isfinite(rd) || std::abs(rd) <= kEps) return false;
    const double sc1 = ra / rd;
    const double sc3 = rb / rd;
    const double sc2 = e2 * (1.0 - nps) / rd;
    if (phi > kEps) tension = std::min(tension, cohesion / std::tan(phi));

    auto frame = principal_frame(stress_pa);
    auto& p1 = frame.values[0];
    auto& p2 = frame.values[1];
    auto& p3 = frame.values[2];
    const double fs = nph * p3 - p1 - csn;
    const double fsd = fs / rc;
    const double ftz = p3 - tension;

    if (fsd > 0.0 && fsd >= ftz) {
        p1 += fs * sc1;
        p2 += fs * sc2;
        p3 += fs * sc3;
        plastic_state |= kShearNow;
    } else if (ftz > 0.0 && ftz >= fsd) {
        const double ftx = p1 - tension;
        const double fty = p2 - tension;
        if (ftx > 0.0) {
            p1 = p2 = p3 = tension;
        } else if (fty > 0.0) {
            p1 -= (fty + ftz) * e2 / (e1 + e2);
            p2 = p3 = tension;
        } else {
            const double tco = ftz * e2 / e1;
            p1 -= tco;
            p2 -= tco;
            p3 = tension;
        }
        plastic_state |= kTensionNow;
    }

    if (phi > kEps) {
        const double apex = cohesion / std::tan(phi);
        if (p1 >= apex || p2 >= apex || p3 >= apex) {
            p1 = p2 = p3 = apex;
            plastic_state |= kTensionNow;
        }
    }

    reconstruct(frame, stress_pa);
    return std::all_of(stress_pa.begin(), stress_pa.end(),
                       [](double value) { return std::isfinite(value); });
}

}  // namespace rocksolver
