#include "physics/resistive_mhd3d.hpp"

namespace fvm3d::physics {

Eigen::VectorXd ResistiveMHD3D::flux_x(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd F(8);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    F(0) = rho * u;

    // Momentum fluxes
    F(1) = rho * u * u + p_total - Bx * Bx / MU0;
    F(2) = rho * u * v - Bx * By / MU0;
    F(3) = rho * u * w - Bx * Bz / MU0;

    // Energy flux
    double u_dot_B = u*Bx + v*By + w*Bz;
    F(4) = (U(4) + p_total) * u - Bx * u_dot_B / MU0;

    // Magnetic field fluxes (Faraday's law)
    F(5) = 0.0;  // div B = 0 constraint: dBx/dx = 0 in flux form
    F(6) = v * Bx - u * By;
    F(7) = w * Bx - u * Bz;

    return F;
}

Eigen::VectorXd ResistiveMHD3D::flux_y(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd G(8);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    G(0) = rho * v;

    // Momentum fluxes
    G(1) = rho * v * u - By * Bx / MU0;
    G(2) = rho * v * v + p_total - By * By / MU0;
    G(3) = rho * v * w - By * Bz / MU0;

    // Energy flux
    double u_dot_B = u*Bx + v*By + w*Bz;
    G(4) = (U(4) + p_total) * v - By * u_dot_B / MU0;

    // Magnetic field fluxes (Faraday's law)
    G(5) = u * By - v * Bx;
    G(6) = 0.0;  // div B = 0 constraint: dBy/dy = 0 in flux form
    G(7) = w * By - v * Bz;

    return G;
}

Eigen::VectorXd ResistiveMHD3D::flux_z(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p, Bx, By, Bz;
    conservative_to_primitive(U, rho, u, v, w, p, Bx, By, Bz);

    Eigen::VectorXd H(8);
    double B_sq = Bx*Bx + By*By + Bz*Bz;
    double p_total = p + 0.5 * B_sq / MU0;

    // Mass flux
    H(0) = rho * w;

    // Momentum fluxes
    H(1) = rho * w * u - Bz * Bx / MU0;
    H(2) = rho * w * v - Bz * By / MU0;
    H(3) = rho * w * w + p_total - Bz * Bz / MU0;

    // Energy flux
    double u_dot_B = u*Bx + v*By + w*Bz;
    H(4) = (U(4) + p_total) * w - Bz * u_dot_B / MU0;

    // Magnetic field fluxes (Faraday's law)
    H(5) = u * Bz - w * Bx;
    H(6) = v * Bz - w * By;
    H(7) = 0.0;  // div B = 0 constraint: dBz/dz = 0 in flux form

    return H;
}

Eigen::VectorXd ResistiveMHD3D::resistive_source(
    const Eigen::VectorXd& U,
    double dx, double dy, double dz,
    const Eigen::VectorXd& U_xm, const Eigen::VectorXd& U_xp,
    const Eigen::VectorXd& U_ym, const Eigen::VectorXd& U_yp,
    const Eigen::VectorXd& U_zm, const Eigen::VectorXd& U_zp
) const {
    Eigen::VectorXd S = Eigen::VectorXd::Zero(8);

    // If no resistivity, return zero source
    if (eta_ <= 0.0) {
        return S;
    }

    // Extract magnetic field components from neighboring cells
    double Bx_xm = U_xm(5), Bx_xp = U_xp(5);
    double By_ym = U_ym(6), By_yp = U_yp(6);
    double Bz_zm = U_zm(7), Bz_zp = U_zp(7);

    double Bx = U(5), By = U(6), Bz = U(7);

    double Bx_xm_c = U_xm(5), Bx_xp_c = U_xp(5);
    double By_xm = U_xm(6), By_xp = U_xp(6);
    double Bz_xm = U_xm(7), Bz_xp = U_xp(7);

    double Bx_ym = U_ym(5), Bx_yp = U_yp(5);
    double By_ym_c = U_ym(6), By_yp_c = U_yp(6);
    double Bz_ym = U_ym(7), Bz_yp = U_yp(7);

    double Bx_zm = U_zm(5), Bx_zp = U_zp(5);
    double By_zm = U_zm(6), By_zp = U_zp(6);
    double Bz_zm_c = U_zm(7), Bz_zp_c = U_zp(7);

    // Compute current density J = curl(B)
    // Jx = dBz/dy - dBy/dz
    double dBz_dy = (Bz_yp - Bz_ym) / (2.0 * dy);
    double dBy_dz = (By_zp - By_zm) / (2.0 * dz);
    double Jx = dBz_dy - dBy_dz;

    // Jy = dBx/dz - dBz/dx
    double dBx_dz = (Bx_zp - Bx_zm) / (2.0 * dz);
    double dBz_dx = (Bz_xp - Bz_xm) / (2.0 * dx);
    double Jy = dBx_dz - dBz_dx;

    // Jz = dBy/dx - dBx/dy
    double dBy_dx = (By_xp - By_xm) / (2.0 * dx);
    double dBx_dy = (Bx_yp - Bx_ym) / (2.0 * dy);
    double Jz = dBy_dx - dBx_dy;

    // Resistive force: F = eta * J (dissipation)
    // This acts on the magnetic field evolution
    S(5) = eta_ * Jx;  // dBx/dt source
    S(6) = eta_ * Jy;  // dBy/dt source
    S(7) = eta_ * Jz;  // dBz/dt source

    // Energy dissipation: eta * J^2
    double J_sq = Jx*Jx + Jy*Jy + Jz*Jz;
    S(4) = eta_ * J_sq;  // Energy dissipation

    return S;
}

} // namespace fvm3d::physics
