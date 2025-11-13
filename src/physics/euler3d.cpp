#include "physics/euler3d.hpp"

namespace fvm3d::physics {

Eigen::VectorXd EulerEquations3D::flux_x(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);
    flux(0) = rho * u;
    flux(1) = rho * u * u + p;
    flux(2) = rho * u * v;
    flux(3) = rho * u * w;
    flux(4) = (U(4) + p) * u;

    return flux;
}

Eigen::VectorXd EulerEquations3D::flux_y(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);
    flux(0) = rho * v;
    flux(1) = rho * v * u;
    flux(2) = rho * v * v + p;
    flux(3) = rho * v * w;
    flux(4) = (U(4) + p) * v;

    return flux;
}

Eigen::VectorXd EulerEquations3D::flux_z(const Eigen::VectorXd& U) const {
    double rho, u, v, w, p;
    conservative_to_primitive(U, rho, u, v, w, p);

    Eigen::VectorXd flux(5);
    flux(0) = rho * w;
    flux(1) = rho * w * u;
    flux(2) = rho * w * v;
    flux(3) = rho * w * w + p;
    flux(4) = (U(4) + p) * w;

    return flux;
}

} // namespace fvm3d::physics
