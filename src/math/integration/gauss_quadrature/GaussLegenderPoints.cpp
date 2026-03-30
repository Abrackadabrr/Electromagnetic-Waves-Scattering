//
// Created by evgen on 18.03.2026.
//

#include "GaussLegenderPoints.hpp"

// Предопределённы квадратуры для быстрой компиляции
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<2, 2, 2, 2, 2, 2>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<3, 3, 3, 3, 3, 3>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<4, 4, 4, 4, 4, 4>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<2, 2, 2, 2>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<3, 3, 3, 3>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<4, 4, 4, 4>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<2, 2, 2>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<3, 3, 3>;
template struct EMW::DefiniteIntegrals::GaussLegendre::Quadrature<4, 4, 4>;
