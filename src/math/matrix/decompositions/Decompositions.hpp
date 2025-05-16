//
// Created by evgen on 18.03.2025.
//

#ifndef DECOMPOSITIONS_HPP
#define DECOMPOSITIONS_HPP

#include "rsvd.hpp"
#include "adaptive_cross.hpp"
#include "types/Types.hpp"

namespace EMW::Math::LinAgl::Decompositions {

using ComplexRSVD = RSVD<Types::MatrixXc, Types::VectorXc, Types::complex_d>;
using RealRSVD = RSVD<Types::MatrixXd, Types::VectorXd, Types::scalar>;

using ComplexACA = ACA<Types::MatrixXc, Types::VectorXc, Types::complex_d>;
using RealACA = ACA<Types::MatrixXd, Types::VectorXd, Types::scalar>;

}


#endif //DECOMPOSITIONS_HPP
