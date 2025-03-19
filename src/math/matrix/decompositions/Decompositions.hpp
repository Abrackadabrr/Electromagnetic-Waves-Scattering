//
// Created by evgen on 18.03.2025.
//

#ifndef DECOMPOSITIONS_HPP
#define DECOMPOSITIONS_HPP

#include "rsvd.hpp"
#include "types/Types.hpp"

namespace EMW::Math::Matrix::Decompositions {

using ComplexRSVD = RSVD<Types::MatrixXc, Types::VectorXc, Types::complex_d>;
using RealRSVD = RSVD<Types::MatrixXd, Types::VectorXd, Types::scalar>;

}


#endif //DECOMPOSITIONS_HPP
