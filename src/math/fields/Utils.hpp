//
// Created by evgen on 11.10.2024.
//

#ifndef UTILS_HPP
#define UTILS_HPP

#include "SurfaceVectorField.hpp"


namespace EMW::Math::FieldUtils {

Math::SurfaceScalarField<Types::complex_d> relativeError(const SurfaceVectorField & v1, const SurfaceVectorField & v2);

}

#endif  //UTILS_HPP
