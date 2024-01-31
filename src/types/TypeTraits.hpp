//
// Created by evgen on 31.01.24.
//

#ifndef ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP
#define ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP

#include "Types.hpp"

namespace EMW::Types::TypeTraits {
    template<typename x_Function>
    struct function_traits;

    template<typename x_Result, typename... x_Args>
    struct function_traits<x_Result (x_Args...)> {
        using arguments = ::std::tuple<x_Args...>;
        using return_type = x_Result;
    };

    template <typename CTy, typename RTy, typename... ATy>
    struct function_traits<RTy (CTy::*)(ATy...) const> : function_traits<RTy(ATy...)> {};
}

#endif //ELECTROMAGNETIC_WAVES_SCATTERING_TYPETRAITS_HPP
