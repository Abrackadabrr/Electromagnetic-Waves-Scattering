//
// Created by evgen on 10.01.2025.
//

#ifndef UTILS_HPP
#define UTILS_HPP

namespace EMW::Utils {

template <typename Container>
void to_csv(const Container &cont1, const Container &cont2, const std::string &name1, const std::string &name2,
            std::ostream &str) {
    str << name1 << "," << name2 << "\n";
    for (int i = 0; i < cont1.size(); i++) {
        str << cont1[i] << ',' << cont2[i] << '\n';
    }
}
}

#endif //UTILS_HPP
