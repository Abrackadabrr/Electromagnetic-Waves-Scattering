/*
 * Copyright 2023-2026 Bulat Valiakhmetov
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
//
// Created by bulat_v on 26.07.23.
//


#include "msk/utility/logging.hpp"


namespace MSk {

namespace log {

Level string_to_level(const std::string &level) {
    if (level == "trace") {
        return Level::trace;
    } else if (level == "debug") {
        return Level::debug;
    } else if (level == "info") {
        return Level::info;
    } else if (level == "warn") {
        return Level::warn;
    } else if (level == "error" || level == "err") {
        return Level::err;
    } else if (level == "critical") {
        return Level::critical;
    } else if (level == "off") {
        return Level::off;
    }
    return Level::off;
}


#ifdef MSK_USE_SPDLOG

std::shared_ptr<Logger> build_logger(const std::string &name,
                                     const std::string &file,
                                     const Level &level) {
    std::shared_ptr<Logger> log{nullptr};
    if (file == "stdout") {
        log = spdlog::stdout_color_mt(name);
    } else if (file == "stderr") {
        log = spdlog::stderr_color_mt(name);
    } else {
        log = spdlog::basic_logger_mt(name, file);
    }
    log->set_level(level);
    log->set_pattern("[%H:%M:%S:%e][%n][%t][%P][%^%l%$] %v");
    // Drop the name from the registry preemptively to avoid loggers duplication
    spdlog::drop(name);
    return log;
}
#else

std::shared_ptr<Logger> build_logger(const std::string &name,
                                     const std::string &file,
                                     const Level &level) {
    return std::make_shared<Logger>(name, file, level);
}

#endif

}
}