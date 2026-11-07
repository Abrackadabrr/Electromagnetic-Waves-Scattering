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
// Created by b.valiakhmetov on 9/11/23.
//


/**
 * @file
 *
 * @brief Containg logging interface utilities to spdlog library.
 */


#ifndef MOSAIC_SKELETON_UTILITY_LOGGING_HPP
#define MOSAIC_SKELETON_UTILITY_LOGGING_HPP


#include <cstddef>
#include <cstdint>
#include <ctime>
#include <chrono>
#include <fstream>
#include <inttypes.h>
#include <iomanip>
#include <ios>
#include <iostream>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>
#include <sstream>
#include <type_traits>

#include "msk_config.h"

#ifdef MSK_USE_SPDLOG
#include "spdlog/spdlog.h"
#include "spdlog/sinks/sink.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/sinks/basic_file_sink.h"
#include <spdlog/common.h>
#include <spdlog/logger.h>
#endif


namespace MSk {

namespace log {



#ifdef MSK_USE_SPDLOG

using Logger = typename spdlog::logger;
using Level = typename spdlog::level::level_enum;

#else


#define MSK_LOG_LEVEL_TRACE 0
#define MSK_LOG_LEVEL_DEBUG 1
#define MSK_LOG_LEVEL_INFO 2
#define MSK_LOG_LEVEL_WARN 3
#define MSK_LOG_LEVEL_ERROR 4
#define MSK_LOG_LEVEL_CRITICAL 5
#define MSK_LOG_LEVEL_OFF 6

struct FormatSpec {
    char fill = ' ';         // Fill symbol
    enum Align { LEFT, RIGHT, CENTER, NUMERIC } align = RIGHT;
    bool sign_plus = false;  // Always show + sign
    bool sign_space = false; // Show space for positive numbers and - for negative
    bool alternate = false;  // Alternative form (#)
    bool zero_pad = false;   // Zero padding for numbers
    int width = 0;           // Minimal window width
    int precision = -1;      // Numbers precision formatting
    char type = '\0';        // Formatting type
    
    // Numeric flags
    bool is_integer = false;
    bool is_floating = false;
    bool is_string = false;
public:
    void specify_type(const char *) {
        type = 's';
        is_string = true;
        is_integer = is_floating = false;
    }

    template<typename T, typename std::enable_if<std::is_same<T, std::string>::value, bool>::type = true>
    void specify_type(const std::string &) {
        type = 's';
        is_string = true;
        is_integer = is_floating = false;
    }
    
    template<typename T, typename std::enable_if<std::is_same<T, char>::value, bool>::type = true>
    void specify_type(const std::string &) {
        type = 'c';
        is_string = true;
        is_integer = is_floating = false;
    }

    template<typename T, typename std::enable_if<std::is_integral<T>::value, bool>::type = true>
    void specify_type(const T &) {
        type = 'd';
        is_integer = true;
        is_string = is_floating = false;
    }

    template<typename T, typename std::enable_if<std::is_floating_point<T>::value, bool>::type = true>
    void specify_type(const T &) {
        type = 'f';
        is_floating = true;
        is_integer = is_string = false;
    }

    template<typename T, typename std::enable_if<
            !std::is_same<T, std::string>::value &&
            !std::is_same<T, char>::value &&
            !std::is_integral<T>::value &&
            !std::is_floating_point<T>::value, 
        bool>::type = true>
    void specify_type(const T &) {
        type = '\0';
        is_string = is_integer = is_floating = false;
    }
};

class FormatParser {
private:
    std::string spec;
    size_t pos = 0;
    
    char peek() const {
        return (pos < spec.size()) ? spec[pos] : '\0';
    }
    
    char get() {
        return (pos < spec.size()) ? spec[pos++] : '\0';
    }
    
    bool match(char c) {
        if (peek() == c) {
            pos++;
            return true;
        }
        return false;
    }
    
    bool is_align_char(char c) {
        return c == '<' || c == '>' || c == '^' || c == '=';
    }
    
    bool is_sign_char(char c) {
        return c == '+' || c == '-' || c == ' ';
    }
    
    bool is_type_char(char c) {
        return std::isalpha(c);
    }
    
    int parse_number() {
        int result = 0;
        while (std::isdigit(peek())) {
            result = result * 10 + (get() - '0');
        }
        return result;
    }

public:
    FormatParser(const std::string& format_spec) : spec(format_spec) {}
    
    FormatSpec parse() {
        FormatSpec fs;
        
        // 1. [[fill]align]
        if (is_align_char(peek())) {
            // Align character only
            char align_char = get();
            switch (align_char) {
                case '<': fs.align = FormatSpec::LEFT; break;
                case '>': fs.align = FormatSpec::RIGHT; break;
                case '^': fs.align = FormatSpec::CENTER; break;
                case '=': fs.align = FormatSpec::NUMERIC; break;
                default: fs.align = FormatSpec::RIGHT;
            }
        } else if (pos + 1 < spec.size() && is_align_char(spec[pos + 1])) {
            // Fill and align characters
            fs.fill = get();
            char align_char = get();
            switch (align_char) {
                case '<': fs.align = FormatSpec::LEFT; break;
                case '>': fs.align = FormatSpec::RIGHT; break;
                case '^': fs.align = FormatSpec::CENTER; break;
                case '=': fs.align = FormatSpec::NUMERIC; break;
            }
        }
        
        // 2. [sign]
        if (match('+')) {
            fs.sign_plus = true;
        } else if (match('-')) {
            // default: only - for negatives
        } else if (match(' ')) {
            fs.sign_space = true;
        }
        
        // 3. [#]
        if (match('#')) {
            fs.alternate = true;
        }
        
        // 4. [0]
        if (match('0')) {
            fs.zero_pad = true;
        }
        
        // 5. [width]
        if (std::isdigit(peek())) {
            fs.width = parse_number();
        }
        
        // 6. [.precision]
        if (match('.')) {
            if (std::isdigit(peek())) {
                fs.precision = parse_number();
            } else {
                throw std::runtime_error("Expected precision after '.'");
            }
        }
        
        // 7. [type]
        if (is_type_char(peek())) {
            fs.type = get();
            
            // Type category
            switch (fs.type) {
                case 'd': case 'b': case 'o': case 'x': case 'X':
                case 'c': case 'n': case 'B': case 'O':
                    fs.is_integer = true;
                    break;
                case 'e': case 'E': case 'f': case 'F': case 'g': case 'G':
                case 'a': case 'A':
                    fs.is_floating = true;
                    break;
                case 's':
                    fs.is_string = true;
                    break;
                default:
                    throw std::runtime_error("Unknown format type: " + std::string(1, fs.type));
            }
        }
        
        return fs;
    }
};



template<typename T>
std::string format_argument(const std::string &format_spec, const T& value) {
    // Parse format string
    FormatSpec spec = FormatParser(format_spec).parse();
    if (spec.type == '\0') {
        spec.specify_type(value);
    }
    if (spec.type == 's') {
        return format_argument(format_spec, value);
    }

    std::ostringstream oss;

    if (spec.type == '\0') {
        oss << value;
    } else if (spec.is_integer) {
        oss << std::noshowbase;
        // Base
        switch (spec.type) {
            case 'b': case 'B':
                if (spec.alternate) oss << "0b";
                oss << std::dec;
                break;
            case 'o': case 'O':
                if (spec.alternate) oss << "0";
                oss << std::oct;
                break;
            case 'x':
                if (spec.alternate) oss << "0x";
                oss << std::hex << std::nouppercase;
                break;
            case 'X':
                if (spec.alternate) oss << "0X";
                oss << std::hex << std::uppercase;
                break;
            case 'd': default:
                oss << std::dec;
                break;
        }
        
        // Sign
        if (spec.sign_plus) {
            oss << std::showpos;
        } else {
            oss << std::noshowpos;
        }
        
        oss << value;
    } else if (spec.is_floating) {
        oss << std::noshowpoint;
        
        // Set precision
        if (spec.precision >= 0) {
            oss << std::setprecision(spec.precision);
            if (spec.type == 'f' || spec.type == 'F' || spec.type == '%') {
                oss << std::fixed;
            } else if (spec.type == 'e' || spec.type == 'E') {
                oss << std::scientific;
                if (spec.type == 'E') oss << std::uppercase;
            } else if (spec.type == 'g' || spec.type == 'G') {
                oss << std::defaultfloat;
                if (spec.type == 'G') oss << std::uppercase;
            } else if (spec.type == 'a' || spec.type == 'A') {
                oss << std::hexfloat;
                if (spec.type == 'A') oss << std::uppercase;
            }
        }
        
        // Sign
        if (spec.sign_plus) {
            oss << std::showpos;
        } else if (spec.sign_space) {
            oss << std::showpos;
        }
        
        // Alternative form
        if (spec.alternate) {
            oss << std::showpoint;
        }
        
        oss << value;
        
    }
    
    std::string result = oss.str();
    
    if ((spec.is_integer || spec.is_floating) && spec.sign_space && 
        !spec.sign_plus && result[0] != '-') {
        result = " " + result;
    }
    
    // Apply width and alignment parameters
    if (spec.width > 0) {
        int len = static_cast<int>(result.length());
        if (len < spec.width) {
            int pad = spec.width - len;
            
            if (spec.zero_pad && spec.align == FormatSpec::NUMERIC) {
                // Numbers zero padding
                size_t sign_pos = result.find_first_of("+- ");
                if (sign_pos != std::string::npos) {
                    result.insert(sign_pos + 1, std::string(pad, '0'));
                } else {
                    result.insert(0, std::string(pad, '0'));
                }
            } else if (spec.zero_pad && spec.align == FormatSpec::RIGHT) {
                // General zero padding
                if (result[0] == '+' || result[0] == '-' || result[0] == ' ') {
                    result[0] = result[0];
                    result.insert(1, std::string(pad, '0'));
                } else {
                    result.insert(0, std::string(pad, '0'));
                }
            } else {
                switch (spec.align) {
                    case FormatSpec::LEFT:
                        result.append(std::string(pad, spec.fill));
                        break;
                    case FormatSpec::RIGHT:
                        result.insert(0, std::string(pad, spec.fill));
                        break;
                    case FormatSpec::CENTER:
                        result.insert(0, std::string(pad / 2, spec.fill));
                        result.append(std::string(pad - pad / 2, spec.fill));
                        break;
                    case FormatSpec::NUMERIC:
                        size_t sign_pos = result.find_first_of("+- ");
                        if (sign_pos != std::string::npos) {
                            result.insert(sign_pos + 1, std::string(pad, spec.fill));
                        } else {
                            result.insert(0, std::string(pad, spec.fill));
                        }
                        break;
                }
            }
        }
    }
    
    return result;
}

// Strings specialization
template<>
inline std::string format_argument<std::string>(const std::string& format_spec, const std::string& value) {
    FormatParser parser(format_spec);
    FormatSpec spec = parser.parse();
    spec.type = 's';
    
    std::string result = value;
    
    // Cutting string
    if (spec.precision >= 0 && spec.precision < static_cast<int>(result.length())) {
        result = result.substr(0, spec.precision);
    }
    
    // Width and alignment
    if (spec.width > 0) {
        int len = static_cast<int>(result.length());
        if (len < spec.width) {
            int pad = spec.width - len;
            switch (spec.align) {
                case FormatSpec::LEFT:
                    result.append(std::string(pad, spec.fill));
                    break;
                case FormatSpec::NUMERIC:
                case FormatSpec::RIGHT:
                    result.insert(0, std::string(pad, spec.fill));
                    break;
                case FormatSpec::CENTER:
                    result.insert(0, std::string(pad / 2, spec.fill));
                    result.append(std::string(pad - pad / 2, spec.fill));
                    break;
            }
        }
    }
    return result;
}

inline std::string format_argument(const std::string& format_spec, const char *value) {
    return format_argument(format_spec, std::string(value));
}


enum class Level: uint8_t {
    off = MSK_LOG_LEVEL_OFF,
    critical = MSK_LOG_LEVEL_CRITICAL,
    err = MSK_LOG_LEVEL_ERROR,
    warn = MSK_LOG_LEVEL_WARN,
    info = MSK_LOG_LEVEL_INFO,
    debug = MSK_LOG_LEVEL_DEBUG,
    trace = MSK_LOG_LEVEL_TRACE,
};


class Stream {
private:
    mutable std::mutex _m;
    mutable std::size_t ref_count{0};
    std::ostream *_stream{nullptr};
    const std::string _file{};

private:
    static std::string format_arg_by_index(const int idx, const std::string &format) {
        return format;
    }

    template<typename T>
    static std::string format_arg_by_index(const int idx, const std::string &format, T arg) {
        if (idx != 0) {
            return "";
        }
        return format_argument(format, arg);
    }

    template<typename T, typename... Args>
    static std::string format_arg_by_index(const int idx, const std::string &format, T arg, Args... args) {
        if (idx > 0) {
            return format_arg_by_index(idx - 1, format, std::forward<Args>(args)...);
        }
        return format_argument(format, arg);

    }

    template<typename... Args>
    std::string parse_format_string(const std::string &format, Args... args) const {
        const std::size_t n_args = sizeof...(args);
        std::stringstream result;   
                
        std::size_t pos = 0;
        std::size_t arg_index = 0;
        
        while (pos < format.length()) {
            // Check for \{ and \}
            if (pos + 1 < format.length() && format[pos] == '\\' && format[pos + 1] == '{') {
                result << '{';
                pos += 2;
            } else if (pos + 1 < format.length() && format[pos] == '\\' && format[pos + 1] == '}') {
                result << '}';
                pos += 2;
            } else if (format[pos] == '{') {
                size_t end_pos = format.find('}', pos);
                if (end_pos == std::string::npos) {
                    result << "<>FURTHER PARSING ERROR<>";
                    break;
                }
                
                std::string token = format.substr(pos + 1, end_pos - pos - 1);
                
                size_t colon_pos = token.find(':');
                std::string index_str = (colon_pos == std::string::npos) ? token : token.substr(0, colon_pos);
                std::string format_str = (colon_pos == std::string::npos) ? "" : token.substr(colon_pos + 1);
                
                size_t idx;
                if (index_str.empty()) {
                    idx = arg_index++;
                } else {
                    idx = std::stoi(index_str);
                }
                if (idx >= n_args) {
                    result << "<>INVALID ARGUMENT INDEX<>";
                } else {
                    result << format_arg_by_index(idx, format_str, args...);
                }
                pos = end_pos + 1;
            } else {
                result << format[pos++];
            }
        }        
        return result.str();
    }

    std::string current_time() const {
        using namespace std::chrono;
        auto now = system_clock::now();
        auto ms = duration_cast<milliseconds>(now.time_since_epoch()) % 1000;
        auto timer = system_clock::to_time_t(now);
        std::tm bt = *std::localtime(&timer);

        std::ostringstream oss;
        oss << std::put_time(&bt, "%F %T");  // Equivalent to %Y-%m-%d %H:%M:%S.%MMM
        oss << '.' << std::setfill('0') << std::setw(3) << ms.count();
        return oss.str();
    }

public:
    Stream(const std::string &file, std::ostream *stream) : _file{file}, _stream{stream} {}

    template<typename... Args>
    void log(const std::string &name, const std::string &format, Args... args) const {
        const std::string timestamp = current_time();
        std::string result = parse_format_string(format, args...);
        result = "[" + timestamp + "] [" + name + "] " + result + "\n";
        {
            std::lock_guard<std::mutex> lock(_m);
            *_stream << result;
        }
    }

    template<typename... Args>
    void log_levelled(const std::string &level, const std::string &name, 
                      const std::string &format, Args... args) const {
        const std::string l_time = "[" + current_time() + "]";
        const std::string l_name = name.empty() ? "" : "[" + name + "]";
        const std::string l_mesg = " " + parse_format_string(format, args...);
        {
            std::lock_guard<std::mutex> lock(_m);
            *_stream << l_time << l_name << l_mesg << "\n";
        }
    }

    std::size_t out() const {
        std::lock_guard<std::mutex> lock(_m);
        if (ref_count > 0) {
            --ref_count;
        }
        if (ref_count == 0 && !_file.empty()) {
            delete _stream;
        }
        return ref_count;
    }

    std::size_t in() const {
        std::lock_guard<std::mutex> lock(_m);
        ++ref_count;
        return ref_count;
    }

    void flush() const {
        std::lock_guard<std::mutex> lock(_m);
        _stream->flush();
    }

};

class Logger {
private:
    const std::string _name{""};
    Level _level{Level::off};
    Stream *_stream{nullptr};

private:
    Logger(const std::string &name = "") : _name{name} {}

public:

    Logger(const std::string &name, std::ostream &stream, const Level level = Level::off) : 
            _name{name}, _level{level}, _stream{new Stream{"", &stream}} {
        _stream->in();
    }

    Logger(const std::string &name, const std::string &file, const Level level = Level::off) : 
            _name{name}, _level{level} {
        if (file == "stdout" || file.empty()) {
            _stream = new Stream{"", &std::cout};
        } else if (file == "stderr" ) {
            _stream = new Stream{"", &std::cerr};
        } else {
            auto stream = new std::ofstream{file, std::ios_base::out | std::ios_base::app};
            _stream = new Stream{file, stream};
        }
        _stream->in();
    }

    ~Logger() {
        if (_stream->out() == 0) {
            // No references left
            delete _stream;
        }
    }

    Logger(const Logger &) = delete;
    
    Logger(Logger &&other) noexcept = delete;
    // Logger(Logger &&other) noexcept : 
    //         _name{other._name}, _level{other._level}, 
    //         _stream{std::move(other._stream)} {
    //     other._level = Level::off;
    //     other._stream = nullptr;
    // }

    Logger &operator=(const Logger &) = delete;

    Logger &operator=(Logger &&other) noexcept = delete;
    //  {
    //     // _name = other.name; // Name cannot be changed
    //     _level = other._level;
    //     _stream = std::move(other._stream);
    //     other._level = Level::off;
    //     other._stream = nullptr;
    //     return *this;
    // }

    std::shared_ptr<Logger> clone(const std::string &name) const {
        auto logger_ptr = new Logger{name};
        logger_ptr->_level = _level;
        logger_ptr->_stream = _stream;
        _stream->in();
        return std::shared_ptr<Logger>{logger_ptr};
    }

    template<typename... Args>
    void log(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::critical) {
            _stream->log(_name, format, args...);
        }
    }

    template<typename... Args>
    void critical(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::critical) {
            _stream->log_levelled("critical", _name, format, args...);
        }
    }

    template<typename... Args>
    void error(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::err) {
            _stream->log_levelled("error", _name, format, args...);
        }
    }

    template<typename... Args>
    void warn(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::warn) {
            _stream->log_levelled("warning", _name, format, args...);
        }
    }

    template<typename... Args>
    void info(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::info) {
            _stream->log_levelled("info", _name, format, args...);
        }      
    }

    template<typename... Args>
    void debug(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::debug) {
            _stream->log_levelled("debug", _name, format, args...);
        }
    }

    template<typename... Args>
    void trace(const std::string &format, Args... args) const {
        if (_stream && _level <= Level::trace) {
            _stream->log_levelled("trace", _name, format, args...);
        }
    }

    // Return true if logging is enabled for the given level
    bool should_log(Level level) const {
        return _level <= level;
    }

    void set_level(Level level) {
        _level = level;
    }

    Level level() const {
        return _level;
    }

    const std::string &name() const {
        return _name;
    }

    void flush() const {
        _stream->flush();
    }
};

#endif

/**
 * @brief Create Logger from its components.
 *
 * @details If filename is empty or equals "stdout" then 
 * logging will be directed to standard output.
 * If "stderr" then to standard errors stream.
 *
 * @param[in] name Logger name.
 * @param[in] file Logger filename.
 * @param[in] level Logging level.
 * @return Logger.
 */
std::shared_ptr<Logger> build_logger(const std::string &name = "",
                                     const std::string &file = "",
                                     const Level &level = Level::info);

/**
 * @brief Convert a string to spdlog level.
 *
 * @param[in] level Level name.
 * @return Level value.
 */
Level string_to_level(const std::string &level);

}

}

#endif //MOSAIC_SKELETON_UTILITY_LOGGING_HPP
