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
// Created by bulat_v on 31.08.22.
//



#include <iostream>
#include <random>

#include "msk/thread_pool/thread_pool.hpp"
#include "msk/utility/logging.hpp"


void rand_num(std::shared_ptr<MSk::log::Logger> logger, std::string pre_msg = "", int id = 0) {
    int num = std::rand() % 100;
    logger->info("{0} {1} {2}", pre_msg, id, num);
}

class Printer {
private:
    std::string name;
    int count{0};
public:
    Printer(const std::string &p_name) : name{p_name} {}

    void run(MSk::exec::ThreadPool &pool) {
        pool.execute(&Printer::print, this, "message #" + std::to_string(count++));
    }

private:
    void print(const std::string &line) {
        std::cout << name << " :: " << line << std::endl;
    }
};


int main() {
    auto log = MSk::log::build_logger("TP", "stdout");
    Printer p("Printer");
    MSk::exec::ThreadPool pool("tp");
    pool.start();
    p.run(pool);
    for (int i = 0; i < 20; ++i) {
        pool.execute(rand_num, log, "before", i);
    }
    pool.stop(true);
    pool.start();
    p.run(pool);
    log->flush();
    return 0;
}