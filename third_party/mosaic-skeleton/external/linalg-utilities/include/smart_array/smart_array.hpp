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


#ifndef SMART_ARRAY_HPP
#define SMART_ARRAY_HPP

#if __cplusplus < 201103L
#define nullptr NULL
#endif


#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <stdexcept>


template <class DataType>
class SmartArray
{
private:
        void *pointer{nullptr};
        const uint64_t &reference_number() const
        {
            return *static_cast<uint64_t *>(pointer);
        }
        uint64_t &reference_number()
        {
            return *static_cast<uint64_t *>(pointer);
        }
        void setsize(uint64_t s)
        {
            static_cast<uint64_t *>(pointer)[1] = s;
        }
        void make_copy()
        {
            void *new_pointer = std::malloc(size() * sizeof(DataType) + 2 * sizeof(uint64_t));
            std::memcpy(new_pointer, pointer, size() * sizeof(DataType) + 2 * sizeof(uint64_t));
            reference_number()--;
            pointer = new_pointer;
            reference_number() = 0;
        }
    public:
        SmartArray() = default;
        ~SmartArray()
        {
            free();
        }
        SmartArray(const SmartArray &other)
        {
            pointer = other.pointer;
            if (pointer != nullptr)
            {
                reference_number()++;
            }
        }
        SmartArray(SmartArray &&other) noexcept {
            pointer = other.pointer;
            other.pointer = nullptr;
        }
        explicit SmartArray(const uint64_t number)
        {
            if (number > 0) {
                pointer = std::malloc(number * sizeof(DataType) + 2 * sizeof(uint64_t));
                reference_number() = 0;
                setsize(number);
            }
        }
        SmartArray &operator=(const SmartArray &other) {
            if (this == &other || pointer == other.pointer) {
                return *this;
            }
            free();
            pointer = other.pointer;
            if (pointer != nullptr)
            {
                reference_number()++;
            }
            return *this;
        }
        SmartArray &operator=(SmartArray &&other) noexcept {
            if (this == &other || pointer == other.pointer) {
                return *this;
            }
            free();
            pointer = other.pointer;
            other.pointer = nullptr;
            return *this;
        }


        uint64_t size() const
        {
            if (pointer == nullptr) return 0;
            return static_cast<uint64_t *>(pointer)[1];
        }
        uint64_t byte_size() const
        {
            return size() * sizeof(DataType);
        }
        const DataType &operator[](uint64_t i) const
        {
            if (pointer == nullptr || i >= size())
            {
                throw std::runtime_error("Trying to access non existing element of SmartArray");
            }
            return reinterpret_cast<DataType *>(static_cast<uint64_t *>(pointer) + 2)[i];
        }
        DataType &operator[](uint64_t i)
        {
            if (pointer == nullptr || i >= size())
            {
                throw std::runtime_error("Trying to access non existing element of SmartArray");
            }
            if (reference_number() > 0)
            {
                make_copy();
            }
            return reinterpret_cast<DataType *>(static_cast<uint64_t *>(pointer) + 2)[i];
        }
        const DataType *data() const
        {
            return reinterpret_cast<const DataType *>(static_cast<uint64_t *>(pointer) + 2);
        }
        DataType *data()
        {
            if (pointer != nullptr && reference_number() > 0)
            {
                make_copy();
            }
            return reinterpret_cast<DataType *>(static_cast<uint64_t *>(pointer) + 2);
        }
        void resize(uint64_t new_size)
        {
            if (pointer == nullptr && new_size != 0)
            {
                pointer = std::malloc(new_size * sizeof(DataType) + 2 * sizeof(uint64_t));
                reference_number() = 0;
                setsize(new_size);
            }
            else if (size() != new_size)
            {
                if (new_size == 0)
                {
                    free();
                }
                else if (reference_number() > 0)
                {
                    void *new_pointer = std::malloc(new_size * sizeof(DataType) + 2 * sizeof(uint64_t));
                    std::memcpy(new_pointer, pointer, std::min(size(), new_size) * sizeof(DataType) + 2 * sizeof(uint64_t));
                    reference_number()--;
                    pointer = new_pointer;
                    reference_number() = 0;
                    setsize(new_size);
                }
                else
                {
                    pointer = std::realloc(pointer, new_size * sizeof(DataType) + 2 * sizeof(uint64_t));
                    setsize(new_size);
                }
            }
        }
        uint64_t reserve(uint64_t need_size)
        {
            if (need_size > size()) {
                resize(need_size);
            }
            return size();
        }
        void free()
        {
            if (pointer != nullptr)
            {
                if (reference_number() > 0)
                {
                    reference_number()--;
                }
                else
                {
                    std::free(pointer);
                    pointer = nullptr;
                }
            }
        }
};

#endif //SMART_ARRAY_HPP
