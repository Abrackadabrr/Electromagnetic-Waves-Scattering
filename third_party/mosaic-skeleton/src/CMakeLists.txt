# Copyright 2023-2026 Bulat Valiakhmetov
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


add_subdirectory(block)
add_subdirectory(computation)
add_subdirectory(evaluation)
add_subdirectory(geometry)
add_subdirectory(index)
add_subdirectory(interface)
add_subdirectory(matrix)
add_subdirectory(oracle)
add_subdirectory(thread_pool)
add_subdirectory(utility)


set(SOURCE_FILES
        mosaic_skeleton.cpp
)

target_sources(zaimsk PRIVATE ${SOURCE_FILES})