# zaiMSk

Matrix approximation in the mosaic-skeleton format.


## Installation

The library can be installed using `CMake`. 
```bash
cmake -S . -B build --prefix=/path/to/install
cmake --build build --target install
```
Some distributions of `CMake` provide `--install-prefix` option instead of `--prefix`.
More orthodox way is to set the path manually and use `make`.
```bash
mkdir -p build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install
make install
```

## Usage

You can use the installed version of MSk or as your CMake subproject.
For more examples see the [examples](https://gitlab.com/bulatral/mosaic-skeleton/tree/main/examples) directory.

### Installed library

The installed library can be found with CMake's tools and linked directly to your executable or library.

```cmake
find_package(zaimsk REQUIRED)

add_executable(my_exe main.cpp)

target_link_libraries(my_exe zaimsk::zaimsk)
```

### CMake subproject

You can also add the library as your CMake project dependency and build them together.

```cmake
add_subdirectory(mosaic-skeleton)

add_executable(my_exe main.cpp)

target_link_libraries(my_exe zaimsk::zaimsk)
```


## Support

If you have questions or bug-cases please open an Issue at this repository.


## Roadmap

In the future we will add test for separate modules and more examples. 
Soon we are aiming to extend our package for another approximation formats.


## Contributing

In case you want to contribute to zaiMSk package please open a merge request at this repository.


## Authors and acknowledgment

Bulat Valiakhmetov, Lomonosov Moscow State University

Dmitry Zheltkov, Marchuk Institute of Numerical Mathematics of Russian Academy of Sciences


## References

### Package
* Valiakhmetov, B., Tyrtyshnikov, E. (2020). MSk - the package for a dense matrix approximation in the mosaic-skeleton format.
  Russian Supercomputing Days : Proceedings of the International Conference, (Moscow, Russia), 20–27.

### Method
* Tyrtyshnikov, E. (1996). Mosaic-skeleton approximations. Calcolo, 33, 47-57.
* Tyrtyshnikov, E. (2000). Incomplete cross approximation in the mosaic-skeleton method. Computing, 64(4), 367-380.

### Applications
* Valiakhmetov, B., Lukyanenko, D., Tyrtyshnikov, E. (2026). Acceleration of Iterative Methods for Solving Linear Inverse Problems Based on Low-Rank Approximation. Computational Mathematics and Mathematical Physics, 66(1), 1-14. https://doi.org/10.1134/S0965542525701696
* Lukyanenko, D., Valiakhmetov, B., Tyrtyshnikov, E. & Yagola, A. (2025). An a priori method for estimating the informativeness of the configuration of sensor placement when solving inverse problems of remote sensing. Journal of Inverse and Ill-posed Problems, 33(1), 15-30. https://doi.org/10.1515/jiip-2024-0017
* Dyachenko, R. R., Matveev, S. A., & Valiakhmetov, B. I. (2025). Mosaic-skeleton approximation is all you need for Smoluchowski equations. [arXiv preprint](https://doi.org/10.48550/arXiv.2501.10206)
* Gladkov, A. O., Valiakhmetov, B. I., Tyrtyshnikov, E. E., & Samokhin, A. B. (2025). Comparison of Interpolation-Based and Mosaic-Skeleton Methods for Solving Integral Equations with a Convolutional Kernel. Computational Mathematics and Mathematical Physics, 65(6), 1206-1219. https://doi.org/10.1134/S0965542525700423


## License

This product is licensed under the terms of the Apache License, Version 2.0.

    Copyright (C) 2023-2026 Bulat Valiakhmetov

All rights reserved.
