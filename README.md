# libifxc

`libifxc` is a C library for evaluating integral features in integral-feature
density functional theory. It is heavily inspired by the
[Libxc project](https://libxc.gitlab.io/).

## Build and install

The build requires a C99 compiler and CMake 3.16 or newer.

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build
cmake --install build --prefix /path/to/prefix
```

Generated C sources are included, so Python and Maple are not required for a
normal build. Python 3 is needed only for generated-source validation when
testing is enabled. Shared libraries are built by default; pass
`-DBUILD_SHARED_LIBS=OFF` for a static library.

Installed CMake packages can be consumed with:

```cmake
find_package(IFXC CONFIG REQUIRED)
target_link_libraries(my_target PRIVATE IFXC::ifxc)
```

A `libifxc.pc` file is also installed for `pkg-config` users.

## Using the library

Include [`ifxc.h`](include/ifxc.h) and follow the standard lifecycle:

1. Create a feature-set handle with `ifxc_init()`.
2. Evaluate one or more requests with `ifxc_eval()`.
3. Release the handle with `ifxc_end()`.

Use `ifxc_output_size()` to size result buffers. Feature-set and feature
metadata are available through `ifxc_feature_set_info()` and
`ifxc_feature_info()`.

| Feature set | Features | Maximum derivative order |
| --- | ---: | ---: |
| ML25 | 66 | 3 |
| ML26 | 69 | 2 |

### C example

This program evaluates ML25 feature integrands and their integrated values for
two spin-unpolarized grid points:

```c
#include <stdio.h>
#include "ifxc.h"

int main(void)
{
  ifxc_func_type func;
  ifxc_status status;

  double rho[2] = {0.30, 0.40};
  double sigma[2] = {0.05, 0.02};
  double tau[2] = {0.10, 0.11};
  double weights[2] = {1.0, 0.5};

  double local[IFXC_ML25_NFEATURES * 2];
  double integral[IFXC_ML25_NFEATURES];

  ifxc_input input = {
    .npoints = 2,
    .rho = rho,
    .sigma = sigma,
    .lapl = NULL,
    .tau = tau,
    .weights = weights
  };

  ifxc_deriv_entry entries[2] = {
    {
      .target = IFXC_TARGET_LOCAL,
      .order = 0,
      .out = local
    },
    {
      .target = IFXC_TARGET_INTEGRAL,
      .order = 0,
      .out = integral
    }
  };

  status = ifxc_init(&func, IFXC_FEATURE_SET_ML25, IFXC_UNPOLARIZED);
  if(status != IFXC_OK){
    fprintf(stderr, "ifxc_init failed: %s\n", ifxc_strerror(status));
    return 1;
  }

  status = ifxc_eval(&func, &input, 2, entries);
  if(status != IFXC_OK){
    fprintf(stderr, "ifxc_eval failed: %s\n", ifxc_strerror(status));
    ifxc_end(&func);
    return 1;
  }

  printf("feature 0 local at point 0: %g\n",
         local[IFXC_ML25_LAK_X]);
  printf("feature 0 integral: %g\n",
         integral[IFXC_ML25_LAK_X]);

  ifxc_end(&func);
  return 0;
}
```

With an installed library, save the program as `example.c` and compile it
with:

```sh
cc example.c $(pkg-config --cflags --libs libifxc) -lm -o example
```

The header documents input components and output layouts. See
[`tests/test_api.c`](tests/test_api.c) for complete local, integral, and
derivative examples.

## Adding an integral-feature functional

During development, add a Maple source such as
`maple/if_mgga/mgga_xc_myfunctional.mpl`, including its `type: if_mgga`,
`feature_set`, `nfeatures`, `max_order`, and `variables` headers. With Maple
available, generate the C source directly:

```sh
python3 scripts/maple2c.py \
  --srcdir . \
  --functional mgga_xc_myfunctional \
  --maxorder 2
```

The generated file is written to
`src/maple2c/if_mgga/mgga_xc_myfunctional.c`. Commit it with the matching
feature metadata and backend registration; do not edit the generated C by
hand.

## License

`libifxc` is licensed under MPL-2.0. See [`LICENSE`](LICENSE) and
[`NOTICE`](NOTICE).

## Recommended citation

**Software:** D. Zhang, Y. Shu, and D. G. Truhlar, *libifxc*, version 0.4.0,
University of Minnesota, Minneapolis, 2026.
[https://github.com/Dayou-Zhang/libifxc](https://github.com/Dayou-Zhang/libifxc).

**ML25:** Zhang, D.; Shu, Y.; Truhlar, D. G. Reinventing Density Functional
Theory with Machine Learning on Integral Features. *J. Chem. Theory Comput.*
**2026**, *22*, 6295–6303.
[https://doi.org/10.1021/acs.jctc.6c00999](https://doi.org/10.1021/acs.jctc.6c00999).

**ML26:** Zhang, D.; Shu, Y.; Truhlar, D. G. Manuscript in preparation.
