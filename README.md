# libifxc

`libifxc` is a small C library for evaluating integral-feature exchange-correlation kernels. The first release focuses on ML25 and exposes a simple IFXC-native API.

This is **not** a Libxc compatibility layer. The codebase borrows useful implementation ideas from Libxc, but the public surface is intentionally different.

## How IFXC differs from Libxc

- One feature set only: ML25.
- ML25 is always evaluated as all 66 features together.
- No functional registry or functional-number lookup.
- No `xc_func_type`, `xc_func_init`, `xc_mgga`, or other Libxc public entry points.
- No external-parameter selection path.
- No one-functional-at-a-time wrappers.
- No public LDA/GGA/MGGA family split.
- The public API is centered on `ifxc_init()`, `ifxc_eval()`, and `ifxc_end()`.
- Feature metadata is fixed and public, but it is not a registry system.
- Generated ML25 formulas come from one combined Maple source, not many per-feature hand-maintained C files.

## Public API summary

The main header is [`include/ifxc.h`](include/ifxc.h).

Typical flow:

1. Initialize a handle with `ifxc_init()`.
2. Prepare input arrays and one or more derivative entries.
3. Call `ifxc_eval()`.
4. Destroy the handle with `ifxc_end()`.

For ML25, the public input uses density-variable arrays laid out by component, then point.

## Minimal usage example

This example evaluates ML25 for two unpolarized grid points and requests both local and integral feature values.

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

  printf("feature 0 local at point 0: %g\n", local[IFXC_ML25_LAK_X * input.npoints + 0]);
  printf("feature 0 integral: %g\n", integral[IFXC_ML25_LAK_X]);

  ifxc_end(&func);
  return 0;
}
```

### Build the example

If the library is already built and installed:

```sh
cc $(pkg-config --cflags libifxc) example.c $(pkg-config --libs libifxc) -lm
```

If you are building inside this repository with CMake, the usual flow is:

```sh
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

## Notes

- `IFXC_ML25_NFEATURES` is 66.
- `ifxc_output_size()` can be used to size result buffers dynamically.
- Higher-order derivative requests are validated against the generated maximum derivative order.
- Unsupported variables, such as `IFXC_VAR_LAPL` for ML25, fail explicitly.
