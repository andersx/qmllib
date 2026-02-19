# CI Workflow Comparison

Both Ubuntu and macOS workflows have been updated to follow the same structure and use modern tooling.

## Workflow Structure (Identical)

Both workflows follow this pattern:

1. **Checkout code** - `actions/checkout@v4`
2. **Install system dependencies** - OS-specific package managers
3. **Set up uv** - `astral-sh/setup-uv@v5` with Python 3.12
4. **Build & install** - Explicit dependency installation + package build
5. **Run tests** - Unit tests only (exclude integration tests)

## Side-by-Side Comparison

| Step | Ubuntu | macOS |
|------|--------|-------|
| **OS** | `ubuntu-latest` | `macos-latest` |
| **Python** | `3.12` | `3.12` |
| **Package Manager** | `apt-get` | `brew` |
| **System Deps** | `gfortran libomp-dev libopenblas-dev` | `gcc libomp llvm` |
| **Python Setup** | `setup-uv@v5` | `setup-uv@v5` |
| **Build Deps** | `scikit-build-core pybind11 setuptools` | `scikit-build-core pybind11 setuptools` |
| **Install** | `uv pip install -e .[test]` | `uv pip install -e .[test]` |
| **Test Command** | `uv run pytest -m "not integration" -v` | `uv run pytest -m "not integration" -v` |

## macOS-Specific Configuration

The macOS workflow includes additional environment variables to help CMake find the Homebrew-installed compilers and libraries:

```yaml
env:
  FC: gfortran-14
  CMAKE_PREFIX_PATH: /opt/homebrew
  OpenMP_ROOT: /opt/homebrew/opt/libomp
  CMAKE_ARGS: >-
    -DCMAKE_Fortran_COMPILER=gfortran-14
    -DOpenMP_ROOT=/opt/homebrew/opt/libomp
```

These are not needed on Ubuntu because the system packages install to standard locations.

## Benefits of Unified Structure

1. **Consistency** - Same tools and commands across platforms
2. **Maintainability** - Easier to update both workflows together
3. **Speed** - `uv` is faster than traditional `pip`
4. **Caching** - Both workflows enable dependency caching
5. **Fast CI** - Only unit tests run (~15s), integration tests excluded

## Testing

To verify both workflows:

1. Create a PR from `feature/macos-compile` to `main`
2. Both workflows will run automatically
3. Check GitHub Actions tab for results

Expected results:
- ✅ Ubuntu: ~2-3 minutes total (build + test)
- ✅ macOS: ~5-10 minutes total (build + test)

## Complete Workflow Files

### Ubuntu: `.github/workflows/test.ubuntu.yml`
```yaml
name: Test Ubuntu

on:
  pull_request:
    branches: [main]

jobs:
  test:
    name: Testing ${{matrix.os}} py-${{matrix.python-version}}
    runs-on: ${{matrix.os}}

    strategy:
      matrix:
        os: ['ubuntu-latest']
        python-version: ['3.12']

    steps:
      - uses: actions/checkout@v4

      - name: Install system dependencies
        run: |
          sudo apt-get update
          sudo apt-get install -y gfortran libomp-dev libopenblas-dev

      - name: Set up uv
        uses: astral-sh/setup-uv@v5
        with:
          python-version: ${{ matrix.python-version }}
          enable-cache: true

      - name: Build & install
        run: |
          uv pip install scikit-build-core pybind11 setuptools
          uv pip install -e .[test] --verbose

      - name: Run unit tests (exclude integration tests)
        run: uv run pytest -m "not integration" -v
```

### macOS: `.github/workflows/test.macos.yml`
```yaml
name: Test MacOS

on:
  pull_request:
    branches: [main]

jobs:
  test:
    name: Testing ${{matrix.os}} py-${{matrix.python-version}}
    runs-on: ${{matrix.os}}

    strategy:
      matrix:
        os: ['macos-latest']
        python-version: ['3.12']

    env:
      HOMEBREW_NO_AUTO_UPDATE: "1"

    steps:
      - uses: actions/checkout@v4

      - name: Install Homebrew dependencies
        run: |
          brew install gcc libomp llvm

      - name: Set up uv
        uses: astral-sh/setup-uv@v5
        with:
          python-version: ${{ matrix.python-version }}
          enable-cache: true

      - name: Build & install (macOS only)
        env:
          FC: gfortran-14
          CMAKE_PREFIX_PATH: /opt/homebrew
          OpenMP_ROOT: /opt/homebrew/opt/libomp
          CMAKE_ARGS: >-
            -DCMAKE_Fortran_COMPILER=gfortran-14
            -DOpenMP_ROOT=/opt/homebrew/opt/libomp
        run: |
          uv pip install scikit-build-core pybind11 setuptools
          uv pip install -e .[test] --verbose

      - name: Run unit tests (exclude integration tests)
        run: uv run pytest -m "not integration" -v
```
