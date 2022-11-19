name: test

on:
  push:
    branches: [ "master" ]
  pull_request:
    branches: [ "master" ]
  workflow_dispatch:

permissions:
  contents: read

jobs:
  build:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        r-version: ['3.6.3', '4.1.1']

    steps:
      - uses: actions/checkout@v3
      - name: Set up R ${{ matrix.r-version }}
        uses: r-lib/actions/setup-r@v2
        with:
          r-version: ${{ matrix.r-version }}
      - name: Install dependencies
        run: |
          install.packages(c("remotes", "rcmdcheck"))
          remotes::install_deps(dependencies = TRUE)
        shell: Rscript {0}
      - name: Check
        run: rcmdcheck::rcmdcheck(args = "--no-manual", error_on = "error")
        shell: Rscript {0}