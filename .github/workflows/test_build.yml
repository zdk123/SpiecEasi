name: test

on:
  push:
    branches: [ "master", "gh-actions" ]
  pull_request:
    branches: [ "master" ]
  workflow_dispatch:

permissions:
  contents: read

jobs:
  build:
    runs-on: ubuntu-latest
    # strategy:
    #   matrix:
    #     r-version: ['3.6.3', '4.1.1']


    steps:
    - uses: actions/checkout@v3
    - uses: r-lib/actions/setup-r@v2
    - uses: r-lib/actions/setup-r-dependencies@v2
      with:
        cache: true
        cache-version: 1
        extra-packages: |
          any::rcmdcheck
          any::remotes
          any::devtools
          any::pulsar
          local::.
        needs: check
    - uses: r-lib/actions/check-r-package@v2
      with:
        args: 'c("--no-manual", "--as-cran")'
        error-on: '"error"'
        check-dir: '"check"'
