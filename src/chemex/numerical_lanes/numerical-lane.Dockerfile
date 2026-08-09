# syntax=docker/dockerfile:1.7
# The lane definition hashes this complete recipe. Build one image per CPython
# source tuple below; no tag is allowed to stand in for the Debian manifest.
ARG PLATFORM_MANIFEST=debian:bookworm-slim@sha256:362e64223cc0da95422b3b13c045186fc0a81250e765d31c025fbddf257f6143
FROM --platform=linux/amd64 ${PLATFORM_MANIFEST}
ARG PLATFORM_MANIFEST

ARG PYTHON_VERSION
ARG PYTHON_SOURCE_SHA256
ARG UV_VERSION=0.11.15

# Canonical build arguments:
# - 3.13.5 / e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385
# - 3.14.5 / 9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc

RUN test -n "${PYTHON_VERSION}" && test -n "${PYTHON_SOURCE_SHA256}"
RUN printf 'deb [check-valid-until=no] http://snapshot.debian.org/archive/debian/20260803T000000Z bookworm main\n' > /etc/apt/sources.list \
    && apt-get -o Acquire::Check-Valid-Until=false update \
    && apt-get install --yes --no-install-recommends \
        build-essential ca-certificates curl libbz2-dev libffi-dev libgdbm-dev \
        liblzma-dev libncursesw5-dev libreadline-dev libsqlite3-dev \
        libssl-dev tk-dev uuid-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

RUN curl --fail --location --output /tmp/Python.tgz \
        "https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz" \
    && echo "${PYTHON_SOURCE_SHA256}  /tmp/Python.tgz" | sha256sum --check --strict \
    && tar --extract --file /tmp/Python.tgz --directory /tmp \
    && cd "/tmp/Python-${PYTHON_VERSION}" \
    && ./configure --prefix=/opt/python --enable-optimizations --with-lto \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf /tmp/Python.tgz "/tmp/Python-${PYTHON_VERSION}"

ENV PATH=/opt/python/bin:${PATH} \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    CHEMEX_NUMERICAL_LANE_WORKERS=1 \
    NPY_DISABLE_CPU_FEATURES=AVX512F,AVX512CD,AVX512_SKX

COPY pyproject.toml uv.lock /workspace/
COPY src /workspace/src
WORKDIR /workspace
RUN printf '%s\n' \
        "uv==${UV_VERSION} --hash=sha256:98edf1bdaf82447014852051d93e3ee95012509c567bf057fd117e6bdbd9a807" \
        > /tmp/uv-requirements.txt \
    && /opt/python/bin/python -m pip install --no-cache-dir --only-binary=:all: \
        --require-hashes --requirement /tmp/uv-requirements.txt \
    && uv sync --locked --no-dev \
    && .venv/bin/python -c "import numpy; c = numpy.__config__.CONFIG['Build Dependencies']['blas']; f = numpy._core._multiarray_umath.__cpu_features__; assert (c['name'], c['version']) == ('scipy-openblas', '0.3.33.112.0'); assert f['X86_V3'] and not f['X86_V4']"

# Provenance is derived from the image's immutable build inputs, never from
# launcher-provided environment variables. Runtime attestation reads this file.
RUN install --directory /opt/chemex-numerical-lane \
    && recipe_hash="$(sha256sum src/chemex/numerical_lanes/numerical-lane.Dockerfile | cut -d ' ' -f 1)" \
    && lockfile_hash="$(sha256sum uv.lock | cut -d ' ' -f 1)" \
    && printf '{"build_recipe_hash":"%s","dependency_lock_hash":"%s","platform_manifest":"%s","python_source_hash":"%s"}\n' \
        "$recipe_hash" "$lockfile_hash" "$PLATFORM_MANIFEST" "$PYTHON_SOURCE_SHA256" \
        > /opt/chemex-numerical-lane/provenance.json \
    && chmod 0444 /opt/chemex-numerical-lane/provenance.json
