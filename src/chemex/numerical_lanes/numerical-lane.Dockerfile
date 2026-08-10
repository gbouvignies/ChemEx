# Build one image per declared CPython source tuple through build-image.sh. The
# final OCI digest is observed outside the image and combined with post-import
# process facts.
ARG PLATFORM_MANIFEST=debian:bookworm-slim@sha256:362e64223cc0da95422b3b13c045186fc0a81250e765d31c025fbddf257f6143
ARG SOURCE_DATE_EPOCH=1785715200
FROM --platform=linux/amd64 ${PLATFORM_MANIFEST}
ARG PLATFORM_MANIFEST
ARG SOURCE_DATE_EPOCH

ARG PYTHON_VERSION
ARG PYTHON_SOURCE_SHA256
ARG UV_VERSION=0.11.15
ARG UV_WHEEL_SHA256=98edf1bdaf82447014852051d93e3ee95012509c567bf057fd117e6bdbd9a807

# Canonical build arguments:
# - 3.13.5 / e6190f52699b534ee203d9f417bdbca05a92f23e35c19c691a50ed2942835385
# - 3.14.5 / 9c22bfe9939a6c5418fc74b289a5f1cc41859ae82ac6b163016b5844bd0a86bc

RUN test -n "${PYTHON_VERSION}" && test -n "${PYTHON_SOURCE_SHA256}"

# Replace every source inherited from the base image. Leaving debian.sources in
# place would let live bookworm/security packages outrank this frozen snapshot.
RUN rm -f /etc/apt/sources.list /etc/apt/sources.list.d/* \
    && printf '%s\n' \
        'Types: deb' \
        'URIs: http://snapshot.debian.org/archive/debian/20260803T000000Z' \
        'Suites: bookworm bookworm-updates' \
        'Components: main' \
        'Check-Valid-Until: no' \
        'Signed-By: /usr/share/keyrings/debian-archive-keyring.gpg' \
        '' \
        'Types: deb' \
        'URIs: http://snapshot.debian.org/archive/debian-security/20260803T000000Z' \
        'Suites: bookworm-security' \
        'Components: main' \
        'Check-Valid-Until: no' \
        'Signed-By: /usr/share/keyrings/debian-archive-keyring.gpg' \
        > /etc/apt/sources.list.d/snapshot.sources \
    && apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential ca-certificates curl libbz2-dev libffi-dev libgdbm-dev \
        liblzma-dev libncursesw5-dev libreadline-dev libsqlite3-dev \
        libssl-dev tk-dev uuid-dev zlib1g-dev \
    && rm -rf /var/cache/* /var/lib/apt/lists/* /var/log/*

RUN curl --fail --location --output /tmp/Python.tgz \
        "https://www.python.org/ftp/python/${PYTHON_VERSION}/Python-${PYTHON_VERSION}.tgz" \
    && echo "${PYTHON_SOURCE_SHA256}  /tmp/Python.tgz" | sha256sum --check --strict \
    && tar --extract --file /tmp/Python.tgz --directory /tmp \
    && cd "/tmp/Python-${PYTHON_VERSION}" \
    && ./configure --prefix=/opt/python \
    && make -j"$(nproc)" \
    && make install \
    && rm -rf /tmp/Python.tgz "/tmp/Python-${PYTHON_VERSION}"

ENV SOURCE_DATE_EPOCH=${SOURCE_DATE_EPOCH} \
    PATH=/opt/python/bin:${PATH} \
    PYTHONPATH=/workspace/src \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    OPENBLAS_CORETYPE=Haswell \
    MKL_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    CHEMEX_NUMERICAL_LANE_WORKERS=1 \
    NPY_DISABLE_CPU_FEATURES=X86_V4,AVX512_ICL,AVX512_SPR

COPY pyproject.toml uv.lock /workspace/
WORKDIR /workspace

# Download the exact runtime wheels selected from uv's locked export, retain
# their hashes as evidence, and install only from that closed wheelhouse.
RUN printf '%s\n' \
        "uv==${UV_VERSION} --hash=sha256:${UV_WHEEL_SHA256}" \
        > /tmp/uv-requirements.txt \
    && /opt/python/bin/python3 -m pip install --no-cache-dir --only-binary=:all: \
        --require-hashes --requirement /tmp/uv-requirements.txt \
    && uv export --locked --no-dev --no-emit-project --format requirements.txt \
        --output-file /tmp/runtime-requirements.txt \
    && install --directory /opt/chemex-numerical-lane/wheels \
    && /opt/python/bin/python3 -m pip download --disable-pip-version-check \
        --only-binary=:all: --require-hashes \
        --requirement /tmp/runtime-requirements.txt \
        --dest /opt/chemex-numerical-lane/wheels \
    && find /opt/chemex-numerical-lane/wheels -type f -name '*.whl' \
        -exec sha256sum {} + \
        | sed 's#  .*/#  #' | LC_ALL=C sort \
        > /opt/chemex-numerical-lane/wheel-manifest.txt \
    && /opt/python/bin/python3 -m venv /workspace/.venv \
    && /workspace/.venv/bin/python -m pip install --disable-pip-version-check \
        --no-index --only-binary=:all: --require-hashes \
        --find-links /opt/chemex-numerical-lane/wheels \
        --requirement /tmp/runtime-requirements.txt \
    && rm -rf /root/.cache /tmp/*

COPY src /workspace/src

# Read the actual x86 floating-point control words rather than treating NumPy's
# error policy as hardware state.
RUN cc -shared -fPIC -O2 \
        src/chemex/numerical_lanes/fpstate.c \
        -o /opt/chemex-numerical-lane/libchemex_fpstate.so

# These manifests are content-addressed build evidence. Absolute paths never
# enter their hashes; the final image digest is deliberately observed outside.
RUN dpkg-query --show --showformat='${Package}\t${Version}\t${Architecture}\n' \
        | LC_ALL=C sort \
        > /opt/chemex-numerical-lane/os-package-manifest.txt \
    && find pyproject.toml uv.lock src -type f \
        ! -path '*/manifests/*.json' ! -name '*.pyc' -print0 \
        | LC_ALL=C sort -z | xargs -0 sha256sum \
        > /opt/chemex-numerical-lane/build-context-manifest.txt \
    && recipe_hash="$(sha256sum \
        src/chemex/numerical_lanes/build-image.sh \
        src/chemex/numerical_lanes/numerical-lane.Dockerfile \
        | sha256sum | cut -d ' ' -f 1)" \
    && lockfile_hash="$(sha256sum uv.lock | cut -d ' ' -f 1)" \
    && context_hash="$(sha256sum /opt/chemex-numerical-lane/build-context-manifest.txt | cut -d ' ' -f 1)" \
    && os_hash="$(sha256sum /opt/chemex-numerical-lane/os-package-manifest.txt | cut -d ' ' -f 1)" \
    && wheel_hash="$(sha256sum /opt/chemex-numerical-lane/wheel-manifest.txt | cut -d ' ' -f 1)" \
    && printf '{"build_context_hash":"%s","build_recipe_hash":"%s","dependency_lock_hash":"%s","os_package_manifest_hash":"%s","platform_manifest":"%s","python_source_hash":"%s","uv_version":"%s","uv_wheel_hash":"%s","wheel_manifest_hash":"%s"}\n' \
        "$context_hash" "$recipe_hash" "$lockfile_hash" "$os_hash" \
        "$PLATFORM_MANIFEST" "$PYTHON_SOURCE_SHA256" "$UV_VERSION" \
        "$UV_WHEEL_SHA256" "$wheel_hash" \
        > /opt/chemex-numerical-lane/provenance.json \
    && chmod 0444 \
        /opt/chemex-numerical-lane/build-context-manifest.txt \
        /opt/chemex-numerical-lane/libchemex_fpstate.so \
        /opt/chemex-numerical-lane/os-package-manifest.txt \
        /opt/chemex-numerical-lane/provenance.json \
        /opt/chemex-numerical-lane/wheel-manifest.txt \
    && chmod -R a-w /opt/chemex-numerical-lane/wheels
