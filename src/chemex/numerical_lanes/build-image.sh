#!/bin/sh
# Canonical #588 OCI build invocation. BuildKit rewrites every layer timestamp
# to SOURCE_DATE_EPOCH; default provenance attestations are excluded because
# their invocation metadata is not part of the numerical lane semantics.
# The selected builder must use:
# moby/buildkit:v0.26.2@sha256:de10faf919fc71ba4eb1dd7bd6449566d012b0c9436b1c61bfee21d621b009aa
set -eu

if [ "$#" -ne 3 ]; then
    echo "usage: build-image.sh IMAGE_TAG PYTHON_VERSION PYTHON_SOURCE_SHA256" >&2
    exit 2
fi

image_tag=$1
python_version=$2
python_source_sha256=$3
buildx_version=$(docker buildx version)
case "$buildx_version" in
    "github.com/docker/buildx v0.36.1 "*) ;;
    *)
        echo "canonical lane build requires Docker Buildx v0.36.1" >&2
        exit 1
        ;;
esac
buildkit_versions=$(docker buildx inspect --bootstrap \
    | awk '/BuildKit version:/ {print $3}' | LC_ALL=C sort -u)
if [ "$buildkit_versions" != "v0.26.2" ]; then
    echo "canonical lane build requires BuildKit v0.26.2" >&2
    exit 1
fi
archive=$(mktemp "${TMPDIR:-/tmp}/chemex-numerical-lane.XXXXXX.oci.tar")
trap 'rm -f "$archive"' EXIT HUP INT TERM

docker buildx build \
    --platform linux/amd64 \
    --provenance=false \
    --output "type=oci,dest=$archive,rewrite-timestamp=true" \
    --tag "$image_tag" \
    --build-arg SOURCE_DATE_EPOCH=1785715200 \
    --build-arg "PYTHON_VERSION=$python_version" \
    --build-arg "PYTHON_SOURCE_SHA256=$python_source_sha256" \
    --file src/chemex/numerical_lanes/numerical-lane.Dockerfile \
    .
docker load --input "$archive"
