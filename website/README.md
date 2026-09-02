# Website

The ChemEx website is built with [Docusaurus](https://docusaurus.io/).

## Local development

Install the locked dependencies and start the development server:

```sh
npm ci
npm run start
```

Build the production site with:

```sh
npm run build
```

## Documentation versions

Docusaurus's standard documentation-versioning layout is used:

- `docs/` contains documentation for current development. Once at least one
  version has been frozen, Docusaurus labels this version `Next` and serves it
  below `/docs/next/`.
- `versioned_docs/` and `versioned_sidebars/` contain immutable documentation
  snapshots for materially distinct released compatibility surfaces.
- `versions.json` lists frozen documentation versions from newest to oldest.
  Docusaurus serves the newest frozen version at `/docs/` and older versions
  below `/docs/<version>/`.

Normal documentation work belongs in `docs/`. Do not synchronize later edits
back into frozen versions. A critical correction to an old version should be an
explicit, reviewed change.

Normal local and pull-request builds include `Next`. The production release
build sets `CHEMEX_DOCS_INCLUDE_CURRENT_VERSION=false`, so the canonical site
exposes only frozen stable documentation versions.

## Freezing documentation for a release

Documentation versions use ChemEx release identifiers verbatim, without a `v`
prefix. During release preparation, after the documentation for a materially
different user-facing interface is ready, run:

```sh
npm ci
npm run docusaurus docs:version YYYY.MM.MICRO
npm run build
```

Review and commit the generated `versioned_docs/`, `versioned_sidebars/`, and
`versions.json` changes before creating the Git tag and GitHub Release. These
files must therefore be present in the exact tagged source.

Do not create a documentation version mechanically for every software release.
A new version is appropriate when configuration or method syntax, supported
options, CLI behavior, documented APIs, workflows, experiments, models, or
examples have materially changed. A bug-fix release with the same documented
compatibility surface should continue to use the existing latest documentation
version.

## Deployment

Production GitHub Pages is release-bound. Publishing a stable GitHub Release
causes `.github/workflows/website-deploy.yml` to check out that exact tag, build
its committed website, and deploy the result. The deployment workflow consumes
immutable tagged source; it must not generate versions, commit changes, or
deploy development `main` directly.
