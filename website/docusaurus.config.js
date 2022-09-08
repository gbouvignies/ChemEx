// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');

const math = require('remark-math');
const katex = require('rehype-katex');

/** @type {import('@docusaurus/types').Config} */
const config = {
    title: 'ChemEx',
    tagline: 'ChemEx is an analysis program for characterizing chemical exchange detected by NMR.',
    url: 'https://gbouvignies.github.io',
    baseUrl: '/chemex/',
    onBrokenLinks: 'throw',
    onBrokenMarkdownLinks: 'warn',
    favicon: 'img/favicon.ico',
    organizationName: 'gbouvignies',
    projectName: 'chemex',
    deploymentBranch: 'gh-pages',

    presets: [
        [
            'classic',
            /** @type {import('@docusaurus/preset-classic').Options} */
            ({
                docs: {
                    sidebarPath: require.resolve('./sidebars.js'),
                    editUrl: 'https://github.com/gbouvignies/chemex/tree/master/',
                    remarkPlugins: [math],
                    rehypePlugins: [katex],
                },
                theme: {
                    customCss: require.resolve('./src/css/custom.css'),
                },
            }),
        ],
    ],

    stylesheets: [
        {
            href: 'https://cdn.jsdelivr.net/npm/katex@0.13.24/dist/katex.min.css',
            type: 'text/css',
            integrity:
                'sha384-odtC+0UGzzFL/6PNoE8rX/SPcQDXBJ+uRepguP4QkPCm2LBxH3FA3y+fKSiJ+AmM',
            crossorigin: 'anonymous',
        },
    ],

    themeConfig:
        /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
        ({
            navbar: {
                title: 'ChemEx',
                logo: {
                    alt: 'ChemEx Logo',
                    src: 'img/logo-black.svg',
                    srcDark: 'img/logo-white.svg',
                },
                items: [
                    {
                        type: 'doc',
                        docId: 'getting_started',
                        position: 'left',
                        label: 'Docs',
                    },
                    {
                        href: 'https://github.com/gbouvignies/chemex',
                        position: 'right',
                        className: 'header-github-link',
                        'aria-label': 'GitHub repository',
                    },
                ],
            },
            docs: {
                sidebar: {
                    autoCollapseCategories: true,
                },
            },
            footer: {
                style: 'light',
                links: [
                    {
                        title: 'Docs',
                        items: [
                            {
                                label: 'Installation',
                                to: '/docs/getting_started#Installation',
                            },
                            {
                                label: 'User Guide',
                                to: '/docs/user_guide',
                            },
                            {
                                label: 'Experiments',
                                to: '/docs/experiments',
                            },
                            {
                                label: 'Examples',
                                to: '/docs/examples',
                            },
                        ],
                    },
                    {
                        title: 'Community',
                        items: [
                            {
                                label: 'Discussion',
                                href: 'https://github.com/gbouvignies/ChemEx/discussions',
                            },
                            {
                                label: 'Report Issues',
                                href: 'https://github.com/gbouvignies/ChemEx/issues',
                            },
                        ],
                    },
                    {
                        title: 'More',
                        items: [
                            {
                                label: 'GitHub',
                                href: 'https://github.com/gbouvignies/chemex',
                            },
                            {
                                label: 'Twitter',
                                href: 'https://twitter.com/gubouvignies',
                            },
                        ],
                    },
                ],
                copyright: `Copyright © ${new Date().getFullYear()} Guillaume Bouvignies, Inc. Built with Docusaurus.`,
            },
            prism: {
                theme: lightCodeTheme,
                darkTheme: darkCodeTheme,
                additionalLanguages: ['toml', "markdown"],
            },
        }),
};

module.exports = config;
