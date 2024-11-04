// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

import { themes as prismThemes } from 'prism-react-renderer';

const math = require('remark-math');
const katex = require('rehype-katex');

/** @type {import('@docusaurus/types').Config} */
const config = {
    title: 'ChemEx',
    tagline: 'ChemEx is an analysis program for characterizing chemical exchange detected by NMR.',
    url: 'https://gbouvignies.github.io',
    baseUrl: '/ChemEx/',
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
                    sidebarPath: './sidebars.js',
                    editUrl: 'https://github.com/gbouvignies/chemex/tree/main/website/',
                    remarkPlugins: [math],
                    rehypePlugins: [katex],
                },
                theme: {
                    customCss: './src/css/custom.css',
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
                        docId: 'welcome_to_chemex',
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
                                to: '/docs/welcome_to_chemex#installation',
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
                theme: prismThemes.github,
                darkTheme: prismThemes.dracula,
                additionalLanguages: ['toml', "markdown"],
            },
            algolia: {
                // The application ID provided by Algolia
                appId: 'THPPA40I1Q',

                // Public API key: it is safe to commit it
                apiKey: '2365f685942b90aec2a3c2dff76774e4',

                indexName: 'chemex',

                // Optional: see doc section below
                contextualSearch: true,

                // Optional: Specify domains where the navigation should occur through window.location instead on history.push. Useful when our Algolia config crawls multiple documentation sites and we want to navigate with window.location.href to them.
                // externalUrlRegex: 'external\\.com|domain\\.com',

                // Optional: Algolia search parameters
                searchParameters: {},

                // Optional: path for search page that enabled by default (`false` to disable it)
                searchPagePath: 'search',

                //... other Algolia params
            },
        }),
    future: {
        experimental_faster: true,
    },
};

export default config;
