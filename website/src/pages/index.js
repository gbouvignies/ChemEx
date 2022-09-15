import React from 'react';
import clsx from 'clsx';
import Layout from '@theme/Layout';
import Link from '@docusaurus/Link';
import useDocusaurusContext from '@docusaurus/useDocusaurusContext';
import styles from './index.module.css';
import HomepageFeatures from '@site/src/components/HomepageFeatures';
import Logo from '/img/logo-black.svg';

function HomepageHeader() {
    const { siteConfig } = useDocusaurusContext();
    return (
        <header className={clsx('hero hero--primary', styles.heroBanner)}>
            <div className="container">
                {/* Left */}
                <div className={styles.heroLeft}>
                    <div className={styles.imageLogo}>
                        <Logo />
                    </div>
                </div>
                {/* Right */}
                <div className={styles.heroRight}>
                    <h1 className="hero__title">{siteConfig.title}</h1>
                    <p className="hero__subtitle">{siteConfig.tagline}</p>
                    <div className={styles.buttons}>
                        <Link
                            className="button button--secondary button--lg"
                            to="/docs/getting_started">
                            Get started
                        </Link>
                    </div>
                    <span className={styles.indexCtasGitHubButtonWrapper}>
                        <iframe
                            className={styles.indexCtasGitHubButton}
                            src="https://ghbtns.com/github-btn.html?user=gbouvignies&amp;repo=chemex&amp;type=star&amp;count=true&amp;size=large"
                            width={160}
                            height={30}
                            title="GitHub Stars"
                        />
                    </span>
                </div>
            </div>
        </header>
    );
}

export default function Home() {
    const { siteConfig } = useDocusaurusContext();
    return (
        <Layout
            title={`Hello from ${siteConfig.title}`}
            description="Description will go into a meta tag in <head />">
            <HomepageHeader />
            <main>
                <HomepageFeatures />
            </main>
        </Layout>
    );
}
