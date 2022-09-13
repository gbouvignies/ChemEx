import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
    {
        title: 'Accurate',
        Svg: require('@site/static/img/pulse_sequence.svg').default,
        description: (
            <>
                ChemEx simulates the spin evolution over the pulse sequence,
                which makes it easy to include experimental details, such as pulse
                imperfections, off-resonance effects, or phase cycling,
                into the analysis.
            </>
        ),
    },
    {
        title: 'Versatile',
        Svg: require('@site/static/img/combined_experiments.svg').default,
        description: (
            <>
                Choose from a large set of experiments and kinetic models. ChemEx
                lets you analyze multiple experiments jointly to extract the
                exchange parameters.
            </>
        ),
    },
    {
        title: 'Open Source',
        Svg: require('@site/static/img/open-source-logos.svg').default,
        description: (
            <>
                ChemEx is an open source python (BSD licensed) application that
                builds upon well established packages, including Numpy, Scipy,
                Matplotlib, Pydantic, Rich.
            </>
        ),
    },
];

function Feature({ Svg, title, description }) {
    return (
        <div className={clsx('col col--4')}>
            <div className="text--center">
                <Svg className={styles.featureSvg} role="img" />
            </div>
            <div className="text--center padding-horiz--md">
                <h3>{title}</h3>
                <p>{description}</p>
            </div>
        </div>
    );
}

export default function HomepageFeatures() {
    return (
        <section className={styles.features}>
            <div className="container">
                <div className="row">
                    {FeatureList.map((props, idx) => (
                        <Feature key={idx} {...props} />
                    ))}
                </div>
            </div>
        </section>
    );
}
