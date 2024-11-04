import React from 'react';
import clsx from 'clsx';
import styles from './styles.module.css';

const FeatureList = [
    {
        title: 'Accurate',
        Svg: require('@site/static/img/pulse_sequence.svg').default,
        description: (
            <>
                ChemEx accurately simulates spin evolution across the entire pulse sequence,
                enabling precise modeling of experimental factors like pulse imperfections,
                off-resonance effects, and phase cycling.
            </>
        ),
    },
    {
        title: 'Versatile',
        Svg: require('@site/static/img/combined_experiments.svg').default,
        description: (
            <>
                Supporting diverse experiments and kinetic models, ChemEx allows joint analysis for high-precision exchange parameter extraction, making it essential for complex kinetics studies.
            </>
        ),
    },
    {
        title: 'Open Source',
        Svg: require('@site/static/img/open-source-logos.svg').default,
        description: (
            <>
                ChemEx is a GPLv3-licensed open-source Python application,
                built on established packages like Numpy, Scipy, Matplotlib, Pydantic, and Rich.
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
