"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[5734],{3901:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>h,contentTitle:()=>d,default:()=>f,frontMatter:()=>c,metadata:()=>u,toc:()=>m});var s=t(4848),a=t(8453),i=t(1470),r=t(9365);const l=t.p+"assets/images/cest_26hz_fit-cd540bd92fd0770c065e108770e4d2e1.png",o=t.p+"assets/images/cpmg_800mhz_fit-ec19160c7c18453507ae0aa05c4e295b.png",c={sidebar_position:8},d="Outputs",u={id:"user_guide/fitting/outputs",title:"Outputs",description:"Location",source:"@site/docs/user_guide/fitting/outputs.mdx",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/outputs",permalink:"/ChemEx/docs/user_guide/fitting/outputs",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/outputs.mdx",tags:[],version:"current",sidebarPosition:8,frontMatter:{sidebar_position:8},sidebar:"tutorialSidebar",previous:{title:"Kinetic models",permalink:"/ChemEx/docs/user_guide/fitting/kinetic_models"},next:{title:"Additional modules",permalink:"/ChemEx/docs/user_guide/additional_modules"}},h={},m=[{value:"Location",id:"location",level:2},{value:"Multi-step fit",id:"multi-step-fit",level:3},{value:"Multi-group fit",id:"multi-group-fit",level:3},{value:"Content",id:"content",level:2},{value:"<code>Parameters/</code>",id:"parameters",level:3},{value:"Example files",id:"example-files",level:4},{value:"<code>Plot/</code>",id:"plot",level:3},{value:"<code>Data/</code>",id:"data",level:3},{value:"<code>statistics.toml</code>",id:"statisticstoml",level:3}];function p(e){const n={a:"a",admonition:"admonition",annotation:"annotation",code:"code",h1:"h1",h2:"h2",h3:"h3",h4:"h4",math:"math",mi:"mi",mn:"mn",mrow:"mrow",msup:"msup",p:"p",pre:"pre",semantics:"semantics",span:"span",...(0,a.R)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(n.h1,{id:"outputs",children:"Outputs"}),"\n",(0,s.jsx)(n.h2,{id:"location",children:"Location"}),"\n",(0,s.jsxs)(n.p,{children:["The results of the fits are written in the directory ",(0,s.jsx)(n.code,{children:"Output/"})," by default.\nHowever, you can change this location using the command-line option ",(0,s.jsx)(n.code,{children:"--output"}),"\n(or ",(0,s.jsx)(n.code,{children:"-o"}),") followed by the desired path name."]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-bash",children:"chemex fit -o <path> [other options]\n"})}),"\n",(0,s.jsx)(n.h3,{id:"multi-step-fit",children:"Multi-step fit"}),"\n",(0,s.jsx)(n.p,{children:"If the fitting method has multiple fitting steps, each step will create its own\noutput subdirectory with the name of the fitting step."}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-shell",children:"Output/\n\u251c\u2500\u2500 STEP1/\n\u2514\u2500\u2500 STEP2/\n"})}),"\n",(0,s.jsx)(n.h3,{id:"multi-group-fit",children:"Multi-group fit"}),"\n",(0,s.jsxs)(n.p,{children:["During any fitting step, if the dataset can be divided in multiple groups\ndepending on distinct sets of fitting parameters, ChemEx will fit each group of\ndata separately. The results of these individual fits are then stored in\nseparate subfolders placed in the the ",(0,s.jsx)(n.code,{children:"Groups/"})," folder. Subfolders are named\nwith a number and an ID, which depends on the parameters that have been\noptimized for the specific group. The results are also put together in the\n",(0,s.jsx)(n.code,{children:"All/"})," folder for convenience."]}),"\n",(0,s.jsxs)(n.p,{children:["For example, if all global parameters (p",(0,s.jsx)("sub",{children:"B"}),", k",(0,s.jsx)("sub",{children:"ex"}),", etc.) are\nfixed or residue-specific fitting model is used (e.g. the model\n",(0,s.jsx)(n.a,{href:"/ChemEx/docs/user_guide/fitting/kinetic_models",children:(0,s.jsx)(n.code,{children:"2st_rs"})}),"), then all parameters are allowed to vary are\nresidue-specific. Residue-specific fits will then be run and two separate\nsubdirectories ",(0,s.jsx)(n.code,{children:"All/"})," and ",(0,s.jsx)(n.code,{children:"Groups/"})," will be created, which contain fitting\nresults for all residues and each individual residue, respectively."]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-shell",children:"Output/\n\u2514\u2500\u2500 STEP2/\n    \u251c\u2500\u2500 All/\n    \u2514\u2500\u2500 Groups/\n        \u251c\u2500\u2500 10_11N/\n        \u251c\u2500\u2500 11_12N/\n        \u251c\u2500\u2500 12_13N/\n"})}),"\n",(0,s.jsx)(n.h2,{id:"content",children:"Content"}),"\n",(0,s.jsx)(n.p,{children:"The fitting output typically contains the following files and directories:"}),"\n",(0,s.jsx)(n.h3,{id:"parameters",children:(0,s.jsx)(n.code,{children:"Parameters/"})}),"\n",(0,s.jsxs)(n.p,{children:["Contains fitting results as three separate files ",(0,s.jsx)(n.code,{children:"fitted.toml"}),", ",(0,s.jsx)(n.code,{children:"fixed.toml"})," and\n",(0,s.jsx)(n.code,{children:"constrained.toml"}),", which contain output parameters that are fitted, fixed and\nconstrained during the fitting process, respectively."]}),"\n",(0,s.jsx)(n.h4,{id:"example-files",children:"Example files"}),"\n","\n",(0,s.jsxs)(i.A,{children:[(0,s.jsxs)(r.A,{value:"fitted",label:"fitted.toml",default:!0,children:[(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:'[GLOBAL]\nKEX_AB =  3.81511e+02 # \xb18.90870e+00\nPB     =  7.02971e-02 # \xb11.14784e-03\n\n[DW_AB]\n15N =  2.00075e+00 # \xb12.30817e-02\n31N =  1.98968e+00 # \xb11.90842e-02\n33N =  1.82003e+00 # \xb12.44821e-02\n34N =  3.63170e+00 # \xb13.62801e-02\n37N =  1.69183e+00 # \xb12.41070e-02\n\n["R2_A, B0->500.0MHZ"]\n15N =  3.98674e+00 # \xb12.55793e-01\n31N =  5.85923e+00 # \xb11.94734e-01\n33N =  4.02099e+00 # \xb12.78003e-01\n34N =  4.16615e+00 # \xb12.05190e-01\n37N =  3.67705e+00 # \xb13.04621e-01\n\n["R2_A, B0->800.0MHZ"]\n15N =  6.33712e+00 # \xb13.86894e-01\n31N =  7.99927e+00 # \xb12.81972e-01\n33N =  6.23967e+00 # \xb15.82366e-01\n34N =  5.99535e+00 # \xb13.17832e-01\n37N =  5.37295e+00 # \xb15.69929e-01\n'})}),(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsx)(n.p,{children:'If uncertainties on fitted parameters have been calculated \u2013 this depends on the\nfitting algorithm used \u2013, they are reported as comments at the end of the line\nand are preceded by the sign "\xb1". Uncertainties on fitted parameters are,\ntypically, estimated through the covariance matrix obtained from the\nLevenberg-Marquardt optimization.'})})]}),(0,s.jsx)(r.A,{value:"fixed",label:"fixed.toml",children:(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:'[CS_A]\n15N =  1.19849e+02 # (fixed)\n31N =  1.26388e+02 # (fixed)\n33N =  1.18762e+02 # (fixed)\n34N =  1.14897e+02 # (fixed)\n37N =  1.21618e+02 # (fixed)\n\n["R1_A, B0->500.0MHZ"]\n15N =  1.50000e+00 # (fixed)\n31N =  1.50000e+00 # (fixed)\n33N =  1.50000e+00 # (fixed)\n34N =  1.50000e+00 # (fixed)\n37N =  1.50000e+00 # (fixed)\n\n["R1_A, B0->800.0MHZ"]\n15N =  1.50000e+00 # (fixed)\n31N =  1.50000e+00 # (fixed)\n33N =  1.50000e+00 # (fixed)\n34N =  1.50000e+00 # (fixed)\n37N =  1.50000e+00 # (fixed)\n'})})}),(0,s.jsxs)(r.A,{value:"constrained",label:"constrained.toml",children:[(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",children:'[GLOBAL]\nKAB =  2.68192e+01 # \xb13.06068e-01 ([KEX_AB] * [PB])\nKBA =  3.54692e+02 # \xb18.28245e+00 ([KEX_AB] * [PA])\nPA  =  9.29703e-01 # \xb11.14784e-03 (1.0 - [PB])\n\n[CS_B]\n15N =  1.21850e+02 # \xb12.30817e-02 ([CS_A, NUC->15N] + [DW_AB, NUC->15N])\n31N =  1.28378e+02 # \xb11.90842e-02 ([CS_A, NUC->31N] + [DW_AB, NUC->31N])\n33N =  1.20582e+02 # \xb12.44821e-02 ([CS_A, NUC->33N] + [DW_AB, NUC->33N])\n34N =  1.18529e+02 # \xb13.62801e-02 ([CS_A, NUC->34N] + [DW_AB, NUC->34N])\n37N =  1.23310e+02 # \xb12.41070e-02 ([CS_A, NUC->37N] + [DW_AB, NUC->37N])\n\n["R1_B, B0->500.0MHZ"]\n15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->500.0MHZ])\n31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->500.0MHZ])\n33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->500.0MHZ])\n34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->500.0MHZ])\n37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->500.0MHZ])\n\n["R1_B, B0->800.0MHZ"]\n15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->800.0MHZ])\n31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->800.0MHZ])\n33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->800.0MHZ])\n34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->800.0MHZ])\n37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->800.0MHZ])\n\n["R2_B, B0->500.0MHZ"]\n15N =  3.98674e+00 # \xb12.55793e-01 ([R2_A, NUC->15N, B0->500.0MHZ])\n31N =  5.85923e+00 # \xb11.94734e-01 ([R2_A, NUC->31N, B0->500.0MHZ])\n33N =  4.02099e+00 # \xb12.78003e-01 ([R2_A, NUC->33N, B0->500.0MHZ])\n34N =  4.16615e+00 # \xb12.05190e-01 ([R2_A, NUC->34N, B0->500.0MHZ])\n37N =  3.67705e+00 # \xb13.04621e-01 ([R2_A, NUC->37N, B0->500.0MHZ])\n\n["R2_B, B0->800.0MHZ"]\n15N =  6.33712e+00 # \xb13.86894e-01 ([R2_A, NUC->15N, B0->800.0MHZ])\n31N =  7.99927e+00 # \xb12.81972e-01 ([R2_A, NUC->31N, B0->800.0MHZ])\n33N =  6.23967e+00 # \xb15.82366e-01 ([R2_A, NUC->33N, B0->800.0MHZ])\n34N =  5.99535e+00 # \xb13.17832e-01 ([R2_A, NUC->34N, B0->800.0MHZ])\n37N =  5.37295e+00 # \xb15.69929e-01 ([R2_A, NUC->37N, B0->800.0MHZ])\n'})}),(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsx)(n.p,{children:"Propagated uncertainties \u2013 when available \u2013 and applied constrained are reported\nat the end of the line as comments."})})]})]}),"\n",(0,s.jsx)(n.h3,{id:"plot",children:(0,s.jsx)(n.code,{children:"Plot/"})}),"\n",(0,s.jsxs)(n.p,{children:["Contains fitting results as plots (in ",(0,s.jsx)(n.code,{children:".pdf"})," format) and also the raw datasets\n(both the original input and fitted data points) for creating the plots. Example\nfitting results for CPMG and CEST experiments are shown below:"]}),"\n","\n",(0,s.jsxs)("figure",{children:[(0,s.jsx)("img",{src:l,alt:"CEST profile",width:"50%"}),(0,s.jsx)("img",{src:o,alt:"CPMG profile",width:"50%"}),(0,s.jsx)("figcaption",{align:"center",children:(0,s.jsx)("b",{children:"Examples of CEST and CPMG fitting results"})})]}),"\n",(0,s.jsx)(n.admonition,{type:"note",children:(0,s.jsx)(n.p,{children:'In plots of (D-/cos-)CEST fitting results, the positions for ground and excited\nstates are indicated by solid and dashed vertical lines respectively. Besides,\ndata points that are filtered out from the fit are shown in lighter color. In\nplots of D-CEST/COS-CEST fitting results, the "folded" positions for ground and\nexcited states are indicated by "*" at the vertical lines.'})}),"\n",(0,s.jsx)(n.h3,{id:"data",children:(0,s.jsx)(n.code,{children:"Data/"})}),"\n",(0,s.jsxs)(n.p,{children:["Contains all the data values used for the fitting along with the back-calculated\nvalues. These files can be used to calculate the ",(0,s.jsxs)(n.span,{className:"katex",children:[(0,s.jsx)(n.span,{className:"katex-mathml",children:(0,s.jsx)(n.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(n.semantics,{children:[(0,s.jsx)(n.mrow,{children:(0,s.jsxs)(n.msup,{children:[(0,s.jsx)(n.mi,{children:"\u03c7"}),(0,s.jsx)(n.mn,{children:"2"})]})}),(0,s.jsx)(n.annotation,{encoding:"application/x-tex",children:"\u03c7^2"})]})})}),(0,s.jsx)(n.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(n.span,{className:"base",children:[(0,s.jsx)(n.span,{className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,s.jsxs)(n.span,{className:"mord",children:[(0,s.jsx)(n.span,{className:"mord mathnormal",children:"\u03c7"}),(0,s.jsx)(n.span,{className:"msupsub",children:(0,s.jsx)(n.span,{className:"vlist-t",children:(0,s.jsx)(n.span,{className:"vlist-r",children:(0,s.jsx)(n.span,{className:"vlist",style:{height:"0.8141em"},children:(0,s.jsxs)(n.span,{style:{top:"-3.063em",marginRight:"0.05em"},children:[(0,s.jsx)(n.span,{className:"pstrut",style:{height:"2.7em"}}),(0,s.jsx)(n.span,{className:"sizing reset-size6 size3 mtight",children:(0,s.jsx)(n.span,{className:"mord mtight",children:"2"})})]})})})})})]})]})})]})," value."]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",metastring:"title=Data/500mhz.dat",children:"[15N]\n#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)\n             0    3.47059800e+04    1.77491406e+02    3.47055362e+04\n            30    3.05930380e+04    1.77491406e+02    3.06111963e+04\n             1    1.81234230e+04    1.77491406e+02    1.80856300e+04\n            28    3.07144730e+04    1.77491406e+02    3.05838767e+04\n             2    2.02155120e+04    1.77491406e+02    2.02110184e+04\n            26    3.05222020e+04    1.77491406e+02    3.05501255e+04\n             3    2.23056070e+04    1.77491406e+02    2.23601656e+04\n            24    3.05381830e+04    1.77491406e+02    3.05077686e+04\n             4    2.43783050e+04    1.77491406e+02    2.44111932e+04\n            22    3.06981570e+04    1.77491406e+02    3.04536377e+04\n             5    2.58673980e+04    1.77491406e+02    2.59884541e+04\n            20    3.06069180e+04    1.77491406e+02    3.03829761e+04\n             6    2.70660870e+04    1.77491406e+02    2.71148924e+04\n            18    3.00982700e+04    1.77491406e+02    3.02883916e+04\n             7    2.81512990e+04    1.77491406e+02    2.79178594e+04\n            16    3.02515700e+04    1.77491406e+02    3.01579211e+04\n             8    2.84045570e+04    1.77491406e+02    2.84985729e+04\n            14    2.97528530e+04    1.77491406e+02    2.99712519e+04\n             9    2.84536650e+04    1.77491406e+02    2.89264158e+04\n            13    2.98219180e+04    1.77491406e+02    2.98461022e+04\n            10    2.96936410e+04    1.77491406e+02    2.92495591e+04\n            12    2.95269900e+04    1.77491406e+02    2.96918706e+04\n            11    2.92540210e+04    1.77491406e+02    2.94972506e+04\n             2    2.04476470e+04    1.77491406e+02    2.02110184e+04\n            28    3.05466550e+04    1.77491406e+02    3.05838767e+04\n             8    2.86183900e+04    1.77491406e+02    2.84985729e+04\n\n[31N]\n#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)\n             0    4.71577550e+04    1.77491406e+02    4.71537007e+04\n            30    4.00023250e+04    1.77491406e+02    3.94823953e+04\n             1    2.30167260e+04    1.77491406e+02    2.34136180e+04\n            28    4.00615990e+04    1.77491406e+02    3.94547546e+04\n             2    2.61934530e+04    1.77491406e+02    2.61653780e+04\n            26    3.96898190e+04    1.77491406e+02    3.94190638e+04\n             3    2.86882570e+04    1.77491406e+02    2.88729551e+04\n            24    3.98238690e+04    1.77491406e+02    3.93722573e+04\n             4    3.21031440e+04    1.77491406e+02    3.13942278e+04\n            22    3.96733310e+04    1.77491406e+02    3.93097577e+04\n             5    3.35871430e+04    1.77491406e+02    3.34628362e+04\n            20    3.86957600e+04    1.77491406e+02    3.92245398e+04\n             6    3.53600830e+04    1.77491406e+02    3.49641428e+04\n            18    3.88115960e+04    1.77491406e+02    3.91054958e+04\n             7    3.58938150e+04    1.77491406e+02    3.59453347e+04\n            16    3.85843960e+04    1.77491406e+02    3.89345050e+04\n             8    3.61450060e+04    1.77491406e+02    3.66900206e+04\n            14    3.85711880e+04    1.77491406e+02    3.86811094e+04\n             9    3.72815770e+04    1.77491406e+02    3.72429094e+04\n            13    3.82358560e+04    1.77491406e+02    3.85092366e+04\n            10    3.76426640e+04    1.77491406e+02    3.76803042e+04\n            12    3.78470460e+04    1.77491406e+02    3.82929922e+04\n            11    3.76000060e+04    1.77491406e+02    3.80220076e+04\n             2    2.62301400e+04    1.77491406e+02    2.61653780e+04\n            28    3.97731200e+04    1.77491406e+02    3.94547546e+04\n             8    3.63510140e+04    1.77491406e+02    3.66900206e+04\n\n[...]\n"})}),"\n",(0,s.jsx)(n.h3,{id:"statisticstoml",children:(0,s.jsx)(n.code,{children:"statistics.toml"})}),"\n",(0,s.jsxs)(n.p,{children:["Contains all goodness-of-fit statistics, such as ",(0,s.jsxs)(n.span,{className:"katex",children:[(0,s.jsx)(n.span,{className:"katex-mathml",children:(0,s.jsx)(n.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,s.jsxs)(n.semantics,{children:[(0,s.jsx)(n.mrow,{children:(0,s.jsxs)(n.msup,{children:[(0,s.jsx)(n.mi,{children:"\u03c7"}),(0,s.jsx)(n.mn,{children:"2"})]})}),(0,s.jsx)(n.annotation,{encoding:"application/x-tex",children:"\u03c7^2"})]})})}),(0,s.jsx)(n.span,{className:"katex-html","aria-hidden":"true",children:(0,s.jsxs)(n.span,{className:"base",children:[(0,s.jsx)(n.span,{className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,s.jsxs)(n.span,{className:"mord",children:[(0,s.jsx)(n.span,{className:"mord mathnormal",children:"\u03c7"}),(0,s.jsx)(n.span,{className:"msupsub",children:(0,s.jsx)(n.span,{className:"vlist-t",children:(0,s.jsx)(n.span,{className:"vlist-r",children:(0,s.jsx)(n.span,{className:"vlist",style:{height:"0.8141em"},children:(0,s.jsxs)(n.span,{style:{top:"-3.063em",marginRight:"0.05em"},children:[(0,s.jsx)(n.span,{className:"pstrut",style:{height:"2.7em"}}),(0,s.jsx)(n.span,{className:"sizing reset-size6 size3 mtight",children:(0,s.jsx)(n.span,{className:"mord mtight",children:"2"})})]})})})})})]})]})})]}),"."]}),"\n",(0,s.jsx)(n.pre,{children:(0,s.jsx)(n.code,{className:"language-toml",metastring:"title='\"statistics.toml\"'",children:'"number of data points"                = 230\n"number of variables"                  = 17\n"chi-square"                           =  4.34824e+02\n"reduced-chi-square"                   =  2.04143e+00\n"chi-squared test"                     =  0.00000e+00\n"Kolmogorov-Smirnov test"              =  6.72236e-02\n"Akaike Information Criterion (AIC)"   =  1.80479e+02\n"Bayesian Information Criterion (BIC)" =  2.38926e+02\n\n'})})]})}function f(e={}){const{wrapper:n}={...(0,a.R)(),...e.components};return n?(0,s.jsx)(n,{...e,children:(0,s.jsx)(p,{...e})}):p(e)}},9365:(e,n,t)=>{t.d(n,{A:()=>r});t(6540);var s=t(4164);const a={tabItem:"tabItem_Ymn6"};var i=t(4848);function r(e){let{children:n,hidden:t,className:r}=e;return(0,i.jsx)("div",{role:"tabpanel",className:(0,s.A)(a.tabItem,r),hidden:t,children:n})}},1470:(e,n,t)=>{t.d(n,{A:()=>C});var s=t(6540),a=t(4164),i=t(3104),r=t(6347),l=t(205),o=t(7485),c=t(1682),d=t(9466);function u(e){return s.Children.toArray(e).filter((e=>"\n"!==e)).map((e=>{if(!e||(0,s.isValidElement)(e)&&function(e){const{props:n}=e;return!!n&&"object"==typeof n&&"value"in n}(e))return e;throw new Error(`Docusaurus error: Bad <Tabs> child <${"string"==typeof e.type?e.type:e.type.name}>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.`)}))?.filter(Boolean)??[]}function h(e){const{values:n,children:t}=e;return(0,s.useMemo)((()=>{const e=n??function(e){return u(e).map((e=>{let{props:{value:n,label:t,attributes:s,default:a}}=e;return{value:n,label:t,attributes:s,default:a}}))}(t);return function(e){const n=(0,c.X)(e,((e,n)=>e.value===n.value));if(n.length>0)throw new Error(`Docusaurus error: Duplicate values "${n.map((e=>e.value)).join(", ")}" found in <Tabs>. Every value needs to be unique.`)}(e),e}),[n,t])}function m(e){let{value:n,tabValues:t}=e;return t.some((e=>e.value===n))}function p(e){let{queryString:n=!1,groupId:t}=e;const a=(0,r.W6)(),i=function(e){let{queryString:n=!1,groupId:t}=e;if("string"==typeof n)return n;if(!1===n)return null;if(!0===n&&!t)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return t??null}({queryString:n,groupId:t});return[(0,o.aZ)(i),(0,s.useCallback)((e=>{if(!i)return;const n=new URLSearchParams(a.location.search);n.set(i,e),a.replace({...a.location,search:n.toString()})}),[i,a])]}function f(e){const{defaultValue:n,queryString:t=!1,groupId:a}=e,i=h(e),[r,o]=(0,s.useState)((()=>function(e){let{defaultValue:n,tabValues:t}=e;if(0===t.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(n){if(!m({value:n,tabValues:t}))throw new Error(`Docusaurus error: The <Tabs> has a defaultValue "${n}" but none of its children has the corresponding value. Available values are: ${t.map((e=>e.value)).join(", ")}. If you intend to show no default tab, use defaultValue={null} instead.`);return n}const s=t.find((e=>e.default))??t[0];if(!s)throw new Error("Unexpected error: 0 tabValues");return s.value}({defaultValue:n,tabValues:i}))),[c,u]=p({queryString:t,groupId:a}),[f,x]=function(e){let{groupId:n}=e;const t=function(e){return e?`docusaurus.tab.${e}`:null}(n),[a,i]=(0,d.Dv)(t);return[a,(0,s.useCallback)((e=>{t&&i.set(e)}),[t,i])]}({groupId:a}),N=(()=>{const e=c??f;return m({value:e,tabValues:i})?e:null})();(0,l.A)((()=>{N&&o(N)}),[N]);return{selectedValue:r,selectValue:(0,s.useCallback)((e=>{if(!m({value:e,tabValues:i}))throw new Error(`Can't select invalid tab value=${e}`);o(e),u(e),x(e)}),[u,x,i]),tabValues:i}}var x=t(2303);const N={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};var g=t(4848);function j(e){let{className:n,block:t,selectedValue:s,selectValue:r,tabValues:l}=e;const o=[],{blockElementScrollPositionUntilNextRender:c}=(0,i.a_)(),d=e=>{const n=e.currentTarget,t=o.indexOf(n),a=l[t].value;a!==s&&(c(n),r(a))},u=e=>{let n=null;switch(e.key){case"Enter":d(e);break;case"ArrowRight":{const t=o.indexOf(e.currentTarget)+1;n=o[t]??o[0];break}case"ArrowLeft":{const t=o.indexOf(e.currentTarget)-1;n=o[t]??o[o.length-1];break}}n?.focus()};return(0,g.jsx)("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,a.A)("tabs",{"tabs--block":t},n),children:l.map((e=>{let{value:n,label:t,attributes:i}=e;return(0,g.jsx)("li",{role:"tab",tabIndex:s===n?0:-1,"aria-selected":s===n,ref:e=>o.push(e),onKeyDown:u,onClick:d,...i,className:(0,a.A)("tabs__item",N.tabItem,i?.className,{"tabs__item--active":s===n}),children:t??n},n)}))})}function b(e){let{lazy:n,children:t,selectedValue:a}=e;const i=(Array.isArray(t)?t:[t]).filter(Boolean);if(n){const e=i.find((e=>e.props.value===a));return e?(0,s.cloneElement)(e,{className:"margin-top--md"}):null}return(0,g.jsx)("div",{className:"margin-top--md",children:i.map(((e,n)=>(0,s.cloneElement)(e,{key:n,hidden:e.props.value!==a})))})}function v(e){const n=f(e);return(0,g.jsxs)("div",{className:(0,a.A)("tabs-container",N.tabList),children:[(0,g.jsx)(j,{...n,...e}),(0,g.jsx)(b,{...n,...e})]})}function C(e){const n=(0,x.A)();return(0,g.jsx)(v,{...e,children:u(e.children)},String(n))}},8453:(e,n,t)=>{t.d(n,{R:()=>r,x:()=>l});var s=t(6540);const a={},i=s.createContext(a);function r(e){const n=s.useContext(i);return s.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function l(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(a):e.components||a:r(e.components),s.createElement(i.Provider,{value:n},e.children)}}}]);