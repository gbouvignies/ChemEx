"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[9763],{4567:(e,n,i)=>{i.r(n),i.d(n,{assets:()=>l,contentTitle:()=>r,default:()=>d,frontMatter:()=>o,metadata:()=>a,toc:()=>h});var t=i(4848),s=i(8453);const o={sidebar_label:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",sidebar_position:1,description:'"shift_15n_sqmq"'},r="\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",a={id:"experiments/shift/shift_15n_sqmq",title:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",description:'"shift_15n_sqmq"',source:"@site/docs/experiments/shift/shift_15n_sqmq.md",sourceDirName:"experiments/shift",slug:"/experiments/shift/shift_15n_sqmq",permalink:"/ChemEx/docs/experiments/shift/shift_15n_sqmq",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/shift/shift_15n_sqmq.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_label:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",sidebar_position:1,description:'"shift_15n_sqmq"'},sidebar:"tutorialSidebar",previous:{title:"Exchange-induced chemical shift experiments",permalink:"/ChemEx/docs/experiments/shift/"},next:{title:"Examples",permalink:"/ChemEx/docs/examples/"}},l={},h=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"References",id:"references",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}];function c(e){const n={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",li:"li",p:"p",pre:"pre",strong:"strong",ul:"ul",...(0,s.R)(),...e.components};return(0,t.jsxs)(t.Fragment,{children:[(0,t.jsx)(n.h1,{id:"n-exchange-induced-shifts-with-nh-hsqchmqc",children:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC"}),"\n",(0,t.jsx)(n.h2,{id:"module-name",children:"Module name"}),"\n",(0,t.jsx)(n.p,{children:(0,t.jsx)(n.code,{children:'"shift_15n_sqmq"'})}),"\n",(0,t.jsx)(n.h2,{id:"description",children:"Description"}),"\n",(0,t.jsx)(n.p,{children:"Analyzes exchange induced \xb9\u2075N chemical shift changes measured in (\xb9\u2075N\u2013\xb9HN) HMQC\nand HSQC data sets."}),"\n",(0,t.jsx)(n.admonition,{type:"note",children:(0,t.jsx)(n.p,{children:"Since this experiment is used for determining the sign of \u0394\u03d6, it is usually\ncombined with other CPMG experiments."})}),"\n",(0,t.jsx)(n.h2,{id:"references",children:"References"}),"\n",(0,t.jsxs)(n.ul,{children:["\n",(0,t.jsxs)(n.li,{children:["N.R. Skrynnikov, F.W. Dahlquist, L.E. Kay. ",(0,t.jsx)(n.em,{children:"J. Am. Chem. Soc."})," ",(0,t.jsx)(n.strong,{children:"124"}),",\n12352-12360 (2002)"]}),"\n",(0,t.jsxs)(n.li,{children:["P. Vallurupalli, G. Bouvignies, and L.E. Kay. ",(0,t.jsx)(n.em,{children:"J. Phys. Chem. B"})," ",(0,t.jsx)(n.strong,{children:"115"}),",\n14891-14900 (2011)"]}),"\n"]}),"\n",(0,t.jsx)(n.h2,{id:"example",children:"Example"}),"\n",(0,t.jsxs)(n.p,{children:["An example use of the module associated with \xb9\u2075N and \xb9H CPMG datasets is given\n",(0,t.jsx)(n.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/Shifts/",children:"here"}),"."]}),"\n",(0,t.jsx)(n.h2,{id:"sample-configuration-file",children:"Sample configuration file"}),"\n",(0,t.jsx)(n.pre,{children:(0,t.jsx)(n.code,{className:"language-toml",metastring:'title="experiment.toml"',children:'## This is a sample configuration file for the module \'shift_15n_sqmq\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "shift_15n_sqmq"\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Filename of the file containing the list of the shifts in ppb.\n## The file should be formatted as follow:\n##\n## #     name   shift  error\n##     G2N-HN    10.9    0.5\n##     H3N-HN    32.1    0.5\n##     K4N-HN   -54.3    1.5\n##     S5N-HN     0.7    0.5\n##     L6N-HN   -15.2    0.5\n##\n## The name of the spin systems should follow the Sparky-NMR\n## conventions.\nshifts = "sqmq.txt"\n'})})]})}function d(e={}){const{wrapper:n}={...(0,s.R)(),...e.components};return n?(0,t.jsx)(n,{...e,children:(0,t.jsx)(c,{...e})}):c(e)}},8453:(e,n,i)=>{i.d(n,{R:()=>r,x:()=>a});var t=i(6540);const s={},o=t.createContext(s);function r(e){const n=t.useContext(o);return t.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function a(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(s):e.components||s:r(e.components),t.createElement(o.Provider,{value:n},e.children)}}}]);