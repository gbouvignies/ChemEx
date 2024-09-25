"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[5087],{5612:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>s,contentTitle:()=>r,default:()=>m,frontMatter:()=>a,metadata:()=>l,toc:()=>c});var i=t(4848),o=t(8453);const a={sidebar_label:"Methyl \xb9H double quantum",sidebar_position:13,description:'"cpmg_ch3_1h_dq"'},r="Methyl \xb9H double quantum CPMG",l={id:"experiments/cpmg/cpmg_ch3_1h_dq",title:"Methyl \xb9H double quantum CPMG",description:'"cpmg_ch3_1h_dq"',source:"@site/docs/experiments/cpmg/cpmg_ch3_1h_dq.md",sourceDirName:"experiments/cpmg",slug:"/experiments/cpmg/cpmg_ch3_1h_dq",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_dq",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/cpmg/cpmg_ch3_1h_dq.md",tags:[],version:"current",sidebarPosition:13,frontMatter:{sidebar_label:"Methyl \xb9H double quantum",sidebar_position:13,description:'"cpmg_ch3_1h_dq"'},sidebar:"tutorialSidebar",previous:{title:"Methyl \xb9H single quantum",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_sq"},next:{title:"Methyl \xb9H triple quantum",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_tq"}},s={},c=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}];function d(e){const n={a:"a",code:"code",em:"em",h1:"h1",h2:"h2",header:"header",li:"li",p:"p",pre:"pre",strong:"strong",ul:"ul",...(0,o.R)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(n.header,{children:(0,i.jsx)(n.h1,{id:"methyl-h-double-quantum-cpmg",children:"Methyl \xb9H double quantum CPMG"})}),"\n",(0,i.jsx)(n.h2,{id:"module-name",children:"Module name"}),"\n",(0,i.jsx)(n.p,{children:(0,i.jsx)(n.code,{children:'"cpmg_ch3_1h_dq"'})}),"\n",(0,i.jsx)(n.h2,{id:"description",children:"Description"}),"\n",(0,i.jsx)(n.p,{children:"Measures methyl proton chemical exchange recorded on site-specifically\n13CH3-labeled proteins in a highly deuterated background. Magnetization is\ninitially anti-phase and is read out as anti-phase. The calculation uses a\nsimplified scheme that includes only (6n)\xd7(6n) basis set, two-spin matrix, where\nn is the number of states:"}),"\n",(0,i.jsx)(n.p,{children:"{ 2HDQxCz(a), 2HDQyCz(a), 2HzCz(a), HDQx(a), HDQy(a), Hz(a),\n2HDQxCz(b), 2HDQyCz(b), 2HzCz(b), HDQx(b), HDQy(b), Hz(b), ... }"}),"\n",(0,i.jsx)(n.h2,{id:"reference",children:"Reference"}),"\n",(0,i.jsxs)(n.ul,{children:["\n",(0,i.jsxs)(n.li,{children:["Gopalan, T. Yuwen, L.E. Kay, and P. Vallurupalli. ",(0,i.jsx)(n.em,{children:"J. Biomol. NMR"})," ",(0,i.jsx)(n.strong,{children:"72"}),",\n79-91 (2018)"]}),"\n"]}),"\n",(0,i.jsx)(n.h2,{id:"example",children:"Example"}),"\n",(0,i.jsxs)(n.p,{children:["An example, where joint fit of\n",(0,i.jsx)(n.a,{href:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_dq",children:"methyl \xb9H double quantum CPMG (cpmg_ch3_1h_dq)"})," and\n",(0,i.jsx)(n.a,{href:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_tq",children:"methyl \xb9H triple quantum CPMG (cpmg_ch3_1h_tq)"})," experiments\nis performed, is available\n",(0,i.jsx)(n.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/CPMG_CH3_1H_DQ_TQ/",children:"here"}),"."]}),"\n",(0,i.jsx)(n.h2,{id:"sample-configuration-file",children:"Sample configuration file"}),"\n",(0,i.jsx)(n.pre,{children:(0,i.jsx)(n.code,{className:"language-toml",metastring:'title="experiment.toml"',children:'## This is a sample configuration file for the module \'cpmg_ch3_1h_dq\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "cpmg_ch3_1h_dq"\n\n## CPMG relaxation delay, in seconds\ntime_t2 = 0.04\n\n## Position of the 1H carrier during the CPMG period, in ppm\ncarrier = 0.5\n\n## 1H 90 degree pulse width of CPMG pulses, in seconds\npw90 = 12.0e-6\n\n## Performs CPMG using 90-180-90 composite pulses [optional, default: true]\n# comp180_flg = true\n\n## Enters the CPMG period with equal amount of initial IP and\n## AP DQ magnetizations [optional, default: false]\n# ipap_flg = false\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  L3HD1-CD1  = "L3CD1-QD1.out"\n  L3HD2-CD2  = "L3CD2-QD2.out"\n  I12HD1-CD1 = "I12CD1-QD1.out"\n  V25HG1-CG1 = "V25CG1-QG1.out"\n  V25HG2-CG2 = "V25CG2-QG2.out"\n  M36HE-CE   = "M36CE-QE.out"\n'})})]})}function m(e={}){const{wrapper:n}={...(0,o.R)(),...e.components};return n?(0,i.jsx)(n,{...e,children:(0,i.jsx)(d,{...e})}):d(e)}},8453:(e,n,t)=>{t.d(n,{R:()=>r,x:()=>l});var i=t(6540);const o={},a=i.createContext(o);function r(e){const n=i.useContext(a);return i.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function l(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(o):e.components||o:r(e.components),i.createElement(a.Provider,{value:n},e.children)}}}]);