"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[5159],{3905:(e,t,n)=>{n.d(t,{Zo:()=>c,kt:()=>u});var i=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function a(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);t&&(i=i.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,i)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?a(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,i,r=function(e,t){if(null==e)return{};var n,i,r={},a=Object.keys(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var p=i.createContext({}),s=function(e){var t=i.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},c=function(e){var t=s(e.components);return i.createElement(p.Provider,{value:t},e.children)},m={inlineCode:"code",wrapper:function(e){var t=e.children;return i.createElement(i.Fragment,{},t)}},d=i.forwardRef((function(e,t){var n=e.components,r=e.mdxType,a=e.originalType,p=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),d=s(n),u=r,f=d["".concat(p,".").concat(u)]||d[u]||m[u]||a;return n?i.createElement(f,o(o({ref:t},c),{},{components:n})):i.createElement(f,o({ref:t},c))}));function u(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var a=n.length,o=new Array(a);o[0]=d;var l={};for(var p in t)hasOwnProperty.call(t,p)&&(l[p]=t[p]);l.originalType=e,l.mdxType="string"==typeof e?e:r,o[1]=l;for(var s=2;s<a;s++)o[s]=n[s];return i.createElement.apply(null,o)}return i.createElement.apply(null,n)}d.displayName="MDXCreateElement"},4669:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>p,contentTitle:()=>o,default:()=>m,frontMatter:()=>a,metadata:()=>l,toc:()=>s});var i=n(7462),r=(n(7294),n(3905));const a={sidebar_label:"Methyl \xb9H triple quantum with gradients",sidebar_position:15,description:'"cpmg_ch3_1h_tq_diff"'},o="Methyl \xb9H triple quantum CPMG with gradients",l={unversionedId:"experiments/cpmg/cpmg_ch3_1h_tq_diff",id:"experiments/cpmg/cpmg_ch3_1h_tq_diff",title:"Methyl \xb9H triple quantum CPMG with gradients",description:'"cpmg_ch3_1h_tq_diff"',source:"@site/docs/experiments/cpmg/cpmg_ch3_1h_tq_diff.md",sourceDirName:"experiments/cpmg",slug:"/experiments/cpmg/cpmg_ch3_1h_tq_diff",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_tq_diff",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/cpmg/cpmg_ch3_1h_tq_diff.md",tags:[],version:"current",sidebarPosition:15,frontMatter:{sidebar_label:"Methyl \xb9H triple quantum with gradients",sidebar_position:15,description:'"cpmg_ch3_1h_tq_diff"'},sidebar:"tutorialSidebar",previous:{title:"Methyl \xb9H triple quantum",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_1h_tq"},next:{title:"Pure anti-phase methyl \xb9H",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_chd2_1h_ap"}},p={},s=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],c={toc:s};function m(e){let{components:t,...n}=e;return(0,r.kt)("wrapper",(0,i.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("h1",{id:"methyl-h-triple-quantum-cpmg-with-gradients"},"Methyl \xb9H triple quantum CPMG with gradients"),(0,r.kt)("h2",{id:"module-name"},"Module name"),(0,r.kt)("p",null,(0,r.kt)("inlineCode",{parentName:"p"},'"cpmg_ch3_1h_tq_diff"')),(0,r.kt)("h2",{id:"description"},"Description"),(0,r.kt)("p",null,"Measures methyl proton chemical exchange recorded on site-specifically\n13CH3-labeled proteins in a highly deuterated background. The experiment is\nspecifically designed for measuring translational diffusion constants of\ninvisible states using a pulsed-field gradient approach that exploits methyl 1H\ntriple-quantum coherences. Note that the diffusion constants are given in \u03bcm\xb2/s.\nMagnetization is initially anti-phase and is read out as anti-phase. The\ncalculation uses a simplified scheme that includes only (6n)\xd7(6n) basis set,\ntwo-spin matrix, where n is the number of states:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"{ 2HTQxCz(a), 2HTQyCz(a), 2HzCz(a), HTQx(a), HTQy(a), Hz(a),\n  2HTQxCz(b), 2HTQyCz(b), 2HzCz(b), HTQx(b), HTQy(b), Hz(b), ... }\n")),(0,r.kt)("h2",{id:"reference"},"Reference"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},"T. Yuwen, A. Sekhar, A.J. Baldwin, P. Vallurupalli, and L.E. Kay. ",(0,r.kt)("em",{parentName:"li"},"Angew.\nChem. Int. Ed.")," (2018) 57, 16777-16780")),(0,r.kt)("h2",{id:"example"},"Example"),(0,r.kt)("p",null,"An example use of the module is given\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_1H_TQ_DIFF/"},"here"),"."),(0,r.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'cpmg_ch3_1h_tq_diff\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "cpmg_ch3_1h_tq_diff"\n\n## CPMG relaxation delay, in seconds\ntime_t2 = 0.04\n\n## Position of the 1H carrier during the CPMG period, in ppm\ncarrier = 0.5\n\n## 1H 90 degree pulse width of CPMG pulses, in seconds\npw90 = 12.0e-6\n\n## Duration of the encode/decode gradient, in seconds\ndelta = 1.0e-3\n\n## Dephaing gradient strength, in T/m\ngradient = 30.0e-2\n\n## Gradient recovery delay, in seconds [optional, default: 500e-6]\n# tau = 500.0e-6\n\n## Performs CPMG using 90-180-90 composite pulses [optional, default: true]\n# comp180_flg = true\n\n## Enters the CPMG period with equal amount of initial IP and\n## AP TQ magnetizations [optional, default: false]\n# ipap_flg = false\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n[data]\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  L3HD1-CD1  = "L3CD1-QD1.out"\n  L3HD2-CD2  = "L3CD2-QD2.out"\n  I12HD1-CD1 = "I12CD1-QD1.out"\n  V25HG1-CG1 = "V25CG1-QG1.out"\n  V25HG2-CG2 = "V25CG2-QG2.out"\n  M36HE-CE   = "M36CE-QE.out"\n')))}m.isMDXComponent=!0}}]);