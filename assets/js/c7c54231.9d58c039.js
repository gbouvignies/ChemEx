"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[4711],{3905:(e,n,t)=>{t.d(n,{Zo:()=>s,kt:()=>u});var r=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function a(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function o(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?a(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):a(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,r,i=function(e,n){if(null==e)return{};var t,r,i={},a=Object.keys(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var p=r.createContext({}),c=function(e){var n=r.useContext(p),t=n;return e&&(t="function"==typeof e?e(n):o(o({},n),e)),t},s=function(e){var n=c(e.components);return r.createElement(p.Provider,{value:n},e.children)},m={inlineCode:"code",wrapper:function(e){var n=e.children;return r.createElement(r.Fragment,{},n)}},d=r.forwardRef((function(e,n){var t=e.components,i=e.mdxType,a=e.originalType,p=e.parentName,s=l(e,["components","mdxType","originalType","parentName"]),d=c(t),u=i,h=d["".concat(p,".").concat(u)]||d[u]||m[u]||a;return t?r.createElement(h,o(o({ref:n},s),{},{components:t})):r.createElement(h,o({ref:n},s))}));function u(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var a=t.length,o=new Array(a);o[0]=d;var l={};for(var p in n)hasOwnProperty.call(n,p)&&(l[p]=n[p]);l.originalType=e,l.mdxType="string"==typeof e?e:i,o[1]=l;for(var c=2;c<a;c++)o[c]=t[c];return r.createElement.apply(null,o)}return r.createElement.apply(null,t)}d.displayName="MDXCreateElement"},2426:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>p,contentTitle:()=>o,default:()=>m,frontMatter:()=>a,metadata:()=>l,toc:()=>c});var r=t(7462),i=(t(7294),t(3905));const a={sidebar_label:"Pure anti-phase methyl \xb9H",sidebar_position:16,description:'"cpmg_chd2_1h_ap"'},o="Pure anti-phase methyl \xb9H CPMG",l={unversionedId:"experiments/cpmg/cpmg_chd2_1h_ap",id:"experiments/cpmg/cpmg_chd2_1h_ap",title:"Pure anti-phase methyl \xb9H CPMG",description:'"cpmg_chd2_1h_ap"',source:"@site/docs/experiments/cpmg/cpmg_chd2_1h_ap.md",sourceDirName:"experiments/cpmg",slug:"/experiments/cpmg/cpmg_chd2_1h_ap",permalink:"/chemex/docs/experiments/cpmg/cpmg_chd2_1h_ap",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/master/docs/experiments/cpmg/cpmg_chd2_1h_ap.md",tags:[],version:"current",sidebarPosition:16,frontMatter:{sidebar_label:"Pure anti-phase methyl \xb9H",sidebar_position:16,description:'"cpmg_chd2_1h_ap"'},sidebar:"tutorialSidebar",previous:{title:"Methyl \xb9H triple quantum with gradients",permalink:"/chemex/docs/experiments/cpmg/cpmg_ch3_1h_tq_diff"},next:{title:"Relaxation",permalink:"/chemex/docs/experiments/relaxation/"}},p={},c=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],s={toc:c};function m(e){let{components:n,...t}=e;return(0,i.kt)("wrapper",(0,r.Z)({},s,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"pure-anti-phase-methyl-h-cpmg"},"Pure anti-phase methyl \xb9H CPMG"),(0,i.kt)("h2",{id:"module-name"},"Module name"),(0,i.kt)("p",null,(0,i.kt)("inlineCode",{parentName:"p"},'"cpmg_chd2_1h_ap"')),(0,i.kt)("h2",{id:"description"},"Description"),(0,i.kt)("p",null,"Measures methyl proton chemical exchange recorded on site-specifically\n13CHD2-labeled proteins in a highly deuterated background. Magnetization is\ninitially anti-phase and is read out as anti-phase prior to \xb9\xb3C evolution.\nResulting magnetization intensity after the CPMG block is calculated using the\n(6n)\xd7(6n), two-spin matrix, where n is the number of states:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),\n  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }\n")),(0,i.kt)("h2",{id:"reference"},"Reference"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"A.J. Baldwin, T.M. Religa, D.F. Hansen, G. Bouvignies, and L.E. Kay. ",(0,i.kt)("em",{parentName:"li"},"J. Am.\nChem. Soc.")," ",(0,i.kt)("strong",{parentName:"li"},"132"),", 10992-10995 (2010)")),(0,i.kt)("h2",{id:"example"},"Example"),(0,i.kt)("p",null,"An example use of the module is given\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CHD2_1H_AP/"},"here"),"."),(0,i.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'cpmg_chd2_1h_ap\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "cpmg_chd2_1h_ap"\n\n## CPMG relaxation delay, in seconds\ntime_t2 = 0.02\n\n## Position of the 1H carrier during the CPMG period, in ppm\ncarrier = 0.5\n\n## 1H 90 degree pulse width of CPMG pulses, in seconds\npw90 = 10.0e-6\n\n## Equilibration delay at the end of the CPMG period, in seconds\n## [optional, default: 0.0]\n# time_equil = 0.0\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n## Labeling scheme of the sample, for deuterated samples "2H" should\n## be used to obtain accurate initial estimates of relaxation rates\n## based on model-free parameters [optional, default: []]\n# label = ["2H"]\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  L3HD1-CD1  = "L3CD1-QD1.out"\n  L3HD2-CD2  = "L3CD2-QD2.out"\n  I12HD1-CD1 = "I12CD1-QD1.out"\n  V25HG1-CG1 = "V25CG1-QG1.out"\n  V25HG2-CG2 = "V25CG2-QG2.out"\n  M36HE-CE   = "M36CE-QE.out"\n')))}m.isMDXComponent=!0}}]);