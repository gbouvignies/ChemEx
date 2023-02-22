"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[4320],{3905:(e,n,t)=>{t.d(n,{Zo:()=>s,kt:()=>d});var r=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function o(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function a(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?o(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):o(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,r,i=function(e,n){if(null==e)return{};var t,r,i={},o=Object.keys(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(r=0;r<o.length;r++)t=o[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var c=r.createContext({}),p=function(e){var n=r.useContext(c),t=n;return e&&(t="function"==typeof e?e(n):a(a({},n),e)),t},s=function(e){var n=p(e.components);return r.createElement(c.Provider,{value:n},e.children)},m={inlineCode:"code",wrapper:function(e){var n=e.children;return r.createElement(r.Fragment,{},n)}},u=r.forwardRef((function(e,n){var t=e.components,i=e.mdxType,o=e.originalType,c=e.parentName,s=l(e,["components","mdxType","originalType","parentName"]),u=p(t),d=i,h=u["".concat(c,".").concat(d)]||u[d]||m[d]||o;return t?r.createElement(h,a(a({ref:n},s),{},{components:t})):r.createElement(h,a({ref:n},s))}));function d(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var o=t.length,a=new Array(o);a[0]=u;var l={};for(var c in n)hasOwnProperty.call(n,c)&&(l[c]=n[c]);l.originalType=e,l.mdxType="string"==typeof e?e:i,a[1]=l;for(var p=2;p<o;p++)a[p]=t[p];return r.createElement.apply(null,a)}return r.createElement.apply(null,t)}u.displayName="MDXCreateElement"},9797:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>c,contentTitle:()=>a,default:()=>m,frontMatter:()=>o,metadata:()=>l,toc:()=>p});var r=t(7462),i=(t(7294),t(3905));const o={sidebar_label:"Methyl \xb9\xb3C (\u201cH-to-C\u201d)",sidebar_position:10,description:'"cpmg_ch3_13c_h2c"'},a="Methyl \xb9\xb3C CPMG (\u201cH-to-C\u201d)",l={unversionedId:"experiments/cpmg/cpmg_ch3_13c_h2c",id:"experiments/cpmg/cpmg_ch3_13c_h2c",title:"Methyl \xb9\xb3C CPMG (\u201cH-to-C\u201d)",description:'"cpmg_ch3_13c_h2c"',source:"@site/docs/experiments/cpmg/cpmg_ch3_13c_h2c.md",sourceDirName:"experiments/cpmg",slug:"/experiments/cpmg/cpmg_ch3_13c_h2c",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_13c_h2c",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/cpmg/cpmg_ch3_13c_h2c.md",tags:[],version:"current",sidebarPosition:10,frontMatter:{sidebar_label:"Methyl \xb9\xb3C (\u201cH-to-C\u201d)",sidebar_position:10,description:'"cpmg_ch3_13c_h2c"'},sidebar:"tutorialSidebar",previous:{title:"Pure anti-phase carbonyl \xb9\xb3C",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_13co_ap"},next:{title:"Methyl \xb9\xb3C\u2013\xb9H multiple-quantum",permalink:"/ChemEx/docs/experiments/cpmg/cpmg_ch3_mq"}},c={},p=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],s={toc:p};function m(e){let{components:n,...t}=e;return(0,i.kt)("wrapper",(0,r.Z)({},s,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"methyl-c-cpmg-h-to-c"},"Methyl \xb9\xb3C CPMG (\u201cH-to-C\u201d)"),(0,i.kt)("h2",{id:"module-name"},"Module name"),(0,i.kt)("p",null,(0,i.kt)("inlineCode",{parentName:"p"},'"cpmg_ch3_13c_h2c"')),(0,i.kt)("h2",{id:"description"},"Description"),(0,i.kt)("p",null,"Measures methyl carbon chemical exchange recorded on site-specifically\n13CH3-labeled proteins in a highly deuterated background. Magnetization is\ninitially anti-phase and is read out as in-phase. Because of the P-element only\neven ncyc should be recorded. Resulting magnetization intensity after the CPMG\nblock is calculated using the (6n)\xd7(6n), two-spin matrix, where n is the number\nof states:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),\n  Ix(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }\n")),(0,i.kt)("h2",{id:"reference"},"Reference"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"P. Lundstr\xf6m, P. Vallurupalli, T.M. Religa, F.W. Dahlquist, and L.E. Kay. ",(0,i.kt)("em",{parentName:"li"},"J.\nBiomol. NMR")," ",(0,i.kt)("strong",{parentName:"li"},"38"),", 79-88 (2007)")),(0,i.kt)("h2",{id:"example"},"Example"),(0,i.kt)("p",null,"An example use of the module is given\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_CH3_13C_H2C/"},"here"),"."),(0,i.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'cpmg_ch3_13c_h2c\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "cpmg_ch3_13c_h2c"\n\n## CPMG relaxation delay, in seconds\ntime_t2 = 0.04\n\n## Position of the \xb9\xb3C carrier during the CPMG period, in ppm\ncarrier = 20.0\n\n## \xb9\xb3C 90 degree pulse width of CPMG pulses, in seconds\npw90 = 15.0e-6\n\n## Equilibration delay at the end of the CPMG period, in seconds\n## [optional, default: 0.0]\n# time_equil = 0.0\n\n## P-element delay = 1/4J, in seconds [optional, default: 2.0e-3]\n# taub = 2.0e-3\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  L3CD1-HD1  = "L3CD1-QD1.out"\n  L3CD2-HD2  = "L3CD2-QD2.out"\n  I12CD1-HD1 = "I12CD1-QD1.out"\n  V25CG1-HG1 = "V25CG1-QG1.out"\n  V25CG2-HG2 = "V25CG2-QG2.out"\n  M36CE-HE   = "M36CE-QE.out"\n')))}m.isMDXComponent=!0}}]);