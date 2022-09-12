"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[3541],{3905:(e,n,t)=>{t.d(n,{Zo:()=>c,kt:()=>d});var a=t(7294);function r(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function i(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);n&&(a=a.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,a)}return t}function o(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?i(Object(t),!0).forEach((function(n){r(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):i(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,a,r=function(e,n){if(null==e)return{};var t,a,r={},i=Object.keys(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||(r[t]=e[t]);return r}(e,n);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(a=0;a<i.length;a++)t=i[a],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(r[t]=e[t])}return r}var p=a.createContext({}),s=function(e){var n=a.useContext(p),t=n;return e&&(t="function"==typeof e?e(n):o(o({},n),e)),t},c=function(e){var n=s(e.components);return a.createElement(p.Provider,{value:n},e.children)},m={inlineCode:"code",wrapper:function(e){var n=e.children;return a.createElement(a.Fragment,{},n)}},u=a.forwardRef((function(e,n){var t=e.components,r=e.mdxType,i=e.originalType,p=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),u=s(t),d=r,f=u["".concat(p,".").concat(d)]||u[d]||m[d]||i;return t?a.createElement(f,o(o({ref:n},c),{},{components:t})):a.createElement(f,o({ref:n},c))}));function d(e,n){var t=arguments,r=n&&n.mdxType;if("string"==typeof e||r){var i=t.length,o=new Array(i);o[0]=u;var l={};for(var p in n)hasOwnProperty.call(n,p)&&(l[p]=n[p]);l.originalType=e,l.mdxType="string"==typeof e?e:r,o[1]=l;for(var s=2;s<i;s++)o[s]=t[s];return a.createElement.apply(null,o)}return a.createElement.apply(null,t)}u.displayName="MDXCreateElement"},1547:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>p,contentTitle:()=>o,default:()=>m,frontMatter:()=>i,metadata:()=>l,toc:()=>s});var a=t(7462),r=(t(7294),t(3905));const i={sidebar_label:"Pure anti-phase amide \xb9H with [0013] phase cycle",sidebar_position:6,description:'"cpmg_1hn_ap_0013"'},o="Pure anti-phase amide \xb9H CPMG with [0013] phase cycle",l={unversionedId:"experiments/cpmg/cpmg_1hn_ap_0013",id:"experiments/cpmg/cpmg_1hn_ap_0013",title:"Pure anti-phase amide \xb9H CPMG with [0013] phase cycle",description:'"cpmg_1hn_ap_0013"',source:"@site/docs/experiments/cpmg/cpmg_1hn_ap_0013.md",sourceDirName:"experiments/cpmg",slug:"/experiments/cpmg/cpmg_1hn_ap_0013",permalink:"/chemex/docs/experiments/cpmg/cpmg_1hn_ap_0013",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/master/docs/experiments/cpmg/cpmg_1hn_ap_0013.md",tags:[],version:"current",sidebarPosition:6,frontMatter:{sidebar_label:"Pure anti-phase amide \xb9H with [0013] phase cycle",sidebar_position:6,description:'"cpmg_1hn_ap_0013"'},sidebar:"tutorialSidebar",previous:{title:"Pure anti-phase amide \xb9H",permalink:"/chemex/docs/experiments/cpmg/cpmg_1hn_ap"},next:{title:"Amide \xb9\u2075N\u2013\xb9H double-quantum/zero-quantum",permalink:"/chemex/docs/experiments/cpmg/cpmg_hn_dq_zq"}},p={},s=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],c={toc:s};function m(e){let{components:n,...t}=e;return(0,r.kt)("wrapper",(0,a.Z)({},c,t,{components:n,mdxType:"MDXLayout"}),(0,r.kt)("h1",{id:"pure-anti-phase-amide-h-cpmg-with-0013-phase-cycle"},"Pure anti-phase amide \xb9H CPMG with ","[0013]"," phase cycle"),(0,r.kt)("h2",{id:"module-name"},"Module name"),(0,r.kt)("p",null,(0,r.kt)("inlineCode",{parentName:"p"},'"cpmg_1hn_ap_0013"')),(0,r.kt)("h2",{id:"description"},"Description"),(0,r.kt)("p",null,"Analyzes amide proton chemical exchange that is maintained as anti-phase\nmagnetization throughout the CPMG block. This results in lower intrinsic\nrelaxation rates and therefore better sensitivity. The calculations use the\n(6n)\xd7(6n), two-spin matrix, where n is the number of states:"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre"},"{ Hx(a), Hy(a), Hz(a), 2HxNz(a), 2HyNz(a), 2HzNz(a),\n  Hx(b), Hy(b), Hz(b), 2HxNz(b), 2HyNz(b), 2HzNz(b), ... }\n")),(0,r.kt)("p",null,"This version is modified such that CPMG pulses are applied with ","[0013]"," phase\ncycle to help better overcome off-resonance effects."),(0,r.kt)("h2",{id:"reference"},"Reference"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},"T. Yuwen, and L.E. Kay. ",(0,r.kt)("em",{parentName:"li"},"J. Biomol. NMR")," ",(0,r.kt)("strong",{parentName:"li"},"73"),", 641-650 (2019)")),(0,r.kt)("h2",{id:"example"},"Example"),(0,r.kt)("p",null,"An example use of the module is given\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CPMG_1HN_AP_0013/"},"here"),"."),(0,r.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'cpmg_1hn_ap_0013\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "cpmg_1hn_ap_0013"\n\n## CPMG relaxation delay, in seconds\ntime_t2 = 0.02\n\n## Position of the 1H carrier during the CPMG period, in ppm\ncarrier = 8.3\n\n## 1H 90 degree pulse width of CPMG pulses, in seconds\npw90 = 10.0e-6\n\n## Maximum number of cycles\nncyc_max = 40\n\n## 1/4J, in seconds [optional, default: 2.38e-3]\n# taua = 2.38e-3\n\n## Apply IPAP scheme for IP/AP differential relaxation suppression\n## [optional, default: false]\n# ipap_flg = false\n\n## Use the EBURP scheme instead of the standard 180 central pulse\n## [optional, default: false]\n# eburp_flg = false\n\n## Use the REBURP scheme instead of the standard 180 central pulse\n## [optional, default: false]\n# reburp_flg = false\n\n## 1H EBURP pulse width of CPMG pulses, in seconds [optional, default: 1.4e-3]\n# pw_eburp = 1.4e-3\n\n## 1H REBURP pulse width of CPMG pulses, in seconds [optional, default: 1.52e-3]\n# pw_reburp = 1.52e-3\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n## Labeling scheme of the sample, for deuterated samples "2H" should\n## be used to obtain accurate initial estimates of relaxation rates\n## based on model-free parameters [optional, default: []]\n# label = ["2H"]\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  G2HN-N = "G2N-HN.out"\n  H3HN-N = "H3N-HN.out"\n  K4HN-N = "K4N-HN.out"\n  S5HN-N = "S5N-HN.out"\n  L6HN-N = "L6N-HN.out"\n')))}m.isMDXComponent=!0}}]);