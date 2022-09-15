"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[4741],{3905:(e,n,t)=>{t.d(n,{Zo:()=>c,kt:()=>d});var r=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function a(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function o(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?a(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):a(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function s(e,n){if(null==e)return{};var t,r,i=function(e,n){if(null==e)return{};var t,r,i={},a=Object.keys(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var l=r.createContext({}),p=function(e){var n=r.useContext(l),t=n;return e&&(t="function"==typeof e?e(n):o(o({},n),e)),t},c=function(e){var n=p(e.components);return r.createElement(l.Provider,{value:n},e.children)},u={inlineCode:"code",wrapper:function(e){var n=e.children;return r.createElement(r.Fragment,{},n)}},m=r.forwardRef((function(e,n){var t=e.components,i=e.mdxType,a=e.originalType,l=e.parentName,c=s(e,["components","mdxType","originalType","parentName"]),m=p(t),d=i,f=m["".concat(l,".").concat(d)]||m[d]||u[d]||a;return t?r.createElement(f,o(o({ref:n},c),{},{components:t})):r.createElement(f,o({ref:n},c))}));function d(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var a=t.length,o=new Array(a);o[0]=m;var s={};for(var l in n)hasOwnProperty.call(n,l)&&(s[l]=n[l]);s.originalType=e,s.mdxType="string"==typeof e?e:i,o[1]=s;for(var p=2;p<a;p++)o[p]=t[p];return r.createElement.apply(null,o)}return r.createElement.apply(null,t)}m.displayName="MDXCreateElement"},7801:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>l,contentTitle:()=>o,default:()=>u,frontMatter:()=>a,metadata:()=>s,toc:()=>p});var r=t(7462),i=(t(7294),t(3905));const a={sidebar_label:"Pure in-phase \xb9\u2075N D-CEST",sidebar_position:1,description:'"dcest_15n"'},o="Pure in-phase \xb9\u2075N DANTE-CEST",s={unversionedId:"experiments/dcest/dcest_15n",id:"experiments/dcest/dcest_15n",title:"Pure in-phase \xb9\u2075N DANTE-CEST",description:'"dcest_15n"',source:"@site/docs/experiments/dcest/dcest_15n.md",sourceDirName:"experiments/dcest",slug:"/experiments/dcest/dcest_15n",permalink:"/ChemEx/docs/experiments/dcest/dcest_15n",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/dcest/dcest_15n.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_label:"Pure in-phase \xb9\u2075N D-CEST",sidebar_position:1,description:'"dcest_15n"'},sidebar:"tutorialSidebar",previous:{title:"D-CEST/cos-CEST",permalink:"/ChemEx/docs/experiments/dcest/"},next:{title:"Pure in-phase \xb9\xb3C D-CEST",permalink:"/ChemEx/docs/experiments/dcest/dcest_13c"}},l={},p=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"References",id:"references",level:2},{value:"Examples",id:"examples",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],c={toc:p};function u(e){let{components:n,...t}=e;return(0,i.kt)("wrapper",(0,r.Z)({},c,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"pure-in-phase-n-dante-cest"},"Pure in-phase \xb9\u2075N DANTE-CEST"),(0,i.kt)("h2",{id:"module-name"},"Module name"),(0,i.kt)("p",null,(0,i.kt)("inlineCode",{parentName:"p"},'"dcest_15n"')),(0,i.kt)("h2",{id:"description"},"Description"),(0,i.kt)("p",null,"Analyzes chemical exchange in the presence of \xb9H composite decoupling during the\nD-CEST block. This keeps the spin system purely in-phase throughout, and is\ncalculated using the (3n)\xd7(3n), single-spin matrix, where n is the number of\nstates:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"{ Ix(a), Iy(a), Iz(a),\n  Ix(b), Iy(b), Iz(b), ... }\n")),(0,i.kt)("h2",{id:"references"},"References"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"T. Yuwen, L.E. Kay, and G. Bouvignies. ",(0,i.kt)("em",{parentName:"li"},"ChemPhysChem")," ",(0,i.kt)("strong",{parentName:"li"},"19"),", 1707-1710 (2018)"),(0,i.kt)("li",{parentName:"ul"},"T. Yuwen, A. Bah, Brady, F. Ferrage, G. Bouvignies, and L.E. Kay. ",(0,i.kt)("em",{parentName:"li"},"J. Phys.\nChem. B")," ",(0,i.kt)("strong",{parentName:"li"},"122"),", 11206-11217 (2018)")),(0,i.kt)("h2",{id:"examples"},"Examples"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"An example for studying a two-state exchange system is given\n",(0,i.kt)("a",{parentName:"li",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N/"},"here"),"."),(0,i.kt)("li",{parentName:"ul"},"An example for studying a three-state exchange system is given\n",(0,i.kt)("a",{parentName:"li",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_3States/"},"here"),"."),(0,i.kt)("li",{parentName:"ul"},"An example for H/D solvent exchange measurement using \xb9\u2075N D-CEST can be found\n",(0,i.kt)("a",{parentName:"li",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_HD_EXCH/"},"here"),".")),(0,i.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'dcest_15n\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "dcest_15n"\n\n## CEST relaxation delay, in seconds\ntime_t1 = 0.5\n\n## Position of the \xb9\u2075N carrier during the CEST period, in ppm\ncarrier = 118.0\n\n## Pulse width of a 90 pulse at the power used during the DANTE, in seconds\npw90 = 45e-6\n\n## DANTE "spectral width", in Hz\nsw = 800.0\n\n## B1 radio-frequency field strength, in Hz\nb1_frq = 25.0\n\n## B1 inhomogeneity expressed as a fraction of \'b1_frq\'. If set to "inf",\n## a faster calculation takes place assuming full dephasing of the\n## magnetization components orthogonal to the effective field.\n## [optional, default: 0.1]\n# b1_inh_scale = 0.1\n\n## Number of points used to simulate B1 inhomogeneity, the larger\n## the longer the calculation. [optional, default: 11]\n# b1_inh_res = 11\n\n## Equilibration delay at the end of the CEST period, in seconds\n## [optional, default: 0.0]\n# time_equil = 0.0\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## \xb9H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n## The amount of D2O in the sample, in mass fraction\n## [optional, only used for HD exchange measurement]\n# d2o = 0.1\n\n## Labeling scheme of the sample, for deuterated samples "2H" should\n## be used to obtain accurate initial estimates of relaxation rates\n## based on model-free parameters, for uniformly \xb9\xb3C-labeled samples "13C"\n## should be used to account for 1JCC properly [optional, default: []]\n# label = ["2H", "13C"]\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "scatter": uncertainties are calculated from the baseline\n## [optional, default: "file"]\n# error = "file"\n\n## List of offsets relative to the main resonance position\n## (nu) and bandwidths (delta_nu) defining regions where\n## points are excluded from the calculation (nu +/- 0.5 * delta_nu),\n## both are in Hz [optional, default: [[0.0, 0.0]] ]\n# filter_offsets = [\n#   [0.0, 0.0],\n# ]\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  G2N = "G2N-HN.out"\n  H3N = "H3N-HN.out"\n  K4N = "K4N-HN.out"\n  S5N = "S5N-HN.out"\n  L6N = "L6N-HN.out"\n')))}u.isMDXComponent=!0}}]);