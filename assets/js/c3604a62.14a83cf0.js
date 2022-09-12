"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[7375],{3905:(e,n,t)=>{t.d(n,{Zo:()=>c,kt:()=>u});var r=t(7294);function i(e,n,t){return n in e?Object.defineProperty(e,n,{value:t,enumerable:!0,configurable:!0,writable:!0}):e[n]=t,e}function a(e,n){var t=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);n&&(r=r.filter((function(n){return Object.getOwnPropertyDescriptor(e,n).enumerable}))),t.push.apply(t,r)}return t}function o(e){for(var n=1;n<arguments.length;n++){var t=null!=arguments[n]?arguments[n]:{};n%2?a(Object(t),!0).forEach((function(n){i(e,n,t[n])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(t)):a(Object(t)).forEach((function(n){Object.defineProperty(e,n,Object.getOwnPropertyDescriptor(t,n))}))}return e}function l(e,n){if(null==e)return{};var t,r,i=function(e,n){if(null==e)return{};var t,r,i={},a=Object.keys(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||(i[t]=e[t]);return i}(e,n);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(r=0;r<a.length;r++)t=a[r],n.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(e,t)&&(i[t]=e[t])}return i}var s=r.createContext({}),p=function(e){var n=r.useContext(s),t=n;return e&&(t="function"==typeof e?e(n):o(o({},n),e)),t},c=function(e){var n=p(e.components);return r.createElement(s.Provider,{value:n},e.children)},m={inlineCode:"code",wrapper:function(e){var n=e.children;return r.createElement(r.Fragment,{},n)}},d=r.forwardRef((function(e,n){var t=e.components,i=e.mdxType,a=e.originalType,s=e.parentName,c=l(e,["components","mdxType","originalType","parentName"]),d=p(t),u=i,f=d["".concat(s,".").concat(u)]||d[u]||m[u]||a;return t?r.createElement(f,o(o({ref:n},c),{},{components:t})):r.createElement(f,o({ref:n},c))}));function u(e,n){var t=arguments,i=n&&n.mdxType;if("string"==typeof e||i){var a=t.length,o=new Array(a);o[0]=d;var l={};for(var s in n)hasOwnProperty.call(n,s)&&(l[s]=n[s]);l.originalType=e,l.mdxType="string"==typeof e?e:i,o[1]=l;for(var p=2;p<a;p++)o[p]=t[p];return r.createElement.apply(null,o)}return r.createElement.apply(null,t)}d.displayName="MDXCreateElement"},752:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>s,contentTitle:()=>o,default:()=>m,frontMatter:()=>a,metadata:()=>l,toc:()=>p});var r=t(7462),i=(t(7294),t(3905));const a={sidebar_label:"Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation",sidebar_position:2,description:'"relaxation_hznz"'},o="Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation",l={unversionedId:"experiments/relaxation/relaxation_hznz",id:"experiments/relaxation/relaxation_hznz",title:"Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation",description:'"relaxation_hznz"',source:"@site/docs/experiments/relaxation/relaxation_hznz.md",sourceDirName:"experiments/relaxation",slug:"/experiments/relaxation/relaxation_hznz",permalink:"/chemex/docs/experiments/relaxation/relaxation_hznz",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/master/docs/experiments/relaxation/relaxation_hznz.md",tags:[],version:"current",sidebarPosition:2,frontMatter:{sidebar_label:"Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation",sidebar_position:2,description:'"relaxation_hznz"'},sidebar:"tutorialSidebar",previous:{title:"\xb9\u2075N longitudinal relaxation",permalink:"/chemex/docs/experiments/relaxation/relaxation_nz"},next:{title:"Exchange-induced chemical shift experiments",permalink:"/chemex/docs/experiments/shift/"}},s={},p=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"Reference",id:"reference",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],c={toc:p};function m(e){let{components:n,...t}=e;return(0,i.kt)("wrapper",(0,r.Z)({},c,t,{components:n,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"amide-h-n-longitudinal-two-spin-order-relaxation"},"Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation"),(0,i.kt)("h2",{id:"module-name"},"Module name"),(0,i.kt)("p",null,(0,i.kt)("inlineCode",{parentName:"p"},'"relaxation_hznz"')),(0,i.kt)("h2",{id:"description"},"Description"),(0,i.kt)("p",null,"Analyzes \xb9H-\xb9\u2075N longitudinal two-spin order relaxation experiments. Decay is\ncalculated using the (2n)\xd7(2n), two-spin matrix, where n is the number of\nstates:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre"},"{ Iz(a), 2IzSz(a),\n  Iz(b), 2IzSz(b), ... }\n")),(0,i.kt)("h2",{id:"reference"},"Reference"),(0,i.kt)("ul",null,(0,i.kt)("li",{parentName:"ul"},"D.F. Hansen, D. Yang, H. Feng, Z. Zhou, S. Wiesner, Y. Bai, and L.E. Kay. ",(0,i.kt)("em",{parentName:"li"},"J.\nAm. Chem. Soc.")," ",(0,i.kt)("strong",{parentName:"li"},"129"),", 11468-11479 (2007)")),(0,i.kt)("h2",{id:"example"},"Example"),(0,i.kt)("p",null,"An example use of the module is given\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/RELAXATION_HZNZ/"},"here"),"."),(0,i.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'relaxation_hznz\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "relaxation_hznz"\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n## Labeling scheme of the sample, for deuterated samples "2H" should\n## be used to obtain accurate initial estimates of relaxation rates\n## based on model-free parameters [optional, default: []]\n# label = ["2H"]\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "duplicates": uncertainties are calculated from the duplicate points\n## [optional, default: "file"]\n# error = "file"\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  G2N-HN = "G2N-HN.out"\n  H3N-HN = "H3N-HN.out"\n  K4N-HN = "K4N-HN.out"\n  S5N-HN = "S5N-HN.out"\n  L6N-HN = "L6N-HN.out"\n')))}m.isMDXComponent=!0}}]);