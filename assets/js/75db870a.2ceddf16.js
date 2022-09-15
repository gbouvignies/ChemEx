"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[788],{3905:(e,t,r)=>{r.d(t,{Zo:()=>m,kt:()=>f});var n=r(7294);function a(e,t,r){return t in e?Object.defineProperty(e,t,{value:r,enumerable:!0,configurable:!0,writable:!0}):e[t]=r,e}function o(e,t){var r=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),r.push.apply(r,n)}return r}function i(e){for(var t=1;t<arguments.length;t++){var r=null!=arguments[t]?arguments[t]:{};t%2?o(Object(r),!0).forEach((function(t){a(e,t,r[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(r)):o(Object(r)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(r,t))}))}return e}function s(e,t){if(null==e)return{};var r,n,a=function(e,t){if(null==e)return{};var r,n,a={},o=Object.keys(e);for(n=0;n<o.length;n++)r=o[n],t.indexOf(r)>=0||(a[r]=e[r]);return a}(e,t);if(Object.getOwnPropertySymbols){var o=Object.getOwnPropertySymbols(e);for(n=0;n<o.length;n++)r=o[n],t.indexOf(r)>=0||Object.prototype.propertyIsEnumerable.call(e,r)&&(a[r]=e[r])}return a}var p=n.createContext({}),c=function(e){var t=n.useContext(p),r=t;return e&&(r="function"==typeof e?e(t):i(i({},t),e)),r},m=function(e){var t=c(e.components);return n.createElement(p.Provider,{value:t},e.children)},l={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},d=n.forwardRef((function(e,t){var r=e.components,a=e.mdxType,o=e.originalType,p=e.parentName,m=s(e,["components","mdxType","originalType","parentName"]),d=c(r),f=a,u=d["".concat(p,".").concat(f)]||d[f]||l[f]||o;return r?n.createElement(u,i(i({ref:t},m),{},{components:r})):n.createElement(u,i({ref:t},m))}));function f(e,t){var r=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var o=r.length,i=new Array(o);i[0]=d;var s={};for(var p in t)hasOwnProperty.call(t,p)&&(s[p]=t[p]);s.originalType=e,s.mdxType="string"==typeof e?e:a,i[1]=s;for(var c=2;c<o;c++)i[c]=r[c];return n.createElement.apply(null,i)}return n.createElement.apply(null,r)}d.displayName="MDXCreateElement"},4019:(e,t,r)=>{r.r(t),r.d(t,{assets:()=>p,contentTitle:()=>i,default:()=>l,frontMatter:()=>o,metadata:()=>s,toc:()=>c});var n=r(7462),a=(r(7294),r(3905));const o={sidebar_position:4,sidebar_label:"Measuring Excited state RDCs with CPMG"},i="RDC measurement for excited state with \xb9\u2075N TROSY CPMG",s={unversionedId:"examples/trosy_cpmg_rdc",id:"examples/trosy_cpmg_rdc",title:"RDC measurement for excited state with \xb9\u2075N TROSY CPMG",description:"You can get",source:"@site/docs/examples/trosy_cpmg_rdc.md",sourceDirName:"examples",slug:"/examples/trosy_cpmg_rdc",permalink:"/ChemEx/docs/examples/trosy_cpmg_rdc",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/examples/trosy_cpmg_rdc.md",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_position:4,sidebar_label:"Measuring Excited state RDCs with CPMG"},sidebar:"tutorialSidebar",previous:{title:"Measuring diffusion constants of invisible protein conformers",permalink:"/ChemEx/docs/examples/diffusion_tqcpmg"},next:{title:"CEST experiments for \xb9\xb3C, \xb9\u2075N-labeled samples",permalink:"/ChemEx/docs/examples/cest_13c_15n"}},p={},c=[],m={toc:c};function l(e){let{components:t,...r}=e;return(0,a.kt)("wrapper",(0,n.Z)({},m,r,{components:t,mdxType:"MDXLayout"}),(0,a.kt)("h1",{id:"rdc-measurement-for-excited-state-with-n-trosy-cpmg"},"RDC measurement for excited state with \xb9\u2075N TROSY CPMG"),(0,a.kt)("admonition",{title:"Files",type:"note"},(0,a.kt)("p",{parentName:"admonition"},"You can get\n",(0,a.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/N15_NH_RDC"},"the example files from the GitHub ChemEx page"),".")),(0,a.kt)("p",null,"This example demonstrates the use of CPMG experiment to measure RDC parameters\nfor the excited state. In \xb9\u2075N TROSY CPMG experiment\n(",(0,a.kt)("a",{parentName:"p",href:"../experiments/cpmg/cpmg_15n_tr"},(0,a.kt)("inlineCode",{parentName:"a"},"cpmg_15n_tr")),"), by analyzing datasets measured\nfor TROSY and anti-TROSY components, the information about \xb9\u2075N \u0394\u03d6 for both\ncomponents can be obtained, therefore RDC parameters for excited state can be\nderived based on RDC parameters for ground state together with \u0394\u03d6 of these two\ncomponents ",(0,a.kt)("sup",{parentName:"p",id:"fnref-1"},(0,a.kt)("a",{parentName:"sup",href:"#fn-1",className:"footnote-ref"},"1")),"."),(0,a.kt)("p",null,"Note that ",(0,a.kt)("inlineCode",{parentName:"p"},"antitrosy")," key should be set properly in experiment files to indicate\nwhether the datasets are measured for TROSY or anti-TROSY component. This\nexample also includes datasets measured from pure in-phase CPMG experiment\n(",(0,a.kt)("a",{parentName:"p",href:"../experiments/cpmg/cpmg_15n_ip"},(0,a.kt)("inlineCode",{parentName:"a"},"cpmg_15n_ip")),'), which is optional and may help to\nobtain even more accurate results. Note that \xb9\u2075N chemical shifts provided in\nparameter files should correspond to the "actual" value in the absence of 1JHN,\ntherefore \xb9\u2075N chemical shifts from pure in-phase experiment should be used\ninstead of the TROSY version.'),(0,a.kt)("div",{className:"footnotes"},(0,a.kt)("hr",{parentName:"div"}),(0,a.kt)("ol",{parentName:"div"},(0,a.kt)("li",{parentName:"ol",id:"fn-1"},"P. Vallurupalli, D. F. Hansen, E. Stollar, E. Meirovitch, and L. E. Kay.\nMeasurement of Bond Vector Orientations in Invisible Excited States of\nProteins ",(0,a.kt)("em",{parentName:"li"},"Proc. Natl. Acad. Sci. USA")," ",(0,a.kt)("strong",{parentName:"li"},"104"),", 18473-18477 (2007).\n",(0,a.kt)("a",{parentName:"li",href:"https://doi.org/10.1073/pnas.0708296104"},"https://doi.org/10.1073/pnas.0708296104"),(0,a.kt)("a",{parentName:"li",href:"#fnref-1",className:"footnote-backref"},"\u21a9")))))}l.isMDXComponent=!0}}]);