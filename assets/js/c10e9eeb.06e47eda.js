"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[7878],{3905:(e,t,a)=>{a.d(t,{Zo:()=>m,kt:()=>h});var n=a(7294);function r(e,t,a){return t in e?Object.defineProperty(e,t,{value:a,enumerable:!0,configurable:!0,writable:!0}):e[t]=a,e}function s(e,t){var a=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),a.push.apply(a,n)}return a}function i(e){for(var t=1;t<arguments.length;t++){var a=null!=arguments[t]?arguments[t]:{};t%2?s(Object(a),!0).forEach((function(t){r(e,t,a[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(a)):s(Object(a)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(a,t))}))}return e}function o(e,t){if(null==e)return{};var a,n,r=function(e,t){if(null==e)return{};var a,n,r={},s=Object.keys(e);for(n=0;n<s.length;n++)a=s[n],t.indexOf(a)>=0||(r[a]=e[a]);return r}(e,t);if(Object.getOwnPropertySymbols){var s=Object.getOwnPropertySymbols(e);for(n=0;n<s.length;n++)a=s[n],t.indexOf(a)>=0||Object.prototype.propertyIsEnumerable.call(e,a)&&(r[a]=e[a])}return r}var p=n.createContext({}),c=function(e){var t=n.useContext(p),a=t;return e&&(a="function"==typeof e?e(t):i(i({},t),e)),a},m=function(e){var t=c(e.components);return n.createElement(p.Provider,{value:t},e.children)},l={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},u=n.forwardRef((function(e,t){var a=e.components,r=e.mdxType,s=e.originalType,p=e.parentName,m=o(e,["components","mdxType","originalType","parentName"]),u=c(a),h=r,d=u["".concat(p,".").concat(h)]||u[h]||l[h]||s;return a?n.createElement(d,i(i({ref:t},m),{},{components:a})):n.createElement(d,i({ref:t},m))}));function h(e,t){var a=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var s=a.length,i=new Array(s);i[0]=u;var o={};for(var p in t)hasOwnProperty.call(t,p)&&(o[p]=t[p]);o.originalType=e,o.mdxType="string"==typeof e?e:r,i[1]=o;for(var c=2;c<s;c++)i[c]=a[c];return n.createElement.apply(null,i)}return n.createElement.apply(null,a)}u.displayName="MDXCreateElement"},1069:(e,t,a)=>{a.r(t),a.d(t,{assets:()=>p,contentTitle:()=>i,default:()=>l,frontMatter:()=>s,metadata:()=>o,toc:()=>c});var n=a(7462),r=(a(7294),a(3905));const s={sidebar_position:1},i="Introduction",o={unversionedId:"user_guide/introduction",id:"user_guide/introduction",title:"Introduction",description:"Context",source:"@site/docs/user_guide/introduction.md",sourceDirName:"user_guide",slug:"/user_guide/introduction",permalink:"/ChemEx/docs/user_guide/introduction",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/introduction.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_position:1},sidebar:"tutorialSidebar",previous:{title:"User Guide",permalink:"/ChemEx/docs/user_guide/"},next:{title:"Running ChemEx",permalink:"/ChemEx/docs/user_guide/running_chemex"}},p={},c=[{value:"Context",id:"context",level:2},{value:"About ChemEx",id:"about-chemex",level:2}],m={toc:c};function l(e){let{components:t,...s}=e;return(0,r.kt)("wrapper",(0,n.Z)({},m,s,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("h1",{id:"introduction"},"Introduction"),(0,r.kt)("h2",{id:"context"},"Context"),(0,r.kt)("p",null,"Biological macromolecules such as proteins and nucleic acids are inherently\nflexible molecules, whose function (or malfunction) critically depends on\ninterconversions between different conformational states. Understanding how\nbiomolecules work therefore requires a quantitative characterization of their\nthermally accessible conformational landscape. The task is challenging, however,\nbecause many of the populated conformers are only transiently formed and\nmarginally populated, so that they remain invisible to most standard biophysical\ntechniques. Over the past two decades, nuclear magnetic resonance (NMR) has\nemerged as an extremely powerful tool to study these elusive states at atomic\nresolution. Central to this success has been the development of chemical\nexchange-based approaches, such as (CPMG, R",(0,r.kt)("span",{parentName:"p",className:"math math-inline"},(0,r.kt)("span",{parentName:"span",className:"katex"},(0,r.kt)("span",{parentName:"span",className:"katex-mathml"},(0,r.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},(0,r.kt)("semantics",{parentName:"math"},(0,r.kt)("mrow",{parentName:"semantics"},(0,r.kt)("msub",{parentName:"mrow"},(0,r.kt)("mrow",{parentName:"msub"}),(0,r.kt)("mrow",{parentName:"msub"},(0,r.kt)("mn",{parentName:"mrow"},"1"),(0,r.kt)("mi",{parentName:"mrow"},"\u03c1")))),(0,r.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"_{1\u03c1}")))),(0,r.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,r.kt)("span",{parentName:"span",className:"base"},(0,r.kt)("span",{parentName:"span",className:"strut",style:{height:"0.5872em",verticalAlign:"-0.2861em"}}),(0,r.kt)("span",{parentName:"span",className:"mord"},(0,r.kt)("span",{parentName:"span"}),(0,r.kt)("span",{parentName:"span",className:"msupsub"},(0,r.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,r.kt)("span",{parentName:"span",className:"vlist-r"},(0,r.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.3011em"}},(0,r.kt)("span",{parentName:"span",style:{top:"-2.55em",marginRight:"0.05em"}},(0,r.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,r.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,r.kt)("span",{parentName:"span",className:"mord mtight"},(0,r.kt)("span",{parentName:"span",className:"mord mtight"},"1"),(0,r.kt)("span",{parentName:"span",className:"mord mathnormal mtight"},"\u03c1"))))),(0,r.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,r.kt)("span",{parentName:"span",className:"vlist-r"},(0,r.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.2861em"}},(0,r.kt)("span",{parentName:"span"})))))))))),') relaxation dispersion and\nchemical exchange saturation transfer (CEST) experiments, with the numerical\ntools needed to extract the kinetic (rates) and thermodynamic (populations)\nparameters associated with the exchange process and the structural information\n(chemical shifts) of the sparsely populated, "invisible" excited states.'),(0,r.kt)("p",null,(0,r.kt)("img",{alt:"exchange_cest_cpmg_figure",src:a(4358).Z,width:"2000",height:"535"})),(0,r.kt)("h2",{id:"about-chemex"},"About ChemEx"),(0,r.kt)("p",null,"ChemEx is an open-source ",(0,r.kt)("a",{parentName:"p",href:"https://www.python.org"},"Python")," application for the\nanalysis of NMR chemical exchange data, whose general idea is to integrate the\nevolution matrix of the spin-system of interest over the pulse sequence and\nextract the best-fit parameters by least-squares optimization. As ChemEx does\nnot rely on analytical equations to fit the data sets, any type of experiments\n(e.g., D-CEST/COS-CEST) or kinetic models (e.g., 3-state exchange model) can be\nsimulated, and most experimental details (e.g., finite pulse width,\noff-resonance effects, etc.) can be taken into account. ChemEx offers a wide\nrange of pulse sequences and kinetic models to choose from. Some of its main\nfeatures include multi-step fits, joint analysis of datasets from different\nexperiments, error analysis, grid search, etc. This documentation provides an\noverview of these different features, along with a description of different\nexamples illustrating their use."),(0,r.kt)("p",null,"ChemEx is a pure ",(0,r.kt)("a",{parentName:"p",href:"https://www.python.org"},"Python")," package that builds upon well\nestablished open-source packages, including ",(0,r.kt)("a",{parentName:"p",href:"https://numpy.org"},"NumPy"),",\n",(0,r.kt)("a",{parentName:"p",href:"https://scipy.org"},"SciPy"),", ",(0,r.kt)("a",{parentName:"p",href:"https://matplotlib.org"},"Matplotlib"),",\n",(0,r.kt)("a",{parentName:"p",href:"https://rich.readthedocs.io/en/stable/"},"Rich"),", and\n",(0,r.kt)("a",{parentName:"p",href:"https://pydantic-docs.helpmanual.io"},"pydantic"),". The minimization procedure is\ncarried out using the ",(0,r.kt)("a",{parentName:"p",href:"https://lmfit.github.io/lmfit-py/"},"LMFIT")," module, which\nsupports many of the optimization algorithms available in\n",(0,r.kt)("a",{parentName:"p",href:"https://scipy.org/"},"SciPy"),"."))}l.isMDXComponent=!0},4358:(e,t,a)=>{a.d(t,{Z:()=>n});const n=a.p+"assets/images/exchange_cest_cpmg_figure-40ebdafcd8b8d21523ec595bbcd00cf5.png"}}]);