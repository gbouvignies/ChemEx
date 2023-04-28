"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[9676],{3905:(e,t,n)=>{n.d(t,{Zo:()=>d,kt:()=>u});var i=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function a(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);t&&(i=i.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,i)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?a(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function s(e,t){if(null==e)return{};var n,i,r=function(e,t){if(null==e)return{};var n,i,r={},a=Object.keys(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var l=i.createContext({}),p=function(e){var t=i.useContext(l),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},d=function(e){var t=p(e.components);return i.createElement(l.Provider,{value:t},e.children)},c={inlineCode:"code",wrapper:function(e){var t=e.children;return i.createElement(i.Fragment,{},t)}},m=i.forwardRef((function(e,t){var n=e.components,r=e.mdxType,a=e.originalType,l=e.parentName,d=s(e,["components","mdxType","originalType","parentName"]),m=p(n),u=r,f=m["".concat(l,".").concat(u)]||m[u]||c[u]||a;return n?i.createElement(f,o(o({ref:t},d),{},{components:n})):i.createElement(f,o({ref:t},d))}));function u(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var a=n.length,o=new Array(a);o[0]=m;var s={};for(var l in t)hasOwnProperty.call(t,l)&&(s[l]=t[l]);s.originalType=e,s.mdxType="string"==typeof e?e:r,o[1]=s;for(var p=2;p<a;p++)o[p]=n[p];return i.createElement.apply(null,o)}return i.createElement.apply(null,n)}m.displayName="MDXCreateElement"},5843:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>l,contentTitle:()=>o,default:()=>c,frontMatter:()=>a,metadata:()=>s,toc:()=>p});var i=n(7462),r=(n(7294),n(3905));const a={sidebar_position:9,sidebar_label:"Protein-ligand binding study"},o="Protein-ligand binding study",s={unversionedId:"examples/binding",id:"examples/binding",title:"Protein-ligand binding study",description:"You can get",source:"@site/docs/examples/binding.md",sourceDirName:"examples",slug:"/examples/binding",permalink:"/ChemEx/docs/examples/binding",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/examples/binding.md",tags:[],version:"current",sidebarPosition:9,frontMatter:{sidebar_position:9,sidebar_label:"Protein-ligand binding study"},sidebar:"tutorialSidebar",previous:{title:"\xb9H\u1d3a COS-CEST experiment",permalink:"/ChemEx/docs/examples/coscest"},next:{title:"3-state exchange kinetic model",permalink:"/ChemEx/docs/examples/state"}},l={},p=[],d={toc:p};function c(e){let{components:t,...n}=e;return(0,r.kt)("wrapper",(0,i.Z)({},d,n,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("h1",{id:"protein-ligand-binding-study"},"Protein-ligand binding study"),(0,r.kt)("admonition",{title:"Files",type:"note"},(0,r.kt)("p",{parentName:"admonition"},"You can get\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding"},"the example files from the GitHub ChemEx page"),".")),(0,r.kt)("p",null,"This example demonstrates the study of protein-ligand binding with CPMG and CEST\nexperiments. In this example, CPMG and CEST experiments have been carried out on\nsamples with different protein/ligand molar ratios. Since the amount of ligand\nis much less than protein in each sample, the protein-ligand bound state is only\nsparsely populated, therefore CPMG and CEST can be used to probe the bound state\neven if the resonances are not visible in the spectra."),(0,r.kt)("p",null,"Here, it is assumed that the same binding constant K",(0,r.kt)("sub",null,"d")," is shared among\ndifferent samples, which are distinguished by different ",(0,r.kt)("inlineCode",{parentName:"p"},"p_total")," and ",(0,r.kt)("inlineCode",{parentName:"p"},"l_total"),"\nkeys in experiment files. The ",(0,r.kt)("inlineCode",{parentName:"p"},"2st_binding")," kinetic model is used to fit all the\ndatasets with a two-step fitting scheme. In the first step only selective\nresidues near the protein-ligand interaction sites are fitted, and the goal is\nto get a good estimate of k",(0,r.kt)("sub",null,"off"),"; in the second step all residues are\nincluded while k",(0,r.kt)("sub",null,"off")," is fixed to the estimated value in the previous\nstep, the major purpose is to obtain residue-specific parameters."),(0,r.kt)("admonition",{type:"info"},(0,r.kt)("p",{parentName:"admonition"},"Further analysis of the fitting results can be carried out with the aid of\nadditional functions in ChemEx (refer to\n",(0,r.kt)("a",{parentName:"p",href:"/ChemEx/docs/user_guide/additional_modules#plotting-best-fit-parameters"},(0,r.kt)("inlineCode",{parentName:"a"},"Plotting best-fit parameters")),"\nsubsection). During the whole fitting process K",(0,r.kt)("sub",null,"d")," is fixed to the\nvalue obtained from ITC experiments, more details about this example can be\nfound in the reference.",(0,r.kt)("sup",{parentName:"p",id:"fnref-1"},(0,r.kt)("a",{parentName:"sup",href:"#fn-1",className:"footnote-ref"},"1")))),(0,r.kt)("div",{className:"footnotes"},(0,r.kt)("hr",{parentName:"div"}),(0,r.kt)("ol",{parentName:"div"},(0,r.kt)("li",{parentName:"ol",id:"fn-1"},"C. Charlier, G. Bouvignies, P. Pelupessy, A. Walrant, R. Marquant, M.\nKozlov, P. De Ioannes, N. Bolik-Coulon, S. Sagan, P. Cortes, A. K. Aggarwal,\nL. Carlier, and F. Ferrage. Structure and Dynamics of an Intrinsically\nDisordered Protein Region That Partially Folds upon Binding by\nChemical-Exchange NMR ",(0,r.kt)("em",{parentName:"li"},"J. Am. Chem. Soc."),", ",(0,r.kt)("em",{parentName:"li"},"139"),", 12219-12227 (2017).\n",(0,r.kt)("a",{parentName:"li",href:"https://doi.org/10.1021/jacs.7b05823"},"https://doi.org/10.1021/jacs.7b05823"),(0,r.kt)("a",{parentName:"li",href:"#fnref-1",className:"footnote-backref"},"\u21a9")))))}c.isMDXComponent=!0}}]);