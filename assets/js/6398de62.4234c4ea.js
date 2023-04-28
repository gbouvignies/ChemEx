"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[2021],{3905:(e,t,n)=>{n.d(t,{Zo:()=>d,kt:()=>c});var i=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function r(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);t&&(i=i.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,i)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?r(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):r(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,i,a=function(e,t){if(null==e)return{};var n,i,a={},r=Object.keys(e);for(i=0;i<r.length;i++)n=r[i],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(i=0;i<r.length;i++)n=r[i],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var p=i.createContext({}),m=function(e){var t=i.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},d=function(e){var t=m(e.components);return i.createElement(p.Provider,{value:t},e.children)},s={inlineCode:"code",wrapper:function(e){var t=e.children;return i.createElement(i.Fragment,{},t)}},u=i.forwardRef((function(e,t){var n=e.components,a=e.mdxType,r=e.originalType,p=e.parentName,d=l(e,["components","mdxType","originalType","parentName"]),u=m(n),c=a,f=u["".concat(p,".").concat(c)]||u[c]||s[c]||r;return n?i.createElement(f,o(o({ref:t},d),{},{components:n})):i.createElement(f,o({ref:t},d))}));function c(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var r=n.length,o=new Array(r);o[0]=u;var l={};for(var p in t)hasOwnProperty.call(t,p)&&(l[p]=t[p]);l.originalType=e,l.mdxType="string"==typeof e?e:a,o[1]=l;for(var m=2;m<r;m++)o[m]=n[m];return i.createElement.apply(null,o)}return i.createElement.apply(null,n)}u.displayName="MDXCreateElement"},5931:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>p,contentTitle:()=>o,default:()=>s,frontMatter:()=>r,metadata:()=>l,toc:()=>m});var i=n(7462),a=(n(7294),n(3905));const r={sidebar_position:1},o="Starting a fit",l={unversionedId:"user_guide/fitting/chemex_fit",id:"user_guide/fitting/chemex_fit",title:"Starting a fit",description:"Command",source:"@site/docs/user_guide/fitting/chemex_fit.md",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/chemex_fit",permalink:"/ChemEx/docs/user_guide/fitting/chemex_fit",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/chemex_fit.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_position:1},sidebar:"tutorialSidebar",previous:{title:"Fitting datasets",permalink:"/ChemEx/docs/user_guide/fitting/"},next:{title:"Experiment files",permalink:"/ChemEx/docs/user_guide/fitting/experiment_files"}},p={},m=[{value:"Command",id:"command",level:2},{value:"Options",id:"options",level:2},{value:"Combining multiple experiments",id:"combining-multiple-experiments",level:2},{value:"Example",id:"example",level:2}],d={toc:m};function s(e){let{components:t,...n}=e;return(0,a.kt)("wrapper",(0,i.Z)({},d,n,{components:t,mdxType:"MDXLayout"}),(0,a.kt)("h1",{id:"starting-a-fit"},"Starting a fit"),(0,a.kt)("h2",{id:"command"},"Command"),(0,a.kt)("p",null,"To initiate a fit of NMR chemical exchange dataset(s), run the command\n",(0,a.kt)("inlineCode",{parentName:"p"},"chemex fit")," from the command-line followed by a set of options to set and\ncontrol the fitting process."),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-shell"},"chemex fit <options>\n")),(0,a.kt)("h2",{id:"options"},"Options"),(0,a.kt)("p",null,"The list of available options is given in the table below. More details about\nthe options and their associated files are given in following sections. Some of\nthe option names of the table are clickable, so that you can reach the related\nsection directly. Let's dive in."),(0,a.kt)("table",null,(0,a.kt)("thead",{parentName:"table"},(0,a.kt)("tr",{parentName:"thead"},(0,a.kt)("th",{parentName:"tr",align:null},"Name"),(0,a.kt)("th",{parentName:"tr",align:null},"Description"))),(0,a.kt)("tbody",{parentName:"table"},(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/experiment_files"},(0,a.kt)("inlineCode",{parentName:"a"},"-e"))," or ",(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/experiment_files"},(0,a.kt)("inlineCode",{parentName:"a"},"--experiments"))),(0,a.kt)("td",{parentName:"tr",align:null},"Specify the files containing the experimental setup and data location.")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/parameter_files"},(0,a.kt)("inlineCode",{parentName:"a"},"-p"))," or ",(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/parameter_files"},(0,a.kt)("inlineCode",{parentName:"a"},"--parameters"))),(0,a.kt)("td",{parentName:"tr",align:null},"Specify the files containing the initial values of fitting parameters.")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/method_files"},(0,a.kt)("inlineCode",{parentName:"a"},"-m"))," or ",(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/method_files"},(0,a.kt)("inlineCode",{parentName:"a"},"--methods"))),(0,a.kt)("td",{parentName:"tr",align:null},"Specify the file containing the fitting method (optional).")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/kinetic_models"},(0,a.kt)("inlineCode",{parentName:"a"},"-d"))," or ",(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/kinetic_models"},(0,a.kt)("inlineCode",{parentName:"a"},"--model"))),(0,a.kt)("td",{parentName:"tr",align:null},"Specify the kinetic model used to fit the datasets (optional, default: ",(0,a.kt)("inlineCode",{parentName:"td"},"2st"),").")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/outputs"},(0,a.kt)("inlineCode",{parentName:"a"},"-o"))," or ",(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/outputs"},(0,a.kt)("inlineCode",{parentName:"a"},"--output"))),(0,a.kt)("td",{parentName:"tr",align:null},"Specify the output directory (optional, default: ",(0,a.kt)("inlineCode",{parentName:"td"},"./Output"),").")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("a",{parentName:"td",href:"/ChemEx/docs/user_guide/fitting/method_files#plotting"},(0,a.kt)("inlineCode",{parentName:"a"},"--plot {nothing,normal,all}"))),(0,a.kt)("td",{parentName:"tr",align:null},"Plotting level (optional, default: ",(0,a.kt)("inlineCode",{parentName:"td"},"normal"),").")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("inlineCode",{parentName:"td"},"--include")),(0,a.kt)("td",{parentName:"tr",align:null},"Residue(s) to include in the fit (optional).")),(0,a.kt)("tr",{parentName:"tbody"},(0,a.kt)("td",{parentName:"tr",align:null},(0,a.kt)("inlineCode",{parentName:"td"},"--exclude")),(0,a.kt)("td",{parentName:"tr",align:null},"Residue(s) to exclude from the fit (optional).")))),(0,a.kt)("admonition",{type:"note"},(0,a.kt)("p",{parentName:"admonition"},(0,a.kt)("inlineCode",{parentName:"p"},"--experiments")," and ",(0,a.kt)("inlineCode",{parentName:"p"},"--parameters")," options are mandatory.")),(0,a.kt)("admonition",{type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"For arguments requiring input files (e.g., ",(0,a.kt)("inlineCode",{parentName:"p"},"-p"),"), it is possible to provide\nmultiple files by using the wildcard character (i.e., ",(0,a.kt)("inlineCode",{parentName:"p"},"*"),"), instead of\nspecifying the name of each file individually.")),(0,a.kt)("admonition",{title:"TOML File formats",type:"important"},(0,a.kt)("p",{parentName:"admonition"},"The input and output files of ChemEx use the\n",(0,a.kt)("a",{parentName:"p",href:"https://en.wikipedia.org/wiki/TOML"},"TOML")," file format. You can find a detailed\ndescription of this file format on the ",(0,a.kt)("a",{parentName:"p",href:"https://toml.io/"},"TOML website"),".")),(0,a.kt)("h2",{id:"combining-multiple-experiments"},"Combining multiple experiments"),(0,a.kt)("p",null,"ChemEx offers a very simple way to jointly analyse multiple experiments. For\neach experiment, simply add the corresponding\n",(0,a.kt)("a",{parentName:"p",href:"/ChemEx/docs/user_guide/fitting/experiment_files"},"experiment file")," right after the ",(0,a.kt)("inlineCode",{parentName:"p"},"--experiments")," (or ",(0,a.kt)("inlineCode",{parentName:"p"},"-e"),")\noption in the command-line."),(0,a.kt)("p",null,"This is useful, for example, when you want to fit CEST experiments with\ndifferent B",(0,a.kt)("sub",null,"1")," field together, or CPMG relaxation dispersion\nexperiments recorded at different B",(0,a.kt)("sub",null,"0")," field together, or even a\ncombination of the two."),(0,a.kt)("p",null,"An example illustrating\n",(0,a.kt)("a",{parentName:"p",href:"/ChemEx/docs/examples/binding"},"the study of protein-ligand binding with CPMG and CEST experiments"),"\nis available\n",(0,a.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding/"},"here"),"."),(0,a.kt)("h2",{id:"example"},"Example"),(0,a.kt)("p",null,"A typical command line invocation would be:"),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-shell"},"chemex fit -e Experiments/*.toml \\\n           -p Parameters/parameters.toml \\\n           -m Methods/method.toml \\\n           -o Output\n")),(0,a.kt)("p",null,"Output files are located in the directory specified with the option ",(0,a.kt)("inlineCode",{parentName:"p"},"-o"),". Plots\nshowing the line of best fit are also produced and saved to the disk by default\n\u2013 you can modulate this behaviour with the option ",(0,a.kt)("inlineCode",{parentName:"p"},"--plot"),"."),(0,a.kt)("admonition",{type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"You can save the command in a shell script (typically called ",(0,a.kt)("inlineCode",{parentName:"p"},"run.sh")," in most\nexamples) to save some typing efforts."),(0,a.kt)("pre",{parentName:"admonition"},(0,a.kt)("code",{parentName:"pre",className:"language-shell",metastring:"title=run.sh",title:"run.sh"},"#!/bin/sh\n\nchemex fit -e Experiments/*.toml \\\n           -p Parameters/parameters.toml \\\n           -m Methods/method.toml \\\n           -o Output\n")),(0,a.kt)("p",{parentName:"admonition"},"Make it executable with the following command:"),(0,a.kt)("pre",{parentName:"admonition"},(0,a.kt)("code",{parentName:"pre",className:"language-shell"},"chmod +x run.sh\n"))))}s.isMDXComponent=!0}}]);