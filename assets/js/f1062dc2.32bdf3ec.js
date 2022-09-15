"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[9572],{3905:(e,t,n)=>{n.d(t,{Zo:()=>d,kt:()=>c});var a=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function r(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function l(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?r(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):r(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function o(e,t){if(null==e)return{};var n,a,i=function(e,t){if(null==e)return{};var n,a,i={},r=Object.keys(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);for(a=0;a<r.length;a++)n=r[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var p=a.createContext({}),s=function(e){var t=a.useContext(p),n=t;return e&&(n="function"==typeof e?e(t):l(l({},t),e)),n},d=function(e){var t=s(e.components);return a.createElement(p.Provider,{value:t},e.children)},m={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},u=a.forwardRef((function(e,t){var n=e.components,i=e.mdxType,r=e.originalType,p=e.parentName,d=o(e,["components","mdxType","originalType","parentName"]),u=s(n),c=i,h=u["".concat(p,".").concat(c)]||u[c]||m[c]||r;return n?a.createElement(h,l(l({ref:t},d),{},{components:n})):a.createElement(h,l({ref:t},d))}));function c(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var r=n.length,l=new Array(r);l[0]=u;var o={};for(var p in t)hasOwnProperty.call(t,p)&&(o[p]=t[p]);o.originalType=e,o.mdxType="string"==typeof e?e:i,l[1]=o;for(var s=2;s<r;s++)l[s]=n[s];return a.createElement.apply(null,l)}return a.createElement.apply(null,n)}u.displayName="MDXCreateElement"},7983:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>d,contentTitle:()=>p,default:()=>c,frontMatter:()=>o,metadata:()=>s,toc:()=>m});var a=n(7462),i=(n(7294),n(3905));const r=n.p+"assets/images/cest_26hz_simu-3c69b4b3a3f8c1bd6950dea040dc6ead.png",l=n.p+"assets/images/cpmg_800mhz_simu-2ef3390c127eef53afdd1384f08323d3.png",o={sidebar_position:4},p="Additional modules",s={unversionedId:"user_guide/additional_modules",id:"user_guide/additional_modules",title:"Additional modules",description:"Simulating CPMG and CEST profiles",source:"@site/docs/user_guide/additional_modules.md",sourceDirName:"user_guide",slug:"/user_guide/additional_modules",permalink:"/ChemEx/docs/user_guide/additional_modules",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/additional_modules.md",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_position:4},sidebar:"tutorialSidebar",previous:{title:"Outputs",permalink:"/ChemEx/docs/user_guide/fitting/outputs"},next:{title:"Experiments",permalink:"/ChemEx/docs/experiments/"}},d={},m=[{value:"Simulating CPMG and CEST profiles",id:"simulating-cpmg-and-cest-profiles",level:2},{value:"Options",id:"options",level:3},{value:"Example",id:"example",level:3},{value:"Plotting best-fit parameters",id:"plotting-best-fit-parameters",level:2},{value:"Options",id:"options-1",level:3},{value:"Example",id:"example-1",level:3},{value:"Getting initial estimates of \u0394\u03d6 for CEST experiments",id:"getting-initial-estimates-of-\u03b4\u03d6-for-cest-experiments",level:2},{value:"Example",id:"example-2",level:3}],u={toc:m};function c(e){let{components:t,...n}=e;return(0,i.kt)("wrapper",(0,a.Z)({},u,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"additional-modules"},"Additional modules"),(0,i.kt)("h2",{id:"simulating-cpmg-and-cest-profiles"},"Simulating CPMG and CEST profiles"),(0,i.kt)("p",null,"ChemEx allows simulating CPMG or CEST profiles based on a given set of input\nparameters. Such simulations may be useful to learn about the effects of each\nindividual parameter on the final results."),(0,i.kt)("p",null,"A typical command for simulation purposes with ChemEx is like this:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-bash"},"\n   chemex simulate -e <FILE> \\\n                   -p <FILE> \\\n                   -d <MODEL> \\\n                   -o <DIR>\n\n")),(0,i.kt)("p",null,"Example simulation results for CPMG and CEST experiments are shown below:"),(0,i.kt)("figure",null,(0,i.kt)("img",{src:r,alt:"CEST profile",width:"50%"}),(0,i.kt)("img",{src:l,alt:"CPMG profile",width:"50%"}),(0,i.kt)("figcaption",{align:"center"},(0,i.kt)("b",null,"Examples of CEST and CPMG simulation results"))),(0,i.kt)("h3",{id:"options"},"Options"),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:null},"Option"),(0,i.kt)("th",{parentName:"tr",align:null},"Description"))),(0,i.kt)("tbody",{parentName:"table"},(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-e"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--experiments")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the files containing experimental setup and data location")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-p"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--parameters")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the files containing the initial values of fitting parameters")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-d"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--model")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the exchange model used to fit the datasets (default: ",(0,i.kt)("inlineCode",{parentName:"td"},"2st"),")")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-o"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--output")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the output directory (default: ",(0,i.kt)("inlineCode",{parentName:"td"},"./Output"),")")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"--plot {nothing,normal}")),(0,i.kt)("td",{parentName:"tr",align:null},"Plotting level (default: ",(0,i.kt)("inlineCode",{parentName:"td"},"normal"),")")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"--include")),(0,i.kt)("td",{parentName:"tr",align:null},"Residue(s) to include in the fit (optional)")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"--exclude")),(0,i.kt)("td",{parentName:"tr",align:null},"Residue(s) to exclude from the fit (optional)")))),(0,i.kt)("h3",{id:"example"},"Example"),(0,i.kt)("p",null,"An example use of the module is available\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N/"},"here"),"\nwith the script ",(0,i.kt)("inlineCode",{parentName:"p"},"simulate.sh"),"."),(0,i.kt)("h2",{id:"plotting-best-fit-parameters"},"Plotting best-fit parameters"),(0,i.kt)("p",null,"ChemEx comes with a module ",(0,i.kt)("inlineCode",{parentName:"p"},"plot_param")," that allows visualizing the fitting\nresults interactively."),(0,i.kt)("h3",{id:"options-1"},"Options"),(0,i.kt)("table",null,(0,i.kt)("thead",{parentName:"table"},(0,i.kt)("tr",{parentName:"thead"},(0,i.kt)("th",{parentName:"tr",align:null},"Option"),(0,i.kt)("th",{parentName:"tr",align:null},"Description"))),(0,i.kt)("tbody",{parentName:"table"},(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-p"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--parameters")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the files containing the fitted parameters to be plotted")),(0,i.kt)("tr",{parentName:"tbody"},(0,i.kt)("td",{parentName:"tr",align:null},(0,i.kt)("inlineCode",{parentName:"td"},"-n"),", ",(0,i.kt)("inlineCode",{parentName:"td"},"--parname")),(0,i.kt)("td",{parentName:"tr",align:null},"Specify the name of the parameter to plot")))),(0,i.kt)("h3",{id:"example-1"},"Example"),(0,i.kt)("p",null,"An example use of the module is given in the\n",(0,i.kt)("a",{parentName:"p",href:"/ChemEx/docs/examples/binding"},"protein-ligand binding example"),". After finish running\n",(0,i.kt)("inlineCode",{parentName:"p"},"run.sh"),", the chemical shift differences between the free and bound states can\nbe displayed with:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n DW_AB\n")),(0,i.kt)("p",null,"and the transverse relaxation rates of both states can be compared with:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-bash"},"chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n R2\n")),(0,i.kt)("p",null,"These two commands are saved in the ",(0,i.kt)("inlineCode",{parentName:"p"},"plot_param.sh")," script in\n",(0,i.kt)("a",{parentName:"p",href:"/ChemEx/docs/examples/binding"},"this example"),". From these two observables, the core\nregion of the interaction site can be clearly located. Aside from the core\nregion, there is also a tail with increased R",(0,i.kt)("sub",null,"2")," rates located at\nC-terminal end of the interaction site and with very little chemical shift\nperturbation. This region is likely involved in the transient interactions with\nthe binding partner, which causes certain degree of steric restriction to this\nregion."),(0,i.kt)("h2",{id:"getting-initial-estimates-of-\u03b4\u03d6-for-cest-experiments"},"Getting initial estimates of \u0394\u03d6 for CEST experiments"),(0,i.kt)("p",null,"In CEST (and also D-CEST/COS-CEST) experiments, it is necessary to choose\nsuitable initial value of \u0394\u03d6 to avoid getting trapped in a local minimum. ChemEx\ncomes with a module ",(0,i.kt)("inlineCode",{parentName:"p"},"pick_cest")," for manually picking the major and minor dips of\nCEST profiles, which correspond to the ground and excited states, respectively.\nA typical command for such purpose is like this:"),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-bash"},"chemex pick_cest -e <FILE> -o <DIR>\n")),(0,i.kt)("p",null,"After typing this command, a window showing all CEST profiles will appear. For\neach profile first click on the major dip and then on the minor dip(s). Note\nthat in certain profiles only one dip could be visible, which indicates the\nminor dip is overlapped with the major dip, therefore the major dip should be\nclicked twice. When done with any profile, click the ",(0,i.kt)("inlineCode",{parentName:"p"},"Next")," or ",(0,i.kt)("inlineCode",{parentName:"p"},"Previous")," button\nto proceed to the next or previous profile. The ",(0,i.kt)("inlineCode",{parentName:"p"},"Swap")," button allows switching\nbetween the major and minor states. The ",(0,i.kt)("inlineCode",{parentName:"p"},"Clear")," button allows cleaning the\nselection in the current profile. Two separate files will be created in\nreal-time during the dip picking process: ",(0,i.kt)("inlineCode",{parentName:"p"},"cs_a.toml")," and ",(0,i.kt)("inlineCode",{parentName:"p"},"dw_ab.toml")," that\ncontain chemical shifts of the major state and chemical shift difference between\nthe major and minor states, respectively."),(0,i.kt)("h3",{id:"example-2"},"Example"),(0,i.kt)("p",null,"Try to run the ",(0,i.kt)("inlineCode",{parentName:"p"},"pick_cest.sh")," script under\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N/"},(0,i.kt)("inlineCode",{parentName:"a"},"CEST_15N/")," example"),"\nand ",(0,i.kt)("inlineCode",{parentName:"p"},"pick_dcest.sh")," script under\n",(0,i.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N/"},(0,i.kt)("inlineCode",{parentName:"a"},"DCEST_15N/")," example"),"\nto learn how to make use of this function for CEST and D-CEST experiments,\nrespectively."))}c.isMDXComponent=!0}}]);