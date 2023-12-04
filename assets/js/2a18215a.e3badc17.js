"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[3549],{5306:(e,i,t)=>{t.r(i),t.d(i,{assets:()=>c,contentTitle:()=>a,default:()=>m,frontMatter:()=>d,metadata:()=>o,toc:()=>h});var n=t(5893),s=t(1151);const r=t.p+"assets/images/cest_26hz_simu-3c69b4b3a3f8c1bd6950dea040dc6ead.png",l=t.p+"assets/images/cpmg_800mhz_simu-2ef3390c127eef53afdd1384f08323d3.png",d={sidebar_position:3},a="Additional modules",o={id:"user_guide/additional_modules",title:"Additional modules",description:"Simulating CPMG and CEST profiles",source:"@site/docs/user_guide/additional_modules.mdx",sourceDirName:"user_guide",slug:"/user_guide/additional_modules",permalink:"/ChemEx/docs/user_guide/additional_modules",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/additional_modules.mdx",tags:[],version:"current",sidebarPosition:3,frontMatter:{sidebar_position:3},sidebar:"tutorialSidebar",previous:{title:"Outputs",permalink:"/ChemEx/docs/user_guide/fitting/outputs"},next:{title:"Experiments",permalink:"/ChemEx/docs/experiments/"}},c={},h=[{value:"Simulating CPMG and CEST profiles",id:"simulating-cpmg-and-cest-profiles",level:2},{value:"Options",id:"options",level:3},{value:"Example",id:"example",level:3},{value:"Plotting best-fit parameters",id:"plotting-best-fit-parameters",level:2},{value:"Options",id:"options-1",level:3},{value:"Example",id:"example-1",level:3},{value:"Getting initial estimates of \u0394\u03d6 for CEST experiments",id:"getting-initial-estimates-of-\u03b4\u03d6-for-cest-experiments",level:2},{value:"Example",id:"example-2",level:3}];function p(e){const i={a:"a",code:"code",h1:"h1",h2:"h2",h3:"h3",p:"p",pre:"pre",table:"table",tbody:"tbody",td:"td",th:"th",thead:"thead",tr:"tr",...(0,s.a)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(i.h1,{id:"additional-modules",children:"Additional modules"}),"\n",(0,n.jsx)(i.h2,{id:"simulating-cpmg-and-cest-profiles",children:"Simulating CPMG and CEST profiles"}),"\n",(0,n.jsx)(i.p,{children:"ChemEx allows simulating CPMG or CEST profiles based on a given set of input\nparameters. Such simulations may be useful to learn about the effects of each\nindividual parameter on the final results."}),"\n",(0,n.jsx)(i.p,{children:"A typical command for simulation purposes with ChemEx is like this:"}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-bash",children:"\n   chemex simulate -e <FILE> \\\n                   -p <FILE> \\\n                   -d <MODEL> \\\n                   -o <DIR>\n\n"})}),"\n",(0,n.jsx)(i.p,{children:"Example simulation results for CPMG and CEST experiments are shown below:"}),"\n","\n","\n",(0,n.jsxs)("figure",{children:[(0,n.jsx)("img",{src:r,alt:"CEST profile",width:"50%"}),(0,n.jsx)("img",{src:l,alt:"CPMG profile",width:"50%"}),(0,n.jsx)("figcaption",{align:"center",children:(0,n.jsx)("b",{children:"Examples of CEST and CPMG simulation results"})})]}),"\n",(0,n.jsx)(i.h3,{id:"options",children:"Options"}),"\n",(0,n.jsxs)(i.table,{children:[(0,n.jsx)(i.thead,{children:(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.th,{children:"Option"}),(0,n.jsx)(i.th,{children:"Description"})]})}),(0,n.jsxs)(i.tbody,{children:[(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-e"}),", ",(0,n.jsx)(i.code,{children:"--experiments"})]}),(0,n.jsx)(i.td,{children:"Specify the files containing experimental setup and data location"})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-p"}),", ",(0,n.jsx)(i.code,{children:"--parameters"})]}),(0,n.jsx)(i.td,{children:"Specify the files containing the initial values of fitting parameters"})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-d"}),", ",(0,n.jsx)(i.code,{children:"--model"})]}),(0,n.jsxs)(i.td,{children:["Specify the exchange model used to fit the datasets (default: ",(0,n.jsx)(i.code,{children:"2st"}),")"]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-o"}),", ",(0,n.jsx)(i.code,{children:"--output"})]}),(0,n.jsxs)(i.td,{children:["Specify the output directory (default: ",(0,n.jsx)(i.code,{children:"./Output"}),")"]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.code,{children:"--plot {nothing,normal}"})}),(0,n.jsxs)(i.td,{children:["Plotting level (default: ",(0,n.jsx)(i.code,{children:"normal"}),")"]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.code,{children:"--include"})}),(0,n.jsx)(i.td,{children:"Residue(s) to include in the fit (optional)"})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.code,{children:"--exclude"})}),(0,n.jsx)(i.td,{children:"Residue(s) to exclude from the fit (optional)"})]})]})]}),"\n",(0,n.jsx)(i.h3,{id:"example",children:"Example"}),"\n",(0,n.jsxs)(i.p,{children:["An example use of the module is available\n",(0,n.jsx)(i.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N/",children:"here"}),"\nwith the script ",(0,n.jsx)(i.code,{children:"simulate.sh"}),"."]}),"\n",(0,n.jsx)(i.h2,{id:"plotting-best-fit-parameters",children:"Plotting best-fit parameters"}),"\n",(0,n.jsxs)(i.p,{children:["ChemEx comes with a module ",(0,n.jsx)(i.code,{children:"plot_param"})," that allows visualizing the fitting\nresults interactively."]}),"\n",(0,n.jsx)(i.h3,{id:"options-1",children:"Options"}),"\n",(0,n.jsxs)(i.table,{children:[(0,n.jsx)(i.thead,{children:(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.th,{children:"Option"}),(0,n.jsx)(i.th,{children:"Description"})]})}),(0,n.jsxs)(i.tbody,{children:[(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-p"}),", ",(0,n.jsx)(i.code,{children:"--parameters"})]}),(0,n.jsx)(i.td,{children:"Specify the files containing the fitted parameters to be plotted"})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.code,{children:"-n"}),", ",(0,n.jsx)(i.code,{children:"--parname"})]}),(0,n.jsx)(i.td,{children:"Specify the name of the parameter to plot"})]})]})]}),"\n",(0,n.jsx)(i.h3,{id:"example-1",children:"Example"}),"\n",(0,n.jsxs)(i.p,{children:["An example use of the module is given in the\n",(0,n.jsx)(i.a,{href:"/ChemEx/docs/examples/binding",children:"protein-ligand binding example"}),". After finish running\n",(0,n.jsx)(i.code,{children:"run.sh"}),", the chemical shift differences between the free and bound states can\nbe displayed with:"]}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-shell",children:"chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n DW_AB\n"})}),"\n",(0,n.jsx)(i.p,{children:"and the transverse relaxation rates of both states can be compared with:"}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-bash",children:"chemex plot_param -p Output/STEP2/All/Parameters/fitted.toml -n R2\n"})}),"\n",(0,n.jsxs)(i.p,{children:["These two commands are saved in the ",(0,n.jsx)(i.code,{children:"plot_param.sh"})," script in\n",(0,n.jsx)(i.a,{href:"/ChemEx/docs/examples/binding",children:"this example"}),". From these two observables, the core\nregion of the interaction site can be clearly located. Aside from the core\nregion, there is also a tail with increased R",(0,n.jsx)("sub",{children:"2"})," rates located at\nC-terminal end of the interaction site and with very little chemical shift\nperturbation. This region is likely involved in the transient interactions with\nthe binding partner, which causes certain degree of steric restriction to this\nregion."]}),"\n",(0,n.jsx)(i.h2,{id:"getting-initial-estimates-of-\u03b4\u03d6-for-cest-experiments",children:"Getting initial estimates of \u0394\u03d6 for CEST experiments"}),"\n",(0,n.jsxs)(i.p,{children:["In CEST (and also D-CEST/COS-CEST) experiments, it is necessary to choose\nsuitable initial value of \u0394\u03d6 to avoid getting trapped in a local minimum. ChemEx\ncomes with a module ",(0,n.jsx)(i.code,{children:"pick_cest"})," for manually picking the major and minor dips of\nCEST profiles, which correspond to the ground and excited states, respectively.\nA typical command for such purpose is like this:"]}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-bash",children:"chemex pick_cest -e <FILE> -o <DIR>\n"})}),"\n",(0,n.jsxs)(i.p,{children:["After typing this command, a window showing all CEST profiles will appear. For\neach profile first click on the major dip and then on the minor dip(s). Note\nthat in certain profiles only one dip could be visible, which indicates the\nminor dip is overlapped with the major dip, therefore the major dip should be\nclicked twice. When done with any profile, click the ",(0,n.jsx)(i.code,{children:"Next"})," or ",(0,n.jsx)(i.code,{children:"Previous"})," button\nto proceed to the next or previous profile. The ",(0,n.jsx)(i.code,{children:"Swap"})," button allows switching\nbetween the major and minor states. The ",(0,n.jsx)(i.code,{children:"Clear"})," button allows cleaning the\nselection in the current profile. Two separate files will be created in\nreal-time during the dip picking process: ",(0,n.jsx)(i.code,{children:"cs_a.toml"})," and ",(0,n.jsx)(i.code,{children:"dw_ab.toml"})," that\ncontain chemical shifts of the major state and chemical shift difference between\nthe major and minor states, respectively."]}),"\n",(0,n.jsx)(i.h3,{id:"example-2",children:"Example"}),"\n",(0,n.jsxs)(i.p,{children:["Try to run the ",(0,n.jsx)(i.code,{children:"pick_cest.sh"})," script under\n",(0,n.jsxs)(i.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/CEST_15N/",children:[(0,n.jsx)(i.code,{children:"CEST_15N/"})," example"]}),"\nand ",(0,n.jsx)(i.code,{children:"pick_dcest.sh"})," script under\n",(0,n.jsxs)(i.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N/",children:[(0,n.jsx)(i.code,{children:"DCEST_15N/"})," example"]}),"\nto learn how to make use of this function for CEST and D-CEST experiments,\nrespectively."]})]})}function m(e={}){const{wrapper:i}={...(0,s.a)(),...e.components};return i?(0,n.jsx)(i,{...e,children:(0,n.jsx)(p,{...e})}):p(e)}},1151:(e,i,t)=>{t.d(i,{Z:()=>d,a:()=>l});var n=t(7294);const s={},r=n.createContext(s);function l(e){const i=n.useContext(r);return n.useMemo((function(){return"function"==typeof e?e(i):{...i,...e}}),[i,e])}function d(e){let i;return i=e.disableParentContext?"function"==typeof e.components?e.components(s):e.components||s:l(e.components),n.createElement(r.Provider,{value:i},e.children)}}}]);