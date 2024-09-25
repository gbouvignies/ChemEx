"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[1140],{1236:(e,i,t)=>{t.r(i),t.d(i,{assets:()=>l,contentTitle:()=>d,default:()=>a,frontMatter:()=>r,metadata:()=>o,toc:()=>c});var n=t(4848),s=t(8453);const r={sidebar_position:1},d="Starting a Fit",o={id:"user_guide/fitting/chemex_fit",title:"Starting a Fit",description:"Command",source:"@site/docs/user_guide/fitting/chemex_fit.md",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/chemex_fit",permalink:"/ChemEx/docs/user_guide/fitting/chemex_fit",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/chemex_fit.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_position:1},sidebar:"tutorialSidebar",previous:{title:"Fitting datasets",permalink:"/ChemEx/docs/user_guide/fitting/"},next:{title:"Experiment files",permalink:"/ChemEx/docs/user_guide/fitting/experiment_files"}},l={},c=[{value:"Command",id:"command",level:2},{value:"Options",id:"options",level:2},{value:"Combining Multiple Experiments",id:"combining-multiple-experiments",level:2},{value:"Example",id:"example",level:2}];function h(e){const i={a:"a",admonition:"admonition",code:"code",h1:"h1",h2:"h2",header:"header",p:"p",pre:"pre",table:"table",tbody:"tbody",td:"td",th:"th",thead:"thead",tr:"tr",...(0,s.R)(),...e.components};return(0,n.jsxs)(n.Fragment,{children:[(0,n.jsx)(i.header,{children:(0,n.jsx)(i.h1,{id:"starting-a-fit",children:"Starting a Fit"})}),"\n",(0,n.jsx)(i.h2,{id:"command",children:"Command"}),"\n",(0,n.jsxs)(i.p,{children:["To initiate a fit of NMR chemical exchange datasets, use the ",(0,n.jsx)(i.code,{children:"chemex fit"})," command in the command line. This command can be accompanied by various options to customize the fitting process:"]}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-shell",children:"chemex fit <options>\n"})}),"\n",(0,n.jsx)(i.h2,{id:"options",children:"Options"}),"\n",(0,n.jsxs)(i.p,{children:["The following table lists the options for ",(0,n.jsx)(i.code,{children:"chemex fit"}),". Each option's name is a link to its detailed explanation, and these options are further elaborated in subsequent sections."]}),"\n",(0,n.jsxs)(i.table,{children:[(0,n.jsx)(i.thead,{children:(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.th,{children:"Name"}),(0,n.jsx)(i.th,{children:"Description"})]})}),(0,n.jsxs)(i.tbody,{children:[(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/experiment_files",children:(0,n.jsx)(i.code,{children:"-e"})})," or ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/experiment_files",children:(0,n.jsx)(i.code,{children:"--experiments"})})]}),(0,n.jsx)(i.td,{children:"Specify files containing experimental setup and data."})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/parameter_files",children:(0,n.jsx)(i.code,{children:"-p"})})," or ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/parameter_files",children:(0,n.jsx)(i.code,{children:"--parameters"})})]}),(0,n.jsx)(i.td,{children:"Specify files containing initial fitting parameters."})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/method_files",children:(0,n.jsx)(i.code,{children:"-m"})})," or ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/method_files",children:(0,n.jsx)(i.code,{children:"--methods"})})]}),(0,n.jsx)(i.td,{children:"Indicate the fitting method file (optional)."})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/kinetic_models",children:(0,n.jsx)(i.code,{children:"-d"})})," or ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/kinetic_models",children:(0,n.jsx)(i.code,{children:"--model"})})]}),(0,n.jsxs)(i.td,{children:["Specify the kinetic model for fitting (optional, default: ",(0,n.jsx)(i.code,{children:"2st"}),")."]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsxs)(i.td,{children:[(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/outputs",children:(0,n.jsx)(i.code,{children:"-o"})})," or ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/outputs",children:(0,n.jsx)(i.code,{children:"--output"})})]}),(0,n.jsxs)(i.td,{children:["Set the output directory (optional, default: ",(0,n.jsx)(i.code,{children:"./Output"}),")."]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/method_files#plotting",children:(0,n.jsx)(i.code,{children:"--plot {nothing,normal,all}"})})}),(0,n.jsxs)(i.td,{children:["Select the plotting level (optional, default: ",(0,n.jsx)(i.code,{children:"normal"}),")."]})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.code,{children:"--include"})}),(0,n.jsx)(i.td,{children:"Define residues to include in the fit (optional)."})]}),(0,n.jsxs)(i.tr,{children:[(0,n.jsx)(i.td,{children:(0,n.jsx)(i.code,{children:"--exclude"})}),(0,n.jsx)(i.td,{children:"Define residues to exclude from the fit (optional)."})]})]})]}),"\n",(0,n.jsx)(i.admonition,{type:"note",children:(0,n.jsxs)(i.p,{children:["The ",(0,n.jsx)(i.code,{children:"--experiments"})," and ",(0,n.jsx)(i.code,{children:"--parameters"})," options are mandatory."]})}),"\n",(0,n.jsx)(i.admonition,{type:"tip",children:(0,n.jsxs)(i.p,{children:["For file-input arguments (e.g., ",(0,n.jsx)(i.code,{children:"-p"}),"), you can use the wildcard character (",(0,n.jsx)(i.code,{children:"*"}),") to collectively input multiple files instead of listing each file individually."]})}),"\n",(0,n.jsx)(i.admonition,{title:"TOML File Formats",type:"important",children:(0,n.jsxs)(i.p,{children:["ChemEx uses the ",(0,n.jsx)(i.a,{href:"https://toml.io/",children:"TOML"})," file format for its input and output files. You can find detailed information about this format on the ",(0,n.jsx)(i.a,{href:"https://en.wikipedia.org/wiki/TOML",children:"TOML website"}),"."]})}),"\n",(0,n.jsx)(i.h2,{id:"combining-multiple-experiments",children:"Combining Multiple Experiments"}),"\n",(0,n.jsxs)(i.p,{children:["ChemEx facilitates the combined analysis of multiple experiments. To include various experiments in a fit, add their corresponding ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/user_guide/fitting/experiment_files",children:"experiment files"})," after the ",(0,n.jsx)(i.code,{children:"--experiments"})," (or ",(0,n.jsx)(i.code,{children:"-e"}),") option in the command line. This is particularly useful for fitting different types of CEST or CPMG experiments, or a combination thereof."]}),"\n",(0,n.jsxs)(i.p,{children:["See an ",(0,n.jsx)(i.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/2stBinding/",children:"example of protein-ligand binding analysis"})," using both CPMG and CEST experiments ",(0,n.jsx)(i.a,{href:"/ChemEx/docs/examples/binding",children:"here"}),"."]}),"\n",(0,n.jsx)(i.h2,{id:"example",children:"Example"}),"\n",(0,n.jsx)(i.p,{children:"A typical command for invoking ChemEx is as follows:"}),"\n",(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-shell",children:"chemex fit -e Experiments/*.toml \\\n           -p Parameters/parameters.toml \\\n           -m Methods/method.toml \\\n           -o Output\n"})}),"\n",(0,n.jsxs)(i.p,{children:["The output files are saved in the directory specified by ",(0,n.jsx)(i.code,{children:"-o"}),". The default setting includes the generation of plots to illustrate the best fit lines, which can be adjusted using the ",(0,n.jsx)(i.code,{children:"--plot"})," option."]}),"\n",(0,n.jsxs)(i.admonition,{type:"tip",children:[(0,n.jsxs)(i.p,{children:["For convenience, consider saving this command in a shell script, commonly named ",(0,n.jsx)(i.code,{children:"run.sh"})," in examples:"]}),(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-shell",metastring:"title=run.sh",children:"#!/bin/sh\n\nchemex fit -e Experiments/*.toml \\\n           -p Parameters/parameters.toml \\\n           -m Methods/method.toml \\\n           -o Output\n"})}),(0,n.jsx)(i.p,{children:"To make the script executable, use:"}),(0,n.jsx)(i.pre,{children:(0,n.jsx)(i.code,{className:"language-shell",children:"chmod +x run.sh\n"})})]})]})}function a(e={}){const{wrapper:i}={...(0,s.R)(),...e.components};return i?(0,n.jsx)(i,{...e,children:(0,n.jsx)(h,{...e})}):h(e)}},8453:(e,i,t)=>{t.d(i,{R:()=>d,x:()=>o});var n=t(6540);const s={},r=n.createContext(s);function d(e){const i=n.useContext(r);return n.useMemo((function(){return"function"==typeof e?e(i):{...i,...e}}),[i,e])}function o(e){let i;return i=e.disableParentContext?"function"==typeof e.components?e.components(s):e.components||s:d(e.components),n.createElement(r.Provider,{value:i},e.children)}}}]);