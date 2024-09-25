"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[3467],{1018:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>i,default:()=>m,frontMatter:()=>r,metadata:()=>a,toc:()=>d});var s=n(4848),o=n(8453);const r={sidebar_position:4,sidebar_label:"Measuring Excited state RDCs with CPMG"},i="RDC measurement for excited state with \xb9\u2075N TROSY CPMG",a={id:"examples/trosy_cpmg_rdc",title:"RDC measurement for excited state with \xb9\u2075N TROSY CPMG",description:"You can get",source:"@site/docs/examples/trosy_cpmg_rdc.md",sourceDirName:"examples",slug:"/examples/trosy_cpmg_rdc",permalink:"/ChemEx/docs/examples/trosy_cpmg_rdc",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/examples/trosy_cpmg_rdc.md",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_position:4,sidebar_label:"Measuring Excited state RDCs with CPMG"},sidebar:"tutorialSidebar",previous:{title:"Measuring diffusion constants of invisible protein conformers",permalink:"/ChemEx/docs/examples/diffusion_tqcpmg"},next:{title:"CEST experiments for \xb9\xb3C, \xb9\u2075N-labeled samples",permalink:"/ChemEx/docs/examples/cest_13c_15n"}},c={},d=[];function l(e){const t={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",header:"header",li:"li",ol:"ol",p:"p",section:"section",strong:"strong",sup:"sup",...(0,o.R)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(t.header,{children:(0,s.jsx)(t.h1,{id:"rdc-measurement-for-excited-state-with-n-trosy-cpmg",children:"RDC measurement for excited state with \xb9\u2075N TROSY CPMG"})}),"\n",(0,s.jsx)(t.admonition,{title:"Files",type:"note",children:(0,s.jsxs)(t.p,{children:["You can get\n",(0,s.jsx)(t.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/N15_NH_RDC",children:"the example files from the GitHub ChemEx page"}),"."]})}),"\n",(0,s.jsxs)(t.p,{children:["This example demonstrates the use of CPMG experiment to measure RDC parameters\nfor the excited state. In \xb9\u2075N TROSY CPMG experiment\n(",(0,s.jsx)(t.a,{href:"../experiments/cpmg/cpmg_15n_tr",children:(0,s.jsx)(t.code,{children:"cpmg_15n_tr"})}),"), by analyzing datasets measured\nfor TROSY and anti-TROSY components, the information about \xb9\u2075N \u0394\u03d6 for both\ncomponents can be obtained, therefore RDC parameters for excited state can be\nderived based on RDC parameters for ground state together with \u0394\u03d6 of these two\ncomponents ",(0,s.jsx)(t.sup,{children:(0,s.jsx)(t.a,{href:"#user-content-fn-1",id:"user-content-fnref-1","data-footnote-ref":!0,"aria-describedby":"footnote-label",children:"1"})}),"."]}),"\n",(0,s.jsxs)(t.p,{children:["Note that ",(0,s.jsx)(t.code,{children:"antitrosy"})," key should be set properly in experiment files to indicate\nwhether the datasets are measured for TROSY or anti-TROSY component. This\nexample also includes datasets measured from pure in-phase CPMG experiment\n(",(0,s.jsx)(t.a,{href:"../experiments/cpmg/cpmg_15n_ip",children:(0,s.jsx)(t.code,{children:"cpmg_15n_ip"})}),'), which is optional and may help to\nobtain even more accurate results. Note that \xb9\u2075N chemical shifts provided in\nparameter files should correspond to the "actual" value in the absence of 1JHN,\ntherefore \xb9\u2075N chemical shifts from pure in-phase experiment should be used\ninstead of the TROSY version.']}),"\n","\n",(0,s.jsxs)(t.section,{"data-footnotes":!0,className:"footnotes",children:[(0,s.jsx)(t.h2,{className:"sr-only",id:"footnote-label",children:"Footnotes"}),"\n",(0,s.jsxs)(t.ol,{children:["\n",(0,s.jsxs)(t.li,{id:"user-content-fn-1",children:["\n",(0,s.jsxs)(t.p,{children:["P. Vallurupalli, D. F. Hansen, E. Stollar, E. Meirovitch, and L. E. Kay.\nMeasurement of Bond Vector Orientations in Invisible Excited States of\nProteins ",(0,s.jsx)(t.em,{children:"Proc. Natl. Acad. Sci. USA"})," ",(0,s.jsx)(t.strong,{children:"104"}),", 18473-18477 (2007).\n",(0,s.jsx)(t.a,{href:"https://doi.org/10.1073/pnas.0708296104",children:"https://doi.org/10.1073/pnas.0708296104"})," ",(0,s.jsx)(t.a,{href:"#user-content-fnref-1","data-footnote-backref":"","aria-label":"Back to reference 1",className:"data-footnote-backref",children:"\u21a9"})]}),"\n"]}),"\n"]}),"\n"]})]})}function m(e={}){const{wrapper:t}={...(0,o.R)(),...e.components};return t?(0,s.jsx)(t,{...e,children:(0,s.jsx)(l,{...e})}):l(e)}},8453:(e,t,n)=>{n.d(t,{R:()=>i,x:()=>a});var s=n(6540);const o={},r=s.createContext(o);function i(e){const t=s.useContext(r);return s.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function a(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(o):e.components||o:i(e.components),s.createElement(r.Provider,{value:t},e.children)}}}]);