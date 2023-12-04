"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[4914],{6344:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>o,default:()=>h,frontMatter:()=>r,metadata:()=>a,toc:()=>d});var s=n(5893),i=n(1151);const r={sidebar_position:7,sidebar_label:"\xb9\u2075N D-CEST experiment"},o="\xb9\u2075N D-CEST experiment",a={id:"examples/dcest",title:"\xb9\u2075N D-CEST experiment",description:"You can get the example files from the GitHub ChemEx page:",source:"@site/docs/examples/dcest.md",sourceDirName:"examples",slug:"/examples/dcest",permalink:"/ChemEx/docs/examples/dcest",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/examples/dcest.md",tags:[],version:"current",sidebarPosition:7,frontMatter:{sidebar_position:7,sidebar_label:"\xb9\u2075N D-CEST experiment"},sidebar:"tutorialSidebar",previous:{title:"\xb9H CEST experiments",permalink:"/ChemEx/docs/examples/cest_ip_ap"},next:{title:"\xb9H\u1d3a COS-CEST experiment",permalink:"/ChemEx/docs/examples/coscest"}},c={},d=[];function l(e){const t={a:"a",admonition:"admonition",code:"code",em:"em",h1:"h1",h2:"h2",li:"li",ol:"ol",p:"p",section:"section",strong:"strong",sup:"sup",...(0,i.a)(),...e.components};return(0,s.jsxs)(s.Fragment,{children:[(0,s.jsx)(t.h1,{id:"n-d-cest-experiment",children:"\xb9\u2075N D-CEST experiment"}),"\n",(0,s.jsx)(t.admonition,{title:"Files",type:"note",children:(0,s.jsxs)(t.p,{children:["You can get the example files from the GitHub ChemEx page:\n",(0,s.jsx)(t.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N",children:"regular"}),"\nand\n",(0,s.jsx)(t.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/DCEST_15N_HD_EXCH",children:"H/D solvent exchange"}),"\nstudies."]})}),"\n",(0,s.jsxs)(t.p,{children:["The two examples demonstrate the application of the \xb9\u2075N D-CEST experiment\n(",(0,s.jsx)(t.a,{href:"/ChemEx/docs/experiments/dcest/dcest_15n",children:(0,s.jsx)(t.code,{children:"dcest_15n"})}),"). D-CEST experiment accelerates\ndata acquisition by using the DANTE scheme to perform multifrequency irradiation\nat a desired effective B",(0,s.jsx)("sub",{children:"1"})," field. In D-CEST experiment, the saturation\nelement is replaced with DANTE pulses to perform multifrequency irradiation,\nwhich results in dramatic savings in measurement time.",(0,s.jsx)(t.sup,{children:(0,s.jsx)(t.a,{href:"#user-content-fn-1",id:"user-content-fnref-1","data-footnote-ref":!0,"aria-describedby":"footnote-label",children:"1"})})]}),"\n",(0,s.jsxs)(t.p,{children:["An important application of D-CEST experiment is to measure H/D solvent exchange\nrates, the results should be more accurate than traditional CLEANEX experiment,\nespecially for studying highly flexible biomolecules such as intrinsically\ndisordered proteins.",(0,s.jsx)(t.sup,{children:(0,s.jsx)(t.a,{href:"#user-content-fn-2",id:"user-content-fnref-2","data-footnote-ref":!0,"aria-describedby":"footnote-label",children:"2"})})]}),"\n",(0,s.jsx)(t.admonition,{type:"info",children:(0,s.jsxs)(t.p,{children:["Use the kinetic model ",(0,s.jsx)(t.code,{children:"2st_hd"})," for H/D solvent exchange study. Besides, you\nshould provide the amount of D",(0,s.jsx)("sub",{children:"2"}),"O in the sample with the ",(0,s.jsx)(t.code,{children:"d2o"})," key in\nthe section ",(0,s.jsx)(t.code,{children:"[conditions]"})," of the ",(0,s.jsx)(t.strong,{children:"experiment.toml"})," file. The ",(0,s.jsx)(t.code,{children:"d2o"})," key is\nincluded in experiment files to differentiate datasets recorded with different\namounts of D",(0,s.jsx)("sub",{children:"2"}),"O. Note that the amount of D",(0,s.jsx)("sub",{children:"2"}),"O is a parameter\nthat can be fitted to account for potential errors during sample preparation."]})}),"\n",(0,s.jsxs)(t.section,{"data-footnotes":!0,className:"footnotes",children:[(0,s.jsx)(t.h2,{className:"sr-only",id:"footnote-label",children:"Footnotes"}),"\n",(0,s.jsxs)(t.ol,{children:["\n",(0,s.jsxs)(t.li,{id:"user-content-fn-1",children:["\n",(0,s.jsxs)(t.p,{children:["T. Yuwen, L. E. Kay, and G. Bouvignies. Dramatic Decrease in CEST\nMeasurement Times Using Multi-Site Excitation. ",(0,s.jsx)(t.em,{children:"ChemPhysChem"})," ",(0,s.jsx)(t.strong,{children:"19"}),",\n1707-1710 (2018). ",(0,s.jsx)(t.a,{href:"https://doi.org/10.1002/cphc.201800249",children:"https://doi.org/10.1002/cphc.201800249"})," ",(0,s.jsx)(t.a,{href:"#user-content-fnref-1","data-footnote-backref":"","aria-label":"Back to reference 1",className:"data-footnote-backref",children:"\u21a9"})]}),"\n"]}),"\n",(0,s.jsxs)(t.li,{id:"user-content-fn-2",children:["\n",(0,s.jsxs)(t.p,{children:["T. Yuwen, A. Bah, J. P. Brady, F. Ferrage, G. Bouvignies, and L. E. Kay.\nMeasuring Solvent Hydrogen Exchange Rates by Multifrequency Excitation \xb9\u2075N\nCEST: Application to Protein Phase Separation. ",(0,s.jsx)(t.em,{children:"J. Phys. Chem. B"})," ",(0,s.jsx)(t.strong,{children:"122"}),",\n11206-11217 (2018). ",(0,s.jsx)(t.a,{href:"https://doi.org/10.1021/acs.jpcb.8b06820",children:"https://doi.org/10.1021/acs.jpcb.8b06820"})," ",(0,s.jsx)(t.a,{href:"#user-content-fnref-2","data-footnote-backref":"","aria-label":"Back to reference 2",className:"data-footnote-backref",children:"\u21a9"})]}),"\n"]}),"\n"]}),"\n"]})]})}function h(e={}){const{wrapper:t}={...(0,i.a)(),...e.components};return t?(0,s.jsx)(t,{...e,children:(0,s.jsx)(l,{...e})}):l(e)}},1151:(e,t,n)=>{n.d(t,{Z:()=>a,a:()=>o});var s=n(7294);const i={},r=s.createContext(i);function o(e){const t=s.useContext(r);return s.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function a(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(i):e.components||i:o(e.components),s.createElement(r.Provider,{value:t},e.children)}}}]);