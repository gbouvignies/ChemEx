"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[2413],{6665:(e,s,n)=>{n.r(s),n.d(s,{assets:()=>c,contentTitle:()=>r,default:()=>d,frontMatter:()=>a,metadata:()=>l,toc:()=>o});var i=n(5893),t=n(1151);const a={sidebar_position:5},r="Parameter files",l={id:"user_guide/fitting/parameter_files",title:"Parameter files",description:"Description",source:"@site/docs/user_guide/fitting/parameter_files.md",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/parameter_files",permalink:"/ChemEx/docs/user_guide/fitting/parameter_files",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/parameter_files.md",tags:[],version:"current",sidebarPosition:5,frontMatter:{sidebar_position:5},sidebar:"tutorialSidebar",previous:{title:"Data files",permalink:"/ChemEx/docs/user_guide/fitting/data_files"},next:{title:"Method files",permalink:"/ChemEx/docs/user_guide/fitting/method_files"}},c={},o=[{value:"Description",id:"description",level:2},{value:"Example",id:"example",level:2},{value:"Setting bounds",id:"setting-bounds",level:2}];function m(e){const s={admonition:"admonition",annotation:"annotation",code:"code",h1:"h1",h2:"h2",li:"li",math:"math",mi:"mi",mn:"mn",mrow:"mrow",msup:"msup",p:"p",pre:"pre",semantics:"semantics",span:"span",ul:"ul",...(0,t.a)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(s.h1,{id:"parameter-files",children:"Parameter files"}),"\n",(0,i.jsx)(s.h2,{id:"description",children:"Description"}),"\n",(0,i.jsxs)(s.p,{children:["The parameter files contain initial estimates of the parameters to be used\nduring the fitting process. These files are provided to ChemEx using the ",(0,i.jsx)(s.code,{children:"-p"})," or\n",(0,i.jsx)(s.code,{children:"--parameters"})," option."]}),"\n",(0,i.jsx)(s.pre,{children:(0,i.jsx)(s.code,{className:"language-shell",children:"chemex fit [...] -p <parameter_file> [...]\n"})}),"\n",(0,i.jsx)(s.p,{children:"Parameter files are divided in multiple sections:"}),"\n",(0,i.jsxs)(s.ul,{children:["\n",(0,i.jsxs)(s.li,{children:["The parameter values under the section ",(0,i.jsx)(s.code,{children:"[GLOBAL]"})," apply to all residues."]}),"\n",(0,i.jsxs)(s.li,{children:["Residue-specific parameters are specified under sections named after the\nparameter name, such as ",(0,i.jsx)(s.code,{children:"[CS_A]"}),". Multiple parameter files can be provided if\nnecessary."]}),"\n"]}),"\n",(0,i.jsx)(s.admonition,{type:"warning",children:(0,i.jsxs)(s.p,{children:["Due to the multidimensional searching feature of ",(0,i.jsxs)(s.span,{className:"katex",children:[(0,i.jsx)(s.span,{className:"katex-mathml",children:(0,i.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,i.jsxs)(s.semantics,{children:[(0,i.jsx)(s.mrow,{children:(0,i.jsxs)(s.msup,{children:[(0,i.jsx)(s.mi,{children:"\u03c7"}),(0,i.jsx)(s.mn,{children:"2"})]})}),(0,i.jsx)(s.annotation,{encoding:"application/x-tex",children:"\u03c7^2"})]})})}),(0,i.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,i.jsxs)(s.span,{className:"base",children:[(0,i.jsx)(s.span,{className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,i.jsxs)(s.span,{className:"mord",children:[(0,i.jsx)(s.span,{className:"mord mathnormal",children:"\u03c7"}),(0,i.jsx)(s.span,{className:"msupsub",children:(0,i.jsx)(s.span,{className:"vlist-t",children:(0,i.jsx)(s.span,{className:"vlist-r",children:(0,i.jsx)(s.span,{className:"vlist",style:{height:"0.8141em"},children:(0,i.jsxs)(s.span,{style:{top:"-3.063em",marginRight:"0.05em"},children:[(0,i.jsx)(s.span,{className:"pstrut",style:{height:"2.7em"}}),(0,i.jsx)(s.span,{className:"sizing reset-size6 size3 mtight",children:(0,i.jsx)(s.span,{className:"mord mtight",children:"2"})})]})})})})})]})]})})]})," minimization process, it\nis essential to set a suitable initial value for each parameter to avoid getting\ntrapped in a local minimum."]})}),"\n",(0,i.jsx)(s.admonition,{type:"info",children:(0,i.jsx)(s.p,{children:"When no starting value is provided in the parameter files, a default value is\nassigned as the initial value."})}),"\n",(0,i.jsx)(s.h2,{id:"example",children:"Example"}),"\n",(0,i.jsx)(s.p,{children:"Here is an example parameter file:"}),"\n",(0,i.jsx)(s.pre,{children:(0,i.jsx)(s.code,{className:"language-toml",metastring:'title="parameters.toml"',children:"[GLOBAL]\nPB     = 0.015\nKEX_AB = 70.0\nTAUC_A = 10.0\n\n[CS_A]\n13N = 108.207\n26N = 115.711\n28N = 113.882\n29N = 115.318\n33N = 115.636\n37N = 116.159\n41N = 114.635\n42N = 113.525\n43N = 108.876\n50N = 107.855\n52N = 111.358\n55N = 128.301\n59N = 116.388\n66N = 119.429\n67N = 114.454\n68N = 120.595\n\n[DW_AB]\n13N = 4.0\n26N = 5.5\n28N = 6.5\n29N = 6.0\n33N = 4.5\n37N = 5.0\n41N = 6.0\n42N = 6.0\n43N = 12.5\n50N = 8.0\n52N = 8.5\n55N = -6.5\n59N = 6.5\n66N = 4.0\n67N = 8.0\n68N = 4.5\n"})}),"\n",(0,i.jsx)(s.admonition,{type:"tip",children:(0,i.jsxs)(s.p,{children:["Setting model-free parameters (e.g. ",(0,i.jsx)(s.code,{children:"TAUC_A"}),") is a simple way to obtain initial\nestimates of relaxation parameters (e.g., ",(0,i.jsx)(s.code,{children:"R1_A"}),", ",(0,i.jsx)(s.code,{children:"R2_A"}),", etc.). For every 2.6\nkDa molecular weight, the overall tumbling time is approximately 1 ns at T = 300\nK for biomolecules in H",(0,i.jsx)("sub",{children:"2"}),"O. Assuming similar molecular structure at\ndifferent conditions, the overall tumbling time is proportional to \u03b7/T, where \u03b7\nis solution viscosity and T is temperature in Kelvin."]})}),"\n",(0,i.jsx)(s.h2,{id:"setting-bounds",children:"Setting bounds"}),"\n",(0,i.jsx)(s.p,{children:"You can set upper and lower bounds to any fitting parameters by replacing the\ninitial value by a list of three elements:"}),"\n",(0,i.jsx)(s.pre,{children:(0,i.jsx)(s.code,{className:"language-toml",metastring:'title="parameters.toml"',children:"PARAMETER_WITH_NO_BOUNDS = <initial_value>\nPARAMETER_WITH_BOUNDS = [<initial_value>, <lower_bound>, <upper_bound>]\n"})}),"\n",(0,i.jsxs)(s.p,{children:["Such boundaries can help to prevent parameters wandering off to unrealistic\nvalues to minimize ",(0,i.jsxs)(s.span,{className:"katex",children:[(0,i.jsx)(s.span,{className:"katex-mathml",children:(0,i.jsx)(s.math,{xmlns:"http://www.w3.org/1998/Math/MathML",children:(0,i.jsxs)(s.semantics,{children:[(0,i.jsx)(s.mrow,{children:(0,i.jsxs)(s.msup,{children:[(0,i.jsx)(s.mi,{children:"\u03c7"}),(0,i.jsx)(s.mn,{children:"2"})]})}),(0,i.jsx)(s.annotation,{encoding:"application/x-tex",children:"\u03c7^2"})]})})}),(0,i.jsx)(s.span,{className:"katex-html","aria-hidden":"true",children:(0,i.jsxs)(s.span,{className:"base",children:[(0,i.jsx)(s.span,{className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,i.jsxs)(s.span,{className:"mord",children:[(0,i.jsx)(s.span,{className:"mord mathnormal",children:"\u03c7"}),(0,i.jsx)(s.span,{className:"msupsub",children:(0,i.jsx)(s.span,{className:"vlist-t",children:(0,i.jsx)(s.span,{className:"vlist-r",children:(0,i.jsx)(s.span,{className:"vlist",style:{height:"0.8141em"},children:(0,i.jsxs)(s.span,{style:{top:"-3.063em",marginRight:"0.05em"},children:[(0,i.jsx)(s.span,{className:"pstrut",style:{height:"2.7em"}}),(0,i.jsx)(s.span,{className:"sizing reset-size6 size3 mtight",children:(0,i.jsx)(s.span,{className:"mord mtight",children:"2"})})]})})})})})]})]})})]}),". However, one should be careful not to set too\nstringent boundaries either, as this can result in convergence problems. Certain\nminimization algorithms (e.g. ",(0,i.jsx)(s.code,{children:"differential_evolution"}),") require finite bounds on\nall fitted parameters."]})]})}function d(e={}){const{wrapper:s}={...(0,t.a)(),...e.components};return s?(0,i.jsx)(s,{...e,children:(0,i.jsx)(m,{...e})}):m(e)}},1151:(e,s,n)=>{n.d(s,{Z:()=>l,a:()=>r});var i=n(7294);const t={},a=i.createContext(t);function r(e){const s=i.useContext(a);return i.useMemo((function(){return"function"==typeof e?e(s):{...s,...e}}),[s,e])}function l(e){let s;return s=e.disableParentContext?"function"==typeof e.components?e.components(t):e.components||t:r(e.components),i.createElement(a.Provider,{value:s},e.children)}}}]);