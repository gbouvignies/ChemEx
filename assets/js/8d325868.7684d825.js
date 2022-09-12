"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[8570],{3905:(e,t,n)=>{n.d(t,{Zo:()=>m,kt:()=>c});var a=n(7294);function i(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function s(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);t&&(a=a.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,a)}return n}function r(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?s(Object(n),!0).forEach((function(t){i(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):s(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,a,i=function(e,t){if(null==e)return{};var n,a,i={},s=Object.keys(e);for(a=0;a<s.length;a++)n=s[a],t.indexOf(n)>=0||(i[n]=e[n]);return i}(e,t);if(Object.getOwnPropertySymbols){var s=Object.getOwnPropertySymbols(e);for(a=0;a<s.length;a++)n=s[a],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(i[n]=e[n])}return i}var o=a.createContext({}),p=function(e){var t=a.useContext(o),n=t;return e&&(n="function"==typeof e?e(t):r(r({},t),e)),n},m=function(e){var t=p(e.components);return a.createElement(o.Provider,{value:t},e.children)},u={inlineCode:"code",wrapper:function(e){var t=e.children;return a.createElement(a.Fragment,{},t)}},d=a.forwardRef((function(e,t){var n=e.components,i=e.mdxType,s=e.originalType,o=e.parentName,m=l(e,["components","mdxType","originalType","parentName"]),d=p(n),c=i,N=d["".concat(o,".").concat(c)]||d[c]||u[c]||s;return n?a.createElement(N,r(r({ref:t},m),{},{components:n})):a.createElement(N,r({ref:t},m))}));function c(e,t){var n=arguments,i=t&&t.mdxType;if("string"==typeof e||i){var s=n.length,r=new Array(s);r[0]=d;var l={};for(var o in t)hasOwnProperty.call(t,o)&&(l[o]=t[o]);l.originalType=e,l.mdxType="string"==typeof e?e:i,r[1]=l;for(var p=2;p<s;p++)r[p]=n[p];return a.createElement.apply(null,r)}return a.createElement.apply(null,n)}d.displayName="MDXCreateElement"},5162:(e,t,n)=>{n.d(t,{Z:()=>r});var a=n(7294),i=n(6010);const s="tabItem_Ymn6";function r(e){let{children:t,hidden:n,className:r}=e;return a.createElement("div",{role:"tabpanel",className:(0,i.Z)(s,r),hidden:n},t)}},5488:(e,t,n)=>{n.d(t,{Z:()=>c});var a=n(7462),i=n(7294),s=n(6010),r=n(2389),l=n(7392),o=n(7094),p=n(2466);const m="tabList__CuJ",u="tabItem_LNqP";function d(e){var t,n;const{lazy:r,block:d,defaultValue:c,values:N,groupId:f,className:h}=e,g=i.Children.map(e.children,(e=>{if((0,i.isValidElement)(e)&&"value"in e.props)return e;throw new Error("Docusaurus error: Bad <Tabs> child <"+("string"==typeof e.type?e.type:e.type.name)+'>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.')})),k=null!=N?N:g.map((e=>{let{props:{value:t,label:n,attributes:a}}=e;return{value:t,label:n,attributes:a}})),b=(0,l.l)(k,((e,t)=>e.value===t.value));if(b.length>0)throw new Error('Docusaurus error: Duplicate values "'+b.map((e=>e.value)).join(", ")+'" found in <Tabs>. Every value needs to be unique.');const v=null===c?c:null!=(t=null!=c?c:null==(n=g.find((e=>e.props.default)))?void 0:n.props.value)?t:g[0].props.value;if(null!==v&&!k.some((e=>e.value===v)))throw new Error('Docusaurus error: The <Tabs> has a defaultValue "'+v+'" but none of its children has the corresponding value. Available values are: '+k.map((e=>e.value)).join(", ")+". If you intend to show no default tab, use defaultValue={null} instead.");const{tabGroupChoices:C,setTabGroupChoices:y}=(0,o.U)(),[_,x]=(0,i.useState)(v),w=[],{blockElementScrollPositionUntilNextRender:A}=(0,p.o5)();if(null!=f){const e=C[f];null!=e&&e!==_&&k.some((t=>t.value===e))&&x(e)}const B=e=>{const t=e.currentTarget,n=w.indexOf(t),a=k[n].value;a!==_&&(A(t),x(a),null!=f&&y(f,String(a)))},E=e=>{var t;let n=null;switch(e.key){case"ArrowRight":{var a;const t=w.indexOf(e.currentTarget)+1;n=null!=(a=w[t])?a:w[0];break}case"ArrowLeft":{var i;const t=w.indexOf(e.currentTarget)-1;n=null!=(i=w[t])?i:w[w.length-1];break}}null==(t=n)||t.focus()};return i.createElement("div",{className:(0,s.Z)("tabs-container",m)},i.createElement("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,s.Z)("tabs",{"tabs--block":d},h)},k.map((e=>{let{value:t,label:n,attributes:r}=e;return i.createElement("li",(0,a.Z)({role:"tab",tabIndex:_===t?0:-1,"aria-selected":_===t,key:t,ref:e=>w.push(e),onKeyDown:E,onFocus:B,onClick:B},r,{className:(0,s.Z)("tabs__item",u,null==r?void 0:r.className,{"tabs__item--active":_===t})}),null!=n?n:t)}))),r?(0,i.cloneElement)(g.filter((e=>e.props.value===_))[0],{className:"margin-top--md"}):i.createElement("div",{className:"margin-top--md"},g.map(((e,t)=>(0,i.cloneElement)(e,{key:t,hidden:e.props.value!==_})))))}function c(e){const t=(0,r.Z)();return i.createElement(d,(0,a.Z)({key:String(t)},e))}},9258:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>d,contentTitle:()=>m,default:()=>f,frontMatter:()=>p,metadata:()=>u,toc:()=>c});var a=n(7462),i=(n(7294),n(3905)),s=n(5488),r=n(5162);const l=n.p+"assets/images/cest_26hz_fit-cd540bd92fd0770c065e108770e4d2e1.png",o=n.p+"assets/images/cpmg_800mhz_fit-ec19160c7c18453507ae0aa05c4e295b.png",p={sidebar_position:8},m="Outputs",u={unversionedId:"user_guide/fitting/outputs",id:"user_guide/fitting/outputs",title:"Outputs",description:"Location",source:"@site/docs/user_guide/fitting/outputs.md",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/outputs",permalink:"/chemex/docs/user_guide/fitting/outputs",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/master/docs/user_guide/fitting/outputs.md",tags:[],version:"current",sidebarPosition:8,frontMatter:{sidebar_position:8},sidebar:"tutorialSidebar",previous:{title:"Kinetic models",permalink:"/chemex/docs/user_guide/fitting/kinetic_models"},next:{title:"Additional modules",permalink:"/chemex/docs/user_guide/additional_modules"}},d={},c=[{value:"Location",id:"location",level:2},{value:"Multi-step fit",id:"multi-step-fit",level:3},{value:"Multi-group fit",id:"multi-group-fit",level:3},{value:"Content",id:"content",level:2},{value:"<code>Parameters/</code>",id:"parameters",level:3},{value:"Example files",id:"example-files",level:4},{value:"<code>Plot/</code>",id:"plot",level:3},{value:"<code>Data/</code>",id:"data",level:3},{value:"<code>statistics.toml</code>",id:"statisticstoml",level:3}],N={toc:c};function f(e){let{components:t,...n}=e;return(0,i.kt)("wrapper",(0,a.Z)({},N,n,{components:t,mdxType:"MDXLayout"}),(0,i.kt)("h1",{id:"outputs"},"Outputs"),(0,i.kt)("h2",{id:"location"},"Location"),(0,i.kt)("p",null,"The results of the fits are written in the directory ",(0,i.kt)("inlineCode",{parentName:"p"},"Output/")," by default.\nHowever, you can change this location using the command-line option ",(0,i.kt)("inlineCode",{parentName:"p"},"--output"),"\n(or ",(0,i.kt)("inlineCode",{parentName:"p"},"-o"),") followed by the desired path name."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-bash"},"chemex fit -o <path> [other options]\n")),(0,i.kt)("h3",{id:"multi-step-fit"},"Multi-step fit"),(0,i.kt)("p",null,"If the fitting method has multiple fitting steps, each step will create its own\noutput subdirectory with the name of the fitting step."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"Output/\n\u251c\u2500\u2500 STEP1/\n\u2514\u2500\u2500 STEP2/\n")),(0,i.kt)("h3",{id:"multi-group-fit"},"Multi-group fit"),(0,i.kt)("p",null,"During any fitting step, if the dataset can be divided in multiple groups\ndepending on distinct sets of fitting parameters, ChemEx will fit each group of\ndata separately. The results of these individual fits are then stored in\nseparate subfolders placed in the the ",(0,i.kt)("inlineCode",{parentName:"p"},"Groups/")," folder. Subfolders are named\nwith a number and an ID, which depends on the parameters that have been\noptimized for the specific group. The results are also put together in the\n",(0,i.kt)("inlineCode",{parentName:"p"},"All/")," folder for convenience."),(0,i.kt)("p",null,"For example, if all global parameters (p",(0,i.kt)("sub",null,"B"),", k",(0,i.kt)("sub",null,"ex"),", etc.) are\nfixed or residue-specific fitting model is used (e.g. the model\n",(0,i.kt)("a",{parentName:"p",href:"/chemex/docs/user_guide/fitting/kinetic_models"},(0,i.kt)("inlineCode",{parentName:"a"},"2st_rs")),"), then all parameters are allowed to vary are\nresidue-specific. Residue-specific fits will then be run and two separate\nsubdirectories ",(0,i.kt)("inlineCode",{parentName:"p"},"All/")," and ",(0,i.kt)("inlineCode",{parentName:"p"},"Groups/")," will be created, which contain fitting\nresults for all residues and each individual residue, respectively."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-shell"},"Output/\n\u2514\u2500\u2500 STEP2/\n    \u251c\u2500\u2500 All/\n    \u2514\u2500\u2500 Groups/\n        \u251c\u2500\u2500 10_11N/\n        \u251c\u2500\u2500 11_12N/\n        \u251c\u2500\u2500 12_13N/\n")),(0,i.kt)("h2",{id:"content"},"Content"),(0,i.kt)("p",null,"The fitting output typically contains the following files and directories:"),(0,i.kt)("h3",{id:"parameters"},(0,i.kt)("inlineCode",{parentName:"h3"},"Parameters/")),(0,i.kt)("p",null,"Contains fitting results as three separate files ",(0,i.kt)("inlineCode",{parentName:"p"},"fitted.toml"),", ",(0,i.kt)("inlineCode",{parentName:"p"},"fixed.toml")," and\n",(0,i.kt)("inlineCode",{parentName:"p"},"constrained.toml"),", which contain output parameters that are fitted, fixed and\nconstrained during the fitting process, respectively."),(0,i.kt)("h4",{id:"example-files"},"Example files"),(0,i.kt)(s.Z,{mdxType:"Tabs"},(0,i.kt)(r.Z,{value:"fitted",label:"fitted.toml",default:!0,mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml"},'[GLOBAL]\nKEX_AB =  3.81511e+02 # \xb18.90870e+00\nPB     =  7.02971e-02 # \xb11.14784e-03\n\n[DW_AB]\n15N =  2.00075e+00 # \xb12.30817e-02\n31N =  1.98968e+00 # \xb11.90842e-02\n33N =  1.82003e+00 # \xb12.44821e-02\n34N =  3.63170e+00 # \xb13.62801e-02\n37N =  1.69183e+00 # \xb12.41070e-02\n\n["R2_A, B0->500.0MHZ"]\n15N =  3.98674e+00 # \xb12.55793e-01\n31N =  5.85923e+00 # \xb11.94734e-01\n33N =  4.02099e+00 # \xb12.78003e-01\n34N =  4.16615e+00 # \xb12.05190e-01\n37N =  3.67705e+00 # \xb13.04621e-01\n\n["R2_A, B0->800.0MHZ"]\n15N =  6.33712e+00 # \xb13.86894e-01\n31N =  7.99927e+00 # \xb12.81972e-01\n33N =  6.23967e+00 # \xb15.82366e-01\n34N =  5.99535e+00 # \xb13.17832e-01\n37N =  5.37295e+00 # \xb15.69929e-01\n')),(0,i.kt)("admonition",{type:"note"},(0,i.kt)("p",{parentName:"admonition"},'If uncertainties on fitted parameters have been calculated \u2013 this depends on the\nfitting algorithm used \u2013, they are reported as comments at the end of the line\nand are preceded by the sign "\xb1". Uncertainties on fitted parameters are,\ntypically, estimated through the covariance matrix obtained from the\nLevenberg-Marquardt optimization.'))),(0,i.kt)(r.Z,{value:"fixed",label:"fixed.toml",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml"},'[CS_A]\n15N =  1.19849e+02 # (fixed)\n31N =  1.26388e+02 # (fixed)\n33N =  1.18762e+02 # (fixed)\n34N =  1.14897e+02 # (fixed)\n37N =  1.21618e+02 # (fixed)\n\n["R1_A, B0->500.0MHZ"]\n15N =  1.50000e+00 # (fixed)\n31N =  1.50000e+00 # (fixed)\n33N =  1.50000e+00 # (fixed)\n34N =  1.50000e+00 # (fixed)\n37N =  1.50000e+00 # (fixed)\n\n["R1_A, B0->800.0MHZ"]\n15N =  1.50000e+00 # (fixed)\n31N =  1.50000e+00 # (fixed)\n33N =  1.50000e+00 # (fixed)\n34N =  1.50000e+00 # (fixed)\n37N =  1.50000e+00 # (fixed)\n'))),(0,i.kt)(r.Z,{value:"constrained",label:"constrained.toml",mdxType:"TabItem"},(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml"},'[GLOBAL]\nKAB =  2.68192e+01 # \xb13.06068e-01 ([KEX_AB] * [PB])\nKBA =  3.54692e+02 # \xb18.28245e+00 ([KEX_AB] * [PA])\nPA  =  9.29703e-01 # \xb11.14784e-03 (1.0 - [PB])\n\n[CS_B]\n15N =  1.21850e+02 # \xb12.30817e-02 ([CS_A, NUC->15N] + [DW_AB, NUC->15N])\n31N =  1.28378e+02 # \xb11.90842e-02 ([CS_A, NUC->31N] + [DW_AB, NUC->31N])\n33N =  1.20582e+02 # \xb12.44821e-02 ([CS_A, NUC->33N] + [DW_AB, NUC->33N])\n34N =  1.18529e+02 # \xb13.62801e-02 ([CS_A, NUC->34N] + [DW_AB, NUC->34N])\n37N =  1.23310e+02 # \xb12.41070e-02 ([CS_A, NUC->37N] + [DW_AB, NUC->37N])\n\n["R1_B, B0->500.0MHZ"]\n15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->500.0MHZ])\n31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->500.0MHZ])\n33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->500.0MHZ])\n34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->500.0MHZ])\n37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->500.0MHZ])\n\n["R1_B, B0->800.0MHZ"]\n15N =  1.50000e+00 # ([R1_A, NUC->15N, B0->800.0MHZ])\n31N =  1.50000e+00 # ([R1_A, NUC->31N, B0->800.0MHZ])\n33N =  1.50000e+00 # ([R1_A, NUC->33N, B0->800.0MHZ])\n34N =  1.50000e+00 # ([R1_A, NUC->34N, B0->800.0MHZ])\n37N =  1.50000e+00 # ([R1_A, NUC->37N, B0->800.0MHZ])\n\n["R2_B, B0->500.0MHZ"]\n15N =  3.98674e+00 # \xb12.55793e-01 ([R2_A, NUC->15N, B0->500.0MHZ])\n31N =  5.85923e+00 # \xb11.94734e-01 ([R2_A, NUC->31N, B0->500.0MHZ])\n33N =  4.02099e+00 # \xb12.78003e-01 ([R2_A, NUC->33N, B0->500.0MHZ])\n34N =  4.16615e+00 # \xb12.05190e-01 ([R2_A, NUC->34N, B0->500.0MHZ])\n37N =  3.67705e+00 # \xb13.04621e-01 ([R2_A, NUC->37N, B0->500.0MHZ])\n\n["R2_B, B0->800.0MHZ"]\n15N =  6.33712e+00 # \xb13.86894e-01 ([R2_A, NUC->15N, B0->800.0MHZ])\n31N =  7.99927e+00 # \xb12.81972e-01 ([R2_A, NUC->31N, B0->800.0MHZ])\n33N =  6.23967e+00 # \xb15.82366e-01 ([R2_A, NUC->33N, B0->800.0MHZ])\n34N =  5.99535e+00 # \xb13.17832e-01 ([R2_A, NUC->34N, B0->800.0MHZ])\n37N =  5.37295e+00 # \xb15.69929e-01 ([R2_A, NUC->37N, B0->800.0MHZ])\n')),(0,i.kt)("admonition",{type:"note"},(0,i.kt)("p",{parentName:"admonition"},"Propagated uncertainties \u2013 when available \u2013 and applied constrained are reported\nat the end of the line as comments.")))),(0,i.kt)("h3",{id:"plot"},(0,i.kt)("inlineCode",{parentName:"h3"},"Plot/")),(0,i.kt)("p",null,"Contains fitting results as plots (in ",(0,i.kt)("inlineCode",{parentName:"p"},".pdf")," format) and also the raw datasets\n(both the original input and fitted data points) for creating the plots. Example\nfitting results for CPMG and CEST experiments are shown below:"),(0,i.kt)("figure",null,(0,i.kt)("img",{src:l,alt:"CEST profile",width:"50%"}),(0,i.kt)("img",{src:o,alt:"CPMG profile",width:"50%"}),(0,i.kt)("figcaption",{align:"center"},(0,i.kt)("b",null,"Examples of CEST and CPMG fitting results"))),(0,i.kt)("admonition",{type:"note"},(0,i.kt)("p",{parentName:"admonition"},'In plots of (D-/cos-)CEST fitting results, the positions for ground and excited\nstates are indicated by solid and dashed vertical lines respectively. Besides,\ndata points that are filtered out from the fit are shown in lighter color. In\nplots of D-CEST/COS-CEST fitting results, the "folded" positions for ground and\nexcited states are indicated by "',"*",'" at the vertical lines.')),(0,i.kt)("h3",{id:"data"},(0,i.kt)("inlineCode",{parentName:"h3"},"Data/")),(0,i.kt)("p",null,"Contains all the data values used for the fitting along with the back-calculated\nvalues. These files can be used to calculate the ",(0,i.kt)("span",{parentName:"p",className:"math math-inline"},(0,i.kt)("span",{parentName:"span",className:"katex"},(0,i.kt)("span",{parentName:"span",className:"katex-mathml"},(0,i.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},(0,i.kt)("semantics",{parentName:"math"},(0,i.kt)("mrow",{parentName:"semantics"},(0,i.kt)("msup",{parentName:"mrow"},(0,i.kt)("mi",{parentName:"msup"},"\u03c7"),(0,i.kt)("mn",{parentName:"msup"},"2"))),(0,i.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\u03c7^2")))),(0,i.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,i.kt)("span",{parentName:"span",className:"base"},(0,i.kt)("span",{parentName:"span",className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,i.kt)("span",{parentName:"span",className:"mord"},(0,i.kt)("span",{parentName:"span",className:"mord mathnormal"},"\u03c7"),(0,i.kt)("span",{parentName:"span",className:"msupsub"},(0,i.kt)("span",{parentName:"span",className:"vlist-t"},(0,i.kt)("span",{parentName:"span",className:"vlist-r"},(0,i.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.8141em"}},(0,i.kt)("span",{parentName:"span",style:{top:"-3.063em",marginRight:"0.05em"}},(0,i.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,i.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,i.kt)("span",{parentName:"span",className:"mord mtight"},"2"))))))))))))," value."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:"title=Data/500mhz.dat",title:"Data/500mhz.dat"},"[15N]\n#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)\n             0    3.47059800e+04    1.77491406e+02    3.47055362e+04\n            30    3.05930380e+04    1.77491406e+02    3.06111963e+04\n             1    1.81234230e+04    1.77491406e+02    1.80856300e+04\n            28    3.07144730e+04    1.77491406e+02    3.05838767e+04\n             2    2.02155120e+04    1.77491406e+02    2.02110184e+04\n            26    3.05222020e+04    1.77491406e+02    3.05501255e+04\n             3    2.23056070e+04    1.77491406e+02    2.23601656e+04\n            24    3.05381830e+04    1.77491406e+02    3.05077686e+04\n             4    2.43783050e+04    1.77491406e+02    2.44111932e+04\n            22    3.06981570e+04    1.77491406e+02    3.04536377e+04\n             5    2.58673980e+04    1.77491406e+02    2.59884541e+04\n            20    3.06069180e+04    1.77491406e+02    3.03829761e+04\n             6    2.70660870e+04    1.77491406e+02    2.71148924e+04\n            18    3.00982700e+04    1.77491406e+02    3.02883916e+04\n             7    2.81512990e+04    1.77491406e+02    2.79178594e+04\n            16    3.02515700e+04    1.77491406e+02    3.01579211e+04\n             8    2.84045570e+04    1.77491406e+02    2.84985729e+04\n            14    2.97528530e+04    1.77491406e+02    2.99712519e+04\n             9    2.84536650e+04    1.77491406e+02    2.89264158e+04\n            13    2.98219180e+04    1.77491406e+02    2.98461022e+04\n            10    2.96936410e+04    1.77491406e+02    2.92495591e+04\n            12    2.95269900e+04    1.77491406e+02    2.96918706e+04\n            11    2.92540210e+04    1.77491406e+02    2.94972506e+04\n             2    2.04476470e+04    1.77491406e+02    2.02110184e+04\n            28    3.05466550e+04    1.77491406e+02    3.05838767e+04\n             8    2.86183900e+04    1.77491406e+02    2.84985729e+04\n\n[31N]\n#         NCYC   INTENSITY (EXP)       ERROR (EXP)  INTENSITY (CALC)\n             0    4.71577550e+04    1.77491406e+02    4.71537007e+04\n            30    4.00023250e+04    1.77491406e+02    3.94823953e+04\n             1    2.30167260e+04    1.77491406e+02    2.34136180e+04\n            28    4.00615990e+04    1.77491406e+02    3.94547546e+04\n             2    2.61934530e+04    1.77491406e+02    2.61653780e+04\n            26    3.96898190e+04    1.77491406e+02    3.94190638e+04\n             3    2.86882570e+04    1.77491406e+02    2.88729551e+04\n            24    3.98238690e+04    1.77491406e+02    3.93722573e+04\n             4    3.21031440e+04    1.77491406e+02    3.13942278e+04\n            22    3.96733310e+04    1.77491406e+02    3.93097577e+04\n             5    3.35871430e+04    1.77491406e+02    3.34628362e+04\n            20    3.86957600e+04    1.77491406e+02    3.92245398e+04\n             6    3.53600830e+04    1.77491406e+02    3.49641428e+04\n            18    3.88115960e+04    1.77491406e+02    3.91054958e+04\n             7    3.58938150e+04    1.77491406e+02    3.59453347e+04\n            16    3.85843960e+04    1.77491406e+02    3.89345050e+04\n             8    3.61450060e+04    1.77491406e+02    3.66900206e+04\n            14    3.85711880e+04    1.77491406e+02    3.86811094e+04\n             9    3.72815770e+04    1.77491406e+02    3.72429094e+04\n            13    3.82358560e+04    1.77491406e+02    3.85092366e+04\n            10    3.76426640e+04    1.77491406e+02    3.76803042e+04\n            12    3.78470460e+04    1.77491406e+02    3.82929922e+04\n            11    3.76000060e+04    1.77491406e+02    3.80220076e+04\n             2    2.62301400e+04    1.77491406e+02    2.61653780e+04\n            28    3.97731200e+04    1.77491406e+02    3.94547546e+04\n             8    3.63510140e+04    1.77491406e+02    3.66900206e+04\n\n[...]\n")),(0,i.kt)("h3",{id:"statisticstoml"},(0,i.kt)("inlineCode",{parentName:"h3"},"statistics.toml")),(0,i.kt)("p",null,"Contains all goodness-of-fit statistics, such as ",(0,i.kt)("span",{parentName:"p",className:"math math-inline"},(0,i.kt)("span",{parentName:"span",className:"katex"},(0,i.kt)("span",{parentName:"span",className:"katex-mathml"},(0,i.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},(0,i.kt)("semantics",{parentName:"math"},(0,i.kt)("mrow",{parentName:"semantics"},(0,i.kt)("msup",{parentName:"mrow"},(0,i.kt)("mi",{parentName:"msup"},"\u03c7"),(0,i.kt)("mn",{parentName:"msup"},"2"))),(0,i.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\u03c7^2")))),(0,i.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,i.kt)("span",{parentName:"span",className:"base"},(0,i.kt)("span",{parentName:"span",className:"strut",style:{height:"1.0085em",verticalAlign:"-0.1944em"}}),(0,i.kt)("span",{parentName:"span",className:"mord"},(0,i.kt)("span",{parentName:"span",className:"mord mathnormal"},"\u03c7"),(0,i.kt)("span",{parentName:"span",className:"msupsub"},(0,i.kt)("span",{parentName:"span",className:"vlist-t"},(0,i.kt)("span",{parentName:"span",className:"vlist-r"},(0,i.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.8141em"}},(0,i.kt)("span",{parentName:"span",style:{top:"-3.063em",marginRight:"0.05em"}},(0,i.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,i.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,i.kt)("span",{parentName:"span",className:"mord mtight"},"2")))))))))))),"."),(0,i.kt)("pre",null,(0,i.kt)("code",{parentName:"pre",className:"language-toml",metastring:"title='\"statistics.toml\"'",title:"'\"statistics.toml\"'"},'"number of data points"                = 230\n"number of variables"                  = 17\n"chi-square"                           =  4.34824e+02\n"reduced-chi-square"                   =  2.04143e+00\n"chi-squared test"                     =  0.00000e+00\n"Kolmogorov-Smirnov test"              =  6.72236e-02\n"Akaike Information Criterion (AIC)"   =  1.80479e+02\n"Bayesian Information Criterion (BIC)" =  2.38926e+02\n\n')))}f.isMDXComponent=!0}}]);