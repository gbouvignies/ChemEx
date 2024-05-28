"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[9227],{5273:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>c,contentTitle:()=>o,default:()=>f,frontMatter:()=>l,metadata:()=>u,toc:()=>d});var r=n(4848),a=n(8453),i=n(1470),s=n(9365);const l={sidebar_position:4},o="Data files",u={id:"user_guide/fitting/data_files",title:"Data files",description:"Description",source:"@site/docs/user_guide/fitting/data_files.mdx",sourceDirName:"user_guide/fitting",slug:"/user_guide/fitting/data_files",permalink:"/ChemEx/docs/user_guide/fitting/data_files",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/fitting/data_files.mdx",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_position:4},sidebar:"tutorialSidebar",previous:{title:"Experiment files",permalink:"/ChemEx/docs/user_guide/fitting/experiment_files"},next:{title:"Parameter files",permalink:"/ChemEx/docs/user_guide/fitting/parameter_files"}},c={},d=[{value:"Description",id:"description",level:2},{value:"Example",id:"example",level:2}];function h(e){const t={a:"a",admonition:"admonition",code:"code",h1:"h1",h2:"h2",p:"p",pre:"pre",table:"table",tbody:"tbody",td:"td",th:"th",thead:"thead",tr:"tr",...(0,a.R)(),...e.components};return(0,r.jsxs)(r.Fragment,{children:[(0,r.jsx)(t.h1,{id:"data-files",children:"Data files"}),"\n",(0,r.jsx)(t.h2,{id:"description",children:"Description"}),"\n",(0,r.jsxs)(t.p,{children:["The location of data files is specified within\n",(0,r.jsx)(t.a,{href:"/ChemEx/docs/user_guide/fitting/experiment_files",children:"experiment files"}),". Data files typically have three columns\ncontaining the following information:"]}),"\n",(0,r.jsxs)(t.table,{children:[(0,r.jsx)(t.thead,{children:(0,r.jsxs)(t.tr,{children:[(0,r.jsx)(t.th,{children:"Experiment type"}),(0,r.jsx)(t.th,{children:"First column"}),(0,r.jsx)(t.th,{children:"Second column"}),(0,r.jsx)(t.th,{children:"Third column"})]})}),(0,r.jsxs)(t.tbody,{children:[(0,r.jsxs)(t.tr,{children:[(0,r.jsx)(t.td,{children:"CEST / D-CEST / cos-CEST"}),(0,r.jsx)(t.td,{children:"Offset (Hz)"}),(0,r.jsx)(t.td,{children:"Intensity"}),(0,r.jsx)(t.td,{children:"Uncertainty"})]}),(0,r.jsxs)(t.tr,{children:[(0,r.jsx)(t.td,{children:"CPMG"}),(0,r.jsx)(t.td,{children:"Number of cycles"}),(0,r.jsx)(t.td,{children:"Intensity"}),(0,r.jsx)(t.td,{children:"Uncertainty"})]}),(0,r.jsxs)(t.tr,{children:[(0,r.jsx)(t.td,{children:"Relaxation"}),(0,r.jsx)(t.td,{children:"Delay (s)"}),(0,r.jsx)(t.td,{children:"Intensity"}),(0,r.jsx)(t.td,{children:"Uncertainty"})]})]})]}),"\n",(0,r.jsxs)(t.p,{children:["Data files are peak intensity tables that can be obtained from the output of\nmany peak fitting programs, such as the ",(0,r.jsx)(t.code,{children:"autoFit"})," subroutine in the\n",(0,r.jsx)(t.a,{href:"https://www.ibbr.umd.edu/nmrpipe/index.html",children:"NMRPipe"})," program."]}),"\n",(0,r.jsx)(t.h2,{id:"example",children:"Example"}),"\n","\n",(0,r.jsxs)(i.A,{children:[(0,r.jsx)(s.A,{value:"cest",label:"CEST",default:!0,children:(0,r.jsx)(t.pre,{children:(0,r.jsx)(t.code,{className:"language-python",metastring:'title="./Data/30Hz/F52N-H.out"',children:"# Offset (Hz)        Intensity    Uncertainty\n   -2.000e+04    2.5352430e+07  7.5500000e+04\n   -4.000e+02    1.7378230e+07  7.5500000e+04\n   -3.500e+02    1.7043830e+07  7.5500000e+04\n   -3.000e+02    1.5819080e+07  7.5500000e+04\n   -2.500e+02    8.8238340e+06  7.5500000e+04\n   -2.000e+02    3.0651600e+06  7.5500000e+04\n   -1.500e+02    1.4991090e+07  7.5500000e+04\n   -1.000e+02    1.6768950e+07  7.5500000e+04\n   -5.000e+01    1.7215280e+07  7.5500000e+04\n    0.000e+00    1.7636170e+07  7.5500000e+04\n    5.000e+01    1.7612980e+07  7.5500000e+04\n    1.000e+02    1.7666720e+07  7.5500000e+04\n    1.500e+02    1.7694720e+07  7.5500000e+04\n    2.000e+02    1.7692340e+07  7.5500000e+04\n    2.500e+02    1.7747680e+07  7.5500000e+04\n    3.000e+02    1.7633460e+07  7.5500000e+04\n    3.500e+02    1.7461530e+07  7.5500000e+04\n    4.000e+02    1.7377560e+07  7.5500000e+04\n"})})}),(0,r.jsx)(s.A,{value:"cpmg",label:"CPMG",children:(0,r.jsx)(t.pre,{children:(0,r.jsx)(t.code,{className:"language-python",metastring:'title="./Data/800MHz/G39N-H.out"',children:"#    ncyc_cp        Intensity      Esd(Int.)\n   0.000e+00    1.2584470e+07  9.8068130e+03\n   4.000e+01    9.1841462e+06  9.8068130e+03\n   1.000e+00    6.7679280e+06  9.8068130e+03\n   3.600e+01    9.0545262e+06  9.8068130e+03\n   2.000e+00    6.7616357e+06  9.8068130e+03\n   3.200e+01    8.8758267e+06  9.8068130e+03\n   3.000e+00    6.7377252e+06  9.8068130e+03\n   3.000e+01    8.7499820e+06  9.8068130e+03\n   4.000e+00    6.8950311e+06  9.8068130e+03\n   2.800e+01    8.6316880e+06  9.8068130e+03\n   5.000e+00    6.8094567e+06  9.8068130e+03\n   2.600e+01    8.5410798e+06  9.8068130e+03\n   6.000e+00    6.6798367e+06  9.8068130e+03\n   2.400e+01    8.4227858e+06  9.8068130e+03\n   7.000e+00    6.6458586e+06  9.8068130e+03\n   2.200e+01    8.2566708e+06  9.8068130e+03\n   8.000e+00    6.7125563e+06  9.8068130e+03\n   2.000e+01    8.0012060e+06  9.8068130e+03\n   9.000e+00    6.6949380e+06  9.8068130e+03\n   1.800e+01    7.8099221e+06  9.8068130e+03\n   1.000e+01    6.7440175e+06  9.8068130e+03\n   1.600e+01    7.5330637e+06  9.8068130e+03\n   1.100e+01    6.8925142e+06  9.8068130e+03\n   1.400e+01    7.2436209e+06  9.8068130e+03\n   1.200e+01    6.9504028e+06  9.8068130e+03\n   1.000e+00    6.8107152e+06  9.8068130e+03\n   4.000e+01    9.1904384e+06  9.8068130e+03\n   2.000e+00    6.7855462e+06  9.8068130e+03\n   3.600e+01    9.0104805e+06  9.8068130e+03\n"})})}),(0,r.jsx)(s.A,{value:"relaxation",label:"Relaxation",children:(0,r.jsx)(t.pre,{children:(0,r.jsx)(t.code,{className:"language-python",metastring:'title="./Data/800MHz/S5N-H.out"',children:"#  Times (s)        Intensity      Esd(Int.)\n   5.000e-02    9.2822624e+05  1.0000000e+04\n   1.000e-01    8.6879557e+05  1.0000000e+04\n   1.500e-01    7.8922229e+05  1.0000000e+04\n   2.000e-01    7.7218645e+05  1.0000000e+04\n   2.500e-01    6.8675373e+05  1.0000000e+04\n   3.000e-01    6.4862935e+05  1.0000000e+04\n   3.500e-01    6.1765521e+05  1.0000000e+04\n   4.000e-01    5.5875618e+05  1.0000000e+04\n   4.500e-01    5.2201228e+05  1.0000000e+04\n   5.000e-01    4.8749523e+05  1.0000000e+04\n"})})})]}),"\n",(0,r.jsx)(t.admonition,{type:"note",children:(0,r.jsxs)(t.p,{children:["The input data files for ",(0,r.jsx)(t.code,{children:"shift_experiments"})," have special format, refer to the\n",(0,r.jsx)(t.code,{children:"Shifts/"})," example under ",(0,r.jsx)(t.code,{children:"examples/Combinations/"})," to learn how to create data\nfiles for ",(0,r.jsx)(t.code,{children:"shift_experiments"}),"."]})})]})}function f(e={}){const{wrapper:t}={...(0,a.R)(),...e.components};return t?(0,r.jsx)(t,{...e,children:(0,r.jsx)(h,{...e})}):h(e)}},9365:(e,t,n)=>{n.d(t,{A:()=>s});n(6540);var r=n(4164);const a={tabItem:"tabItem_Ymn6"};var i=n(4848);function s(e){let{children:t,hidden:n,className:s}=e;return(0,i.jsx)("div",{role:"tabpanel",className:(0,r.A)(a.tabItem,s),hidden:n,children:t})}},1470:(e,t,n)=>{n.d(t,{A:()=>w});var r=n(6540),a=n(4164),i=n(3104),s=n(6347),l=n(205),o=n(7485),u=n(1682),c=n(9466);function d(e){return r.Children.toArray(e).filter((e=>"\n"!==e)).map((e=>{if(!e||(0,r.isValidElement)(e)&&function(e){const{props:t}=e;return!!t&&"object"==typeof t&&"value"in t}(e))return e;throw new Error(`Docusaurus error: Bad <Tabs> child <${"string"==typeof e.type?e.type:e.type.name}>: all children of the <Tabs> component should be <TabItem>, and every <TabItem> should have a unique "value" prop.`)}))?.filter(Boolean)??[]}function h(e){const{values:t,children:n}=e;return(0,r.useMemo)((()=>{const e=t??function(e){return d(e).map((e=>{let{props:{value:t,label:n,attributes:r,default:a}}=e;return{value:t,label:n,attributes:r,default:a}}))}(n);return function(e){const t=(0,u.X)(e,((e,t)=>e.value===t.value));if(t.length>0)throw new Error(`Docusaurus error: Duplicate values "${t.map((e=>e.value)).join(", ")}" found in <Tabs>. Every value needs to be unique.`)}(e),e}),[t,n])}function f(e){let{value:t,tabValues:n}=e;return n.some((e=>e.value===t))}function p(e){let{queryString:t=!1,groupId:n}=e;const a=(0,s.W6)(),i=function(e){let{queryString:t=!1,groupId:n}=e;if("string"==typeof t)return t;if(!1===t)return null;if(!0===t&&!n)throw new Error('Docusaurus error: The <Tabs> component groupId prop is required if queryString=true, because this value is used as the search param name. You can also provide an explicit value such as queryString="my-search-param".');return n??null}({queryString:t,groupId:n});return[(0,o.aZ)(i),(0,r.useCallback)((e=>{if(!i)return;const t=new URLSearchParams(a.location.search);t.set(i,e),a.replace({...a.location,search:t.toString()})}),[i,a])]}function m(e){const{defaultValue:t,queryString:n=!1,groupId:a}=e,i=h(e),[s,o]=(0,r.useState)((()=>function(e){let{defaultValue:t,tabValues:n}=e;if(0===n.length)throw new Error("Docusaurus error: the <Tabs> component requires at least one <TabItem> children component");if(t){if(!f({value:t,tabValues:n}))throw new Error(`Docusaurus error: The <Tabs> has a defaultValue "${t}" but none of its children has the corresponding value. Available values are: ${n.map((e=>e.value)).join(", ")}. If you intend to show no default tab, use defaultValue={null} instead.`);return t}const r=n.find((e=>e.default))??n[0];if(!r)throw new Error("Unexpected error: 0 tabValues");return r.value}({defaultValue:t,tabValues:i}))),[u,d]=p({queryString:n,groupId:a}),[m,x]=function(e){let{groupId:t}=e;const n=function(e){return e?`docusaurus.tab.${e}`:null}(t),[a,i]=(0,c.Dv)(n);return[a,(0,r.useCallback)((e=>{n&&i.set(e)}),[n,i])]}({groupId:a}),b=(()=>{const e=u??m;return f({value:e,tabValues:i})?e:null})();(0,l.A)((()=>{b&&o(b)}),[b]);return{selectedValue:s,selectValue:(0,r.useCallback)((e=>{if(!f({value:e,tabValues:i}))throw new Error(`Can't select invalid tab value=${e}`);o(e),d(e),x(e)}),[d,x,i]),tabValues:i}}var x=n(2303);const b={tabList:"tabList__CuJ",tabItem:"tabItem_LNqP"};var g=n(4848);function j(e){let{className:t,block:n,selectedValue:r,selectValue:s,tabValues:l}=e;const o=[],{blockElementScrollPositionUntilNextRender:u}=(0,i.a_)(),c=e=>{const t=e.currentTarget,n=o.indexOf(t),a=l[n].value;a!==r&&(u(t),s(a))},d=e=>{let t=null;switch(e.key){case"Enter":c(e);break;case"ArrowRight":{const n=o.indexOf(e.currentTarget)+1;t=o[n]??o[0];break}case"ArrowLeft":{const n=o.indexOf(e.currentTarget)-1;t=o[n]??o[o.length-1];break}}t?.focus()};return(0,g.jsx)("ul",{role:"tablist","aria-orientation":"horizontal",className:(0,a.A)("tabs",{"tabs--block":n},t),children:l.map((e=>{let{value:t,label:n,attributes:i}=e;return(0,g.jsx)("li",{role:"tab",tabIndex:r===t?0:-1,"aria-selected":r===t,ref:e=>o.push(e),onKeyDown:d,onClick:c,...i,className:(0,a.A)("tabs__item",b.tabItem,i?.className,{"tabs__item--active":r===t}),children:n??t},t)}))})}function v(e){let{lazy:t,children:n,selectedValue:a}=e;const i=(Array.isArray(n)?n:[n]).filter(Boolean);if(t){const e=i.find((e=>e.props.value===a));return e?(0,r.cloneElement)(e,{className:"margin-top--md"}):null}return(0,g.jsx)("div",{className:"margin-top--md",children:i.map(((e,t)=>(0,r.cloneElement)(e,{key:t,hidden:e.props.value!==a})))})}function y(e){const t=m(e);return(0,g.jsxs)("div",{className:(0,a.A)("tabs-container",b.tabList),children:[(0,g.jsx)(j,{...t,...e}),(0,g.jsx)(v,{...t,...e})]})}function w(e){const t=(0,x.A)();return(0,g.jsx)(y,{...e,children:d(e.children)},String(t))}},8453:(e,t,n)=>{n.d(t,{R:()=>s,x:()=>l});var r=n(6540);const a={},i=r.createContext(a);function s(e){const t=r.useContext(i);return r.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function l(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(a):e.components||a:s(e.components),r.createElement(i.Provider,{value:t},e.children)}}}]);