"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[1802],{3905:(e,t,r)=>{r.d(t,{Zo:()=>m,kt:()=>d});var n=r(7294);function a(e,t,r){return t in e?Object.defineProperty(e,t,{value:r,enumerable:!0,configurable:!0,writable:!0}):e[t]=r,e}function i(e,t){var r=Object.keys(e);if(Object.getOwnPropertySymbols){var n=Object.getOwnPropertySymbols(e);t&&(n=n.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),r.push.apply(r,n)}return r}function c(e){for(var t=1;t<arguments.length;t++){var r=null!=arguments[t]?arguments[t]:{};t%2?i(Object(r),!0).forEach((function(t){a(e,t,r[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(r)):i(Object(r)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(r,t))}))}return e}function s(e,t){if(null==e)return{};var r,n,a=function(e,t){if(null==e)return{};var r,n,a={},i=Object.keys(e);for(n=0;n<i.length;n++)r=i[n],t.indexOf(r)>=0||(a[r]=e[r]);return a}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(n=0;n<i.length;n++)r=i[n],t.indexOf(r)>=0||Object.prototype.propertyIsEnumerable.call(e,r)&&(a[r]=e[r])}return a}var o=n.createContext({}),l=function(e){var t=n.useContext(o),r=t;return e&&(r="function"==typeof e?e(t):c(c({},t),e)),r},m=function(e){var t=l(e.components);return n.createElement(o.Provider,{value:t},e.children)},p={inlineCode:"code",wrapper:function(e){var t=e.children;return n.createElement(n.Fragment,{},t)}},u=n.forwardRef((function(e,t){var r=e.components,a=e.mdxType,i=e.originalType,o=e.parentName,m=s(e,["components","mdxType","originalType","parentName"]),u=l(r),d=a,f=u["".concat(o,".").concat(d)]||u[d]||p[d]||i;return r?n.createElement(f,c(c({ref:t},m),{},{components:r})):n.createElement(f,c({ref:t},m))}));function d(e,t){var r=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var i=r.length,c=new Array(i);c[0]=u;var s={};for(var o in t)hasOwnProperty.call(t,o)&&(s[o]=t[o]);s.originalType=e,s.mdxType="string"==typeof e?e:a,c[1]=s;for(var l=2;l<i;l++)c[l]=r[l];return n.createElement.apply(null,c)}return n.createElement.apply(null,r)}u.displayName="MDXCreateElement"},2991:(e,t,r)=>{r.d(t,{Z:()=>b});var n=r(7294),a=r(6010),i=r(3438),c=r(9960),s=r(3919),o=r(5999);const l="cardContainer_fWXF",m="cardTitle_rnsV",p="cardDescription_PWke";function u(e){let{href:t,children:r}=e;return n.createElement(c.Z,{href:t,className:(0,a.Z)("card padding--lg",l)},r)}function d(e){let{href:t,icon:r,title:i,description:c}=e;return n.createElement(u,{href:t},n.createElement("h2",{className:(0,a.Z)("text--truncate",m),title:i},r," ",i),c&&n.createElement("p",{className:(0,a.Z)("text--truncate",p),title:c},c))}function f(e){let{item:t}=e;const r=(0,i.Wl)(t);return r?n.createElement(d,{href:r,icon:"\ud83d\uddc3\ufe0f",title:t.label,description:(0,o.I)({message:"{count} items",id:"theme.docs.DocCard.categoryDescription",description:"The default description for a category card in the generated index about how many items this category includes"},{count:t.items.length})}):null}function h(e){var t;let{item:r}=e;const a=(0,s.Z)(r.href)?"\ud83d\udcc4\ufe0f":"\ud83d\udd17",c=(0,i.xz)(null!=(t=r.docId)?t:void 0);return n.createElement(d,{href:r.href,icon:a,title:r.label,description:null==c?void 0:c.description})}function y(e){let{item:t}=e;switch(t.type){case"link":return n.createElement(h,{item:t});case"category":return n.createElement(f,{item:t});default:throw new Error("unknown item type "+JSON.stringify(t))}}function x(e){let{className:t}=e;const r=(0,i.jA)();return n.createElement(b,{items:r.items,className:t})}function b(e){const{items:t,className:r}=e;if(!t)return n.createElement(x,e);const c=(0,i.MN)(t);return n.createElement("section",{className:(0,a.Z)("row",r)},c.map(((e,t)=>n.createElement("article",{key:t,className:"col col--6 margin-bottom--lg"},n.createElement(y,{item:e})))))}},5849:(e,t,r)=>{r.r(t),r.d(t,{assets:()=>m,contentTitle:()=>o,default:()=>d,frontMatter:()=>s,metadata:()=>l,toc:()=>p});var n=r(7462),a=(r(7294),r(3905)),i=r(2991),c=r(3438);const s={sidebar_label:"Examples",sidebar_position:4},o="Representative Examples",l={unversionedId:"examples/index",id:"examples/index",title:"Representative Examples",description:"This section describes several representative examples in ChemEx. All examples",source:"@site/docs/examples/index.md",sourceDirName:"examples",slug:"/examples/",permalink:"/chemex/docs/examples/",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/master/docs/examples/index.md",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_label:"Examples",sidebar_position:4},sidebar:"tutorialSidebar",previous:{title:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",permalink:"/chemex/docs/experiments/shift/shift_15n_sqmq"},next:{title:"Methyl \xb9H SQ-CPMG",permalink:"/chemex/docs/examples/methyl_1h_sqcpmg"}},m={},p=[],u={toc:p};function d(e){let{components:t,...r}=e;return(0,a.kt)("wrapper",(0,n.Z)({},u,r,{components:t,mdxType:"MDXLayout"}),(0,a.kt)("h1",{id:"representative-examples"},"Representative Examples"),(0,a.kt)("p",null,"This section describes several representative examples in ChemEx. All examples\nare placed under ",(0,a.kt)("inlineCode",{parentName:"p"},"examples/")," directory, each experiment module (see\n",(0,a.kt)("a",{parentName:"p",href:"../experiments/"},(0,a.kt)("inlineCode",{parentName:"a"},"Experiments"))," section) typically has at least one\nexample to demonstrate its application."),(0,a.kt)(i.Z,{items:(0,c.jA)().items,mdxType:"DocCardList"}))}d.isMDXComponent=!0}}]);