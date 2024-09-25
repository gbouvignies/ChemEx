"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[3397],{2211:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>d,contentTitle:()=>l,default:()=>p,frontMatter:()=>c,metadata:()=>a,toc:()=>m});var r=n(4848),s=n(8453),i=n(3514),o=n(5068);const c={sidebar_label:"Examples",sidebar_position:4},l="Representative Examples",a={id:"examples/index",title:"Representative Examples",description:"This section describes several representative examples in ChemEx. All examples",source:"@site/docs/examples/index.mdx",sourceDirName:"examples",slug:"/examples/",permalink:"/ChemEx/docs/examples/",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/examples/index.mdx",tags:[],version:"current",sidebarPosition:4,frontMatter:{sidebar_label:"Examples",sidebar_position:4},sidebar:"tutorialSidebar",previous:{title:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",permalink:"/ChemEx/docs/experiments/shift/shift_15n_sqmq"},next:{title:"Methyl \xb9H SQ-CPMG",permalink:"/ChemEx/docs/examples/methyl_1h_sqcpmg"}},d={},m=[];function u(e){const t={a:"a",code:"code",h1:"h1",header:"header",p:"p",...(0,s.R)(),...e.components};return(0,r.jsxs)(r.Fragment,{children:[(0,r.jsx)(t.header,{children:(0,r.jsx)(t.h1,{id:"representative-examples",children:"Representative Examples"})}),"\n",(0,r.jsxs)(t.p,{children:["This section describes several representative examples in ChemEx. All examples\nare placed under ",(0,r.jsx)(t.code,{children:"examples/"})," directory, each experiment module (see\n",(0,r.jsx)(t.a,{href:"../experiments/",children:(0,r.jsx)(t.code,{children:"Experiments"})})," section) typically has at least one\nexample to demonstrate its application."]}),"\n","\n","\n",(0,r.jsx)(i.A,{items:(0,o.$S)().items})]})}function p(e={}){const{wrapper:t}={...(0,s.R)(),...e.components};return t?(0,r.jsx)(t,{...e,children:(0,r.jsx)(u,{...e})}):u(e)}},3514:(e,t,n)=>{n.d(t,{A:()=>b});n(6540);var r=n(4164),s=n(6972),i=n(8774),o=n(5846),c=n(6654),l=n(1312),a=n(1107);const d={cardContainer:"cardContainer_fWXF",cardTitle:"cardTitle_rnsV",cardDescription:"cardDescription_PWke"};var m=n(4848);function u(e){let{href:t,children:n}=e;return(0,m.jsx)(i.A,{href:t,className:(0,r.A)("card padding--lg",d.cardContainer),children:n})}function p(e){let{href:t,icon:n,title:s,description:i}=e;return(0,m.jsxs)(u,{href:t,children:[(0,m.jsxs)(a.A,{as:"h2",className:(0,r.A)("text--truncate",d.cardTitle),title:s,children:[n," ",s]}),i&&(0,m.jsx)("p",{className:(0,r.A)("text--truncate",d.cardDescription),title:i,children:i})]})}function h(e){let{item:t}=e;const n=(0,s.Nr)(t),r=function(){const{selectMessage:e}=(0,o.W)();return t=>e(t,(0,l.T)({message:"1 item|{count} items",id:"theme.docs.DocCard.categoryDescription.plurals",description:"The default description for a category card in the generated index about how many items this category includes"},{count:t}))}();return n?(0,m.jsx)(p,{href:n,icon:"\ud83d\uddc3\ufe0f",title:t.label,description:t.description??r(t.items.length)}):null}function x(e){let{item:t}=e;const n=(0,c.A)(t.href)?"\ud83d\udcc4\ufe0f":"\ud83d\udd17",r=(0,s.cC)(t.docId??void 0);return(0,m.jsx)(p,{href:t.href,icon:n,title:t.label,description:t.description??r?.description})}function f(e){let{item:t}=e;switch(t.type){case"link":return(0,m.jsx)(x,{item:t});case"category":return(0,m.jsx)(h,{item:t});default:throw new Error(`unknown item type ${JSON.stringify(t)}`)}}function g(e){let{className:t}=e;const n=(0,s.$S)();return(0,m.jsx)(b,{items:n.items,className:t})}function b(e){const{items:t,className:n}=e;if(!t)return(0,m.jsx)(g,{...e});const i=(0,s.d1)(t);return(0,m.jsx)("section",{className:(0,r.A)("row",n),children:i.map(((e,t)=>(0,m.jsx)("article",{className:"col col--6 margin-bottom--lg",children:(0,m.jsx)(f,{item:e})},t)))})}},5068:(e,t,n)=>{n.d(t,{$S:()=>r});n(4586);function r(){return n(4070).$S(...arguments)}},5846:(e,t,n)=>{n.d(t,{W:()=>a});var r=n(6540),s=n(4586);const i=["zero","one","two","few","many","other"];function o(e){return i.filter((t=>e.includes(t)))}const c={locale:"en",pluralForms:o(["one","other"]),select:e=>1===e?"one":"other"};function l(){const{i18n:{currentLocale:e}}=(0,s.A)();return(0,r.useMemo)((()=>{try{return function(e){const t=new Intl.PluralRules(e);return{locale:e,pluralForms:o(t.resolvedOptions().pluralCategories),select:e=>t.select(e)}}(e)}catch(t){return console.error(`Failed to use Intl.PluralRules for locale "${e}".\nDocusaurus will fallback to the default (English) implementation.\nError: ${t.message}\n`),c}}),[e])}function a(){const e=l();return{selectMessage:(t,n)=>function(e,t,n){const r=e.split("|");if(1===r.length)return r[0];r.length>n.pluralForms.length&&console.error(`For locale=${n.locale}, a maximum of ${n.pluralForms.length} plural forms are expected (${n.pluralForms.join(",")}), but the message contains ${r.length}: ${e}`);const s=n.select(t),i=n.pluralForms.indexOf(s);return r[Math.min(i,r.length-1)]}(n,t,e)}}},8453:(e,t,n)=>{n.d(t,{R:()=>o,x:()=>c});var r=n(6540);const s={},i=r.createContext(s);function o(e){const t=r.useContext(i);return r.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function c(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(s):e.components||s:o(e.components),r.createElement(i.Provider,{value:t},e.children)}}}]);