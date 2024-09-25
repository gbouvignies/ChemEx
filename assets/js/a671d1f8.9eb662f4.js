"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[6209],{4259:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>d,contentTitle:()=>a,default:()=>m,frontMatter:()=>o,metadata:()=>l,toc:()=>u});var i=n(4848),r=n(8453),s=n(3514),c=n(5068);const o={sidebar_label:"Exchange-induced chemical shift experiments",sidebar_position:5,description:"Shift"},a="Exchange-induced chemical shift experiments",l={id:"experiments/shift/index",title:"Exchange-induced chemical shift experiments",description:"Shift",source:"@site/docs/experiments/shift/index.mdx",sourceDirName:"experiments/shift",slug:"/experiments/shift/",permalink:"/ChemEx/docs/experiments/shift/",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/shift/index.mdx",tags:[],version:"current",sidebarPosition:5,frontMatter:{sidebar_label:"Exchange-induced chemical shift experiments",sidebar_position:5,description:"Shift"},sidebar:"tutorialSidebar",previous:{title:"Amide \xb9H-\xb9\u2075N longitudinal two-spin order relaxation",permalink:"/ChemEx/docs/experiments/relaxation/relaxation_hznz"},next:{title:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",permalink:"/ChemEx/docs/experiments/shift/shift_15n_sqmq"}},d={},u=[];function h(e){const t={h1:"h1",header:"header",p:"p",...(0,r.R)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(t.header,{children:(0,i.jsx)(t.h1,{id:"exchange-induced-chemical-shift-experiments",children:"Exchange-induced chemical shift experiments"})}),"\n",(0,i.jsx)(t.p,{children:"This section contains information about the exchange-induced chemical shift\nexperiments that are available in ChemEx."}),"\n","\n","\n",(0,i.jsx)(s.A,{items:(0,c.$S)().items})]})}function m(e={}){const{wrapper:t}={...(0,r.R)(),...e.components};return t?(0,i.jsx)(t,{...e,children:(0,i.jsx)(h,{...e})}):h(e)}},3514:(e,t,n)=>{n.d(t,{A:()=>b});n(6540);var i=n(4164),r=n(6972),s=n(8774),c=n(5846),o=n(6654),a=n(1312),l=n(1107);const d={cardContainer:"cardContainer_fWXF",cardTitle:"cardTitle_rnsV",cardDescription:"cardDescription_PWke"};var u=n(4848);function h(e){let{href:t,children:n}=e;return(0,u.jsx)(s.A,{href:t,className:(0,i.A)("card padding--lg",d.cardContainer),children:n})}function m(e){let{href:t,icon:n,title:r,description:s}=e;return(0,u.jsxs)(h,{href:t,children:[(0,u.jsxs)(l.A,{as:"h2",className:(0,i.A)("text--truncate",d.cardTitle),title:r,children:[n," ",r]}),s&&(0,u.jsx)("p",{className:(0,i.A)("text--truncate",d.cardDescription),title:s,children:s})]})}function f(e){let{item:t}=e;const n=(0,r.Nr)(t),i=function(){const{selectMessage:e}=(0,c.W)();return t=>e(t,(0,a.T)({message:"1 item|{count} items",id:"theme.docs.DocCard.categoryDescription.plurals",description:"The default description for a category card in the generated index about how many items this category includes"},{count:t}))}();return n?(0,u.jsx)(m,{href:n,icon:"\ud83d\uddc3\ufe0f",title:t.label,description:t.description??i(t.items.length)}):null}function p(e){let{item:t}=e;const n=(0,o.A)(t.href)?"\ud83d\udcc4\ufe0f":"\ud83d\udd17",i=(0,r.cC)(t.docId??void 0);return(0,u.jsx)(m,{href:t.href,icon:n,title:t.label,description:t.description??i?.description})}function x(e){let{item:t}=e;switch(t.type){case"link":return(0,u.jsx)(p,{item:t});case"category":return(0,u.jsx)(f,{item:t});default:throw new Error(`unknown item type ${JSON.stringify(t)}`)}}function g(e){let{className:t}=e;const n=(0,r.$S)();return(0,u.jsx)(b,{items:n.items,className:t})}function b(e){const{items:t,className:n}=e;if(!t)return(0,u.jsx)(g,{...e});const s=(0,r.d1)(t);return(0,u.jsx)("section",{className:(0,i.A)("row",n),children:s.map(((e,t)=>(0,u.jsx)("article",{className:"col col--6 margin-bottom--lg",children:(0,u.jsx)(x,{item:e})},t)))})}},5068:(e,t,n)=>{n.d(t,{$S:()=>i});n(4586);function i(){return n(4070).$S(...arguments)}},5846:(e,t,n)=>{n.d(t,{W:()=>l});var i=n(6540),r=n(4586);const s=["zero","one","two","few","many","other"];function c(e){return s.filter((t=>e.includes(t)))}const o={locale:"en",pluralForms:c(["one","other"]),select:e=>1===e?"one":"other"};function a(){const{i18n:{currentLocale:e}}=(0,r.A)();return(0,i.useMemo)((()=>{try{return function(e){const t=new Intl.PluralRules(e);return{locale:e,pluralForms:c(t.resolvedOptions().pluralCategories),select:e=>t.select(e)}}(e)}catch(t){return console.error(`Failed to use Intl.PluralRules for locale "${e}".\nDocusaurus will fallback to the default (English) implementation.\nError: ${t.message}\n`),o}}),[e])}function l(){const e=a();return{selectMessage:(t,n)=>function(e,t,n){const i=e.split("|");if(1===i.length)return i[0];i.length>n.pluralForms.length&&console.error(`For locale=${n.locale}, a maximum of ${n.pluralForms.length} plural forms are expected (${n.pluralForms.join(",")}), but the message contains ${i.length}: ${e}`);const r=n.select(t),s=n.pluralForms.indexOf(r);return i[Math.min(s,i.length-1)]}(n,t,e)}}},8453:(e,t,n)=>{n.d(t,{R:()=>c,x:()=>o});var i=n(6540);const r={},s=i.createContext(r);function c(e){const t=i.useContext(s);return i.useMemo((function(){return"function"==typeof e?e(t):{...t,...e}}),[t,e])}function o(e){let t;return t=e.disableParentContext?"function"==typeof e.components?e.components(r):e.components||r:c(e.components),i.createElement(s.Provider,{value:t},e.children)}}}]);