"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[4544],{3905:(e,t,n)=>{n.d(t,{Zo:()=>p,kt:()=>d});var r=n(7294);function a(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function i(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var r=Object.getOwnPropertySymbols(e);t&&(r=r.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,r)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?i(Object(n),!0).forEach((function(t){a(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):i(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function l(e,t){if(null==e)return{};var n,r,a=function(e,t){if(null==e)return{};var n,r,a={},i=Object.keys(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||(a[n]=e[n]);return a}(e,t);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);for(r=0;r<i.length;r++)n=i[r],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(a[n]=e[n])}return a}var s=r.createContext({}),c=function(e){var t=r.useContext(s),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},p=function(e){var t=c(e.components);return r.createElement(s.Provider,{value:t},e.children)},u={inlineCode:"code",wrapper:function(e){var t=e.children;return r.createElement(r.Fragment,{},t)}},m=r.forwardRef((function(e,t){var n=e.components,a=e.mdxType,i=e.originalType,s=e.parentName,p=l(e,["components","mdxType","originalType","parentName"]),m=c(n),d=a,h=m["".concat(s,".").concat(d)]||m[d]||u[d]||i;return n?r.createElement(h,o(o({ref:t},p),{},{components:n})):r.createElement(h,o({ref:t},p))}));function d(e,t){var n=arguments,a=t&&t.mdxType;if("string"==typeof e||a){var i=n.length,o=new Array(i);o[0]=m;var l={};for(var s in t)hasOwnProperty.call(t,s)&&(l[s]=t[s]);l.originalType=e,l.mdxType="string"==typeof e?e:a,o[1]=l;for(var c=2;c<i;c++)o[c]=n[c];return r.createElement.apply(null,o)}return r.createElement.apply(null,n)}m.displayName="MDXCreateElement"},4578:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>s,contentTitle:()=>o,default:()=>u,frontMatter:()=>i,metadata:()=>l,toc:()=>c});var r=n(7462),a=(n(7294),n(3905));const i={sidebar_position:2,lastmod:new Date("2022-04-19T20:14:49.735Z")},o="Running ChemEx",l={unversionedId:"user_guide/running_chemex",id:"user_guide/running_chemex",title:"Running ChemEx",description:"ChemEx is a command-line application. After the installation, you can run the",source:"@site/docs/user_guide/running_chemex.md",sourceDirName:"user_guide",slug:"/user_guide/running_chemex",permalink:"/ChemEx/docs/user_guide/running_chemex",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/user_guide/running_chemex.md",tags:[],version:"current",sidebarPosition:2,frontMatter:{sidebar_position:2,lastmod:"2022-04-19T20:14:49.735Z"},sidebar:"tutorialSidebar",previous:{title:"Introduction",permalink:"/ChemEx/docs/user_guide/introduction"},next:{title:"Fitting datasets",permalink:"/ChemEx/docs/user_guide/fitting/"}},s={},c=[],p={toc:c};function u(e){let{components:t,...n}=e;return(0,a.kt)("wrapper",(0,r.Z)({},p,n,{components:t,mdxType:"MDXLayout"}),(0,a.kt)("h1",{id:"running-chemex"},"Running ChemEx"),(0,a.kt)("p",null,"ChemEx is a command-line application. After the installation, you can run the\nprogram directly from the shell prompt using the command ",(0,a.kt)("inlineCode",{parentName:"p"},"chemex"),":"),(0,a.kt)("pre",null,(0,a.kt)("code",{parentName:"pre",className:"language-shell"},"chemex <subcommand> <options>\n")),(0,a.kt)("p",null,"Four ",(0,a.kt)("strong",{parentName:"p"},"sub-commands")," are available, each corresponding to a specific task:"),(0,a.kt)("ul",null,(0,a.kt)("li",{parentName:"ul"},(0,a.kt)("a",{parentName:"li",href:"/ChemEx/docs/user_guide/fitting/chemex_fit"},(0,a.kt)("inlineCode",{parentName:"a"},"fit"))," starts the fits of experimental datasets: if\nyou are interested in ChemEx, this is likely the sub-command you are looking\nfor."),(0,a.kt)("li",{parentName:"ul"},(0,a.kt)("inlineCode",{parentName:"li"},"simulate")," executes the module, which allows to calculate synthetic profiles\nfrom a given experiment and a given model."),(0,a.kt)("li",{parentName:"ul"},(0,a.kt)("inlineCode",{parentName:"li"},"pick_cest")," is a small GUI application, which plots CEST, cos-CEST and D-CEST\nprofiles for interactive dip picking and obtaining starting minor state\nposition values."),(0,a.kt)("li",{parentName:"ul"},(0,a.kt)("inlineCode",{parentName:"li"},"plot_param")," produces plots of selected parameters based on fitting output\nresults.")),(0,a.kt)("admonition",{type:"tip"},(0,a.kt)("p",{parentName:"admonition"},"You can obtain the list of available sub-commands using the ",(0,a.kt)("inlineCode",{parentName:"p"},"--help")," option."),(0,a.kt)("pre",{parentName:"admonition"},(0,a.kt)("code",{parentName:"pre",className:"language-bash"},"chemex --help\n")),(0,a.kt)("p",{parentName:"admonition"},"You can display the list of options associated with each sub-command using the\nsub-command name followed by the ",(0,a.kt)("inlineCode",{parentName:"p"},"--help")," option."),(0,a.kt)("pre",{parentName:"admonition"},(0,a.kt)("code",{parentName:"pre",className:"language-bash"},"chemex <subcommand> --help\n"))))}u.isMDXComponent=!0}}]);