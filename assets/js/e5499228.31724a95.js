"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[7487],{3905:(e,t,n)=>{n.d(t,{Zo:()=>c,kt:()=>f});var i=n(7294);function r(e,t,n){return t in e?Object.defineProperty(e,t,{value:n,enumerable:!0,configurable:!0,writable:!0}):e[t]=n,e}function a(e,t){var n=Object.keys(e);if(Object.getOwnPropertySymbols){var i=Object.getOwnPropertySymbols(e);t&&(i=i.filter((function(t){return Object.getOwnPropertyDescriptor(e,t).enumerable}))),n.push.apply(n,i)}return n}function o(e){for(var t=1;t<arguments.length;t++){var n=null!=arguments[t]?arguments[t]:{};t%2?a(Object(n),!0).forEach((function(t){r(e,t,n[t])})):Object.getOwnPropertyDescriptors?Object.defineProperties(e,Object.getOwnPropertyDescriptors(n)):a(Object(n)).forEach((function(t){Object.defineProperty(e,t,Object.getOwnPropertyDescriptor(n,t))}))}return e}function s(e,t){if(null==e)return{};var n,i,r=function(e,t){if(null==e)return{};var n,i,r={},a=Object.keys(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||(r[n]=e[n]);return r}(e,t);if(Object.getOwnPropertySymbols){var a=Object.getOwnPropertySymbols(e);for(i=0;i<a.length;i++)n=a[i],t.indexOf(n)>=0||Object.prototype.propertyIsEnumerable.call(e,n)&&(r[n]=e[n])}return r}var l=i.createContext({}),m=function(e){var t=i.useContext(l),n=t;return e&&(n="function"==typeof e?e(t):o(o({},t),e)),n},c=function(e){var t=m(e.components);return i.createElement(l.Provider,{value:t},e.children)},p={inlineCode:"code",wrapper:function(e){var t=e.children;return i.createElement(i.Fragment,{},t)}},d=i.forwardRef((function(e,t){var n=e.components,r=e.mdxType,a=e.originalType,l=e.parentName,c=s(e,["components","mdxType","originalType","parentName"]),d=m(n),f=r,h=d["".concat(l,".").concat(f)]||d[f]||p[f]||a;return n?i.createElement(h,o(o({ref:t},c),{},{components:n})):i.createElement(h,o({ref:t},c))}));function f(e,t){var n=arguments,r=t&&t.mdxType;if("string"==typeof e||r){var a=n.length,o=new Array(a);o[0]=d;var s={};for(var l in t)hasOwnProperty.call(t,l)&&(s[l]=t[l]);s.originalType=e,s.mdxType="string"==typeof e?e:r,o[1]=s;for(var m=2;m<a;m++)o[m]=n[m];return i.createElement.apply(null,o)}return i.createElement.apply(null,n)}d.displayName="MDXCreateElement"},1429:(e,t,n)=>{n.r(t),n.d(t,{assets:()=>l,contentTitle:()=>o,default:()=>p,frontMatter:()=>a,metadata:()=>s,toc:()=>m});var i=n(7462),r=(n(7294),n(3905));const a={sidebar_label:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",sidebar_position:1,description:'"shift_15n_sqmq"'},o="\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",s={unversionedId:"experiments/shift/shift_15n_sqmq",id:"experiments/shift/shift_15n_sqmq",title:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",description:'"shift_15n_sqmq"',source:"@site/docs/experiments/shift/shift_15n_sqmq.md",sourceDirName:"experiments/shift",slug:"/experiments/shift/shift_15n_sqmq",permalink:"/ChemEx/docs/experiments/shift/shift_15n_sqmq",draft:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/shift/shift_15n_sqmq.md",tags:[],version:"current",sidebarPosition:1,frontMatter:{sidebar_label:"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC",sidebar_position:1,description:'"shift_15n_sqmq"'},sidebar:"tutorialSidebar",previous:{title:"Exchange-induced chemical shift experiments",permalink:"/ChemEx/docs/experiments/shift/"},next:{title:"Examples",permalink:"/ChemEx/docs/examples/"}},l={},m=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"References",id:"references",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}],c={toc:m};function p(e){let{components:t,...n}=e;return(0,r.kt)("wrapper",(0,i.Z)({},c,n,{components:t,mdxType:"MDXLayout"}),(0,r.kt)("h1",{id:"n-exchange-induced-shifts-with-nh-hsqchmqc"},"\xb9\u2075N exchange induced shifts with \xb9\u2075N\u2013\xb9H HSQC/HMQC"),(0,r.kt)("h2",{id:"module-name"},"Module name"),(0,r.kt)("p",null,(0,r.kt)("inlineCode",{parentName:"p"},'"shift_15n_sqmq"')),(0,r.kt)("h2",{id:"description"},"Description"),(0,r.kt)("p",null,"Analyzes exchange induced \xb9\u2075N chemical shift changes measured in (\xb9\u2075N\u2013\xb9HN) HMQC\nand HSQC data sets."),(0,r.kt)("admonition",{type:"note"},(0,r.kt)("p",{parentName:"admonition"},"Since this experiment is used for determining the sign of \u0394\u03d6, it is usually\ncombined with other CPMG experiments.")),(0,r.kt)("h2",{id:"references"},"References"),(0,r.kt)("ul",null,(0,r.kt)("li",{parentName:"ul"},"N.R. Skrynnikov, F.W. Dahlquist, L.E. Kay. ",(0,r.kt)("em",{parentName:"li"},"J. Am. Chem. Soc.")," ",(0,r.kt)("strong",{parentName:"li"},"124"),",\n12352-12360 (2002)"),(0,r.kt)("li",{parentName:"ul"},"P. Vallurupalli, G. Bouvignies, and L.E. Kay. ",(0,r.kt)("em",{parentName:"li"},"J. Phys. Chem. B")," ",(0,r.kt)("strong",{parentName:"li"},"115"),",\n14891-14900 (2011)")),(0,r.kt)("h2",{id:"example"},"Example"),(0,r.kt)("p",null,"An example use of the module associated with \xb9\u2075N and \xb9H CPMG datasets is given\n",(0,r.kt)("a",{parentName:"p",href:"https://github.com/gbouvignies/chemex/tree/master/examples/Combinations/Shifts/"},"here"),"."),(0,r.kt)("h2",{id:"sample-configuration-file"},"Sample configuration file"),(0,r.kt)("pre",null,(0,r.kt)("code",{parentName:"pre",className:"language-toml",metastring:'title="experiment.toml"',title:'"experiment.toml"'},'## This is a sample configuration file for the module \'shift_15n_sqmq\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "shift_15n_sqmq"\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## 1H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Filename of the file containing the list of the shifts in ppb.\n## The file should be formatted as follow:\n##\n## #     name   shift  error\n##     G2N-HN    10.9    0.5\n##     H3N-HN    32.1    0.5\n##     K4N-HN   -54.3    1.5\n##     S5N-HN     0.7    0.5\n##     L6N-HN   -15.2    0.5\n##\n## The name of the spin systems should follow the Sparky-NMR\n## conventions.\nshifts = "sqmq.txt"\n')))}p.isMDXComponent=!0}}]);