"use strict";(self.webpackChunkwebsite=self.webpackChunkwebsite||[]).push([[3250],{2664:(e,n,t)=>{t.r(n),t.d(n,{assets:()=>l,contentTitle:()=>o,default:()=>p,frontMatter:()=>s,metadata:()=>r,toc:()=>c});var i=t(5893),a=t(1151);const s={sidebar_label:"In-phase/anti-phase amide \xb9H cos-CEST",sidebar_position:3,description:'"coscest_1hn_ip_ap"'},o="In-phase/anti-phase amide \xb9H cos-CEST",r={id:"experiments/dcest/coscest_1hn_ip_ap",title:"In-phase/anti-phase amide \xb9H cos-CEST",description:'"coscest_1hn_ip_ap"',source:"@site/docs/experiments/dcest/coscest_1hn_ip_ap.md",sourceDirName:"experiments/dcest",slug:"/experiments/dcest/coscest_1hn_ip_ap",permalink:"/ChemEx/docs/experiments/dcest/coscest_1hn_ip_ap",draft:!1,unlisted:!1,editUrl:"https://github.com/gbouvignies/chemex/tree/main/website/docs/experiments/dcest/coscest_1hn_ip_ap.md",tags:[],version:"current",sidebarPosition:3,frontMatter:{sidebar_label:"In-phase/anti-phase amide \xb9H cos-CEST",sidebar_position:3,description:'"coscest_1hn_ip_ap"'},sidebar:"tutorialSidebar",previous:{title:"Pure in-phase \xb9\xb3C D-CEST",permalink:"/ChemEx/docs/experiments/dcest/dcest_13c"},next:{title:"CPMG",permalink:"/ChemEx/docs/experiments/cpmg/"}},l={},c=[{value:"Module name",id:"module-name",level:2},{value:"Description",id:"description",level:2},{value:"References",id:"references",level:2},{value:"Example",id:"example",level:2},{value:"Sample configuration file",id:"sample-configuration-file",level:2}];function d(e){const n={a:"a",code:"code",em:"em",h1:"h1",h2:"h2",li:"li",p:"p",pre:"pre",strong:"strong",ul:"ul",...(0,a.a)(),...e.components};return(0,i.jsxs)(i.Fragment,{children:[(0,i.jsx)(n.h1,{id:"in-phaseanti-phase-amide-h-cos-cest",children:"In-phase/anti-phase amide \xb9H cos-CEST"}),"\n",(0,i.jsx)(n.h2,{id:"module-name",children:"Module name"}),"\n",(0,i.jsx)(n.p,{children:(0,i.jsx)(n.code,{children:'"coscest_1hn_ip_ap"'})}),"\n",(0,i.jsx)(n.h2,{id:"description",children:"Description"}),"\n",(0,i.jsx)(n.p,{children:"Analyzes chemical exchange during the COS-CEST block. Magnetization evolution is\ncalculated using the (6n)\xd7(6n), two-spin matrix, where n is the number of\nstates:"}),"\n",(0,i.jsx)(n.p,{children:"{ Ix(a), Iy(a), Iz(a), IxSz(a), IySz(a), IzSz(a),\nIx(b), Iy(b), Iz(b), IxSz(b), IySz(b), IzSz(b), ... }"}),"\n",(0,i.jsx)(n.h2,{id:"references",children:"References"}),"\n",(0,i.jsxs)(n.ul,{children:["\n",(0,i.jsxs)(n.li,{children:["T. Yuwen, G. Bouvignies, and L.E. Kay. ",(0,i.jsx)(n.em,{children:"J. Mag. Reson."})," ",(0,i.jsx)(n.strong,{children:"292"}),", 1-7 (2018)"]}),"\n"]}),"\n",(0,i.jsx)(n.h2,{id:"example",children:"Example"}),"\n",(0,i.jsxs)(n.p,{children:["An example for studying \xb9\u2075N-labeled sample is given\n",(0,i.jsx)(n.a,{href:"https://github.com/gbouvignies/chemex/tree/master/examples/Experiments/COSCEST_1HN_IP_AP/",children:"here"}),"."]}),"\n",(0,i.jsx)(n.h2,{id:"sample-configuration-file",children:"Sample configuration file"}),"\n",(0,i.jsx)(n.pre,{children:(0,i.jsx)(n.code,{className:"language-toml",metastring:'title="experiment.toml"',children:'## This is a sample configuration file for the module \'coscest_1hn_ip_ap\'\n\n[experiment]\n\n## Name of the chemex module corresponding to the experiment\nname = "coscest_1hn_ip_ap"\n\n## Recycle delay, in seconds\nd1 = 0.1\n\n## CEST relaxation delay, in seconds\ntime_t1 = 0.5\n\n## Position of the carrier during the CEST period, in ppm\ncarrier = 8.3\n\n## Cosine "spectral width", in Hz\nsw = 800.0\n\n## Number of excitation frequencies\ncos_n = 3\n\n## Number of points used to simulate the cosine-modulated shape\n# cos_res = 10\n\n## B1 radio-frequency field strength, in Hz\nb1_frq = 25.0\n\n## B1 inhomogeneity expressed as a fraction of \'b1_frq\'. If set to "inf",\n## a faster calculation takes place assuming full dephasing of the\n## magnetization components orthogonal to the effective field. The "inf" value\n## should not be used with an "eta_block" value larger than 0.\n## [optional, default: 0.1]\n# b1_inh_scale = 0.1\n\n## Number of points used to simulate B1 inhomogeneity, the larger\n## the longer the calculation. [optional, default: 11]\n# b1_inh_res = 11\n\n## State of the observed resonance [optional, default: "a"]\n# observed_state = "a"\n\n[conditions]\n\n## \xb9H Larmor frequency, in MHz\nh_larmor_frq = 800.0\n\n## Sample temperature, in Celsius [optional, depending on the kinetic model]\n# temperature = 25.0\n\n## Protein concentration, in M [optional, depending on the kinetic model]\n# p_total = 500.0e-6\n\n## Ligand concentration, in M [optional, depending on the kinetic model]\n# l_total = 50.0e-6\n\n## Labeling scheme of the sample, for deuterated samples "2H" should\n## be used to obtain accurate initial estimates of relaxation rates\n## based on model-free parameters [optional, default: []]\n# label = ["2H"]\n\n[data]\n\n## Directory containing the profiles [optional, default: "./"]\n# path = "./"\n\n## Option defining how intensity uncertainties are estimated.\n## "file": uncertainties are taken from the profile files\n## "scatter": uncertainties are calculated from the baseline\n## [optional, default: "file"]\n# error = "file"\n\n## List of offsets relative to the main resonance position\n## (nu) and bandwidths (delta_nu) defining regions where\n## points are excluded from the calculation (nu +/- 0.5 * delta_nu),\n## both are in Hz [optional, default: [[0.0, 0.0]] ]\n# filter_offsets = [\n#   [0.0, 0.0],\n# ]\n\n## List of plane indices to be excluded from the calculation.\n## The first plane has index 0 and is usually the reference plane,\n## this plane is always excluded by default [optional, default: [] ]\n# filter_planes = []\n\n  ## List of the profile names and their associated filenames.\n  ## The name of the spin systems should follow the Sparky-NMR\n  ## conventions.\n  [data.profiles]\n  G2HN-N = "G2N-HN.out"\n  H3HN-N = "H3N-HN.out"\n  K4HN-N = "K4N-HN.out"\n  S5HN-N = "S5N-HN.out"\n  L6HN-N = "L6N-HN.out"\n'})})]})}function p(e={}){const{wrapper:n}={...(0,a.a)(),...e.components};return n?(0,i.jsx)(n,{...e,children:(0,i.jsx)(d,{...e})}):d(e)}},1151:(e,n,t)=>{t.d(n,{Z:()=>r,a:()=>o});var i=t(7294);const a={},s=i.createContext(a);function o(e){const n=i.useContext(s);return i.useMemo((function(){return"function"==typeof e?e(n):{...n,...e}}),[n,e])}function r(e){let n;return n=e.disableParentContext?"function"==typeof e.components?e.components(a):e.components||a:o(e.components),i.createElement(s.Provider,{value:n},e.children)}}}]);